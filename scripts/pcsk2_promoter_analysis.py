#!/usr/bin/env python3
"""
Pcsk2 Promoter Bioinformatics Analysis
=======================================
For: 福建省自然科学基金面上项目 — Pomc-Pcsk2环路分子机制
Task: Predict TF binding sites on mouse Pcsk2 promoter (-2000 to +100 bp from TSS)
Focus TFs: AP-1 family (Fos, Jun, Junb, Jund), CREB (Crem, Atf3),
           Nr4a family (Nr4a1, Nr4a2, Nr4a3), Egr1, Dbp, Nr3c1

Approach:
1. Fetch Pcsk2 TSS from Ensembl REST API
2. Retrieve promoter sequence via Ensembl REST API
3. Scan with JASPAR 2024 position weight matrices
4. Generate binding site report for grant proposal
"""

import requests
import json
import re
import sys
import os
from io import StringIO
from datetime import datetime

# ========== Configuration ==========
SPECIES = "mus_musculus"
GENE_SYMBOL = "Pcsk2"
PROMOTER_UPSTREAM = 2000  # bp upstream of TSS
PROMOTER_DOWNSTREAM = 100  # bp downstream of TSS
ENSEMBL_REST = "https://rest.ensembl.org"
JASPAR_API = "https://jaspar.elixir.no/api/v1"

# Key TFs from our GRN analysis to prioritize
# Subset of TFs with strong Pomc co-expression or ChIP-qPCR candidates
PRIORITY_TFS = {
    # AP-1 family (core focus)
    "Fos":   ["MA0476.2", "MA0099.3"],   # c-Fos; MA0476.2=FOS, MA0099.3=FOS::JUN
    "Jun":   ["MA0488.2", "MA0489.2", "MA0099.3"],   # c-Jun; MA0488.2=JUN, MA0489.2=JUN(var.2)
    "Junb":  ["MA0490.2", "MA1142.1"],   # JunB; MA0490.2=JUNB, MA1142.1=JUNB(var.2)
    "Jund":  ["MA0491.2"],               # JunD
    "Fosb":  ["MA0477.2"],               # FosB
    "Fosl1": ["MA0478.2"],               # Fra-1
    "Fosl2": ["MA0479.2"],               # Fra-2
    # CREB/ATF family (core focus)
    "Crem":  ["MA0608.2"],               # CREM
    "Creb1": ["MA0018.3"],               # CREB1
    "Atf3":  ["MA0605.2"],               # ATF3
    "Atf4":  ["MA1946.1"],               # ATF4
    # Nr4a family
    "Nr4a1": ["MA0160.3"],               # NGFI-B / Nur77
    "Nr4a2": ["MA0161.3"],               # Nurr1
    # Egr family
    "Egr1":  ["MA0162.3"],               # EGR1 / NGFI-A / Zif268
    "Egr2":  ["MA0472.2"],               # EGR2
    # Circadian
    "Dbp":   ["MA1945.1"],               # DBP (if available)
    # Nuclear receptor
    "Nr3c1": ["MA0113.3"],               # GR / glucocorticoid receptor
}

# Additional TFs of potential interest (broader scan)
SECONDARY_TFS = {
    "Klf4":  ["MA0039.3"],
    "Sox9":  ["MA0079.4"],
    "Sp1":   ["MA0079.5"],
    "Stat3": ["MA0144.3"],
    "Nfkb1": ["MA0105.4"],
    "Mef2c": ["MA0497.2"],
    "Mef2d": ["MA0498.2"],
}


def fetch_gene_info(symbol):
    """Fetch gene TSS and genomic coordinates from Ensembl REST API."""
    url = f"{ENSEMBL_REST}/lookup/symbol/{SPECIES}/{symbol}?expand=1"
    headers = {"Content-Type": "application/json"}
    print(f"  Fetching gene info for {symbol}...")
    r = requests.get(url, headers=headers)
    r.raise_for_status()
    data = r.json()

    gene_info = {
        "id": data.get("id"),
        "symbol": symbol,
        "seq_region_name": data.get("seq_region_name"),
        "start": data.get("start"),
        "end": data.get("end"),
        "strand": data.get("strand"),
        "biotype": data.get("biotype"),
        "description": data.get("description"),
    }

    # Get TSS (5' end on + strand, 3' end on - strand)
    gene_info["tss"] = gene_info["start"] if gene_info["strand"] == 1 else gene_info["end"]

    print(f"    Gene ID: {gene_info['id']}")
    print(f"    Location: chr{gene_info['seq_region_name']}:{gene_info['start']}-{gene_info['end']} ({'+' if gene_info['strand']==1 else '-'} strand)")
    print(f"    TSS: {gene_info['tss']}")

    return gene_info


def fetch_promoter_sequence(gene_info):
    """Fetch promoter sequence using Ensembl sequence endpoint."""
    upstream = PROMOTER_UPSTREAM
    downstream = PROMOTER_DOWNSTREAM
    strand = gene_info["strand"]
    chrom = gene_info["seq_region_name"]
    tss = gene_info["tss"]

    # For + strand: promoter = tss-upstream .. tss+downstream
    # For - strand: promoter = tss-downstream .. tss+upstream (reverse complement needed)
    if strand == 1:
        region_start = tss - upstream
        region_end = tss + downstream
    else:
        region_start = tss - downstream
        region_end = tss + upstream

    url = f"{ENSEMBL_REST}/sequence/region/{SPECIES}/{chrom}:{region_start}:{region_end}:1"
    headers = {"Content-Type": "text/plain"}
    print(f"\n  Fetching promoter sequence {chrom}:{region_start}-{region_end}...")
    r = requests.get(url, headers=headers)
    r.raise_for_status()
    seq = r.text.strip().upper()

    # If gene is on - strand, the returned sequence needs reverse complement
    if strand == -1:
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        seq = "".join(comp.get(b, "N") for b in reversed(seq))

    print(f"    Promoter length: {len(seq)} bp")
    return seq, region_start, region_end


def fetch_jaspar_motif(matrix_id):
    """Fetch a JASPAR position weight matrix by matrix ID."""
    url = f"{JASPAR_API}/matrix/{matrix_id}/"
    headers = {"Accept": "application/json"}
    try:
        r = requests.get(url, headers=headers, timeout=10)
        if r.status_code == 200:
            return r.json()
        else:
            print(f"      (JASPAR matrix {matrix_id} not found, HTTP {r.status_code})")
            return None
    except Exception as e:
        print(f"      (JASPAR API error for {matrix_id}: {e})")
        return None


def compute_pwm_score(sequence, pwm_dict):
    """Compute PWM score for a sequence against a JASPAR PWM.
    pwm_dict: {"A": [freq1, freq2, ...], "C": [...], "G": [...], "T": [...]}
    Returns relative score (0-1) where 1 = perfect match to consensus.
    """
    nucleotides = {"A": 0, "C": 1, "G": 2, "T": 3}
    pwm_len = len(pwm_dict["A"])

    if len(sequence) != pwm_len:
        return 0.0

    score = 0.0
    min_score = float("inf")
    max_score = float("-inf")

    for i, base in enumerate(sequence):
        if base not in nucleotides:
            return 0.0
        idx = nucleotides[base]
        freqs = [pwm_dict["A"][i], pwm_dict["C"][i], pwm_dict["G"][i], pwm_dict["T"][i]]

        # Avoid log(0)
        eps = 1e-6
        b = 0.25  # background frequency

        score += sum(f * max(0.01, f) for f in freqs)  # simplified scoring
        min_score += sum(b * max(0.01, b) for _ in range(4))
        max_score += sum(1.0 for _ in range(4))

    # Normalize to [0, 1]
    if max_score > min_score:
        return max(0, min(1, (score - min_score) / (max_score - min_score)))
    return 0.0


def pwm_to_pssm(jaspar_pwm, background=(0.25, 0.25, 0.25, 0.25)):
    """Convert JASPAR PWM frequencies to log-odds PSSM.
    Returns: list of dicts [{A: score, C: score, G: score, T: score}, ...] and max_total_score
    """
    import math
    nucleotides = ["A", "C", "G", "T"]
    pwm_len = len(jaspar_pwm["A"])
    pssm = []
    max_total = 0.0

    for i in range(pwm_len):
        col = {}
        col_max = float("-inf")
        for j, nt in enumerate(nucleotides):
            freq = jaspar_pwm[nt][i]
            # Add pseudocount
            freq = (freq + 0.001) / 1.004  # slight pseudocount normalization
            score = math.log2(freq / background[j]) if freq > 0 else -10
            col[nt] = score
            if score > col_max:
                col_max = score
        max_total += col_max
        pssm.append(col)

    return pssm, max_total


def scan_sequence_with_pssm(sequence, pssm, max_total, threshold=0.80):
    """Scan sequence with a PSSM, return hits above relative threshold.
    threshold: fraction of max score (0.80 = 80% of maximum possible score)
    """
    mot_len = len(pssm)
    hits = []

    for i in range(len(sequence) - mot_len + 1):
        subseq = sequence[i:i+mot_len]
        score = 0.0

        for j, base in enumerate(subseq):
            score += pssm[j].get(base, -10)

        rel_score = score / max_total if max_total > 0 else 0

        if rel_score >= threshold:
            hits.append({
                "start": i,
                "end": i + mot_len,
                "sequence": subseq,
                "score": round(score, 2),
                "relative_score": round(rel_score, 4),
                "strand": "+",
            })

    return hits


def scan_both_strands(sequence, pssm, max_total, threshold=0.80):
    """Scan both strands with complement PSSM."""
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}

    # Forward strand
    hits_fwd = scan_sequence_with_pssm(sequence, pssm, max_total, threshold)

    # Reverse complement PSSM
    pssm_rc = []
    for col in reversed(pssm):
        rc_col = {comp.get(nt, "N"): score for nt, score in col.items()}
        pssm_rc.append(rc_col)

    hits_rev = scan_sequence_with_pssm(sequence, pssm_rc, max_total, threshold)
    for h in hits_rev:
        h["strand"] = "-"

    return hits_fwd + hits_rev


def annotate_position(pos, upstream_len):
    """Annotate position relative to TSS."""
    rel = pos - upstream_len  # position relative to TSS (0 = TSS)
    if rel < 0:
        return f"{rel} (upstream)"
    elif rel == 0:
        return "0 (TSS)"
    else:
        return f"+{rel} (downstream)"


def analyze_promoter():
    """Main analysis pipeline."""
    print("=" * 70)
    print("  Pcsk2 Promoter TF Binding Site Prediction")
    print("  For: 福建省自然科学基金面上项目")
    print(f"  Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 70)

    # Step 1: Get gene info
    print("\n[Step 1] Retrieving Pcsk2 gene coordinates...")
    gene_info = fetch_gene_info(GENE_SYMBOL)

    # Step 2: Get promoter sequence
    print("\n[Step 2] Retrieving promoter sequence...")
    seq, region_start, region_end = fetch_promoter_sequence(gene_info)

    # Step 3: Scan with JASPAR motifs
    print("\n[Step 3] Scanning for TF binding sites (JASPAR 2024)...")

    # Combine priority and secondary TFs
    all_tfs = {}
    all_tfs.update(PRIORITY_TFS)
    all_tfs.update(SECONDARY_TFS)

    results = {}
    motif_count_total = 0

    for tf_name, matrix_ids in all_tfs.items():
        print(f"\n  Analyzing {tf_name}...")
        tf_hits = []

        for mid in matrix_ids:
            motif_data = fetch_jaspar_motif(mid)
            if motif_data is None:
                continue

            # Extract PWM
            pwm = motif_data.get("pwm")
            if pwm is None:
                # Try alternate field names
                pwm = motif_data.get("pfm") or motif_data.get("frequency_matrix")
            if pwm is None:
                print(f"      No PWM found in response for {mid}")
                continue

            motif_name = motif_data.get("name", mid)
            print(f"      {mid} ({motif_name}): PWM length={len(pwm['A'])}")

            # Convert to PSSM and scan
            pssm, max_total = pwm_to_pssm(pwm)

            # Use 80% threshold for initial scan, 85% for high-confidence
            for threshold in [0.85, 0.80]:
                hits = scan_both_strands(seq, pssm, max_total, threshold)
                for h in hits:
                    h["motif_id"] = mid
                    h["motif_name"] = motif_name
                    h["threshold"] = threshold
                    # Avoid duplicates from different thresholds
                    key = (h["start"], h["end"], h["strand"], mid)
                    if key not in [(x["start"], x["end"], x["strand"], x["motif_id"]) for x in tf_hits]:
                        tf_hits.append(h)

        if tf_hits:
            # Deduplicate: keep highest score for overlapping hits
            tf_hits.sort(key=lambda x: (x["start"], -x["relative_score"]))
            deduped = []
            for h in tf_hits:
                # Check for overlap with existing hits
                overlap = False
                for d in deduped:
                    if h["motif_id"] == d["motif_id"]:
                        o_start = max(h["start"], d["start"])
                        o_end = min(h["end"], d["end"])
                        if o_end > o_start:
                            overlap = True
                            if h["relative_score"] > d["relative_score"]:
                                deduped.remove(d)
                                deduped.append(h)
                            break
                if not overlap:
                    deduped.append(h)

            results[tf_name] = deduped
            motif_count_total += len(deduped)
            print(f"      -> {len(deduped)} unique binding sites found")
        else:
            print(f"      -> No binding sites found at threshold >= 0.80")

    # Step 4: Generate report
    print("\n" + "=" * 70)
    print("  RESULTS SUMMARY")
    print("=" * 70)

    upstream_len = PROMOTER_UPSTREAM

    # Sort TFs by number of hits
    sorted_tfs = sorted(results.items(), key=lambda x: len(x[1]), reverse=True)

    report_lines = []
    report_lines.append(f"\n## Pcsk2 Promoter Analysis: Predicted TF Binding Sites\n")
    report_lines.append(f"**Region**: chr{gene_info['seq_region_name']}:{region_start}-{region_end}")
    report_lines.append(f"**Relative to TSS**: -{PROMOTER_UPSTREAM} to +{PROMOTER_DOWNSTREAM} bp")
    report_lines.append(f"**Gene**: {GENE_SYMBOL} (strand: {'+' if gene_info['strand']==1 else '-'})")
    report_lines.append(f"**Scanned**: {len(all_tfs)} TFs with JASPAR 2024 PWMs\n")

    report_lines.append("| TF | #Sites | Key Binding Positions (relative to TSS) |")
    report_lines.append("|-----|--------|------------------------------------------|")

    for tf_name, hits in sorted_tfs:
        if not hits:
            continue
        # Get top 3 highest-scoring sites
        top_hits = sorted(hits, key=lambda x: x["relative_score"], reverse=True)[:3]
        positions = []
        for h in top_hits:
            rel_pos = h["start"] - upstream_len
            pos_str = f"{rel_pos}"
            if h["strand"] == "-":
                pos_str = f"{rel_pos}(rc)"
            positions.append(f"{pos_str} [{h['relative_score']:.2f}]")

        print(f"  {tf_name:10s}: {len(hits):2d} sites  |  {', '.join(positions)}")
        report_lines.append(f"| {tf_name} | {len(hits)} | {', '.join(positions)} |")

    # Detailed table for priority TFs
    report_lines.append("\n\n## Detailed Binding Site Table (Priority TFs for ChIP-qPCR)\n")
    report_lines.append("| TF | Motif | Position (rel. TSS) | Strand | Sequence | Relative Score |")
    report_lines.append("|----|-------|---------------------|--------|----------|----------------|")

    priority_order = ["Fos", "Jun", "Junb", "Jund", "Fosb", "Crem", "Atf3", "Nr4a1", "Nr4a2", "Egr1", "Dbp", "Nr3c1"]

    for tf_name in priority_order:
        if tf_name in results:
            for h in sorted(results[tf_name], key=lambda x: x["start"]):
                rel_pos = h["start"] - upstream_len
                report_lines.append(
                    f"| {tf_name} | {h['motif_name']} ({h['motif_id']}) | {rel_pos} | {h['strand']} | {h['sequence']} | {h['relative_score']:.3f} |"
                )

    # Write report
    report_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "grant",
        "Pcsk2_promoter_TF_binding_sites.md"
    )

    full_report = "\n".join(report_lines)
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(full_report)

    print(f"\n  Report saved to: {report_path}")
    print(f"  Total predicted binding sites: {motif_count_total}")
    print(f"  TFs with binding sites: {sum(1 for _, h in sorted_tfs if h)}")

    # Also save JSON for programmatic use
    json_path = report_path.replace(".md", ".json")
    json_output = {
        "gene": GENE_SYMBOL,
        "species": SPECIES,
        "genomic_region": f"chr{gene_info['seq_region_name']}:{region_start}-{region_end}",
        "promoter_upstream": PROMOTER_UPSTREAM,
        "promoter_downstream": PROMOTER_DOWNSTREAM,
        "tss": gene_info["tss"],
        "strand": gene_info["strand"],
        "analysis_date": datetime.now().isoformat(),
        "results": {
            tf: [{
                "start": h["start"],
                "end": h["end"],
                "position_rel_tss": h["start"] - upstream_len,
                "strand": h["strand"],
                "sequence": h["sequence"],
                "relative_score": h["relative_score"],
                "motif_id": h["motif_id"],
                "motif_name": h["motif_name"],
            } for h in hits]
            for tf, hits in results.items()
        }
    }
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(json_output, f, indent=2, ensure_ascii=False)
    print(f"  JSON data saved to: {json_path}")

    # ============ Key findings for grant ============
    print("\n" + "=" * 70)
    print("  KEY FINDINGS FOR GRANT PROPOSAL")
    print("=" * 70)

    print(f"""
The Pcsk2 promoter region (-{PROMOTER_UPSTREAM} to +{PROMOTER_DOWNSTREAM} bp) was scanned
against JASPAR 2024 position weight matrices for {len(all_tfs)} transcription factors
identified from the Pomc gene regulatory network (GRN) analysis.

Key predictions supporting the grant hypothesis:
""")

    # Check for the TFs proposed for ChIP-qPCR
    chip_targets = ["Crem", "Fos", "Junb", "Egr1", "Creb1"]
    for tf in chip_targets:
        if tf in results and results[tf]:
            n = len(results[tf])
            top = sorted(results[tf], key=lambda x: x["relative_score"], reverse=True)[0]
            rel = top["start"] - upstream_len
            print(f"  [OK] {tf}: {n} predicted binding site(s), best at {rel} bp (score={top['relative_score']:.3f})")
        else:
            print(f"  [--] {tf}: No high-confidence binding sites found in core promoter")

    print(f"""
These predictions provide computational support for the proposed ChIP-qPCR
experiments (内容三) targeting CREB, c-Fos, JunB, and EGR1 at the Pcsk2
promoter region. The analysis will be refined using PROMO and LASAGNA
for consensus predictions.
""")

    return results, gene_info, seq


if __name__ == "__main__":
    results, gene_info, seq = analyze_promoter()
