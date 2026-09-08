"""Analyses requested during peer review of the commentary manuscript.

Produces, under results/tables/:
  1. pomc_burden.csv          - implied per-cell Pomc burden under a range of prevalences
  2. randomised_block.txt     - 3 region x 3 condition randomised-block analysis
  3. wilson_ci.txt            - Wilson interval for 5/10
  4. geo_scan_45_exclusions.csv - full 45-row candidate table with exclusion reasons
  5. adjudication_gsm.csv     - per-GSM adjudication table for included series
"""

import csv
import json
import math
import os
import re

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results", "tables")
os.makedirs(RES, exist_ok=True)

PB = os.path.join(RES, "pseudobulk_unified", "pseudobulk_per_library.csv")

# ---------------------------------------------------------------- 1. burden
def burden():
    rows = list(csv.DictReader(open(PB, encoding="utf-8")))
    hypo = [r for r in rows if r["brain_region"] == "Hypothalamus"]
    out = []
    for r in hypo:
        n_cells = int(r["n_cells"])
        tot = int(r["lib_size_total_counts"])
        pomc = int(r["Pomc_total_counts"])
        depth = tot / n_cells                      # mean total UMI per cell
        per_cell_all = pomc / n_cells              # if spread over every cell
        rec = {
            "library": r["sample"],
            "condition": r["condition_code"],
            "n_cells": n_cells,
            "total_UMI": tot,
            "Pomc_UMI": pomc,
            "Pomc_CPM": float(r["Pomc_CPM"]),
            "mean_UMI_per_cell": round(depth, 1),
            "Pomc_pct_of_library": round(100 * pomc / tot, 3),
            "Pomc_UMI_per_cell_if_uniform": round(per_cell_all, 2),
        }
        for p in (0.005, 0.01, 0.075, 0.306, 0.50):
            k = max(1, round(p * n_cells))
            rec["implied_UMI_at_%s" % p] = round(pomc / k, 1)
            rec["pct_of_transcriptome_at_%s" % p] = round(100 * (pomc / k) / depth, 1)
        out.append(rec)
    with open(os.path.join(RES, "pomc_burden.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print("== 1. implied per-cell Pomc burden (hypothalamus) ==")
    for r in out:
        print("  %s (%s): %d cells, %d total UMI, %d Pomc UMI = %.2f%% of library; "
              "mean depth %.0f UMI/cell"
              % (r["library"], r["condition"], r["n_cells"], r["total_UMI"],
                 r["Pomc_UMI"], r["Pomc_pct_of_library"], r["mean_UMI_per_cell"]))
        for p in (0.005, 0.075, 0.306):
            print("     prevalence %5.1f%% -> %8.1f Pomc UMI/cell = %6.1f%% of mean transcriptome"
                  % (100 * p, r["implied_UMI_at_%s" % p],
                     r["pct_of_transcriptome_at_%s" % p]))
    return out


# ------------------------------------------------- 2. randomised block ANOVA
def randomised_block():
    rows = list(csv.DictReader(open(PB, encoding="utf-8")))
    regions = ["Hypothalamus", "Brainstem", "Cortex"]
    conds = ["A1", "A2", "A3"]
    y = np.zeros((3, 3))
    for r in rows:
        i = regions.index(r["brain_region"])
        j = conds.index(r["condition_code"])
        y[i, j] = math.log2(float(r["Pomc_CPM"]) + 1)
    grand = y.mean()
    row_m = y.mean(axis=1)
    col_m = y.mean(axis=0)
    ss_region = 3 * ((row_m - grand) ** 2).sum()
    ss_cond = 3 * ((col_m - grand) ** 2).sum()
    fit = row_m[:, None] + col_m[None, :] - grand
    ss_inter = ((y - fit) ** 2).sum()
    df_inter = (3 - 1) * (3 - 1)
    ms_inter = ss_inter / df_inter
    se_diff = math.sqrt(2 * ms_inter / 3)
    tcrit = 2.776  # t_{4, 0.975}
    lines = []
    lines.append("log2(CPM + 1) layout (rows = region, cols = A1/A2/A3):")
    for i, rg in enumerate(regions):
        lines.append("  %-13s %s" % (rg, "  ".join("%7.3f" % v for v in y[i])))
    lines.append("")
    lines.append("SS region      = %.3f (df 2)" % ss_region)
    lines.append("SS condition   = %.3f (df 2)" % ss_cond)
    lines.append("SS interaction = %.3f (df %d) -> used as error" % (ss_inter, df_inter))
    lines.append("MS error       = %.4f  -> residual SD = %.4f log2 units" %
                 (ms_inter, math.sqrt(ms_inter)))
    lines.append("SE of a condition difference = %.4f" % se_diff)
    lines.append("")
    lines.append("Region x condition interaction is the error term in an unreplicated")
    lines.append("randomised block. It is grossly inflated here because the region effect")
    lines.append("is ~10 log2 units, so this is a deliberately conservative bound.")
    lines.append("")
    for (a, b) in ((1, 0), (2, 0), (2, 1)):
        d = col_m[b] - col_m[a]
        lo, hi = d - tcrit * se_diff, d + tcrit * se_diff
        lines.append("A%s - A%s = %+.3f log2   95%% interval [%+.3f, %+.3f]  (width %.2f)"
                     % (conds[b][1], conds[a][1], d, lo, hi, hi - lo))
    txt = "\n".join(lines)
    open(os.path.join(RES, "randomised_block.txt"), "w", encoding="utf-8").write(txt)
    print("\n== 2. randomised-block analysis ==")
    print(txt)
    return ms_inter


# ------------------------------------------------------------ 3. Wilson CI
def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def wilson_report():
    lo, hi = wilson(5, 10)
    lines = ["5 of 10 series: point 50.0%%, Wilson 95%% CI [%.1f%%, %.1f%%]" % (100 * lo, 100 * hi)]
    for k, n, lab in ((4, 8, "published CNS-relevant only"),
                      (3, 7, "dropping GSE213496 (tissue unresolved)")):
        a, b = wilson(k, n)
        lines.append("  sensitivity %-34s %d/%d = %.0f%%  CI [%.0f%%, %.0f%%]"
                     % (lab, k, n, 100 * k / n, 100 * a, 100 * b))
    txt = "\n".join(lines)
    open(os.path.join(RES, "wilson_ci.txt"), "w", encoding="utf-8").write(txt)
    print("\n== 3. Wilson interval ==")
    print(txt)


# ----------------------------------------------------- 4. 45-row exclusions
SC_TERMS = ["single cell", "single-cell", "single nucleus", "single-nucleus",
            "scrna", "snrna", "10x", "10 x", "drop-seq", "dropseq", "smart-seq",
            "smartseq", "chromium"]
SLEEP_TERMS = ["sleep", "wake", "circadian", "zeitgeber", "zt0", "zt6", "zt12",
               "hypocretin", "orexin"]
SD_TERMS = ["sleep deprivation", "sleep-deprivation", "sleep deprived",
            "sleep-deprived", "sleep restriction", "sleep-restricted",
            "recovery sleep", "rebound sleep", "sleep recovery", "wakefulness",
            "sleep loss", "sleep fragmentation", "sleepiness"]


def classify(title, summary):
    blob = (title + " " + summary).lower()
    has_sc = any(t in blob for t in SC_TERMS)
    has_sleep = any(t in blob for t in SLEEP_TERMS)
    has_sd = any(t in blob for t in SD_TERMS)
    if not has_sc:
        return "Not single-cell / single-nucleus transcriptomics"
    if not has_sleep:
        return "Not sleep- or circadian-related"
    if not has_sd:
        return "Sleep-related but no acute sleep-deprivation or recovery-sleep contrast"
    return "Eligible on title/summary screen"


def exclusion_table():
    src = os.path.join(RES, "geo_scan_datasets.csv")
    rows = list(csv.DictReader(open(src, encoding="utf-8")))
    included = {"GSE137665", "GSE214337", "GSE213496", "GSE218624", "GSE280145",
                "GSE289089", "GSE211088", "GSE245537", "GSE256140", "GSE243489"}
    out = []
    counts = {}
    for r in rows:
        acc = r["accession"]
        reason = classify(r.get("title", ""), r.get("summary_excerpt", ""))
        if acc in included:
            reason = "INCLUDED"
        else:
            if reason == "Eligible on title/summary screen":
                reason = ("Eligible on title/summary but excluded at full-record screen: "
                          "no sample-level metadata resolving libraries per condition")
        counts[reason] = counts.get(reason, 0) + 1
        out.append({
            "accession": acc,
            "year": r.get("year", ""),
            "n_samples": r.get("n_samples", ""),
            "title": (r.get("title", "") or "")[:110],
            "auto_class": r.get("replication_class", ""),
            "included": "yes" if acc in included else "no",
            "exclusion_reason": reason,
        })
    path = os.path.join(RES, "geo_scan_45_exclusions.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    print("\n== 4. 45-row candidate table ==")
    for k in sorted(counts, key=lambda x: -counts[x]):
        print("  %3d  %s" % (counts[k], k))
    print("  total rows:", len(out))
    return out


# ------------------------------------------------------- 5. GSM adjudication
def adjudication():
    src = os.path.join(RES, "geo_scan_curated_samples.txt")
    if not os.path.exists(src):
        print("\n[5] sample file missing")
        return
    txt = open(src, encoding="utf-8").read()
    blocks = re.split(r"={20,}\n", txt)
    rows = []
    for b in blocks:
        m = re.match(r"(GSE\d+)\s*\|\s*([^|]+)\|\s*n=(\d+)\s*\|\s*PMID=(\S+)", b.strip())
        if not m:
            continue
        acc, date, n, pmid = m.group(1), m.group(2).strip(), m.group(3), m.group(4)
        if "SAMPLES:" not in b:
            continue
        samples = b.split("SAMPLES:", 1)[1]
        for line in samples.splitlines():
            line = line.strip()
            if not line.startswith("-"):
                continue
            rows.append({
                "accession": acc,
                "release_date": date,
                "pmid": pmid,
                "sample_title": line.lstrip("- ").strip(),
            })
    path = os.path.join(RES, "adjudication_gsm.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["accession", "release_date", "pmid", "sample_title"])
        w.writeheader()
        w.writerows(rows)
    print("\n== 5. per-sample adjudication rows ==")
    print("  %d sample titles across %d series"
          % (len(rows), len({r["accession"] for r in rows})))


if __name__ == "__main__":
    burden()
    randomised_block()
    wilson_report()
    exclusion_table()
    adjudication()
