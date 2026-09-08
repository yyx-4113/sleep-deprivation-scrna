#!/usr/bin/env python3
"""
build_table1_from_geo.py
=======================
Regenerate Table 1 of the data commentary directly from NCBI GEO `esummary`
records (cached as JSON), joining human-adjudicated series-level columns from
`table1_adjudication.csv`.

Why this script exists
---------------------
The round-3 editorial decision (NBSCR, 8 Sep 2026) required that Table 1 be
*generated* from the records at build time rather than "written by hand", so
that table<->record drift is structurally impossible and any reader can re-run
the inventory on a later date. The automated columns (GEO series title, deposit
year, linked PMID, number of deposited samples, sample titles) come straight
from GEO; the human-judged columns (tissue/modality as interpreted, contrast,
libraries per condition, replication status, orthogonal validation, venue) are
adjudicated by the author from the deposited sample titles and joined by
accession.

Outputs
-------
  results/tables/table1_generated.csv   - machine-readable table
  results/tables/table1_generated.md    - markdown table for the manuscript
  (drift-check report printed to stdout)

Re-run: `python scripts/build_table1_from_geo.py`
Requires network only if a cached JSON is missing; otherwise uses cache.
"""
import csv, json, os, sys, ssl, urllib.request, urllib.parse, datetime

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
GEO_DIR = os.path.join(ROOT, "results", "tables", "geo_esummary")
ADJ_CSV = os.path.join(ROOT, "results", "tables", "table1_adjudication.csv")
OUT_CSV = os.path.join(ROOT, "results", "tables", "table1_generated.csv")
OUT_MD = os.path.join(ROOT, "results", "tables", "table1_generated.md")
RETRIEVAL_DATE = "2026-09-08"

ctx = ssl.create_default_context()
ctx.check_hostname = False
ctx.verify_mode = ssl.CERT_NONE


def fetch_geo(acc):
    uid = "200" + acc[3:]
    cache = os.path.join(GEO_DIR, f"{acc}.json")
    if os.path.exists(cache):
        with open(cache) as f:
            return json.load(f)
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?"
           + urllib.parse.urlencode({"db": "gds", "id": uid, "retmode": "json"}))
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(req, context=ctx) as r:
        data = json.load(r)
    os.makedirs(GEO_DIR, exist_ok=True)
    with open(cache, "w") as f:
        json.dump(data, f, indent=2)
    return data


def geo_fields(acc):
    d = fetch_geo(acc)
    r = d["result"][f"200{acc[3:]}"]
    samples = r.get("samples") or []
    titles = [s.get("title", "") if isinstance(s, dict) else str(s) for s in samples]
    pmid = (r.get("pubmedids") or [])
    return {
        "geo_series_title": r.get("title", ""),
        "geo_pdat": (r.get("pdat") or "")[:4],
        "geo_pmid": ",".join(pmid) if pmid else "(none)",
        "geo_n_samples": len(samples),
        "geo_sample_titles": "; ".join(titles),
    }


def load_adjudication():
    rows = []
    with open(ADJ_CSV, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def main():
    adj = load_adjudication()
    acute = [r for r in adj if r["acute_included"].strip().upper() == "TRUE"]
    chronic = [r for r in adj if r["acute_included"].strip().upper() == "FALSE"]

    table_rows = []
    drift_notes = []
    for r in acute:
        acc = r["accession"]
        g = geo_fields(acc)
        # Drift checks (informational; the human columns are authoritative for the
        # manuscript, but GEO n_samples must be compatible with libs_per_cond).
        table_rows.append({
            "Accession": acc,
            "Year": r["year_display"],
            "Tissue / Modality": r["tissue_modality"],
            "Contrast (as deposited)": r["contrast"],
            "Libraries / cond.": r["libs_per_cond"],
            "Repl. at lib. level?": r["repl_at_lib"],
            "Orthogonal val. reported?": r["ortho_val"],
            "Venue": r["venue"],
            # automated provenance (not printed in the paper table, kept for audit)
            "_geo_n_samples": g["geo_n_samples"],
            "_geo_pmid": g["geo_pmid"],
            "_geo_series_title": g["geo_series_title"],
        })
        # Sanity: if a condition is "replicated" we expect >=2 libs/cond; if not,
        # 1 lib/cond. Flag only gross mismatches vs deposited sample count.
        try:
            expected = int(str(r["libs_per_cond"]).split()[0])
        except Exception:
            expected = None
        if expected is not None and g["geo_n_samples"] < expected:
            drift_notes.append(f"  [WARN] {acc}: GEO deposited {g['geo_n_samples']} samples but "
                               f"adjudicated libs/cond = {r['libs_per_cond']}")

    # ---- write CSV ----
    cols = ["Accession", "Year", "Tissue / Modality", "Contrast (as deposited)",
            "Libraries / cond.", "Repl. at lib. level?", "Orthogonal val. reported?",
            "Venue", "_geo_n_samples", "_geo_pmid", "_geo_series_title"]
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for row in table_rows:
            w.writerow({c: row.get(c, "") for c in cols})

    # ---- write markdown ----
    print_table_md(table_rows, chronic)
    with open(OUT_MD, "w") as f:
        f.write(render_md(table_rows, chronic))

    # ---- drift report ----
    print("=" * 70)
    print(f"Table 1 build complete. Retrieval date: {RETRIEVAL_DATE}")
    print(f"Acute-eligible rows: {len(table_rows)} | Chronic-excluded: {len(chronic)}")
    if drift_notes:
        print("Drift checks:")
        for n in drift_notes:
            print(n)
    else:
        print("Drift checks: no gross GEO<->adjudication mismatches.")
    print("Replication tally (acute): "
          f"{sum(1 for r in table_rows if r['Repl. at lib. level?'].strip().lower().startswith('y'))} replicated, "
          f"{sum(1 for r in table_rows if r['Repl. at lib. level?'].strip().lower().startswith('n'))} not replicated.")
    # CNS subset
    cns = [r for r in table_rows if not r["Tissue / Modality"].lower().startswith(("small intestinal", "lin-kit", "bone marrow", "lin-"))]
    cns_unrep = sum(1 for r in cns if r["Repl. at lib. level?"].strip().lower().startswith("n"))
    print(f"CNS subset: {cns_unrep}/{len(cns)} not replicated at library level.")


def render_md(table_rows, chronic):
    lines = []
    lines.append("| Accession | Year | Tissue / Modality | Contrast (as deposited) | "
                 "Libraries / cond. | Repl. at lib. level? | Orthogonal val. reported? | Venue |")
    lines.append("|-----------|------|-------------------|--------------------------|"
                 "-------------------|----------------------|---------------------------|-------|")
    for r in table_rows:
        lines.append("| {} | {} | {} | {} | {} | {} | {} | {} |".format(
            r["Accession"], r["Year"], r["Tissue / Modality"], r["Contrast (as deposited)"],
            r["Libraries / cond."], r["Repl. at lib. level?"], r["Orthogonal val. reported?"], r["Venue"]))
    for r in chronic:
        lines.append("")
        lines.append("* {} ({}; {}; {}) met title screening but was excluded by criterion (b) "
                     "(acute <=48 h); it also had one library per condition, reinforcing the "
                     "replication gap while lying outside the acute scope.".format(
                         r["accession"], r["venue"], r["tissue_modality"], r["contrast"]))
    return "\n".join(lines) + "\n"


def print_table_md(table_rows, chronic):
    sys.stdout.write(render_md(table_rows, chronic))


if __name__ == "__main__":
    main()
