# sleep-deprivation-scrna — reproducibility bundle

Version-controlled companion repository for the data commentary:

> Yang Y. *Unreplicated by design: pseudoreplication, data leakage and the
> limits of inference in single-cell transcriptomics of sleep deprivation — a
> data commentary.* Neurobiology of Sleep and Circadian Rhythms (Research Report).

This bundle contains the complete analysis pipeline behind the manuscript: the
GEO-driven inventory, the Table 1 build script, the reference audit, the
randomised-block and magnitude-arithmetic computations, and the figure
sources. It is released under the MIT licence. A persistent DOI (Zenodo or
GigaDB / Science Data Bank) will be minted prior to acceptance.

## What is in here

```
scripts/
  build_table1_from_geo.py   regenerate Table 1 from cached GEO esummary
  build_reference_audit.py   verify all references against PubMed
  generate_commentary_figures.py   figure sources (>=600 dpi)
  revision_analyses.py       randomised-block, Wilson CI, Pomc burden
results/tables/
  geo_esummary/*.json        cached GEO esummary (retrieved 2026-09-08)
  pubmed/*.json              cached PubMed esummary (retrieved 2026-09-08)
  table1_adjudication.csv    series-level human adjudication (joined by script)
  table1_generated.csv/.md   Table 1 as produced by the build script
  reference_audit_final.csv  per-reference PubMed verification
  geo_scan_45_exclusions.csv full 45-row candidate list + exclusion reasons
  adjudication_gsm_classified.csv  per-GSM adjudication (107 rows)
  randomised_block.txt, wilson_ci.txt, pomc_burden.csv
supplementary_s1_manifest.md   manifest of Supplementary File S1
author_verification_statement.md  what was checked, by what method, on what date
LICENSE                        MIT
```

## Reproduce Table 1

```bash
python scripts/build_table1_from_geo.py
```

The script reads the cached `geo_esummary/*.json` (or fetches live if missing),
joins the human-adjudicated columns from `table1_adjudication.csv`, runs a
table–record drift check, and writes `results/tables/table1_generated.csv` and
`.md`. Re-running on a later date with fresh GEO records updates the prevalence
estimate automatically.

## Verify references

```bash
python scripts/build_reference_audit.py
```

Outputs `results/tables/reference_audit_final.csv` (13/13 verified).

## Retrieval date

All live records were retrieved from NCBI on **8 September 2026**.
