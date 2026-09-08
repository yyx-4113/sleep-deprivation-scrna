# Supplementary File S1 — manifest

This file lists every material deposited as Supplementary File S1 with the
manuscript "Unreplicated by design: pseudoreplication, data leakage and the
limits of inference in single-cell transcriptomics of sleep deprivation — a
data commentary". The same materials are also released in the version-controlled
repository bundle (MIT licence; persistent DOI to be minted prior to acceptance).

| # | File | Description | Source / how generated |
|---|------|-------------|------------------------|
| 1 | `geo_scan_45_exclusions.csv` | Full 45-row candidate list with exclusion reasons (14 not sc/sn; 17 sleep-related but no acute SD/recovery contrast; 3 unrelated; 1 insufficient metadata; 1 chronic) | GEO inventory + manual curation |
| 2 | `adjudication_gsm_classified.csv` | Per-sample (per-GSM) adjudication: assigned condition, assignment rule, libraries per condition, in-primary-inventory flag (107 rows) | GEO sample titles + keyword rules |
| 3 | `table1_adjudication.csv` | Series-level human adjudication joined to GEO by the Table 1 build script (tissue/modality, contrast, libs/cond, replication status, orthogonal validation, venue) | author adjudication |
| 4 | `geo_esummary/*.json` | Cached NCBI GEO `esummary` records for the 10 inventory series (GSE137665, GSE214337, GSE218624, GSE280145, GSE211088, GSE243489, GSE256140, GSE245537, GSE289089, GSE213496) | E-utilities, retrieved 2026-09-08 |
| 5 | `build_table1_from_geo.py` | Script that regenerates Table 1 from the cached GEO records at build time (P0-G) | this work |
| 6 | `table1_generated.csv` / `table1_generated.md` | Machine- and human-readable Table 1 produced by the build script | build_table1_from_geo.py |
| 7 | `reference_audit_final.csv` | Per-reference PubMed verification (first author, journal, year, volume, pages, PMID; verdict) for all 13 references | build_reference_audit.py |
| 8 | `pubmed/*.json` | Cached NCBI PubMed `esummary` records for the 13 cited PMIDs | E-utilities, retrieved 2026-09-08 |
| 9 | `randomised_block.txt` | Re-derived two-way ANOVA without replication (SS/df/MS, CIs, Tukey non-additivity) | revision script |
| 10 | `wilson_ci.txt` | Wilson 95% CIs for 3/7, 4/9, 5/10 | revision script |
| 11 | `pomc_burden.csv` | *Pomc* UMI-burden arithmetic across prevalence assumptions | revision script |
| 12 | `geo_scan_curated.csv` | Curated inventory summary (10 acute-eligible + GSE213496 chronic) | GEO inventory |
| 13 | `author_verification_statement.md` | Author verification statement (what was checked, by what method, on what date) | this work |
| 14 | `supplementary_s1_manifest.md` | This manifest | this work |
| 15 | `figures/fig1.png` … `fig6.png` (+ `fig3b_burden.png`) | All six printed figures at ≥600 dpi | generate_commentary_figures.py |
| 16 | `scripts/*.py` | All analysis and figure-generation scripts (reproducibility) | this work |

**Retrieval date for all live records:** 8 September 2026.
**Licence:** MIT (code and tables); figures and text CC-BY (upon acceptance).
