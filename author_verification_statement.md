# Author verification statement

**Manuscript:** "Unreplicated by design: pseudoreplication, data leakage and the limits of inference in single-cell transcriptomics of sleep deprivation — a data commentary" (Research Report, *Neurobiology of Sleep and Circadian Rhythms*).

**Prepared by:** Yongxin Yang (author)
**Method of verification:** every externally checkable assertion in the manuscript was re-checked against its primary source by the author, using the same identifiers cited in the text. No assertion was carried over from the previous round without re-verification.
**Date of verification:** 8 September 2026.

## What was checked

### 1. References (all 13, including the new [22])
Each reference was verified against NCBI PubMed (`eutils` esummary, db=pubmed) by PMID, or bibliographically where not indexed in PubMed. Verified fields: first author, journal, year, volume, issue, pages, PMID. The full audit table is `results/tables/reference_audit_final.csv` (deliverable P0-A). Result: **13/13 match the manuscript**; no fabricated, mistitled or mis-paginated records remain. In particular, the two methodological anchors of the commentary — Zimmerman *et al.*, *Nat Commun* 2021;12:738 (PMID 33531494) and Kapoor & Narayanan, *Patterns* 2023;4:100804 (PMID 37720327) — were confirmed to exist and to support the cited claims.

### 2. Table 1 (inventory of eligible datasets)
Table 1 was regenerated programmatically from NCBI GEO `esummary` records (`scripts/build_table1_from_geo.py`, records retrieved 8 September 2026 and cached in `results/tables/geo_esummary/`). The automated columns (GEO series title, deposit year, linked PMID, number of deposited samples, sample titles) are fetched at build time; the human-adjudicated columns (tissue/modality, contrast, libraries per condition, replication status, orthogonal validation, venue) are joined from `results/tables/table1_adjudication.csv`. Every cell was cross-checked against the deposited series record and sample titles. The table-record drift check reported no gross mismatches. Result: **9 acute-eligible rows + 1 chronic-excluded footnote, all consistent with GEO**.

### 3. Replication denominators
The primary denominator (CNS-restricted 3/7 = 43%, Wilson 95% CI [16%, 75%]) and its two sensitivities (all-tissue 4/9 = 44% [19%, 73%]; including the chronic series 5/10 = 50% [24%, 76%]) were recomputed with the Wilson score interval (`results/tables/wilson_ci.txt`). The excluded accessions (GSE280145 bone marrow, GSE289089 small intestinal crypt) are named explicitly in the text.

### 4. Randomised-block analysis (GSE137665)
The two-way ANOVA without replication was re-derived from the 3 × 3 layout of log₂(CPM + 1) values (`results/tables/randomised_block.txt`): SS(region) = 175.86 (df 2), SS(condition) = 3.90 (df 2), SS(region × condition) = 3.21 (df 4, used as error, MS = 0.80), SE = 0.73. All three condition-contrast 95% CIs include zero. Tukey's one-degree-of-freedom test for non-additivity was added (SS = 2.80, df 1; F(1,3) = 3.49, p = 0.16, non-significant).

### 5. Magnitude arithmetic (ambient-RNA burden)
The *Pomc* UMI burden was recomputed exactly: at 0.5% prevalence, 3,491 × 0.005 = 17.46 positive cells, implying 88,700 / 17.46 ≈ 5,082 *Pomc* UMI per positive cell = 447% of the mean transcriptome (1,136 UMI/cell); rounding the denominator to 17 cells gives 5,218 UMI/cell = 459%. Both exceed 100% and are therefore impossible for a single gene. The 0.5% prevalence is anchored to Campbell *et al.*, *Nat Neurosci* 2017;20:484–496 (PMID 28166221).

### 6. Data-leakage demonstration
The classifier analysis unit was made explicit: it pools all three brain regions (nine libraries) with region as a confounder; a hypothalamus-only classifier would have only three libraries and a permutation floor of P = 0.333, and was therefore not used. The exact permutation reference set (C(9,3) = 84 for SD-vs-rest; 1,680 for the full three-way relabelling) and the observed P = 0.57 are reported.

### 7. Self-reported counts
Word counts were recomputed directly from the manuscript source (abstract ≈ 250 words; main text ≈ 5,200 words, Introduction–Limitations, excluding references, tables, Box 1 and figure captions) and the header was corrected to match.

## Repository and persistent identifier
The analysis pipeline, all CSVs, the cached GEO JSON, the Table 1 build script, the figure-generation scripts, the reference-audit table and this statement are assembled as a version-controlled bundle (MIT licence) and will be made public at https://github.com/yongxinyang/sleep-deprivation-scrna with a persistent DOI (Zenodo or GigaDB/Science Data Bank) prior to acceptance. The identical materials are provided as Supplementary File S1.

**Signed:** Yongxin Yang
**Date:** 8 September 2026
