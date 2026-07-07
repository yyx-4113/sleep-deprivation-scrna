# Data and Code Availability Statement

## Raw Data

All raw sequencing data analyzed in this study are publicly available from the NCBI Gene Expression Omnibus (GEO):

| Dataset | Accession | Description | Samples |
|---------|-----------|-------------|---------|
| Primary | [GSE137665](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137665) | Mouse scRNA-seq: hypothalamus, brainstem, cortex under SD and control conditions | 29,571 cells |
| Validation #1 | [GSE211088](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211088) | Mouse prefrontal cortex bulk RNA-seq: 5h SD vs home-cage | 10 samples |
| Validation #2 | [GSE237419](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE237419) | Mouse cerebral cortex bulk RNA-seq: multi-timepoint SD | 42 samples |

## Processed Data

Processed single-cell data (h5ad format), differential expression tables, gene regulatory network edges, in silico knockout predictions, machine learning results, and validation data are available at:

**GitHub Repository:** [https://github.com/example/sleep-deprivation-scrnaseq](https://github.com/example/sleep-deprivation-scrnaseq)

The repository contains:
- `data/processed/` — Processed AnnData objects (h5ad)
- `results/tables/` — All differential expression, GRN, KO prediction, ML, and validation tables
- `results/figures/` — All 30 publication-quality figures (300 dpi PNG)
- `data/external/` — Ensembl-to-Symbol mapping file derived from NCBI gene_info

## Code Availability

All analysis scripts are available in the GitHub repository under `scripts/`:

| Script | Function |
|--------|----------|
| `step1_preprocessing.py` | Single-cell QC, normalization, PCA, UMAP, clustering |
| `step2_pomc_grn_ko.py` | Pomc co-expression network + in silico knockout |
| `step3_de_analysis.py` | Differential expression across brain regions |
| `step4_ml_biomarkers.py` | Machine learning biomarker discovery (LASSO, RF, SVM-RFE) |
| `step5_external_validation.py` | Cross-dataset validation against GSE211088, GSE237419 |
| `step6_tf_activity.py` | TF activity inference from co-expression regulons |
| `step7_pyscenic_local.py` | Gradient boosting regression (GRNBoost2-equivalent) GRN inference |

**Software requirements:** Python 3.12, Scanpy v1.10, scikit-learn v1.8, NumPy, SciPy, Pandas, Matplotlib, Seaborn, NetworkX, dask, arboreto.

A complete conda/pip environment specification is provided in `environment.yml` in the GitHub repository.

## Contact

For questions regarding data or code, please contact:

Yang Yongxin
Fujian Second People's Hospital
Email: 960856791@qq.com
