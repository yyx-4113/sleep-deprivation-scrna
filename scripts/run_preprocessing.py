# Sleep Deprivation scRNA-seq: Full Preprocessing
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=True)

# Load
print("Loading data...")
adata = sc.read_h5ad('/c/Users/1/sleep-deprivation-project/data/GSE137665_raw.h5ad')
print(adata)

# Add QC metrics
adata.var['mt'] = adata.var_names.str.startswith('mt-')
adata.var['ribo'] = adata.var_names.str.startswith(('Rps', 'Rpl'))
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, inplace=True)

# Filter
print(f'\nBefore QC: {adata.n_obs} cells')
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata = adata[adata.obs.n_genes_by_counts < 6000, :]
adata = adata[adata.obs.pct_counts_mt < 15, :]
print(f'After QC: {adata.n_obs} cells, {adata.n_vars} genes')

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
print(f'HVG: {adata.n_vars} genes')

# PCA
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
print(f'PCA: {adata.obsm["X_pca"].shape}')

# Neighbors + UMAP
sc.pp.neighbors(adata, n_pcs=20, n_neighbors=15)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# Clustering
sc.tl.leiden(adata, resolution=0.5)
print(f'Clusters: {adata.obs["leiden"].nunique()}')

# Cell type annotation using brain marker genes
BRAIN_MARKERS = {
    'Excitatory_Neuron': ['Slc17a7', 'Neurod6', 'Tbr1', 'Camk2a', 'Snap25'],
    'Inhibitory_Neuron': ['Gad1', 'Gad2', 'Slc32a1', 'Lhx6', 'Sst', 'Pvalb'],
    'Astrocyte': ['Gfap', 'Aqp4', 'Slc1a3', 'Aldh1l1'],
    'Microglia': ['Cx3cr1', 'Tmem119', 'P2ry12', 'Aif1', 'C1qa'],
    'Oligodendrocyte': ['Mog', 'Mbp', 'Plp1', 'Mag'],
    'OPC': ['Pdgfra', 'Cspg4', 'Sox10'],
    'Endothelial': ['Cldn5', 'Pecam1', 'Flt1'],
    'Ependymal': ['Foxj1', 'Tmem212'],
}

print('\nCell type annotation...')
for ct, genes in BRAIN_MARKERS.items():
    present = [g for g in genes if g in adata.raw.var_names]
    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name=f'score_{ct}', use_raw=True)

score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
adata.obs['cell_type'] = adata.obs[score_cols].idxmax(axis=1).str.replace('score_', '')
max_scores = adata.obs[score_cols].max(axis=1)
adata.obs.loc[max_scores < 0.5, 'cell_type'] = 'Unassigned'

print('Cell type distribution:')
print(adata.obs['cell_type'].value_counts())

# Marker genes per cluster
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=50)
print('\nTop markers per cluster identified.')

# Save
adata.write('/c/Users/1/sleep-deprivation-project/data/GSE137665_processed.h5ad', compression='gzip')
print(f'\nSaved processed data.')

# Quick summary
print(f'\n=== Analysis Summary ===')
print(f'Cells: {adata.n_obs}')
print(f'Clusters: {adata.obs["leiden"].nunique()}')
print(f'Cell types: {adata.obs["cell_type"].nunique() - 1} (+ Unassigned)')
for region in ['Brainstem', 'Hypothalamus', 'Cortex']:
    sub = adata[adata.obs['brain_region'] == region]
    print(f'{region}: {sub.n_obs} cells, {sub.obs["leiden"].nunique()} clusters')
