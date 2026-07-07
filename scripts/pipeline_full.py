"""
Full Pipeline: TSV -> AnnData -> Preprocess -> DE -> GRN -> Perturbation
All paths use raw Windows format to avoid /c/ issues
"""
import sys, os
BASE = r'C:\Users\1\sleep-deprivation-project'
DATA = os.path.join(BASE, 'data')
FIG = os.path.join(BASE, 'results', 'figures')
TBL = os.path.join(BASE, 'results', 'tables')
os.makedirs(DATA, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
os.makedirs(TBL, exist_ok=True)

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import warnings
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from scipy import sparse
import gzip

sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=True)

# ============================================================
# STEP 1: Build AnnData from TSV
# ============================================================
tsv_path = os.path.join(DATA, 'GSE137665_exp_matrix.all.tsv.gz')
h5ad_raw = os.path.join(DATA, 'GSE137665_raw.h5ad')
h5ad_proc = os.path.join(DATA, 'GSE137665_processed.h5ad')

if os.path.exists(h5ad_proc):
    print("Loading processed data...")
    adata = sc.read_h5ad(h5ad_proc)
else:
    if not os.path.exists(h5ad_raw):
        print("Building sparse AnnData from TSV...")
        with gzip.open(tsv_path, 'rt') as f:
            header = f.readline().strip().split('\t')
            barcodes = header
            genes, rows, cols, vals = [], [], [], []
            for gene_idx, line in enumerate(f):
                parts = line.strip().split('\t')
                genes.append(parts[0])
                counts = np.array([float(x) for x in parts[1:]])
                nz = np.where(counts > 0)[0]
                if len(nz) > 0:
                    rows.extend([gene_idx] * len(nz))
                    cols.extend(nz.tolist())
                    vals.extend(counts[nz].tolist())
                if gene_idx % 5000 == 0:
                    print(f'  {gene_idx} genes...')

        X = sparse.coo_matrix((vals, (rows, cols)), shape=(len(genes), len(barcodes))).tocsc()
        adata_raw = sc.AnnData(X=X.T, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
        adata_raw.var_names_make_unique()
        adata_raw.obs['sample'] = [bc.split('-')[-1] if '-' in bc else bc for bc in adata_raw.obs.index]
        region_map = {'VAL630': 'Cortex', 'JHA710': 'Hypothalamus', 'VAL173': 'Brainstem'}
        adata_raw.obs['brain_region'] = adata_raw.obs['sample'].str[:6].map(region_map)
        adata_raw.obs['condition_code'] = adata_raw.obs['sample'].str[-2:]
        adata_raw.write(h5ad_raw, compression='gzip')
        print(f'Raw AnnData saved: {adata_raw.n_obs} cells x {adata_raw.n_vars} genes')
        adata = adata_raw
    else:
        print("Loading raw AnnData...")
        adata = sc.read_h5ad(h5ad_raw)

    # ============================================================
    # STEP 2: Preprocessing
    # ============================================================
    print(f'\nPreprocessing: {adata.n_obs} cells')
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribo'] = adata.var_names.str.startswith(('Rps','Rpl'))
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, inplace=True)

    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.n_genes_by_counts < 6000, :]
    adata = adata[adata.obs.pct_counts_mt < 15, :]
    print(f'After QC: {adata.n_obs} cells ({adata.n_obs/n_before*100:.1f}% kept)')

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=30)
    sc.pp.neighbors(adata, n_pcs=20, n_neighbors=15)
    sc.tl.umap(adata, min_dist=0.3, spread=1.0)
    sc.tl.leiden(adata, resolution=0.5)

    # Cell annotation
    BRAIN_MARKERS = {
        'Excitatory_Neuron': ['Slc17a7','Neurod6','Tbr1','Camk2a','Snap25'],
        'Inhibitory_Neuron': ['Gad1','Gad2','Slc32a1','Lhx6','Sst','Pvalb'],
        'Astrocyte': ['Gfap','Aqp4','Slc1a3','Aldh1l1'],
        'Microglia': ['Cx3cr1','Tmem119','P2ry12','Aif1','C1qa'],
        'Oligodendrocyte': ['Mog','Mbp','Plp1','Mag'],
        'OPC': ['Pdgfra','Cspg4','Sox10'],
        'Endothelial': ['Cldn5','Pecam1','Flt1'],
        'Ependymal': ['Foxj1','Tmem212'],
    }
    for ct, genes in BRAIN_MARKERS.items():
        present = [g for g in genes if g in adata.raw.var_names]
        if present:
            sc.tl.score_genes(adata, gene_list=present, score_name=f'score_{ct}', use_raw=True)
    score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
    adata.obs['cell_type'] = adata.obs[score_cols].idxmax(axis=1).str.replace('score_','')
    max_scores = adata.obs[score_cols].max(axis=1)
    adata.obs.loc[max_scores < 0.5, 'cell_type'] = 'Unassigned'

    print(f'Clusters: {adata.obs["leiden"].nunique()}')
    print(f'Cell types: {adata.obs["cell_type"].value_counts().to_dict()}')
    adata.write(h5ad_proc, compression='gzip')
    print('Processed data saved.')

# ============================================================
# STEP 3: UMAP Figures
# ============================================================
print('\n=== Generating UMAPs ===')
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='brain_region', ax=ax, show=False,
           palette={'Cortex':'#E64B35','Hypothalamus':'#4DBBD5','Brainstem':'#00A087'},
           title='Brain Regions')
plt.savefig(os.path.join(FIG, 'UMAP_brain_regions.png'), dpi=300, bbox_inches='tight')
plt.close()

fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='cell_type', ax=ax, show=False, title='Cell Types')
plt.savefig(os.path.join(FIG, 'UMAP_cell_types.png'), dpi=300, bbox_inches='tight')
plt.close()

# ============================================================
# STEP 4: Load raw for Pomc analysis
# ============================================================
print('\n=== Pomc GRN Analysis ===')
raw = sc.read_h5ad(h5ad_raw)
raw = raw[raw.obs_names.isin(adata.obs_names)]
hypo = raw[raw.obs['brain_region'] == 'Hypothalamus'].copy()
sc.pp.normalize_total(hypo, target_sum=1e4)
sc.pp.log1p(hypo)
print(f'Hypothalamus: {hypo.n_obs} cells')

# Pomc expression
pomc_expr = hypo[:, 'Pomc'].X.toarray().flatten() if sparse.issparse(hypo[:, 'Pomc'].X) else hypo[:, 'Pomc'].X
for code in ['A1','A2','A3']:
    mask = hypo.obs['condition_code'] == code
    print(f'  Pomc {code}: mean={pomc_expr[mask].mean():.3f}')

# ============================================================
# STEP 5: Co-expression network
# ============================================================
expr = hypo.X.toarray() if sparse.issparse(hypo.X) else hypo.X
pomc_idx = list(hypo.var_names).index('Pomc')
pomc_vec = expr[:, pomc_idx]
gene_means = expr.mean(axis=0)
n_top = min(5000, hypo.n_vars)
top_idx = np.argsort(gene_means)[-n_top:]

correlations = []
for i in top_idx:
    if i == pomc_idx: continue
    gv = expr[:, i]
    if np.std(gv) == 0: continue
    r, p = spearmanr(pomc_vec, gv)
    correlations.append({'gene': hypo.var_names[i], 'spearman_r': r, 'pvalue': p})

corr_df = pd.DataFrame(correlations).sort_values('spearman_r', key=abs, ascending=False)
corr_df['padj'] = corr_df['pvalue'].clip(lower=1e-300) * len(corr_df)
corr_df['padj'] = corr_df['padj'].clip(upper=1.0)
corr_df.to_csv(os.path.join(TBL, 'Pomc_coexpression.csv'), index=False)

print('Top 15 Pomc-correlated genes:')
for _, row in corr_df.head(15).iterrows():
    print(f'  {row["gene"]:15s}  r={row["spearman_r"]:+.3f}  padj={row["padj"]:.2e}')

# ============================================================
# STEP 6: TF identification + GRN
# ============================================================
MOUSE_TFS = {
    'Nr3c1','Crem','Fos','Fosb','Jun','Junb','Jund','Atf3','Atf4','Mef2c','Mef2d',
    'Srebf1','Srebf2','Nr4a1','Nr4a2','Nr4a3','Egr1','Egr2','Egr3','Klf4','Klf9',
    'Sox2','Sox9','Mafb','Cebpb','Hlf','Zic1','Npas4','Per1','Per2','Cry1','Cry2',
    'Clock','Arntl','Nfil3','Dbp','Tef','Stat3','Stat5a','Nfkb1','Rela','Sp1','Sp3',
    'Yy1','Ctcf','Pparg','Ppargc1a','Foxo1','Foxo3','Rest','Rcor1','Hdac1','Hdac2',
    'Mecp2','Olig1','Olig2','Sox10','Sox8','Tbr1','Neurod1','Neurod2','Neurod6',
    'Isl1','Lhx6','Dlx1','Dlx2','Dlx5','Nkx2-1','Otp','Sim1','Pou3f2','Pou3f4',
    'Bcl11b','Satb2','Cux1','Cux2','Fezf2','Rxrg','Rorb','Rora','Rorc','Nr1d1',
    'Nr1d2','Thrb','Ppara','Ppard','Nr5a1','Nr5a2','Esr1','Esr2','Ar',
}

pomc_net = corr_df[abs(corr_df['spearman_r']) > 0.3]
pomc_tfs = pomc_net[pomc_net['gene'].isin(MOUSE_TFS)]
print(f'\nTFs in Pomc network (|r|>0.3): {len(pomc_tfs)}')
for _, row in pomc_tfs.iterrows():
    print(f'  {row["gene"]:15s}  r={row["spearman_r"]:+.3f}')
pomc_tfs.to_csv(os.path.join(TBL, 'Pomc_network_TFs.csv'), index=False)

# ============================================================
# STEP 7: Build local GRN
# ============================================================
core_genes = ['Pomc'] + pomc_tfs['gene'].tolist()[:15] + corr_df['gene'].tolist()[:20]
core_genes = list(dict.fromkeys(core_genes))[:35]
core_idx = [list(hypo.var_names).index(g) for g in core_genes if g in hypo.var_names]
core_names = [hypo.var_names[i] for i in core_idx]
core_expr = expr[:, core_idx]

adj = np.zeros((len(core_idx), len(core_idx)))
for i in range(len(core_idx)):
    for j in range(i+1, len(core_idx)):
        if np.std(core_expr[:,i]) > 0 and np.std(core_expr[:,j]) > 0:
            r, _ = spearmanr(core_expr[:,i], core_expr[:,j])
            adj[i,j] = adj[j,i] = r

pomc_pos = core_names.index('Pomc')
pomc_mean = pomc_expr.mean()

perturb = []
for j, gene in enumerate(core_names):
    if gene == 'Pomc': continue
    r = adj[pomc_pos, j]
    pred_fc = -r * pomc_mean
    perturb.append({
        'gene': gene, 'correlation_with_Pomc': r,
        'predicted_log2FC': pred_fc,
        'direction': 'Down' if pred_fc < 0 else 'Up',
        'is_TF': gene in MOUSE_TFS
    })

perturb_df = pd.DataFrame(perturb).sort_values('predicted_log2FC', key=abs, ascending=False)
perturb_df.to_csv(os.path.join(TBL, 'Pomc_knockout_prediction.csv'), index=False)

print('\nPredicted effects of Pomc KO:')
for _, row in perturb_df.head(15).iterrows():
    tf_tag = ' [TF]' if row['is_TF'] else ''
    print(f'  {row["gene"]:15s}  log2FC={row["predicted_log2FC"]:+.3f}{tf_tag}')

# ============================================================
# STEP 8: Figures
# ============================================================
print('\n=== Generating GRN Figures ===')

# Pomc UMAP
fig, ax = plt.subplots(figsize=(8, 7))
hypo_adata = adata[adata.obs['brain_region']=='Hypothalamus']
sc.pl.umap(hypo_adata, color='score_Astrocyte', ax=ax, show=False,
           title='Cell Types (Hypothalamus)', cmap='viridis')
plt.savefig(os.path.join(FIG, 'Hypothalamus_UMAP.png'), dpi=300, bbox_inches='tight')
plt.close()

# GRN network
top_net = ['Pomc'] + perturb_df.head(15)['gene'].tolist()
net_idx_sub = [core_names.index(g) for g in top_net if g in core_names]
net_names_sub = [core_names[i] for i in net_idx_sub]

G = nx.Graph()
for name in net_names_sub:
    G.add_node(name, is_tf=name in MOUSE_TFS, is_pomc=(name=='Pomc'))
for i in range(len(net_idx_sub)):
    for j in range(i+1, len(net_idx_sub)):
        r = adj[net_idx_sub[i], net_idx_sub[j]]
        if abs(r) > 0.2:
            G.add_edge(net_names_sub[i], net_names_sub[j], weight=abs(r), sign='+' if r>0 else '-')

fig, ax = plt.subplots(figsize=(16, 14))
pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
node_colors = ['#E64B35' if n=='Pomc' else '#4DBBD5' if G.nodes[n].get('is_tf') else '#CCCCCC' for n in G.nodes()]
node_sizes = [2000 if n=='Pomc' else 1200 if G.nodes[n].get('is_tf') else 800 for n in G.nodes()]
edge_colors = ['#E64B35' if G[u][v]['sign']=='+' else '#4DBBD5' for u,v in G.edges()]
edge_widths = [abs(G[u][v]['weight'])*3 for u,v in G.edges()]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.9, ax=ax)
nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax)

from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#E64B35', markersize=15, label='Pomc'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#4DBBD5', markersize=12, label='Transcription Factor'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#CCCCCC', markersize=10, label='Target Gene'),
    Line2D([0],[0], color='#E64B35', lw=2, label='Positive correlation'),
    Line2D([0],[0], color='#4DBBD5', lw=2, label='Negative correlation'),
], loc='upper left', frameon=False)
ax.set_title('Pomc Regulatory Network in Hypothalamus', fontsize=14, fontweight='bold')
ax.axis('off')
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'Pomc_GRN_Network.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Network saved')

# KO prediction barplot
fig, ax = plt.subplots(figsize=(10, 6))
top20 = perturb_df.head(20)
colors = ['#E64B35' if x < 0 else '#4DBBD5' for x in top20['predicted_log2FC']]
ax.barh(range(len(top20)), top20['predicted_log2FC'].values[::-1], color=colors[::-1])
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['gene'].values[::-1], fontsize=9)
ax.set_xlabel('Predicted log2 Fold Change (Pomc KO)')
ax.axvline(0, color='black', linewidth=0.5)
ax.set_title('In Silico Pomc Knockout: Predicted Downstream Effects', fontsize=13, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'Pomc_KO_Prediction.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  KO prediction saved')

# Pomc boxplot
fig, ax = plt.subplots(figsize=(6, 5))
data = [pomc_expr[hypo.obs['condition_code']==c] for c in ['A1','A2','A3']]
bp = ax.boxplot(data, labels=['A1','A2','A3'], patch_artist=True)
for patch, c in zip(bp['boxes'], ['#4DBBD5','#E64B35','#00A087']):
    patch.set_facecolor(c); patch.set_alpha(0.7)
ax.set_ylabel('Pomc Expression (log-normalized)')
ax.set_title('Pomc Expression Across Conditions (Hypothalamus)', fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'Pomc_Conditions_Boxplot.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Boxplot saved')

print(f'\n=== Pipeline Complete ===')
print(f'Figures: {len(os.listdir(FIG))} in {FIG}')
print(f'Tables: {len(os.listdir(TBL))} in {TBL}')
