"""
Step 7: GRNBoost2-style GRN inference (sklearn GradientBoostingRegressor) + Pomc perturbation
Equivalent to pySCENIC step 1 without dask dependency (arboreto incompatible with dask version).
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from scipy.stats import spearmanr
from scipy.sparse import issparse
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import os, warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=150, facecolor='white')
RESULTS_DIR = 'results/figures'
TABLES_DIR = 'results/tables'
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# ============================================================
# 1. Load data
# ============================================================
print('Loading data...')
adata = sc.read_h5ad('data/GSE137665_processed.h5ad')
raw = sc.read_h5ad('data/GSE137665_raw.h5ad')

hypo_processed = adata[adata.obs['brain_region'] == 'Hypothalamus'].copy()
hypo_raw = raw[raw.obs_names.isin(hypo_processed.obs_names)].copy()
sc.pp.normalize_total(hypo_raw, target_sum=1e4)
sc.pp.log1p(hypo_raw)
print(f'Hypothalamus: {hypo_raw.n_obs} cells, {hypo_raw.n_vars} genes')

# ============================================================
# 2. Define TFs and target genes of interest
# ============================================================
print('\nSelecting genes for GRN inference...')

MOUSE_TFS = [
    'Nr3c1','Crem','Fos','Fosb','Jun','Junb','Jund','Atf3','Atf4',
    'Mef2c','Mef2d','Srebf1','Srebf2','Nr4a1','Nr4a2','Nr4a3',
    'Egr1','Egr2','Egr3','Klf4','Klf9','Sox2','Sox9','Mafb','Cebpb',
    'Hlf','Zic1','Npas4','Per1','Per2','Cry1','Cry2','Clock','Arntl',
    'Nfil3','Dbp','Tef','Stat3','Stat5a','Nfkb1','Rela','Sp1','Sp3',
    'Yy1','Ctcf','Pparg','Ppargc1a','Foxo1','Foxo3','Rest','Mecp2',
    'Olig1','Olig2','Sox10','Tbr1','Neurod1','Neurod2','Neurod6',
    'Lhx6','Dlx1','Dlx2','Nkx2-1','Sim1','Bcl11b','Satb2',
    'Rxrg','Rorb','Rora','Rorc','Nr1d1','Nr1d2','Thrb','Ppara','Ppard',
    'Nr5a1','Nr5a2','Esr1','Ar','Hdac1','Hdac2','Rbm3','Rcor1','Fezf2',
    'Cux1','Cux2','Pomc','Apoe','Tsc22d3','Bdnf','Arc','Homer1','Pcsk2','Scg2'
]
MOUSE_TFS = list(dict.fromkeys(MOUSE_TFS))

# TFs present in data
tfs_present = sorted([t for t in MOUSE_TFS if t in hypo_raw.var_names])
print(f'TFs available: {len(tfs_present)}')

# Known Pomc-related genes from step2 co-expression analysis + key DEGs
pomc_coexpressed = [
    'Pomc','Fos','Jun','Junb','Jund','Fosb','Egr1','Egr2','Egr3',
    'Nr4a1','Nr4a2','Nr4a3','Per1','Per2','Nfil3','Dbp','Cry2',
    'Nr3c1','Tsc22d3','Srebf1','Atf3','Crem','Rbm3','Klf4','Sox9',
    'Nr1d1','Clock','Arntl','Ppargc1a','Mafb','Npas4','Arc','Bdnf',
    'Homer1','Pcsk2','Scg2','Cga','Apoe','Malat1','Meg3','Ccnd2',
    'Rnaset2a','Mt1','Mt3','Gnas','Cebpb','Yy1','Sp1','Hdac1',
]
# Filter to present genes
target_genes = sorted([g for g in pomc_coexpressed if g in hypo_raw.var_names])
# Add TFs as targets too (for TF→TF regulatory relationships)
for tf in tfs_present:
    if tf not in target_genes:
        target_genes.append(tf)
print(f'Target genes for GRN: {len(target_genes)}')

# ============================================================
# 3. Run GBR-based GRN inference (GRNBoost2 equivalent)
# ============================================================
print('\n=== Running Gradient Boosting Regression for GRN (GRNBoost2 equivalent) ===')

# Get expression data
all_genes = sorted(set(tfs_present + target_genes))
expr_for_grn = hypo_raw[:, all_genes].copy()
if issparse(expr_for_grn.X):
    X_all = expr_for_grn.X.toarray()
else:
    X_all = expr_for_grn.X

# TF expression matrix (features)
tf_indices = [all_genes.index(tf) for tf in tfs_present]
X_tf = X_all[:, tf_indices]  # cells x TFs

# Results collection
all_links = []

print(f'Training regressors: {len(target_genes)} target genes x {len(tfs_present)} TFs...')
n_done = 0
for target_idx, target_gene in enumerate(target_genes):
    if target_gene in tfs_present:
        # For TF targets, exclude self from predictors
        tf_mask = np.array([tf != target_gene for tf in tfs_present])
        X_sub = X_tf[:, tf_mask]
        tf_sub = [tf for tf, m in zip(tfs_present, tf_mask) if m]
    else:
        X_sub = X_tf
        tf_sub = tfs_present

    if X_sub.shape[1] == 0:
        continue

    y = X_all[:, all_genes.index(target_gene)]

    # Skip genes with no variance
    if np.std(y) < 1e-10:
        continue

    try:
        model = GradientBoostingRegressor(
            n_estimators=500, max_depth=3, learning_rate=0.01,
            subsample=0.8, min_samples_leaf=5, random_state=42, verbose=0
        )
        model.fit(X_sub, y)

        # Extract importance for TF→target links
        importances = model.feature_importances_
        for tf_name, imp in zip(tf_sub, importances):
            if imp > 0:
                all_links.append({
                    'TF': tf_name,
                    'target': target_gene,
                    'importance': imp
                })
    except Exception:
        pass

    n_done += 1
    if n_done % 20 == 0:
        print(f'  Progress: {n_done}/{len(target_genes)} targets done, {len(all_links)} links found')

links = pd.DataFrame(all_links, columns=['TF', 'target', 'importance'])
links = links.sort_values('importance', ascending=False).reset_index(drop=True)
print(f'\nGRNBoost2-style (sklearn GBR): {len(links)} regulatory links')
print(f'TFs with predicted targets: {links["TF"].nunique()}')
print(f'Genes with predicted regulators: {links["target"].nunique()}')

links.to_csv(f'{TABLES_DIR}/GRNBoost2_links.csv', index=False)

# ============================================================
# 4. Pomc-centered analysis
# ============================================================
print('\n=== Pomc-Centered Network Analysis ===')

pomc_targets = links[links['TF'] == 'Pomc'].sort_values('importance', ascending=False)
pomc_regulators = links[links['target'] == 'Pomc'].sort_values('importance', ascending=False)

print(f'Pomc -> targets: {len(pomc_targets)}')
if len(pomc_targets) > 0:
    print('Top targets:')
    for _, row in pomc_targets.head(10).iterrows():
        print(f'  {row["target"]:15s}  importance={row["importance"]:.6f}')

print(f'Regulators -> Pomc: {len(pomc_regulators)}')
if len(pomc_regulators) > 0:
    print('Top regulators:')
    for _, row in pomc_regulators.head(10).iterrows():
        print(f'  {row["TF"]:15s}  importance={row["importance"]:.6f}')

# Global TF connectivity
tf_connectivity = links.groupby('TF').agg(
    n_targets=('target', 'nunique'),
    mean_importance=('importance', 'mean'),
    max_importance=('importance', 'max')
).sort_values('n_targets', ascending=False)
print(f'\nTop 15 most connected TFs:')
print(tf_connectivity.head(15))

# ============================================================
# 5. Build Pomc-centered GRN
# ============================================================
print('\n=== Building Pomc GRN ===')

pomc_network_genes = {'Pomc'}
if len(pomc_regulators) > 0:
    pomc_network_genes.update(pomc_regulators.head(10)['TF'].tolist())
if len(pomc_targets) > 0:
    pomc_network_genes.update(pomc_targets.head(15)['target'].tolist())
pomc_network_genes.update(tf_connectivity.head(8).index.tolist())
print(f'Network nodes: {len(pomc_network_genes)}')

G = nx.Graph()
G.add_node('Pomc', type='hub')

for _, row in pomc_regulators.head(10).iterrows():
    if row['TF'] in pomc_network_genes:
        G.add_node(row['TF'], type='tf')
        G.add_edge(row['TF'], 'Pomc', weight=row['importance']*100, edge_type='regulates')

for _, row in pomc_targets.head(15).iterrows():
    if row['target'] in pomc_network_genes:
        G.add_node(row['target'], type='target')
        G.add_edge('Pomc', row['target'], weight=row['importance']*100, edge_type='targets')

tf_nodes_list = [n for n in G.nodes() if G.nodes[n].get('type') == 'tf']
for i, tf1 in enumerate(tf_nodes_list):
    for tf2 in tf_nodes_list[i+1:]:
        shared = set(links[links['TF']==tf1]['target']) & set(links[links['TF']==tf2]['target'])
        if len(shared) >= 2:
            G.add_edge(tf1, tf2, weight=len(shared)*0.5, edge_type='co-regulation')

print(f'Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges')

# ============================================================
# 6. In silico Pomc KO
# ============================================================
print('\n=== In Silico Pomc Knockout ===')

# Use core genes from the GRN + key co-expressed genes
core_genes_all = list(pomc_network_genes)
core_genes = sorted([g for g in core_genes_all if g in hypo_raw.var_names])

# Get expression for core genes
if issparse(hypo_raw[:, core_genes].X):
    core_expr = hypo_raw[:, core_genes].X.toarray().T  # genes x cells
else:
    core_expr = hypo_raw[:, core_genes].X.T

n_genes = len(core_genes)
print(f'Core genes: {n_genes}')

# Spearman adjacency
adj = np.zeros((n_genes, n_genes))
for i in range(n_genes):
    for j in range(i+1, n_genes):
        if np.std(core_expr[i]) > 0 and np.std(core_expr[j]) > 0:
            r, _ = spearmanr(core_expr[i], core_expr[j])
            adj[i, j] = r
            adj[j, i] = r

pomc_idx = core_genes.index('Pomc')
pomc_mean = core_expr[pomc_idx].mean()

# Direct correlation-based prediction
direct_effects = []
for j, gene in enumerate(core_genes):
    if gene == 'Pomc':
        continue
    r = adj[pomc_idx, j]
    pred_fc = -r * pomc_mean
    direct_effects.append({
        'gene': gene,
        'correlation_with_Pomc': r,
        'predicted_log2FC': pred_fc,
        'direction': 'Downregulated' if pred_fc < -0.05 else ('Upregulated' if pred_fc > 0.05 else 'Unchanged')
    })
perturb_df = pd.DataFrame(direct_effects).sort_values('predicted_log2FC', key=abs, ascending=False)

# Network propagation
np.fill_diagonal(adj, 0)
D_inv = np.diag(1.0 / (np.sum(np.abs(adj), axis=1) + 1e-8))
L = D_inv @ adj

delta = np.zeros(n_genes)
delta[pomc_idx] = -pomc_mean

propagated = delta.copy()
for step in range(3):
    propagated = 0.7 * L @ propagated + 0.3 * delta

# Add GBR importance
for j, gene in enumerate(core_genes):
    if gene == 'Pomc':
        continue
    pomc_target_imp = pomc_targets[pomc_targets['target']==gene]['importance'].values
    pomc_reg_imp = pomc_regulators[pomc_regulators['TF']==gene]['importance'].values
    grn_importance = 0
    if len(pomc_target_imp) > 0:
        grn_importance = pomc_target_imp[0]
    if len(pomc_reg_imp) > 0:
        grn_importance = max(grn_importance, pomc_reg_imp[0])
    mask = perturb_df['gene'] == gene
    perturb_df.loc[mask, 'propagated_effect'] = propagated[j]
    perturb_df.loc[mask, 'grnboost2_importance'] = grn_importance

perturb_df['combined_score'] = (
    abs(perturb_df['predicted_log2FC'].fillna(0)) * 0.4 +
    abs(perturb_df['propagated_effect'].fillna(0)) * 0.3 +
    (perturb_df['grnboost2_importance'].fillna(0) / max(perturb_df['grnboost2_importance'].max(), 1e-10)) * 0.3
)

perturb_df.to_csv(f'{TABLES_DIR}/Pomc_KO_GRNBoost2.csv', index=False)
print('\nTop 15 predicted effects of Pomc KO:')
for _, row in perturb_df.head(15).iterrows():
    print(f'  {row["gene"]:15s}  direct={row["predicted_log2FC"]:+.3f}  '
          f'propagated={row["propagated_effect"]:+.4f}  GBR_imp={row["grnboost2_importance"]:.6f}')

# ============================================================
# 7. Figures
# ============================================================

# 7.1 Pomc Regulatory Network
fig, ax = plt.subplots(figsize=(16, 14))
pos = nx.spring_layout(G, k=2.5, iterations=200, seed=42)

node_colors, node_sizes = [], []
for node in G.nodes():
    ntype = G.nodes[node].get('type', 'target')
    if node == 'Pomc':
        node_colors.append('#E64B35')
        node_sizes.append(2500)
    elif ntype == 'tf':
        node_colors.append('#4DBBD5')
        node_sizes.append(1400)
    else:
        node_colors.append('#AAAAAA')
        node_sizes.append(900)

edge_colors, edge_widths = [], []
for u, v, d in G.edges(data=True):
    if d.get('edge_type') == 'regulates':
        edge_colors.append('#E64B35')
        edge_widths.append(d['weight'] * 2)
    elif d.get('edge_type') == 'targets':
        edge_colors.append('#4DBBD5')
        edge_widths.append(d['weight'] * 2)
    else:
        edge_colors.append('#CCCCCC')
        edge_widths.append(d['weight'] * 1.5)
edge_widths = [max(0.5, min(6, w)) for w in edge_widths]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes,
                       alpha=0.9, edgecolors='white', linewidths=1.5, ax=ax)
nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.6, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#E64B35', markersize=15, label='Pomc (Hub)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#4DBBD5', markersize=12, label='Transcription Factor'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#AAAAAA', markersize=10, label='Target Gene'),
    Line2D([0],[0], color='#E64B35', linewidth=2, label='Regulates Pomc'),
    Line2D([0],[0], color='#4DBBD5', linewidth=2, label='Pomc Targets'),
    Line2D([0],[0], color='#CCCCCC', linewidth=1, label='TF Co-regulation'),
]
ax.legend(handles=legend_elements, loc='upper left', frameon=True, fontsize=8)
ax.set_title('Pomc Gene Regulatory Network\n(Gradient Boosting Regression, GRNBoost2-equivalent)', fontsize=16, fontweight='bold')
ax.axis('off')
plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/GRNBoost2_Pomc_Network.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure saved: GRNBoost2_Pomc_Network.png')

# 7.2 Pomc KO Prediction
fig, ax = plt.subplots(figsize=(12, 8))
top_ko = perturb_df.head(25).sort_values('predicted_log2FC')
colors = ['#E64B35' if x < 0 else '#4DBBD5' for x in top_ko['predicted_log2FC']]
ax.barh(range(len(top_ko)), top_ko['predicted_log2FC'].values, color=colors, alpha=0.85)
ax.set_yticks(range(len(top_ko)))
ax.set_yticklabels(top_ko['gene'].values, fontsize=9)
ax.set_xlabel('Predicted log2 Fold Change (Pomc KO)', fontsize=12)
ax.axvline(0, color='black', linewidth=0.8)
for i, (_, row) in enumerate(top_ko.iterrows()):
    if abs(row['propagated_effect']) > 0.001:
        ax.scatter(row['propagated_effect'], i, color='#FFD700', s=60, zorder=5,
                   edgecolors='black', linewidths=0.5)
ax.set_title('In Silico Pomc Knockout: Predicted Downstream Effects\n(Gradient Boosting + Network Propagation)', fontsize=14, fontweight='bold')
from matplotlib.lines import Line2D
legend_elements2 = [
    Line2D([0],[0], color='#E64B35', linewidth=3, label='Predicted Downregulation'),
    Line2D([0],[0], color='#4DBBD5', linewidth=3, label='Predicted Upregulation'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#FFD700', markersize=8,
           markeredgecolor='black', label='Propagated Effect'),
]
ax.legend(handles=legend_elements2, loc='lower left', fontsize=9, frameon=True)
sns.despine()
plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/GRNBoost2_Pomc_KO_Prediction.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure saved: GRNBoost2_Pomc_KO_Prediction.png')

# 7.3 TF Connectivity
fig, ax = plt.subplots(figsize=(10, 8))
top_tfs = tf_connectivity.head(20)
colors_tf = ['#E64B35' if idx == 'Pomc' else '#4DBBD5' for idx in top_tfs.index]
ax.barh(range(len(top_tfs)), top_tfs['n_targets'].values[::-1], color=colors_tf[::-1], alpha=0.85)
ax.set_yticks(range(len(top_tfs)))
ax.set_yticklabels(top_tfs.index[::-1], fontsize=10)
ax.set_xlabel('Number of Target Genes', fontsize=12)
ax.set_title('TF Connectivity in Hypothalamus\n(Gradient Boosting Regression)', fontsize=14, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/GRNBoost2_TF_Connectivity.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure saved: GRNBoost2_TF_Connectivity.png')

# 7.4 Combined: GBR importance vs Spearman correlation
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Left: Top TFs by mean importance
top15_tfs = tf_connectivity.head(15)
axes[0].barh(range(len(top15_tfs)), top15_tfs['mean_importance'].values[::-1],
             color=['#E64B35' if idx == 'Pomc' else '#4DBBD5' for idx in top15_tfs.index[::-1]])
axes[0].set_yticks(range(len(top15_tfs)))
axes[0].set_yticklabels(top15_tfs.index[::-1])
axes[0].set_xlabel('Mean GBR Importance', fontsize=11)
axes[0].set_title('TF Regulatory Importance (GBR)', fontsize=12, fontweight='bold')
sns.despine()

# Right: Importance distribution
axes[1].hist(np.log10(links['importance'].values + 1e-10), bins=60, color='#4DBBD5', alpha=0.7, edgecolor='white')
if len(pomc_targets) > 0:
    pomc_mean_imp = pomc_targets['importance'].mean()
    axes[1].axvline(np.log10(pomc_mean_imp + 1e-10), color='#E64B35', linewidth=2, linestyle='--',
                    label=f'Pomc targets (mean={pomc_mean_imp:.2e})')
    axes[1].legend(fontsize=9)
axes[1].set_xlabel('log10(GBR Importance)', fontsize=11)
axes[1].set_ylabel('Number of Links', fontsize=11)
axes[1].set_title('Importance Distribution', fontsize=12, fontweight='bold')
sns.despine()

plt.tight_layout()
plt.savefig(f'{RESULTS_DIR}/GRNBoost2_Importance.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure saved: GRNBoost2_Importance.png')

# ============================================================
# 8. Summary
# ============================================================
print('\n' + '='*60)
print('GRNBoost2-style (GBR) ANALYSIS SUMMARY')
print('='*60)
print(f'Total regulatory links: {len(links)}')
print(f'TFs with targets: {links["TF"].nunique()}')
print(f'Unique target genes: {links["target"].nunique()}')
print(f'Mean links per TF: {len(links)/links["TF"].nunique():.1f}')
print(f'Top 5 connected TFs: {", ".join(tf_connectivity.head(5).index.tolist())}')
print(f'Pomc targets: {len(pomc_targets)}')
print(f'Pomc regulators: {len(pomc_regulators)}')
if len(pomc_targets) > 0:
    print(f'  Top 5 Pomc targets: {", ".join(pomc_targets.head(5)["target"].tolist())}')
if len(pomc_regulators) > 0:
    print(f'  Top 5 Pomc regulators: {", ".join(pomc_regulators.head(5)["TF"].tolist())}')
print(f'Top 5 predicted KO effects: {", ".join(perturb_df.head(5)["gene"].tolist())}')
print(f'\nAll results saved to {RESULTS_DIR}/ and {TABLES_DIR}/')
print('Done!')
