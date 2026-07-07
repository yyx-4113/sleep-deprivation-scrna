"""
Pomc-centered GRN + In Silico Perturbation Analysis
From GSE137665 sleep deprivation scRNA-seq data
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr, pearsonr
from scipy.sparse import issparse

FIG_DIR = 'C:/Users/1/sleep-deprivation-project/results/figures'
TBL_DIR = 'C:/Users/1/sleep-deprivation-project/results/tables'
import os
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TBL_DIR, exist_ok=True)

# ============= 1. Load Data =============
print("Loading data...")
adata = sc.read_h5ad('C:/Users/1/sleep-deprivation-project/data/GSE137665_processed.h5ad')

# Use the full expression (not just HVGs) for co-expression analysis
raw = sc.read_h5ad('C:/Users/1/sleep-deprivation-project/data/GSE137665_raw.h5ad')
# QC filter to match processed data
raw = raw[raw.obs_names.isin(adata.obs_names)]

# Focus on hypothalamus where Pomc is expressed
hypo = raw[raw.obs['brain_region'] == 'Hypothalamus'].copy()
sc.pp.normalize_total(hypo, target_sum=1e4)
sc.pp.log1p(hypo)

print(f'Hypothalamus subset: {hypo.n_obs} cells, {hypo.n_vars} genes')

# ============= 2. Pomc Expression Pattern =============
print("\n=== Pomc Expression Analysis ===")
pomc_expr = hypo[:, 'Pomc'].X
if issparse(pomc_expr):
    pomc_expr = pomc_expr.toarray().flatten()
hypo.obs['Pomc_expr'] = pomc_expr

print(f'Pomc expression: mean={pomc_expr.mean():.3f}, '
      f'cells with Pomc>0: {(pomc_expr>0).sum()}/{len(pomc_expr)}')

# Pomc across conditions
for code in ['A1','A2','A3']:
    mask = hypo.obs['condition_code'] == code
    print(f'  {code}: mean Pomc={pomc_expr[mask].mean():.3f}')

# ============= 3. Co-expression Analysis =============
print("\n=== Co-expression Network ===")

# Get Pomc-correlated genes
expr_matrix = hypo.X
if issparse(expr_matrix):
    expr_matrix = expr_matrix.toarray()

pomc_idx = list(hypo.var_names).index('Pomc')
pomc_vec = expr_matrix[:, pomc_idx]

# Calculate Spearman correlation (faster, more robust for scRNA)
n_genes = min(5000, hypo.n_vars)  # Top 5000 expressed genes
gene_means = np.array(expr_matrix.mean(axis=0)).flatten()
top_genes_idx = np.argsort(gene_means)[-n_genes:]

correlations = []
for i in top_genes_idx:
    if i == pomc_idx:
        continue
    gene_vec = expr_matrix[:, i]
    if np.std(gene_vec) == 0:
        continue
    r, p = spearmanr(pomc_vec, gene_vec)
    correlations.append({
        'gene': hypo.var_names[i],
        'spearman_r': r,
        'pvalue': p,
        'mean_expr': gene_means[i]
    })

corr_df = pd.DataFrame(correlations)
corr_df['abs_r'] = abs(corr_df['spearman_r'])
corr_df['padj'] = corr_df['pvalue'].apply(lambda x: min(x * n_genes, 1.0))  # Bonferroni
corr_df = corr_df.sort_values('abs_r', ascending=False)

print(f'Top 20 Pomc-correlated genes:')
for _, row in corr_df.head(20).iterrows():
    direction = 'pos' if row['spearman_r'] > 0 else 'neg'
    print(f'  {row["gene"]:15s}  r={row["spearman_r"]:+.3f}  padj={row["padj"]:.2e}  [{direction}]')

corr_df.to_csv(f'{TBL_DIR}/Pomc_coexpression.csv', index=False)

# ============= 4. Identify TFs in Pomc Network =============
print("\n=== TF Identification ===")

# Known mouse TFs from AnimalTFDB / literature
MOUSE_TFS = {
    'Nr3c1', 'Crem', 'Fos', 'Fosb', 'Jun', 'Junb', 'Jund',
    'Atf3', 'Atf4', 'Mef2c', 'Mef2d', 'Srebf1', 'Srebf2',
    'Nr4a1', 'Nr4a2', 'Nr4a3', 'Egr1', 'Egr2', 'Egr3',
    'Klf4', 'Klf9', 'Sox2', 'Sox9', 'Mafb', 'Cebpb',
    'Hlf', 'Zic1', 'Npas4', 'Per1', 'Per2', 'Cry1', 'Cry2',
    'Clock', 'Bmal1', 'Arntl', 'Nfil3', 'Dbp', 'Tef',
    'Stat3', 'Stat5a', 'Nfkb1', 'Rela', 'Sp1', 'Sp3',
    'Yy1', 'Ctcf', 'Pparg', 'Ppargc1a', 'Foxo1', 'Foxo3',
    'Rest', 'Rcor1', 'Hdac1', 'Hdac2', 'Mecp2',
    'Olig1', 'Olig2', 'Sox10', 'Sox8',
    'Tbr1', 'Neurod1', 'Neurod2', 'Neurod6',
    'Isl1', 'Lhx6', 'Dlx1', 'Dlx2', 'Dlx5',
    'Nkx2-1', 'Otp', 'Sim1', 'Pou3f2', 'Pou3f4',
    'Bcl11b', 'Satb2', 'Cux1', 'Cux2', 'Fezf2',
    'Nr5a1', 'Nr5a2', 'Esr1', 'Esr2', 'Ar',
    'Rxrg', 'Rorb', 'Rora', 'Rorc', 'Nr1d1', 'Nr1d2',
    'Thrb', 'Ppara', 'Ppard',
}

# Check which of our correlated genes are TFs
pomc_network = corr_df[corr_df['abs_r'] > 0.3].copy()  # |r| > 0.3
pomc_tfs = pomc_network[pomc_network['gene'].isin(MOUSE_TFS)].copy()
print(f'\nTFs in Pomc co-expression network (|r|>0.3): {len(pomc_tfs)}')
for _, row in pomc_tfs.iterrows():
    print(f'  {row["gene"]:15s}  r={row["spearman_r"]:+.3f}')

pomc_tfs.to_csv(f'{TBL_DIR}/Pomc_network_TFs.csv', index=False)

# ============= 5. Build Local GRN for Pomc =============
print("\n=== Building Local GRN ===")

# Core network genes: Pomc + top correlated TFs + top correlated genes
core_tfs = pomc_tfs.nlargest(15, 'abs_r')['gene'].tolist()
core_genes_list = ['Pomc'] + core_tfs + corr_df.nlargest(20, 'abs_r')['gene'].tolist()
core_genes_list = list(dict.fromkeys(core_genes_list))  # unique, preserve order
core_genes_list = core_genes_list[:40]  # keep manageable

print(f'Core network: {len(core_genes_list)} genes')

# Build partial correlation network among core genes
core_idx = [list(hypo.var_names).index(g) for g in core_genes_list if g in hypo.var_names]
core_names = [hypo.var_names[i] for i in core_idx]
core_expr = expr_matrix[:, core_idx]

# Create adjacency matrix (Spearman correlation + threshold)
adj = np.zeros((len(core_idx), len(core_idx)))
for i in range(len(core_idx)):
    for j in range(i+1, len(core_idx)):
        if np.std(core_expr[:, i]) > 0 and np.std(core_expr[:, j]) > 0:
            r, _ = spearmanr(core_expr[:, i], core_expr[:, j])
            adj[i, j] = r
            adj[j, i] = r

# ============= 6. Simulate Pomc Knockout =============
print("\n=== In Silico Pomc Knockout ===")

# Method: Set Pomc to 0, propagate through network via correlation structure
# This is a simplified version of the CellOracle approach

# Get Pomc's position in the network
pomc_pos = core_names.index('Pomc')

# Perturbation: set Pomc to 0
# Predicted effect on gene j = r(Pomc, j) * mean(Pomc) * sign
perturbed_effects = []
for j, gene in enumerate(core_names):
    if gene == 'Pomc':
        continue
    r = adj[pomc_pos, j]
    pomc_mean = pomc_expr.mean()
    # Simple linear model: delta_expr_j = r_ij * delta_expr_i
    predicted_fc = -r * pomc_mean  # negative because we're removing Pomc
    perturbed_effects.append({
        'gene': gene,
        'correlation_with_Pomc': r,
        'predicted_log2FC': predicted_fc,
        'direction': 'Down' if predicted_fc < 0 else 'Up',
        'is_TF': gene in MOUSE_TFS
    })

perturb_df = pd.DataFrame(perturbed_effects)
perturb_df['abs_predicted_FC'] = abs(perturb_df['predicted_log2FC'])
perturb_df = perturb_df.sort_values('abs_predicted_FC', ascending=False)

print('Top 15 predicted effects of Pomc knockout:')
for _, row in perturb_df.head(15).iterrows():
    tf_tag = ' [TF]' if row['is_TF'] else ''
    print(f'  {row["gene"]:15s}  predicted log2FC={row["predicted_log2FC"]:+.3f}  '
          f'r={row["correlation_with_Pomc"]:+.3f}{tf_tag}')

perturb_df.to_csv(f'{TBL_DIR}/Pomc_knockout_prediction.csv', index=False)

# ============= 7. Visualization =============
print("\n=== Generating Figures ===")

# Fig A: Pomc expression UMAP
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(adata[adata.obs['brain_region']=='Hypothalamus'], color='Pomc',
           ax=ax, show=False, title='Pomc Expression in Hypothalamus',
           cmap='Reds', vmax='p99')
plt.savefig(f'{FIG_DIR}/Pomc_Expression_UMAP.png', dpi=300, bbox_inches='tight')
plt.close()
print('  Pomc UMAP saved')

# Fig B: Pomc co-expression network
# Use top 15 genes for a clean network visualization
import networkx as nx

top_net_genes = ['Pomc'] + perturb_df.nlargest(15, 'abs_predicted_FC')['gene'].tolist()
net_idx = [core_names.index(g) for g in top_net_genes if g in core_names]
net_names = [core_names[i] for i in net_idx]

G = nx.Graph()
for name in net_names:
    is_tf = name in MOUSE_TFS
    G.add_node(name, is_tf=is_tf, is_pomc=(name=='Pomc'))

for i in range(len(net_idx)):
    for j in range(i+1, len(net_idx)):
        r = adj[net_idx[i], net_idx[j]]
        if abs(r) > 0.2:
            G.add_edge(net_names[i], net_names[j], weight=abs(r), sign='+' if r>0 else '-')

fig, ax = plt.subplots(figsize=(16, 14))
pos = nx.circular_layout(G)

# Node colors
node_colors = []
node_sizes = []
for node in G.nodes():
    if node == 'Pomc':
        node_colors.append('#E64B35')  # Red for Pomc
        node_sizes.append(2000)
    elif G.nodes[node].get('is_tf', False):
        node_colors.append('#4DBBD5')  # Blue for TFs
        node_sizes.append(1200)
    else:
        node_colors.append('#CCCCCC')  # Grey for other genes
        node_sizes.append(800)

# Edge colors
edge_colors = ['#E64B35' if G[u][v]['sign']=='+' else '#4DBBD5' for u, v in G.edges()]
edge_widths = [abs(G[u][v]['weight'])*3 for u, v in G.edges()]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.9, ax=ax)
nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax)

# Legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#E64B35', markersize=15, label='Pomc (target)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#4DBBD5', markersize=12, label='Transcription Factor'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#CCCCCC', markersize=10, label='Target Gene'),
    Line2D([0],[0], color='#E64B35', lw=2, label='Positive correlation'),
    Line2D([0],[0], color='#4DBBD5', lw=2, label='Negative correlation'),
]
ax.legend(handles=legend_elements, loc='upper left', frameon=False)
ax.set_title('Pomc Regulatory Network in Hypothalamus\n(Sleep Deprivation scRNA-seq)', fontsize=14, fontweight='bold')
ax.axis('off')
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/Pomc_GRN_Network.png', dpi=300, bbox_inches='tight')
plt.close()
print('  Network saved')

# Fig C: Perturbation bar plot
fig, ax = plt.subplots(figsize=(10, 6))
top20 = perturb_df.head(20).copy()
colors = ['#E64B35' if x < 0 else '#4DBBD5' for x in top20['predicted_log2FC']]
ax.barh(range(len(top20)), top20['predicted_log2FC'].values[::-1], color=colors[::-1])
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['gene'].values[::-1], fontsize=9)
ax.set_xlabel('Predicted log2 Fold Change (Pomc KO)')
ax.axvline(0, color='black', linewidth=0.5)
ax.set_title('Predicted Effects of In Silico Pomc Knockout', fontsize=13, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/Pomc_KO_Prediction.png', dpi=300, bbox_inches='tight')
plt.close()
print('  KO prediction saved')

# Fig D: Pomc expression across conditions
fig, ax = plt.subplots(figsize=(6, 5))
condition_data = [pomc_expr[hypo.obs['condition_code']==c] for c in ['A1','A2','A3']]
bp = ax.boxplot(condition_data, labels=['A1','A2','A3'], patch_artist=True)
for patch, color in zip(bp['boxes'], ['#4DBBD5','#E64B35','#00A087']):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
ax.set_ylabel('Pomc Expression (log-normalized)')
ax.set_title('Pomc Expression Across Sleep Conditions\n(Hypothalamus)', fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/Pomc_Conditions_Boxplot.png', dpi=300, bbox_inches='tight')
plt.close()
print('  Boxplot saved')

print(f'\n=== GRN Analysis Complete ===')
print(f'All figures saved to {FIG_DIR}')
print(f'All tables saved to {TBL_DIR}')
print(f'Key TF targets for Pomc: {core_tfs[:10]}')
