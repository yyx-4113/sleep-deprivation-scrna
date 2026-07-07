"""
Step 3: Pomc GRN + In Silico Knockout
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import warnings, os
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from scipy import sparse
from matplotlib.lines import Line2D

FIG = r'C:\Users\1\sleep-deprivation-project\results\figures'
TBL = r'C:\Users\1\sleep-deprivation-project\results\tables'
os.makedirs(FIG, exist_ok=True)
os.makedirs(TBL, exist_ok=True)

print('Loading data...')
proc = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_processed.h5ad')
raw = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_raw.h5ad')
raw = raw[raw.obs_names.isin(proc.obs_names)]

hypo = raw[raw.obs['brain_region'] == 'Hypothalamus'].copy()
sc.pp.normalize_total(hypo, target_sum=1e4)
sc.pp.log1p(hypo)
print(f'Hypothalamus: {hypo.n_obs} cells')

# Pomc expression
pomc_expr = hypo[:, 'Pomc'].X
if sparse.issparse(pomc_expr):
    pomc_expr = pomc_expr.toarray().flatten()
print(f'Pomc: mean={pomc_expr.mean():.3f}, nonzero={(pomc_expr>0).sum()}/{len(pomc_expr)}')
for code in ['A1', 'A2', 'A3']:
    m = hypo.obs['condition_code'] == code
    print(f'  {code}: Pomc mean={pomc_expr[m].mean():.3f}')

# Co-expression
print('Computing co-expression...')
expr = hypo.X
if sparse.issparse(expr):
    expr = expr.toarray()
pomc_idx = list(hypo.var_names).index('Pomc')
pomc_vec = expr[:, pomc_idx]
gene_means = expr.mean(axis=0)
top_idx = np.argsort(gene_means)[-5000:]

correlations = []
for i in top_idx:
    if i == pomc_idx:
        continue
    gv = expr[:, i]
    if np.std(gv) == 0:
        continue
    r, p = spearmanr(pomc_vec, gv)
    correlations.append({'gene': hypo.var_names[i], 'spearman_r': r, 'pvalue': p})

corr_df = pd.DataFrame(correlations).sort_values('spearman_r', key=abs, ascending=False)
corr_df['padj'] = (corr_df['pvalue'].clip(lower=1e-300) * len(corr_df)).clip(upper=1.0)
corr_df.to_csv(os.path.join(TBL, 'Pomc_coexpression.csv'), index=False)

print('Top 20 Pomc-correlated:')
for _, row in corr_df.head(20).iterrows():
    print(f'  {row["gene"]:15s}  r={row["spearman_r"]:+.3f}  padj={row["padj"]:.2e}')

# TFs
MOUSE_TFS = {
    'Nr3c1', 'Crem', 'Fos', 'Fosb', 'Jun', 'Junb', 'Jund', 'Atf3', 'Atf4',
    'Mef2c', 'Mef2d', 'Srebf1', 'Srebf2', 'Nr4a1', 'Nr4a2', 'Nr4a3',
    'Egr1', 'Egr2', 'Egr3', 'Klf4', 'Klf9', 'Sox2', 'Sox9', 'Mafb', 'Cebpb',
    'Hlf', 'Zic1', 'Npas4', 'Per1', 'Per2', 'Cry1', 'Cry2', 'Clock', 'Arntl',
    'Nfil3', 'Dbp', 'Tef', 'Stat3', 'Nfkb1', 'Rela', 'Sp1', 'Sp3',
    'Yy1', 'Ctcf', 'Pparg', 'Ppargc1a', 'Foxo1', 'Foxo3', 'Rest',
    'Hdac1', 'Hdac2', 'Mecp2', 'Olig1', 'Olig2', 'Sox10',
    'Tbr1', 'Neurod1', 'Neurod2', 'Neurod6', 'Lhx6', 'Dlx1', 'Dlx2',
    'Nkx2-1', 'Sim1', 'Bcl11b', 'Satb2', 'Rxrg', 'Rorb', 'Rora', 'Rorc',
    'Nr1d1', 'Nr1d2', 'Thrb', 'Ppara', 'Nr5a1', 'Esr1',
}

pomc_net = corr_df[abs(corr_df['spearman_r']) > 0.3]
pomc_tfs = pomc_net[pomc_net['gene'].isin(MOUSE_TFS)]
pomc_tfs.to_csv(os.path.join(TBL, 'Pomc_network_TFs.csv'), index=False)
print(f'\nTFs in Pomc network (|r|>0.3): {len(pomc_tfs)}')
for _, row in pomc_tfs.iterrows():
    print(f'  {row["gene"]:15s}  r={row["spearman_r"]:+.3f}')

# Build GRN
core_names = ['Pomc'] + pomc_tfs.head(15)['gene'].tolist() + corr_df.head(20)['gene'].tolist()
core_names = list(dict.fromkeys(core_names))[:35]
core_idx = [list(hypo.var_names).index(g) for g in core_names if g in hypo.var_names]
core_names = [hypo.var_names[i] for i in core_idx]
core_expr = expr[:, core_idx]

adj = np.zeros((len(core_idx), len(core_idx)))
for i in range(len(core_idx)):
    for j in range(i + 1, len(core_idx)):
        if np.std(core_expr[:, i]) > 0 and np.std(core_expr[:, j]) > 0:
            r, p = spearmanr(core_expr[:, i], core_expr[:, j])
            adj[i, j] = r
            adj[j, i] = r

pomc_pos = core_names.index('Pomc')
pomc_mean = pomc_expr.mean()

# Perturbation simulation
perturb = []
for j, gene in enumerate(core_names):
    if gene == 'Pomc':
        continue
    r = adj[pomc_pos, j]
    pred_fc = -r * pomc_mean
    perturb.append({
        'gene': gene,
        'correlation': r,
        'predicted_log2FC': pred_fc,
        'direction': 'Down' if pred_fc < 0 else 'Up',
        'is_TF': gene in MOUSE_TFS
    })

perturb_df = pd.DataFrame(perturb).sort_values('predicted_log2FC', key=abs, ascending=False)
perturb_df.to_csv(os.path.join(TBL, 'Pomc_knockout_prediction.csv'), index=False)

print('\n=== Predicted Pomc KO Effects ===')
for _, row in perturb_df.head(15).iterrows():
    tf_tag = ' [TF]' if row['is_TF'] else ''
    print(f'  {row["gene"]:15s}  log2FC={row["predicted_log2FC"]:+.3f}{tf_tag}')

# ---- Figures ----
print('\nGenerating figures...')

# 1. Hypothalamus UMAP
fig, ax = plt.subplots(figsize=(8, 7))
sc.pl.umap(proc[proc.obs['brain_region'] == 'Hypothalamus'], color='cell_type',
           ax=ax, show=False, title='Hypothalamus Cell Types', legend_loc='right margin')
plt.savefig(os.path.join(FIG, 'Hypothalamus_UMAP.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  UMAP saved')

# 2. Network
top_net = ['Pomc'] + perturb_df.head(15)['gene'].tolist()
net_idx_sub = [core_names.index(g) for g in top_net if g in core_names]
net_names_sub = [core_names[i] for i in net_idx_sub]

G = nx.Graph()
for name in net_names_sub:
    G.add_node(name)
for i in range(len(net_idx_sub)):
    for j in range(i + 1, len(net_idx_sub)):
        r = adj[net_idx_sub[i], net_idx_sub[j]]
        if abs(r) > 0.2:
            G.add_edge(net_names_sub[i], net_names_sub[j], weight=abs(r), sign='+' if r > 0 else '-')

fig, ax = plt.subplots(figsize=(16, 14))
pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
node_colors = ['#E64B35' if n == 'Pomc' else '#4DBBD5' if n in MOUSE_TFS else '#CCCCCC' for n in G.nodes()]
node_sizes = [2000 if n == 'Pomc' else 1200 if n in MOUSE_TFS else 800 for n in G.nodes()]
edge_colors = ['#E64B35' if G[u][v]['sign'] == '+' else '#4DBBD5' for u, v in G.edges()]
edge_widths = [abs(G[u][v]['weight']) * 3 for u, v in G.edges()]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.9, ax=ax)
nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=edge_widths, alpha=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', ax=ax)
ax.legend(handles=[
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#E64B35', markersize=15, label='Pomc (Target)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#4DBBD5', markersize=12, label='Transcription Factor'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#CCCCCC', markersize=10, label='Target Gene'),
    Line2D([0], [0], color='#E64B35', lw=2, label='Positive r'),
    Line2D([0], [0], color='#4DBBD5', lw=2, label='Negative r'),
], loc='upper left', frameon=False)
ax.set_title('Pomc Regulatory Network - Hypothalamus (Sleep Deprivation scRNA-seq)', fontsize=14, fontweight='bold')
ax.axis('off')
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'Pomc_GRN_Network.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Network saved')

# 3. KO prediction barplot
fig, ax = plt.subplots(figsize=(10, 6))
top20 = perturb_df.head(20)
colors = ['#E64B35' if x < 0 else '#4DBBD5' for x in top20['predicted_log2FC']]
ax.barh(range(len(top20)), top20['predicted_log2FC'].values[::-1], color=colors[::-1])
ax.set_yticks(range(len(top20)))
ax.set_yticklabels(top20['gene'].values[::-1], fontsize=9)
ax.set_xlabel('Predicted log2 Fold Change')
ax.axvline(0, color='black', linewidth=0.5)
ax.set_title('In Silico Pomc Knockout: Predicted Downstream Effects', fontsize=13, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'Pomc_KO_Prediction.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  KO prediction saved')

# 4. Pomc boxplot
fig, ax = plt.subplots(figsize=(6, 5))
data = [pomc_expr[hypo.obs['condition_code'] == c] for c in ['A1', 'A2', 'A3']]
bp = ax.boxplot(data, labels=['A1', 'A2', 'A3'], patch_artist=True)
for patch, c in zip(bp['boxes'], ['#4DBBD5', '#E64B35', '#00A087']):
    patch.set_facecolor(c)
    patch.set_alpha(0.7)
ax.set_ylabel('Pomc Expression (log-norm)')
ax.set_title('Pomc Across Sleep Conditions (Hypothalamus)', fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'Pomc_Conditions_Boxplot.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Boxplot saved')

print(f'\n=== GRN Analysis Complete ===')
print(f'Figures in {FIG}: {os.listdir(FIG)}')
print(f'Tables in {TBL}: {os.listdir(TBL)}')
