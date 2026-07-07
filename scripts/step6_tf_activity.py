"""
Step 6: TF Activity Inference — local regulon approach
Builds TF regulons from our co-expression data, infers TF activity per cell
Equivalent to pySCENIC GRNBoost2 step without needing external motif databases
"""
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings, os
warnings.filterwarnings('ignore')
from scipy.stats import spearmanr
from scipy.sparse import issparse

FIG = r'C:\Users\1\sleep-deprivation-project\results\figures'
TBL = r'C:\Users\1\sleep-deprivation-project\results\tables'
os.makedirs(FIG, exist_ok=True)
os.makedirs(TBL, exist_ok=True)

# ============= 1. Load data =============
print("Loading data...")
raw = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_raw.h5ad')
proc = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_processed.h5ad')
raw = raw[raw.obs_names.isin(proc.obs_names)]

# Focus on hypothalamus
hypo = raw[raw.obs['brain_region'] == 'Hypothalamus'].copy()
sc.pp.normalize_total(hypo, target_sum=1e4)
sc.pp.log1p(hypo)
print(f"Hypothalamus: {hypo.n_obs} cells, {hypo.n_vars} genes")

expr = hypo.X
if issparse(expr):
    expr = expr.toarray()

# Load our Pomc network TFs
pomc_tfs_df = pd.read_csv(os.path.join(TBL, 'Pomc_network_TFs.csv'))
pomc_tfs = pomc_tfs_df['gene'].tolist() if len(pomc_tfs_df) > 0 else []

# Extended TF list from literature
MOUSE_TFS = {
    'Nr3c1','Crem','Fos','Fosb','Jun','Junb','Jund','Atf3','Atf4',
    'Mef2c','Mef2d','Srebf1','Srebf2','Nr4a1','Nr4a2','Nr4a3',
    'Egr1','Egr2','Egr3','Klf4','Klf9','Sox2','Sox9','Mafb','Cebpb',
    'Hlf','Zic1','Npas4','Per1','Per2','Cry1','Cry2','Clock','Arntl',
    'Nfil3','Dbp','Tef','Stat3','Nfkb1','Rela','Sp1','Sp3',
    'Yy1','Ctcf','Pparg','Ppargc1a','Foxo1','Foxo3','Rest',
    'Hdac1','Hdac2','Mecp2','Olig1','Olig2','Sox10',
    'Tbr1','Neurod1','Neurod2','Neurod6','Lhx6','Dlx1','Dlx2',
    'Nkx2-1','Sim1','Bcl11b','Satb2','Rxrg','Rorb','Rora','Rorc',
    'Nr1d1','Nr1d2','Thrb','Ppara','Nr5a1','Esr1','Pomc',
    'Rbm3','Tsc22d3', 'Stat5a', 'Rcor1', 'Cux1', 'Cux2', 'Fezf2',
    'Nr5a2', 'Esr2', 'Ar', 'Ppard'
}

# TFs present in data
tfs_present = sorted([t for t in MOUSE_TFS if t in hypo.var_names])
print(f"TFs in data: {len(tfs_present)}")

# ============= 2. Build TF regulons from co-expression =============
print("\n=== Building TF Regulons ===")
# For each TF, find correlated genes and build a signed regulon
# This is equivalent to SCENIC step 1 (GRNBoost2)

# Use top 5000 expressed genes for efficiency
gene_means = expr.mean(axis=0)
top_gene_idx = np.argsort(gene_means)[-5000:]
top_genes = [hypo.var_names[i] for i in top_gene_idx]

regulons = []
for tf in tfs_present:
    if tf not in hypo.var_names:
        continue
    tf_idx = list(hypo.var_names).index(tf)
    tf_vec = expr[:, tf_idx]
    if np.std(tf_vec) == 0:
        continue

    # Find targets among top expressed genes
    targets_pos = []
    targets_neg = []
    for gi_idx in top_gene_idx:
        gi = hypo.var_names[gi_idx]
        if gi == tf:
            continue
        gv = expr[:, gi_idx]
        if np.std(gv) == 0:
            continue
        r, p = spearmanr(tf_vec, gv)
        if abs(r) > 0.25 and p < 0.01:
            if r > 0:
                targets_pos.append((gi, r))
            else:
                targets_neg.append((gi, r))

    # Keep top 50 positive + 50 negative targets per TF
    targets_pos = sorted(targets_pos, key=lambda x: -x[1])[:50]
    targets_neg = sorted(targets_neg, key=lambda x: x[1])[:50]

    for target, r in targets_pos + targets_neg:
        regulons.append({'TF': tf, 'target': target, 'weight': r})

reg_df = pd.DataFrame(regulons)
print(f"Regulons: {len(reg_df)} TF-target relationships")
print(f"TFs with regulons: {reg_df['TF'].nunique()}")

# Save regulon database
reg_df.to_csv(os.path.join(TBL, 'TF_regulons.csv'), index=False)

# Top TFs by regulon size
tf_sizes = reg_df.groupby('TF').size().sort_values(ascending=False)
print(f"\nTop 10 TFs by regulon size:")
for tf, n in tf_sizes.head(10).items():
    print(f"  {tf:15s}: {n} targets")

# ============= 3. TF Activity Scoring (AUCell-like) =============
print("\n=== Computing TF Activity Scores ===")

# For each cell, compute TF activity as weighted mean of regulon gene expression
n_cells = hypo.n_obs
tf_activity = np.zeros((n_cells, len(tfs_present)))
tf_list_sorted = sorted(tfs_present)

for tf_idx, tf in enumerate(tf_list_sorted):
    tf_reg = reg_df[reg_df['TF'] == tf]
    if len(tf_reg) == 0:
        continue

    # Get target gene indices and weights
    target_idx = []
    weights = []
    for _, row in tf_reg.iterrows():
        if row['target'] in hypo.var_names:
            target_idx.append(list(hypo.var_names).index(row['target']))
            weights.append(row['weight'])

    if len(target_idx) == 0:
        continue

    weights = np.array(weights)
    target_expr = expr[:, target_idx]

    # Weighted mean (normalize weights)
    weights_norm = weights / (np.abs(weights).sum() + 1e-8)
    tf_activity[:, tf_idx] = target_expr @ weights_norm

print(f"TF activity matrix: {tf_activity.shape}")

# Store in AnnData
acts = sc.AnnData(X=tf_activity, obs=hypo.obs.copy(), var=pd.DataFrame(index=tf_list_sorted))
print(f"Activity AnnData: {acts.n_obs} cells x {acts.n_vars} TFs")

# ============= 4. Differential TF Activity =============
print("\n=== Differential TF Activity ===")

acts.obs['condition_code'] = hypo.obs['condition_code'].values

# Compare A2 (SD) vs A1 (Normal) for each TF
diff_tf = []
for tf in tf_list_sorted:
    a1 = acts[acts.obs['condition_code'] == 'A1', tf].X.toarray().flatten() if issparse(acts[:, tf].X) else acts[acts.obs['condition_code'] == 'A1', tf].X.flatten()
    a2 = acts[acts.obs['condition_code'] == 'A2', tf].X.toarray().flatten() if issparse(acts[:, tf].X) else acts[acts.obs['condition_code'] == 'A2', tf].X.flatten()
    a3 = acts[acts.obs['condition_code'] == 'A3', tf].X.toarray().flatten() if issparse(acts[:, tf].X) else acts[acts.obs['condition_code'] == 'A3', tf].X.flatten()

    from scipy.stats import mannwhitneyu
    stat, p = mannwhitneyu(a2, a1, alternative='two-sided')

    diff_tf.append({
        'TF': tf,
        'mean_A1': a1.mean(),
        'mean_A2': a2.mean(),
        'mean_A3': a3.mean(),
        'activity_log2FC': np.log2(a2.mean() + 1e-6) - np.log2(a1.mean() + 1e-6),
        'pvalue': p,
        'regulon_size': tf_sizes.get(tf, 0)
    })

diff_tf_df = pd.DataFrame(diff_tf)
diff_tf_df['padj'] = (diff_tf_df['pvalue'].clip(lower=1e-300) * len(diff_tf_df)).clip(upper=1.0)
diff_tf_df = diff_tf_df.sort_values('padj')
diff_tf_df.to_csv(os.path.join(TBL, 'Differential_TF_Activity.csv'), index=False)

print(f"Differentially active TFs (padj<0.05): {(diff_tf_df['padj']<0.05).sum()}")

print("\nTop 20 differentially active TFs:")
for _, row in diff_tf_df.head(20).iterrows():
    print(f"  {row['TF']:15s}  activity_FC={row['activity_log2FC']:+.4f}  padj={row['padj']:.2e}  regulon={row['regulon_size']}")

# ============= 5. Pomc-TF relationships =============
print("\n=== Pomc Upstream Regulators ===")

# Which TFs' regulons include Pomc?
pomc_regulators = reg_df[reg_df['target'] == 'Pomc']
print(f"TFs regulating Pomc (co-expression |r|>0.25): {len(pomc_regulators)}")
if len(pomc_regulators) > 0:
    for _, row in pomc_regulators.sort_values('weight', key=abs, ascending=False).iterrows():
        tf_act = diff_tf_df[diff_tf_df['TF'] == row['TF']]
        act_fc = tf_act['activity_log2FC'].values[0] if len(tf_act) > 0 else np.nan
        print(f"  {row['TF']:15s}  r={row['weight']:+.3f}  TF_activity_FC={act_fc:+.3f}")

# ============= 6. Figures =============
print("\nGenerating figures...")

# 6.1: Top differentially active TFs
fig, ax = plt.subplots(figsize=(9, 7))
top25 = diff_tf_df.nsmallest(25, 'padj')
colors = ['#E64B35' if x < 0 else '#4DBBD5' for x in top25['activity_log2FC']]
ax.barh(range(len(top25)), top25['activity_log2FC'].values[::-1], color=colors[::-1])
ax.set_yticks(range(len(top25)))
ax.set_yticklabels(top25['TF'].values[::-1], fontsize=9)
ax.set_xlabel('TF Activity log2FC (SD vs Normal)', fontsize=11)
ax.axvline(0, color='black', lw=0.5)
ax.set_title('Differential TF Activity in Hypothalamus\n(Sleep Deprivation vs Normal)', fontsize=13, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'TF_Activity_Differential.png'), dpi=300, bbox_inches='tight')
plt.close()

# 6.2: Pomc and its regulators across conditions
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: Pomc expression
ax = axes[0]
for i, (code, color, label) in enumerate(zip(['A1','A2','A3'],
    ['#4DBBD5','#E64B35','#00A087'],
    ['Normal','Sleep Dep','Recovery'])):
    mask = hypo.obs['condition_code'] == code
    pomc_expr = hypo[mask, 'Pomc'].X.toarray().flatten() if issparse(hypo[:, 'Pomc'].X) else hypo[mask, 'Pomc'].X.flatten()
    ax.bar(i, pomc_expr.mean(), color=color, alpha=0.8,
           yerr=pomc_expr.std()/np.sqrt(len(pomc_expr)), capsize=5)
ax.set_xticks(range(3))
ax.set_xticklabels(['Normal (A1)', 'SD (A2)', 'Recovery (A3)'], fontsize=9)
ax.set_title('Pomc Expression', fontsize=11, fontweight='bold')
ax.set_ylabel('log-norm Expression')

# Panel B: Top Pomc regulators
ax = axes[1]
top_regs = pomc_regulators.sort_values('weight', key=abs, ascending=False).head(6)
if 'Pomc' in tf_list_sorted:
    top_regs = pd.concat([pd.DataFrame([{'TF':'Pomc','weight':1.0}]), top_regs], ignore_index=True)

for i, (_, row) in enumerate(top_regs.iterrows()):
    tf = row['TF']
    if tf in tf_list_sorted:
        a1 = acts[acts.obs['condition_code']=='A1', tf].X.toarray().flatten() if issparse(acts[:, tf].X) else acts[acts.obs['condition_code']=='A1', tf].X.flatten()
        a2 = acts[acts.obs['condition_code']=='A2', tf].X.toarray().flatten() if issparse(acts[:, tf].X) else acts[acts.obs['condition_code']=='A2', tf].X.flatten()
        a3 = acts[acts.obs['condition_code']=='A3', tf].X.toarray().flatten() if issparse(acts[:, tf].X) else acts[acts.obs['condition_code']=='A3', tf].X.flatten()
        for j, (vals, color) in enumerate(zip([a1,a2,a3], ['#4DBBD5','#E64B35','#00A087'])):
            ax.bar(i*3 + j, vals.mean(), color=color, alpha=0.8, width=0.8)
ax.set_xticks([i*3+1 for i in range(len(top_regs))])
ax.set_xticklabels(top_regs['TF'].values, fontsize=9, rotation=30)
ax.set_title('TF Activity Across Conditions', fontsize=11, fontweight='bold')
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], color='#4DBBD5', lw=4, label='Normal'),
    Line2D([0],[0], color='#E64B35', lw=4, label='Sleep Dep'),
    Line2D([0],[0], color='#00A087', lw=4, label='Recovery'),
], loc='upper left', frameon=False, fontsize=8)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'TF_Activity_Pomc_Regulators.png'), dpi=300, bbox_inches='tight')
plt.close()

# 6.3: TF Activity Network
import networkx as nx

fig, ax = plt.subplots(figsize=(14, 12))
G = nx.Graph()
G.add_node('Pomc', type='target_gene', size=2000)

# Add top regulators
for _, row in pomc_regulators.sort_values('weight', key=abs, ascending=False).head(10).iterrows():
    tf = row['TF']
    tf_info = diff_tf_df[diff_tf_df['TF'] == tf]
    act_fc = tf_info['activity_log2FC'].values[0] if len(tf_info) > 0 else 0
    G.add_node(tf, type='tf', activity_fc=act_fc, size=1200 + abs(act_fc)*300)
    G.add_edge(tf, 'Pomc', weight=abs(row['weight']))

# Add Pomc targets from co-expression
pomc_targets = reg_df[reg_df['TF'] == 'Pomc'].sort_values('weight', key=abs, ascending=False).head(8)
for _, row in pomc_targets.iterrows():
    G.add_node(row['target'], type='target', size=600)
    G.add_edge('Pomc', row['target'], weight=abs(row['weight']))

pos = nx.spring_layout(G, k=2, iterations=100, seed=42)
node_colors = []
node_sizes = []
for node in G.nodes():
    nd = G.nodes[node]
    if nd.get('type') == 'target_gene':
        node_colors.append('#E64B35')
        node_sizes.append(2000)
    elif nd.get('type') == 'tf':
        fc = nd.get('activity_fc', 0)
        node_colors.append('#E64B35' if fc < 0 else '#4DBBD5')
        node_sizes.append(nd.get('size', 1000))
    else:
        node_colors.append('#CCCCCC')
        node_sizes.append(600)

edge_colors = ['#E64B35' if G[u][v]['weight'] > 0.3 else '#4DBBD5' for u, v in G.edges()]

nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, alpha=0.9, ax=ax)
nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=[G[u][v]['weight']*3 for u,v in G.edges()], alpha=0.5, ax=ax)
nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)

from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#E64B35', markersize=15, label='Pomc'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#4DBBD5', markersize=12, label='TF (activated in SD)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#E64B35', markersize=12, label='TF (repressed in SD)'),
    Line2D([0],[0], marker='o', color='w', markerfacecolor='#CCCCCC', markersize=8, label='Target Gene'),
], loc='upper left', frameon=False)
ax.set_title('Pomc Regulatory Network with TF Activity\n(Hypothalamus, Sleep Deprivation)', fontsize=14, fontweight='bold')
ax.axis('off')
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'TF_Activity_Network.png'), dpi=300, bbox_inches='tight')
plt.close()

print(f"\n=== TF Activity Analysis Complete ===")
print(f"Figures: TF_Activity_Differential.png, TF_Activity_Pomc_Regulators.png, TF_Activity_Network.png")
print(f"Tables: Differential_TF_Activity.csv, TF_regulons.csv")
print(f"\nKey TFs with altered activity in SD:")
for _, row in diff_tf_df[abs(diff_tf_df['activity_log2FC']) > 0.3].nsmallest(10, 'padj').iterrows():
    direction = 'activated' if row['activity_log2FC'] > 0 else 'repressed'
    print(f"  {row['TF']:15s} {direction:10s}  FC={row['activity_log2FC']:+.3f}  padj={row['padj']:.2e}")
