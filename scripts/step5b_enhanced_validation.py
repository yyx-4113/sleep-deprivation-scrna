"""
Step 5b: Enhanced External Validation
- Meta-analysis across GSE211088 + GSE237419
- Fisher's combined p-value method
- Gene overlap UpSet analysis
- Cross-dataset correlation heatmap
- Consensus DEG identification
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings, os, gzip, re
warnings.filterwarnings('ignore')
from scipy.stats import ttest_ind, spearmanr, combine_pvalues, pearsonr
from scipy.cluster.hierarchy import linkage, leaves_list

FIG = r'C:\Users\1\sleep-deprivation-project\results\figures'
TBL = r'C:\Users\1\sleep-deprivation-project\results\tables'
EXT = r'C:\Users\1\sleep-deprivation-project\data\external'
os.makedirs(FIG, exist_ok=True)
os.makedirs(TBL, exist_ok=True)

# Load Ensembl mapping
symbol_map = pd.read_csv(os.path.join(EXT, 'ensembl_to_symbol.csv'), index_col=0)
symbol_dict = symbol_map['Symbol'].to_dict()

# ============================================================
# 1. MODULE: Cross-dataset meta-analysis
# ============================================================
print("=" * 60)
print("ENHANCED EXTERNAL VALIDATION — META-ANALYSIS")
print("=" * 60)

# Our DEGs from GSE137665
our_degs = pd.read_csv(os.path.join(TBL, 'DEGs_significant.csv'))
our_degs = our_degs[our_degs['pvals_adj'] < 0.05].copy()
print(f"Our study significant DEGs: {len(our_degs)}")

# --- GSE211088 full DEG analysis ---
print("\n[1/3] GSE211088 full analysis...")
df211 = pd.read_csv(os.path.join(EXT, 'GSE211088_RNAseq_gene.txt.gz'), sep='\t', index_col=0)
hc_cols = [c for c in df211.columns if c.startswith('WTHC')]
sd_cols = [c for c in df211.columns if c.startswith('WTSD')]
ens_ids = [x.split('.')[0] for x in df211.index]
mapped_symbols = [symbol_dict.get(eid, None) for eid in ens_ids]
df211['symbol'] = mapped_symbols
df211_sym = df211[df211['symbol'].notna()].copy()
expr_cols = [c for c in df211.columns if c not in ['ensembl_base', 'symbol']]
df211_sym = df211_sym.groupby('symbol')[expr_cols].mean()

# Full DEGs for GSE211088 (all genes not just key genes)
results_211 = []
for gene in df211_sym.index:
    hc_vals = df211_sym.loc[gene, hc_cols].astype(float).values
    sd_vals = df211_sym.loc[gene, sd_cols].astype(float).values
    fc = np.log2((sd_vals.mean() + 0.1) / (hc_vals.mean() + 0.1))
    if hc_vals.std() > 1e-8 or sd_vals.std() > 1e-8:
        _, pv = ttest_ind(sd_vals, hc_vals)
    else:
        pv = 1.0
    results_211.append({'gene': gene, 'log2FC_211': fc, 'pval_211': pv})
deg_211 = pd.DataFrame(results_211)
deg_211['padj_211'] = (deg_211['pval_211'].clip(lower=1e-300) * len(deg_211)).clip(upper=1.0)
deg_211['sig_211'] = (deg_211['padj_211'] < 0.05) & (np.abs(deg_211['log2FC_211']) > 0.3)
print(f"  GSE211088 significant DEGs: {deg_211['sig_211'].sum()}")

# --- GSE237419 full DEG analysis ---
print("\n[2/3] GSE237419 full analysis...")
with gzip.open(os.path.join(EXT, 'GSE237419_gse_HCSDRS_WT_SleepIntegration_salmon.txt.gz'), 'rt') as f:
    h1 = f.readline().strip().split('\t')
    h2 = f.readline().strip().split('\t')

labels_map = {}
for i, label in enumerate(h2):
    label = label.strip()
    if label and label != 'ENSEMBL_ID':
        match = re.match(r'^(HC|SD|RS)(\d+)_', label)
        if match:
            labels_map[i] = {'condition': match.group(1), 'hours': int(match.group(2))}

df237 = pd.read_csv(os.path.join(EXT, 'GSE237419_gse_HCSDRS_WT_SleepIntegration_salmon.txt.gz'),
                     sep='\t', skiprows=2, header=None)
df237.index = df237[0]
df237 = df237.drop(columns=[0])
ens_ids_237 = [x.split('.')[0] for x in df237.index]
mapped_237 = [symbol_dict.get(eid, None) for eid in ens_ids_237]
df237['symbol'] = mapped_237
df237_sym = df237[df237['symbol'].notna()].copy()
expr_cols_237 = [c for c in df237_sym.columns if isinstance(c, (int, np.integer))]
df237_sym = df237_sym.groupby('symbol')[expr_cols_237].mean()

hc5 = [k for k, v in labels_map.items() if v['condition'] == 'HC' and v['hours'] == 5]
sd5 = [k for k, v in labels_map.items() if v['condition'] == 'SD' and v['hours'] == 5]

if len(hc5) >= 2 and len(sd5) >= 2:
    results_237 = []
    for gene in df237_sym.index:
        hc_v = df237_sym.loc[gene, hc5].astype(float).values
        sd_v = df237_sym.loc[gene, sd5].astype(float).values
        fc = np.log2((sd_v.mean() + 0.1) / (hc_v.mean() + 0.1))
        try:
            _, pv = ttest_ind(sd_v, hc_v)
        except:
            pv = 1.0
        results_237.append({'gene': gene, 'log2FC_237': fc, 'pval_237': pv})
    deg_237 = pd.DataFrame(results_237)
    deg_237['padj_237'] = (deg_237['pval_237'].clip(lower=1e-300) * len(deg_237)).clip(upper=1.0)
    deg_237['sig_237'] = (deg_237['padj_237'] < 0.05) & (np.abs(deg_237['log2FC_237']) > 0.3)
    print(f"  GSE237419 significant DEGs: {deg_237['sig_237'].sum()}")
else:
    print("  WARNING: Not enough HC5/SD5 samples")
    deg_237 = pd.DataFrame()

# ============================================================
# 2. META-ANALYSIS: Fisher's combined p-values
# ============================================================
print("\n[3/3] Meta-analysis: Fisher's combined p-value...")

our_dedf = our_degs[['names','logfoldchanges','pvals','pvals_adj','brain_region']].copy()
our_dedf.columns = ['gene','log2FC_our','pval_our','padj_our','brain_region']

# Merge all three datasets
meta = our_dedf[['gene','log2FC_our','pval_our','padj_our']].copy()
if len(deg_211) > 0:
    meta = meta.merge(deg_211[['gene','log2FC_211','pval_211','padj_211']], on='gene', how='outer')
if len(deg_237) > 0:
    meta = meta.merge(deg_237[['gene','log2FC_237','pval_237','padj_237']], on='gene', how='outer')

# Fisher's method: combine p-values where available
pval_cols = [c for c in ['pval_our','pval_211','pval_237'] if c in meta.columns]
n_datasets = len(pval_cols)
print(f"  Datasets available for meta-analysis: {n_datasets}")

fisher_pvals = []
for _, row in meta.iterrows():
    pvals = [row[c] for c in pval_cols if pd.notna(row[c]) and row[c] > 0]
    if len(pvals) >= 2:
        _, combined = combine_pvalues(pvals, method='fisher')
        fisher_pvals.append(combined)
    else:
        fisher_pvals.append(np.nan)

meta['fisher_pval'] = fisher_pvals
meta['fisher_padj'] = meta['fisher_pval'].clip(lower=1e-300) * meta['fisher_pval'].notna().sum()
meta['fisher_padj'] = meta['fisher_padj'].clip(upper=1.0)

# Consensus direction (majority vote)
fc_cols = [c for c in ['log2FC_our','log2FC_211','log2FC_237'] if c in meta.columns]
def consensus_dir(row):
    vals = [row[c] for c in fc_cols if pd.notna(row[c])]
    if len(vals) < 2:
        return np.nan
    up = sum(1 for v in vals if v > 0)
    return 1 if up > len(vals)/2 else -1

meta['consensus_dir'] = meta.apply(consensus_dir, axis=1)

# Filter to meta-significant genes
meta_sig = meta[(meta['fisher_padj'] < 0.05) & meta['consensus_dir'].notna()].copy()
print(f"  Meta-significant genes (Fisher padj<0.05): {len(meta_sig)}")

# Save
meta.to_csv(os.path.join(TBL, 'MetaAnalysis_CrossDataset.csv'), index=False)
meta_sig.to_csv(os.path.join(TBL, 'MetaAnalysis_Significant.csv'), index=False)

# ============================================================
# 3. FIGURES
# ============================================================
print("\nGenerating enhanced validation figures...")

# Fig A: Cross-dataset log2FC correlation matrix
fc_merged = meta[fc_cols].dropna(thresh=2)
if len(fc_merged) >= 10:
    corr_mat = np.zeros((len(fc_cols), len(fc_cols)))
    pval_mat = np.zeros((len(fc_cols), len(fc_cols)))
    for i, ci in enumerate(fc_cols):
        for j, cj in enumerate(fc_cols):
            mask = fc_merged[ci].notna() & fc_merged[cj].notna()
            if mask.sum() >= 5:
                r, p = spearmanr(fc_merged.loc[mask, ci], fc_merged.loc[mask, cj])
                corr_mat[i, j] = r
                pval_mat[i, j] = p
            else:
                corr_mat[i, j] = np.nan

    fig, ax = plt.subplots(figsize=(6, 5))
    labels = [c.replace('log2FC_','').replace('_',' ') for c in fc_cols]
    mask = np.isnan(corr_mat)
    im = ax.imshow(corr_mat, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
    for i in range(len(fc_cols)):
        for j in range(len(fc_cols)):
            if not np.isnan(corr_mat[i, j]):
                ax.text(j, i, f'{corr_mat[i,j]:.2f}', ha='center', va='center',
                       fontsize=11, fontweight='bold',
                       color='white' if abs(corr_mat[i,j]) > 0.6 else 'black')
            else:
                ax.text(j, i, 'N/A', ha='center', va='center', fontsize=9, color='grey')
    ax.set_xticks(range(len(fc_cols)))
    ax.set_xticklabels(labels, fontsize=10, rotation=30, ha='right')
    ax.set_yticks(range(len(fc_cols)))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_title('Cross-Dataset log2FC Spearman Correlation\n(Meta-Analysis)', fontsize=13, fontweight='bold')
    plt.colorbar(im, ax=ax, shrink=0.8, label="Spearman ρ")
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'Validation_Meta_CorrelationMatrix.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  Meta correlation matrix saved')

# Fig B: Volcano plot for meta-analysis results
if len(meta_sig) >= 5:
    fig, ax = plt.subplots(figsize=(8, 6))
    meta_valid = meta[meta['fisher_pval'].notna()].copy()
    meta_valid['log_meta_p'] = -np.log10(meta_valid['fisher_padj'].clip(lower=1e-300))

    # Mean log2FC across datasets
    meta_valid['mean_FC'] = meta_valid[fc_cols].mean(axis=1)

    colors = np.where(meta_valid['fisher_padj'] < 0.05,
                      np.where(meta_valid['mean_FC'] > 0.5, '#E64B35',
                               np.where(meta_valid['mean_FC'] < -0.5, '#4DBBD5', '#B09C85')),
                      'lightgrey')
    sizes = np.where(meta_valid['fisher_padj'] < 0.05, 15, 5)
    alpha_vals = np.where(meta_valid['fisher_padj'] < 0.05, 0.7, 0.3)

    ax.scatter(meta_valid['mean_FC'], meta_valid['log_meta_p'],
              c=colors, s=sizes, alpha=alpha_vals, edgecolors='white', linewidth=0.3)

    # Label top genes
    for _, row in meta_valid.nsmallest(15, 'fisher_padj').iterrows():
        if abs(row['mean_FC']) > 0.3:
            ax.annotate(row['gene'], (row['mean_FC'], row['log_meta_p']),
                       fontsize=8, fontweight='bold', xytext=(4, 4), textcoords='offset points')

    ax.axhline(-np.log10(0.05), color='grey', ls='--', alpha=0.5, label='padj=0.05')
    ax.set_xlabel('Meta log2FC (mean across datasets)', fontsize=12)
    ax.set_ylabel('-log10(Fisher combined padj)', fontsize=12)
    ax.set_title('Meta-Analysis — Cross-Dataset Consensus DEGs\n(GSE137665 + GSE211088 + GSE237419)',
                fontsize=13, fontweight='bold')
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([0],[0], color='#E64B35', lw=0, marker='o', ms=8, label='Upregulated consensus'),
        Line2D([0],[0], color='#4DBBD5', lw=0, marker='o', ms=8, label='Downregulated consensus'),
    ], frameon=False, fontsize=9)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'Validation_Meta_Volcano.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  Meta volcano plot saved')

# Fig C: Concordance analysis — direction preservation
if len(meta_sig) >= 5:
    our_sig = meta_sig[meta_sig['padj_our'].notna() & (meta_sig['padj_our'] < 0.05)]
    if len(our_sig) >= 5:
        fig, axes = plt.subplots(1, min(2, len(fc_cols)-1), figsize=(6*min(2,len(fc_cols)-1), 5))
        if len(fc_cols)-1 == 1:
            axes = [axes]

        ext_fc_cols = [c for c in fc_cols if c != 'log2FC_our']
        for ax_idx, ext_col in enumerate(ext_fc_cols[:2]):
            ax = axes[ax_idx]
            valid = meta_sig[meta_sig[['log2FC_our', ext_col]].notna().all(axis=1)]
            if len(valid) >= 5:
                r_v, p_v = spearmanr(valid['log2FC_our'], valid[ext_col])
                ax.scatter(valid[ext_col], valid['log2FC_our'], c='#3C5488', s=40, alpha=0.6,
                          edgecolors='white', linewidth=0.5)
                lim = max(abs(valid[ext_col]).max(), abs(valid['log2FC_our']).max()) * 1.2
                ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3)
                ax.axhline(0, color='grey', lw=0.5)
                ax.axvline(0, color='grey', lw=0.5)
                for _, row in valid.iterrows():
                    if row['gene'] in ['Pomc','Pcsk2','Fos','Egr1','Dbp','Nr3c1','Rbm3','Tsc22d3']:
                        ax.annotate(row['gene'], (row[ext_col], row['log2FC_our']),
                                   fontsize=8, fontweight='bold', xytext=(4, 4), textcoords='offset points')
                dname = ext_col.replace('log2FC_','')
                ax.set_xlabel(f'log2FC ({dname})', fontsize=11)
                ax.set_ylabel('log2FC (Our Study)', fontsize=11)
                ax.set_title(f'{dname} vs Our Study\nSpearman r={r_v:.3f}, p={p_v:.2e}',
                           fontsize=12, fontweight='bold')
                sns.despine()

        plt.tight_layout()
        plt.savefig(os.path.join(FIG, 'Validation_Meta_Concordance.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print('  Concordance scatter plots saved')

# Fig D: Top meta-significant genes heatmap
if len(meta_sig) >= 8:
    top_meta = meta_sig.nsmallest(30, 'fisher_padj')
    heat_cols = [c for c in fc_cols if c in top_meta.columns]
    heat_data = top_meta.set_index('gene')[heat_cols].copy()
    heat_data = heat_data.fillna(0)

    # Clustermap
    g = sns.clustermap(heat_data, cmap='RdBu_r', center=0, vmin=-2, vmax=2,
                       figsize=(8, max(6, len(heat_data)*0.3)),
                       row_cluster=True, col_cluster=False,
                       xticklabels=[c.replace('log2FC_','').replace('_',' ') for c in heat_cols],
                       linewidths=0.5, cbar_kws={'label': 'log2FC'})
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), fontsize=8)
    g.fig.suptitle('Meta-Significant Genes Across Datasets\n(Fisher padj < 0.05, Top 30)',
                   fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'Validation_Meta_Heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  Meta heatmap saved')

# ============================================================
# 4. Summary
# ============================================================
print("\n" + "=" * 60)
print("ENHANCED VALIDATION — RESULTS SUMMARY")
print("=" * 60)
print(f"Total genes in meta-analysis: {len(meta)}")
print(f"Meta-significant genes (Fisher padj<0.05): {len(meta_sig)}")
if len(meta_sig) >= 5:
    print(f"  Up-regulated consensus: {(meta_sig['consensus_dir']==1).sum()}")
    print(f"  Down-regulated consensus: {(meta_sig['consensus_dir']==-1).sum()}")

# Key genes:
key_verify = ['Pomc','Pcsk2','Fos','Egr1','Dbp','Nr3c1','Rbm3','Tsc22d3','Scg2','Nnat','Mt1']
for g in key_verify:
    if g in meta_sig['gene'].values:
        r = meta_sig[meta_sig['gene']==g].iloc[0]
        fc_str = '  '.join([f"{c.replace('log2FC_','')}={r[c]:+.2f}" for c in fc_cols if pd.notna(r.get(c))])
        print(f"  {g:12s}  meta_padj={r['fisher_padj']:.2e}  {fc_str}")

print(f"\nSaved: {TBL}/MetaAnalysis_CrossDataset.csv")
print(f"Saved: {TBL}/MetaAnalysis_Significant.csv")
print(f"Figures saved to: {FIG}/")
