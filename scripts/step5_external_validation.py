"""
Step 5: External Dataset Validation — streamlined version
Validates our key genes in GSE211088 (PFC RNA-seq) and GSE237419 (Cortex RNA-seq)
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings, os, gzip, re
warnings.filterwarnings('ignore')
from scipy.stats import ttest_ind, spearmanr

FIG = r'C:\Users\1\sleep-deprivation-project\results\figures'
TBL = r'C:\Users\1\sleep-deprivation-project\results\tables'
EXT = r'C:\Users\1\sleep-deprivation-project\data\external'
os.makedirs(FIG, exist_ok=True)
os.makedirs(TBL, exist_ok=True)

print("Loading local Ensembl -> Symbol mapping...")
import scanpy as sc
raw = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_raw.h5ad')

# Use locally built Ensembl-to-Symbol mapping (from NCBI gene_info)
symbol_map = pd.read_csv(os.path.join(EXT, 'ensembl_to_symbol.csv'), index_col=0)
symbol_dict = symbol_map['Symbol'].to_dict()
print(f"Loaded mapping: {len(symbol_dict)} Ensembl IDs -> Symbols")

KEY_GENES_OF_INTEREST = ['Pomc','Fos','Fosb','Jun','Junb','Jund','Egr1','Egr2','Egr3',
    'Nr4a1','Nr4a2','Nr4a3','Nr3c1','Crem','Clock','Per1','Per2','Per3',
    'Cry1','Cry2','Arntl','Bmal1','Nr1d1','Nr1d2','Rora','Rorb','Rorc',
    'Dbp','Tef','Nfil3','Stat3','Nfkb1','Rela','Sp1','Yy1','Ctcf',
    'Pparg','Ppargc1a','Foxo1','Foxo3','Rest','Mecp2','Hdac1','Hdac2',
    'Bdnf','Arc','Homer1','Npas4','Srebf1','Srebf2',
    'Rbm3','Tsc22d3','Rnaset2a','Mt1','Mt3','Nnat',
    'Pcsk2','Malat1','Meg3','Gnas','Cidea','Scg2',
    'Gfap','Aif1','Trem2','Slc17a7','Gad1','Gad2','Snap25']

# ============= 1. Parse GSE211088 =============
print("\n=== GSE211088: RNA-seq, PFC ===")
df211 = pd.read_csv(os.path.join(EXT, 'GSE211088_RNAseq_gene.txt.gz'), sep='\t', index_col=0)
print(f"Shape: {df211.shape}")

hc_cols = [c for c in df211.columns if c.startswith('WTHC')]
sd_cols = [c for c in df211.columns if c.startswith('WTSD')]
print(f"HC: {len(hc_cols)}, SD: {len(sd_cols)}")

# Map Ensembl IDs to symbols using local mapping
ens_ids = [x.split('.')[0] for x in df211.index]
mapped_symbols = [symbol_dict.get(eid, None) for eid in ens_ids]
df211['symbol'] = mapped_symbols
print(f"  Mapped via local file: {(pd.Series(mapped_symbols).notna()).sum()} genes")
df211_sym = df211[df211['symbol'].notna()].copy()
# Take mean for duplicate symbols
expr_cols_211 = [c for c in df211.columns if c not in ['ensembl_base', 'symbol']]
df211_sym = df211_sym.groupby('symbol')[expr_cols_211].mean()
print(f"After symbol mapping: {df211_sym.shape}")

# Compute log2FC
df211_sym['mean_HC'] = df211_sym[hc_cols].mean(axis=1)
df211_sym['mean_SD'] = df211_sym[sd_cols].mean(axis=1)
df211_sym['log2FC'] = np.log2((df211_sym['mean_SD'] + 0.1) / (df211_sym['mean_HC'] + 0.1))
df211_sym['pval'] = 1.0
for g in df211_sym.index:
    vals = df211_sym.loc[g, hc_cols + sd_cols].values
    if vals.std() > 1e-8 and df211_sym.loc[g, hc_cols].std() > 1e-8:
        df211_sym.loc[g, 'pval'] = ttest_ind(
            df211_sym.loc[g, sd_cols].values, df211_sym.loc[g, hc_cols].values).pvalue
df211_sym['padj'] = (df211_sym['pval'].clip(lower=1e-300) * len(df211_sym)).clip(upper=1.0)

print(f"Genes with padj<0.05: {(df211_sym['padj']<0.05).sum()}")

# Extract key genes
print("Key gene validation (GSE211088):")
val_rows = []
for gene in KEY_GENES_OF_INTEREST:
    if gene in df211_sym.index:
        r = df211_sym.loc[gene]
        val_rows.append({'gene': gene, 'log2FC_GSE211088': r['log2FC'], 'padj_GSE211088': r['padj']})
        print(f"  {gene:12s}  log2FC={r['log2FC']:+.3f}  padj={r['padj']:.2e}")
val_211 = pd.DataFrame(val_rows)

# ============= 2. Parse GSE237419 =============
print("\n=== GSE237419: Bulk RNA-seq, Cortex ===")
# This file has a two-row header: row0=GSM IDs, row1=condition labels
with gzip.open(os.path.join(EXT, 'GSE237419_gse_HCSDRS_WT_SleepIntegration_salmon.txt.gz'), 'rt') as f:
    header1 = f.readline().strip().split('\t')
    header2 = f.readline().strip().split('\t')

print(f"Header: {len(header1)} cols")
print(f"Sample labels (first 10): {header2[:10]}")
print(f"Sample labels (last 10): {header2[-10:]}")

# Parse condition info from header2
labels_map = {}
for i, label in enumerate(header2):
    label = label.strip()
    if label and label != 'ENSEMBL_ID':
        # Parse prefix: HC=HomeCage, SD=SleepDep, RS=Recovery
        match = re.match(r'^(HC|SD|RS)(\d+)_', label)
        if match:
            labels_map[i] = {'condition': match.group(1), 'hours': int(match.group(2)),
                           'label': label, 'col_idx': i}

print(f"Parsed {len(labels_map)} sample columns")
cond_counts = {}
for v in labels_map.values():
    key = f"{v['condition']}_{v['hours']}h"
    cond_counts[key] = cond_counts.get(key, 0) + 1
print(f"Condition counts: {cond_counts}")

# Re-read with proper indexing
df237 = pd.read_csv(os.path.join(EXT, 'GSE237419_gse_HCSDRS_WT_SleepIntegration_salmon.txt.gz'),
                     sep='\t', skiprows=2, header=None)
df237.index = df237[0]
df237 = df237.drop(columns=[0])

# Map Ensembl IDs using local mapping
ens_ids_237 = [x.split('.')[0] for x in df237.index]
mapped_237 = [symbol_dict.get(eid, None) for eid in ens_ids_237]
df237['symbol'] = mapped_237

# Filter and group
df237_sym = df237[df237['symbol'].notna()].copy()
# Get integer columns for expression data
expr_cols_237 = [c for c in df237_sym.columns if isinstance(c, (int, np.integer))]
df237_sym = df237_sym.groupby('symbol')[expr_cols_237].mean()
print(f"After symbol mapping: {df237_sym.shape}")

# Extract HC5 vs SD5 comparison
hc5_indices = [v['col_idx'] for v in labels_map.values() if v['condition'] == 'HC' and v['hours'] == 5]
sd5_indices = [v['col_idx'] for v in labels_map.values() if v['condition'] == 'SD' and v['hours'] == 5]
print(f"HC 5h: {len(hc5_indices)} samples, SD 5h: {len(sd5_indices)} samples")

if len(hc5_indices) >= 2 and len(sd5_indices) >= 2:
    df237_sym['mean_HC5'] = df237_sym[hc5_indices].mean(axis=1)
    df237_sym['mean_SD5'] = df237_sym[sd5_indices].mean(axis=1)
    df237_sym['log2FC'] = np.log2((df237_sym['mean_SD5'] + 0.1) / (df237_sym['mean_HC5'] + 0.1))

    print("Key gene validation (GSE237419):")
    val_rows_237 = []
    for gene in KEY_GENES_OF_INTEREST:
        if gene in df237_sym.index:
            r = df237_sym.loc[gene]
            try:
                pv = ttest_ind(df237_sym.loc[gene, sd5_indices].values,
                              df237_sym.loc[gene, hc5_indices].values).pvalue
            except:
                pv = np.nan
            val_rows_237.append({'gene': gene, 'log2FC_GSE237419': r['log2FC'], 'pval_GSE237419': pv})
            print(f"  {gene:12s}  log2FC={r['log2FC']:+.3f}  p={pv:.2e}")
    val_237 = pd.DataFrame(val_rows_237)
else:
    print("Not enough samples for HC5 vs SD5 comparison")
    val_237 = pd.DataFrame()

# ============= 3. Merge and compare =============
print("\n=== Cross-Dataset Consistency ===")

# Our results
our_degs = pd.read_csv(os.path.join(TBL, 'DEGs_significant.csv'))
our_degs = our_degs[['names','logfoldchanges','pvals_adj','brain_region']].copy()
our_degs.columns = ['gene','log2FC_ourStudy','padj_ourStudy','brain_region']

pomc_corr = pd.read_csv(os.path.join(TBL, 'Pomc_coexpression.csv'))
pomc_corr = pomc_corr[['gene','spearman_r']].copy()
pomc_corr.columns = ['gene','r_with_Pomc']

val_merged = val_211.copy()
if len(val_237) > 0:
    val_merged = val_merged.merge(val_237, on='gene', how='outer')
val_merged = val_merged.merge(our_degs, on='gene', how='left')
val_merged = val_merged.merge(pomc_corr, on='gene', how='left')

# Consistency check
fc_cols = [c for c in ['log2FC_GSE211088','log2FC_GSE237419','log2FC_ourStudy'] if c in val_merged.columns]
valid_all = val_merged[val_merged[fc_cols].notna().all(axis=1)] if len(fc_cols) >= 2 else pd.DataFrame()

if len(valid_all) >= 3:
    for i in range(len(fc_cols)):
        for j in range(i+1, len(fc_cols)):
            mask = valid_all[fc_cols[i]].notna() & valid_all[fc_cols[j]].notna()
            if mask.sum() >= 3:
                r, p = spearmanr(valid_all.loc[mask, fc_cols[i]], valid_all.loc[mask, fc_cols[j]])
                print(f"Spearman r ({fc_cols[i]} vs {fc_cols[j]}): {r:+.3f}, p={p:.3f}")

# Save
val_merged.to_csv(os.path.join(TBL, 'External_Validation_Results.csv'), index=False)
print(f"Saved: {TBL}/External_Validation_Results.csv")
print(f"Genes validated: {len(val_merged)}")

# ============= 4. Figures =============
print("\nGenerating validation figures...")

# Fig 1: Scatter plot our log2FC vs GSE211088
if 'log2FC_GSE211088' in val_merged.columns:
    valid = val_merged[val_merged[['log2FC_ourStudy','log2FC_GSE211088']].notna().all(axis=1)]
    if len(valid) >= 5:
        fig, ax = plt.subplots(figsize=(7, 7))
        r_val, p_val = spearmanr(valid['log2FC_ourStudy'], valid['log2FC_GSE211088'])

        ax.scatter(valid['log2FC_GSE211088'], valid['log2FC_ourStudy'],
                  c='#4DBBD5', s=50, alpha=0.7, edgecolors='white', linewidth=0.5, zorder=3)
        # Label key genes
        for _, row in valid.iterrows():
            if row['gene'] in ['Pomc','Fos','Egr1','Nr4a1','Dbp','Clock','Per1','Bdnf','Rbm3']:
                ax.annotate(row['gene'], (row['log2FC_GSE211088'], row['log2FC_ourStudy']),
                          fontsize=8, fontweight='bold', xytext=(5, 5), textcoords='offset points')
        # Reference lines
        all_vals = pd.concat([valid['log2FC_GSE211088'], valid['log2FC_ourStudy']])
        lim = max(abs(all_vals.min()), abs(all_vals.max())) * 1.2
        ax.plot([-lim, lim], [-lim, lim], 'k--', alpha=0.3, lw=0.8)
        ax.axhline(0, color='grey', lw=0.5)
        ax.axvline(0, color='grey', lw=0.5)
        ax.set_xlabel('log2FC GSE211088 (PFC RNA-seq)', fontsize=11)
        ax.set_ylabel('log2FC Our Study (scRNA-seq)', fontsize=11)
        ax.set_title(f'Cross-Dataset Validation\nSpearman r={r_val:.3f}, p={p_val:.2e}, n={len(valid)}',
                    fontsize=12, fontweight='bold')
        sns.despine()
        plt.tight_layout()
        plt.savefig(os.path.join(FIG, 'Validation_CrossDataset.png'), dpi=300, bbox_inches='tight')
        plt.close()
        print('  Cross-dataset scatter saved')

# Fig 2: Barplot of key genes across datasets
fc_cols_bar = [c for c in fc_cols if c in val_merged.columns]
if len(fc_cols_bar) >= 1:
    for fc_col in fc_cols_bar:
        df_plot = val_merged[val_merged[fc_col].notna()].sort_values(fc_col, key=abs).tail(25)
        if len(df_plot) > 0:
            fig, ax = plt.subplots(figsize=(8, max(5, len(df_plot)*0.3)))
            colors = ['#E64B35' if x < 0 else '#4DBBD5' for x in df_plot[fc_col]]
            ax.barh(range(len(df_plot)), df_plot[fc_col].values, color=colors, alpha=0.8)
            ax.set_yticks(range(len(df_plot)))
            ax.set_yticklabels(df_plot['gene'].values, fontsize=9)
            ax.axvline(0, color='black', lw=0.5)
            dataset_name = fc_col.replace('log2FC_','')
            ax.set_title(f'Gene Expression Changes in {dataset_name}\n(SD vs Control)', fontsize=12, fontweight='bold')
            sns.despine()
            plt.tight_layout()
            plt.savefig(os.path.join(FIG, f'Validation_{dataset_name}_Barplot.png'), dpi=300, bbox_inches='tight')
            plt.close()
    print('  Barplots saved')

# Fig 3: Pomc expression across datasets
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# GSE211088 Pomc
if 'Pomc' in df211_sym.index:
    ax = axes[0]
    hc_vals = df211_sym.loc['Pomc', hc_cols].astype(float).values
    sd_vals = df211_sym.loc['Pomc', sd_cols].astype(float).values
    x_pos = [0, 1]
    means = [hc_vals.mean(), sd_vals.mean()]
    sems = [hc_vals.std()/np.sqrt(len(hc_vals)), sd_vals.std()/np.sqrt(len(sd_vals))]
    ax.bar(x_pos, means, color=['#4DBBD5','#E64B35'], alpha=0.7, yerr=sems, capsize=5)
    ax.scatter([0]*len(hc_vals), hc_vals, c='#4DBBD5', s=30, alpha=0.6, edgecolors='white')
    ax.scatter([1]*len(sd_vals), sd_vals, c='#E64B35', s=30, alpha=0.6, edgecolors='white')
    pv = ttest_ind(sd_vals, hc_vals).pvalue
    ax.set_xticks(x_pos)
    ax.set_xticklabels(['HomeCage', 'SD 5h'])
    ax.set_title(f'Pomc in GSE211088 (PFC RNA-seq)\np={pv:.2e}', fontsize=11, fontweight='bold')
    ax.set_ylabel('Expression (TPM)', fontsize=10)

# Our data Pomc
ax = axes[1]
hypo = raw[raw.obs['brain_region'] == 'Hypothalamus'].copy()
sc.pp.normalize_total(hypo, target_sum=1e4)
sc.pp.log1p(hypo)
pomc_expr = hypo[:, 'Pomc'].X.toarray().flatten() if hasattr(hypo[:, 'Pomc'].X, 'toarray') else hypo[:, 'Pomc'].X
conds = ['A1', 'A2', 'A3']
c_means = [pomc_expr[hypo.obs['condition_code']==c].mean() for c in conds]
c_sems = [pomc_expr[hypo.obs['condition_code']==c].std()/np.sqrt((hypo.obs['condition_code']==c).sum()) for c in conds]
c_colors = ['#4DBBD5', '#E64B35', '#00A087']
ax.bar(range(3), c_means, color=c_colors, alpha=0.7, yerr=c_sems, capsize=5)
ax.set_xticks(range(3))
ax.set_xticklabels(['Normal (A1)', 'SD (A2)', 'Recovery (A3)'])
ax.set_title('Pomc in Our Study (Hypothalamus scRNA-seq)', fontsize=11, fontweight='bold')
ax.set_ylabel('Expression (log-norm)', fontsize=10)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'Validation_Pomc_External.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Pomc validation figure saved')

print(f"\n=== External Validation Complete ===")
