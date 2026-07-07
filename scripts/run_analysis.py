# Sleep Deprivation: Differential Expression + Visualization
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=True)
sc.logging.print_header()

FIG_DIR = '/c/Users/1/sleep-deprivation-project/results/figures'
TBL_DIR = '/c/Users/1/sleep-deprivation-project/results/tables'
import os
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TBL_DIR, exist_ok=True)

# ============= Load =============
print("Loading processed data...")
adata = sc.read_h5ad('/c/Users/1/sleep-deprivation-project/data/GSE137665_processed.h5ad')
print(adata)

# Set raw for expression lookups
adata = adata.raw.to_adata() if adata.raw is not None else adata

# ============= UMAP Figures =============
print("\nGenerating UMAP plots...")

# Fig 1a: UMAP by brain region
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='brain_region', ax=ax, show=False,
           palette={'Cortex':'#E64B35','Hypothalamus':'#4DBBD5','Brainstem':'#00A087'},
           title='Brain Regions')
plt.savefig(f'{FIG_DIR}/UMAP_brain_regions.png', dpi=300, bbox_inches='tight')
plt.close()

# Fig 1b: UMAP by cell type
fig, ax = plt.subplots(figsize=(10, 8))
sc.pl.umap(adata, color='cell_type', ax=ax, show=False, title='Cell Types')
plt.savefig(f'{FIG_DIR}/UMAP_cell_types.png', dpi=300, bbox_inches='tight')
plt.close()

# Fig 1c: UMAP by condition_code (A1/A2/A3)
fig, axes = plt.subplots(1, 3, figsize=(24, 7))
for i, code in enumerate(['A1','A2','A3']):
    mask = adata.obs['condition_code'] == code
    sc.pl.umap(adata, ax=axes[i], show=False, title=f'Condition {code}')
plt.savefig(f'{FIG_DIR}/UMAP_conditions.png', dpi=300, bbox_inches='tight')
plt.close()

# ============= Cell Type Composition =============
print("Cell type composition...")
comp = pd.crosstab(adata.obs['cell_type'], adata.obs['brain_region'])
comp_pct = comp.div(comp.sum(axis=0), axis=1) * 100

fig, axes = plt.subplots(1, 2, figsize=(16, 6))
comp_pct.T.plot(kind='bar', stacked=True, ax=axes[0], colormap='Set3')
axes[0].set_title('Cell Type Composition by Brain Region (%)')
axes[0].set_ylabel('% of cells')
axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)

comp.T.plot(kind='bar', stacked=True, ax=axes[1], colormap='Set3')
axes[1].set_title('Brain Region Composition by Cell Type')
axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=7)
plt.tight_layout()
plt.savefig(f'{FIG_DIR}/CellType_Composition.png', dpi=300, bbox_inches='tight')
plt.close()

# ============= Marker Gene Dotplot =============
print("Marker gene dotplot...")
marker_genes = ['Slc17a7','Gad1','Gfap','Cx3cr1','Mog','Pdgfra','Cldn5','Snap25','Sst','Aif1','Aqp4','Plp1']
present = [g for g in marker_genes if g in adata.var_names]
fig, ax = plt.subplots(figsize=(12, 5))
sc.pl.dotplot(adata, present, groupby='cell_type', ax=ax, show=False)
plt.savefig(f'{FIG_DIR}/Marker_Dotplot.png', dpi=300, bbox_inches='tight')
plt.close()

# ============= Differential Expression: SD (A2) vs Normal (A1) =============
print("\nDifferential expression: condition A2 vs A1...")
# We don't know exact condition mapping yet, compare A2 vs A1 for all regions
de_results = []

for region in ['Cortex','Hypothalamus','Brainstem']:
    sub = adata[adata.obs['brain_region'] == region].copy()
    a1_cells = sub.obs['condition_code'] == 'A1'
    a2_cells = sub.obs['condition_code'] == 'A2'
    if a1_cells.sum() < 10 or a2_cells.sum() < 10:
        continue
    sc.tl.rank_genes_groups(sub, groupby='condition_code', groups=['A2'], reference='A1',
                            method='wilcoxon')
    de = sc.get.rank_genes_groups_df(sub, group='A2')
    de['brain_region'] = region
    de_results.append(de)
    print(f'  {region}: {len(de[de["pvals_adj"] < 0.05])} DEGs (padj<0.05)')

if de_results:
    all_de = pd.concat(de_results, ignore_index=True)
    all_de.to_csv(f'{TBL_DIR}/DE_A2_vs_A1.csv', index=False)
    print(f'Total DE results saved.')

# ============= Volcano Plots =============
print("Volcano plots...")
for region in ['Cortex','Hypothalamus','Brainstem']:
    ct_de = all_de[all_de['brain_region']==region].copy()
    if len(ct_de) == 0:
        continue
    ct_de['-log10padj'] = -np.log10(ct_de['pvals_adj'].clip(lower=1e-300))
    ct_de['sig'] = 'NS'
    ct_de.loc[(ct_de['logfoldchanges']>0.25)&(ct_de['pvals_adj']<0.05),'sig']='Up'
    ct_de.loc[(ct_de['logfoldchanges']<-0.25)&(ct_de['pvals_adj']<0.05),'sig']='Down'

    n_up = (ct_de['sig']=='Up').sum()
    n_down = (ct_de['sig']=='Down').sum()

    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot NS as small faint background dots
    ns = ct_de[ct_de['sig']=='NS']
    ax.scatter(ns['logfoldchanges'], ns['-log10padj'], c='#B0B0B0', s=4, alpha=0.25,
               label=f'NS ({len(ns)})', rasterized=True)

    # Plot significant points larger and more visible on top
    colors_sig = {'Down': '#1A73E8', 'Up': '#E41E26'}
    for s in ['Down', 'Up']:
        sub = ct_de[ct_de['sig']==s]
        if len(sub) == 0:
            continue
        ax.scatter(sub['logfoldchanges'], sub['-log10padj'],
                   c=colors_sig[s], s=28, alpha=0.75,
                   edgecolors='white', linewidth=0.3,
                   label=f'{s} ({len(sub)})')

    # Label top significant genes (up to 10 per direction)
    top_up = ct_de[ct_de['sig']=='Up'].nsmallest(8, 'pvals_adj')
    top_down = ct_de[ct_de['sig']=='Down'].nsmallest(8, 'pvals_adj')
    top_genes = pd.concat([top_up, top_down])
    for _, row in top_genes.iterrows():
        ax.annotate(row['names'],
                    (row['logfoldchanges'], row['-log10padj']),
                    fontsize=6.5, fontweight='bold',
                    alpha=0.9,
                    xytext=(4, 4), textcoords='offset points',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.7, edgecolor='none'))

    # Threshold lines
    ax.axvline(-0.25, ls='--', c='grey', alpha=0.4, lw=1)
    ax.axvline(0.25, ls='--', c='grey', alpha=0.4, lw=1)
    ax.axhline(-np.log10(0.05), ls='--', c='grey', alpha=0.4, lw=1)

    ax.set_xlabel('log2 Fold Change (A2 vs A1)', fontsize=13, fontweight='bold')
    ax.set_ylabel('−log₁₀(adjusted p-value)', fontsize=13, fontweight='bold')
    ax.set_title(f'{region}: Differential Expression (A2 vs A1)', fontsize=15, fontweight='bold', pad=12)

    ax.legend(loc='upper left', frameon=True, fancybox=True, framealpha=0.8, fontsize=10,
              markerscale=0.8)

    # Add count text
    ax.text(0.98, 0.95, f'↑ {n_up}  ↓ {n_down}',
            transform=ax.transAxes, fontsize=11, fontweight='bold',
            ha='right', va='top',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.7, edgecolor='#CCCCCC'))

    sns.despine()
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/Volcano_{region}.png', dpi=300, bbox_inches='tight')
    plt.close()

# ============= Heatmap of Top DEGs =============
print("Top DEGs heatmap...")
top_genes = all_de.nsmallest(30, 'pvals_adj')['names'].unique()
top_genes_present = [g for g in top_genes if g in adata.var_names][:20]

if top_genes_present:
    fig, ax = plt.subplots(figsize=(12, 8))
    sc.pl.heatmap(adata[adata.obs['cell_type']!='Unassigned'], top_genes_present,
                  groupby='cell_type', ax=ax, show=False)
    plt.savefig(f'{FIG_DIR}/Top_DEGs_Heatmap.png', dpi=300, bbox_inches='tight')
    plt.close()

# ============= Response Score =============
print("Sleep deprivation response score by cell type...")
response = all_de[all_de['pvals_adj']<0.05].groupby('brain_region').agg(
    n_up=('logfoldchanges', lambda x: (x>0.25).sum()),
    n_down=('logfoldchanges', lambda x: (x<-0.25).sum())
).reset_index()
response['total'] = response['n_up'] + response['n_down']

fig, ax = plt.subplots(figsize=(8, 4))
x = range(len(response))
ax.barh(x, response['n_up'], color='#E64B35', label='Up in A2', height=0.5)
ax.barh(x, -response['n_down'], color='#4DBBD5', label='Down in A2', height=0.5)
ax.set_yticks(x)
ax.set_yticklabels(response['brain_region'])
ax.set_xlabel('Number of DEGs (padj < 0.05)')
ax.axvline(0, color='black', linewidth=0.5)
ax.legend(frameon=False)
ax.set_title('Transcriptional Response by Brain Region')
sns.despine()
plt.savefig(f'{FIG_DIR}/Response_by_Region.png', dpi=300, bbox_inches='tight')
plt.close()

print(f'\n=== All figures saved to {FIG_DIR} ===')
print(f'=== All tables saved to {TBL_DIR} ===')
print('Done!')
