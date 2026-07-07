"""
Step 4: ML Biomarker Discovery — LASSO + Random Forest + SVM-RFE
Pseudo-bulk at Leiden cluster level for sufficient observations
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

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.model_selection import StratifiedKFold, cross_val_score, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc

FIG = r'C:\Users\1\sleep-deprivation-project\results\figures'
TBL = r'C:\Users\1\sleep-deprivation-project\results\tables'
os.makedirs(FIG, exist_ok=True)
os.makedirs(TBL, exist_ok=True)

# ============= 1. Build pseudo-bulk at Leiden cluster level =============
print("Building pseudo-bulk (Leiden cluster level)...")
raw = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_raw.h5ad')
proc = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_processed.h5ad')
raw = raw[raw.obs_names.isin(proc.obs_names)]

# Attach leiden clusters to raw data
raw.obs['leiden'] = proc.obs['leiden'].values
raw.obs['cell_type'] = proc.obs['cell_type'].values

# Filter to A1/A2 only (remove recovery A3 for clean ML)
raw_ml = raw[raw.obs['condition_code'].isin(['A1', 'A2'])].copy()
sc.pp.normalize_total(raw_ml, target_sum=1e4)
sc.pp.log1p(raw_ml)

# Load DEG results for feature pre-filtering
deg_all = pd.read_csv(os.path.join(TBL, 'DE_A2_vs_A1.csv'))
deg_sig = deg_all[deg_all['pvals_adj'] < 0.05].copy()
print(f"Significant DEGs: {len(deg_sig)}")

# Load Pomc network TFs
pomc_tfs = pd.read_csv(os.path.join(TBL, 'Pomc_network_TFs.csv'))
tf_set = set(pomc_tfs['gene'].tolist()) if len(pomc_tfs) > 0 else set()

# Feature genes: top DEGs + Pomc TFs + known sleep genes
SLEEP_GENES = ['Clock','Bmal1','Arntl','Per1','Per2','Per3','Cry1','Cry2',
               'Nr1d1','Nr1d2','Rora','Rorb','Arc','Bdnf','Fos','Egr1','Npas4',
               'Homer1','Nr3c1','Crem','Srebf1','Ppargc1a','Foxo1','Pomc']
deg_genes = deg_sig.nsmallest(200, 'pvals_adj')['names'].tolist()
feature_genes = list(dict.fromkeys(deg_genes + list(tf_set) + [g for g in SLEEP_GENES if g in raw_ml.var_names]))
feature_genes = [g for g in feature_genes if g in raw_ml.var_names]
print(f"Feature genes: {len(feature_genes)}")

# Create pseudo-bulk: for each leiden cluster × condition combination
# Only include combinations with >= 20 cells
pseudo_rows = []
for leiden in sorted(raw_ml.obs['leiden'].unique()):
    for cond in ['A1', 'A2']:
        mask = (raw_ml.obs['leiden'] == leiden) & (raw_ml.obs['condition_code'] == cond)
        n_cells = mask.sum()
        if n_cells >= 20:
            expr = raw_ml[mask, feature_genes].X
            if hasattr(expr, 'toarray'):
                expr = expr.toarray()
            pseudo_rows.append({
                'leiden': leiden,
                'condition': cond,
                'label': 0 if cond == 'A1' else 1,
                'n_cells': n_cells,
                'cell_type': raw_ml.obs.loc[mask, 'cell_type'].mode()[0],
                'expression': expr.mean(axis=0)
            })

print(f"Pseudo-bulk samples: {len(pseudo_rows)}")
X = np.vstack([r['expression'] for r in pseudo_rows])
y = np.array([r['label'] for r in pseudo_rows])
print(f"  A1 (Normal): {(y==0).sum()}, A2 (SD): {(y==1).sum()}")

if len(pseudo_rows) < 10:
    print("WARNING: Too few pseudo-bulk samples. ML may be unreliable.")

# Scale
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
n_cv = min(5, min((y==0).sum(), (y==1).sum()))  # adjust CV folds to minority class size
if n_cv < 2:
    n_cv = 2
print(f"Using {n_cv}-fold CV")

# ============= 2. LASSO =============
print("\n=== LASSO Feature Selection ===")
lasso = LogisticRegression(penalty='l1', solver='saga', max_iter=10000, random_state=42)
cv_lasso = GridSearchCV(lasso, {'C': np.logspace(-3, 1, 15)}, cv=StratifiedKFold(n_cv, shuffle=True, random_state=42),
                        scoring='roc_auc')
cv_lasso.fit(X_scaled, y)
print(f"Best C: {cv_lasso.best_params_['C']:.4f}, Best AUC: {cv_lasso.best_score_:.3f}")

lasso_best = cv_lasso.best_estimator_
lasso_coefs = pd.DataFrame({'gene': feature_genes, 'coef': lasso_best.coef_[0]})
lasso_coefs = lasso_coefs.query('coef != 0').sort_values('coef', key=abs, ascending=False)
lasso_set = set(lasso_coefs['gene'].tolist())
print(f"LASSO selected: {len(lasso_set)} genes")

# ============= 3. Random Forest =============
print("\n=== Random Forest ===")
rf = RandomForestClassifier(n_estimators=1000, max_depth=5, random_state=42, n_jobs=-1)
rf.fit(X_scaled, y)
rf_importances = pd.DataFrame({
    'gene': feature_genes,
    'importance': rf.feature_importances_
}).sort_values('importance', ascending=False)
rf_set = set(rf_importances.head(30)['gene'].tolist())

rf_cv_score = cross_val_score(rf, X_scaled, y, cv=StratifiedKFold(n_cv, shuffle=True, random_state=42), scoring='roc_auc')
print(f"RF {n_cv}-fold CV AUC: {rf_cv_score.mean():.3f} +/- {rf_cv_score.std():.3f}")
print(f"RF top 5: {rf_importances.head(5)['gene'].tolist()}")

# ============= 4. SVM-RFE =============
print("\n=== SVM-RFE ===")
n_select = min(30, X_scaled.shape[1] - 5)
svm = SVC(kernel='linear', random_state=42)
rfe = RFE(svm, n_features_to_select=n_select)
rfe.fit(X_scaled, y)
svm_selected = pd.DataFrame({
    'gene': feature_genes,
    'selected': rfe.support_,
    'rank': rfe.ranking_
}).query('selected').sort_values('rank')
svm_set = set(svm_selected['gene'].tolist())
print(f"SVM-RFE selected: {len(svm_set)} genes")

# ============= 5. Consensus =============
print("\n=== Consensus Biomarkers ===")
consensus_2of3 = (lasso_set & rf_set) | (lasso_set & svm_set) | (rf_set & svm_set)
consensus_3of3 = lasso_set & rf_set & svm_set
print(f"Lasso: {len(lasso_set)}, RF: {len(rf_set)}, SVM-RFE: {len(svm_set)}")
print(f"Consensus (>=2 methods): {len(consensus_2of3)} genes")
print(f"Consensus (all 3): {len(consensus_3of3)} genes")

# Use the best available consensus
if len(consensus_3of3) >= 5:
    consensus = consensus_3of3
elif len(consensus_2of3) >= 5:
    consensus = consensus_2of3
elif len(lasso_set & rf_set) >= 5:
    consensus = lasso_set & rf_set
else:
    consensus = lasso_set if len(lasso_set) >= 3 else rf_set

consensus_df = pd.DataFrame({
    'gene': sorted(consensus),
    'in_LASSO': [g in lasso_set for g in sorted(consensus)],
    'in_RF': [g in rf_set for g in sorted(consensus)],
    'in_SVM': [g in svm_set for g in sorted(consensus)],
    'is_TF': [g in tf_set for g in sorted(consensus)]
})
consensus_df.to_csv(os.path.join(TBL, 'ML_consensus_biomarkers.csv'), index=False)

# ============= 6. ROC Evaluation =============
# Use logistic regression on consensus genes for proper multi-gene scoring
consensus_list = sorted(consensus)
consensus_idx = [feature_genes.index(g) for g in consensus_list if g in feature_genes]
X_consensus = X_scaled[:, consensus_idx]

# Cross-validated ROC for consensus panel
from sklearn.linear_model import LogisticRegression as LR
cv_skf = StratifiedKFold(n_cv, shuffle=True, random_state=42)
panel_preds = np.zeros(len(y))
panel_true = np.zeros(len(y))
for train, test in cv_skf.split(X_consensus, y):
    panel_lr = LR(max_iter=5000)
    panel_lr.fit(X_consensus[train], y[train])
    panel_preds[test] = panel_lr.predict_proba(X_consensus[test])[:, 1]
    panel_true[test] = y[test]

fpr_full, tpr_full, _ = roc_curve(panel_true, panel_preds)
auc_full = auc(fpr_full, tpr_full)

# Also fit on all data for full ROC
panel_lr_full = LR(max_iter=5000)
panel_lr_full.fit(X_consensus, y)
panel_scores = panel_lr_full.predict_proba(X_consensus)[:, 1]

# Per-gene AUC
gene_aucs = {}
for g in consensus_list:
    if g in feature_genes:
        idx = feature_genes.index(g)
        scores = X_scaled[:, idx]
        fpr, tpr, _ = roc_curve(y, scores)
        gene_aucs[g] = auc(fpr, tpr)

print(f"\nConsensus biomarker AUC: {auc_full:.3f}")
print(f"Individual gene AUCs:")
for g in sorted(gene_aucs, key=gene_aucs.get, reverse=True)[:10]:
    print(f"  {g:15s} AUC={gene_aucs[g]:.3f}")

# ============= 7. Figures =============
print("\nGenerating ML figures...")

# ROC curves
fig, ax = plt.subplots(figsize=(7, 6))
colors = plt.cm.tab10(np.linspace(0, 1, max(1, min(5, len(consensus_list)))))
for i, g in enumerate(sorted(gene_aucs, key=gene_aucs.get, reverse=True)[:5]):
    idx = feature_genes.index(g)
    scores = X_scaled[:, idx]
    fpr, tpr, _ = roc_curve(y, scores)
    ax.plot(fpr, tpr, color=colors[i], lw=1.2, alpha=0.6, label=f'{g} (AUC={gene_aucs[g]:.2f})')
ax.plot(fpr_full, tpr_full, color='#E64B35', lw=3, label=f'Consensus panel ({len(consensus_idx)} genes, AUC={auc_full:.3f})')
ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5)
ax.set_xlabel('False Positive Rate', fontsize=12)
ax.set_ylabel('True Positive Rate', fontsize=12)
ax.set_title('ROC Curves — Sleep Deprivation Biomarkers', fontsize=13, fontweight='bold')
ax.legend(loc='lower right', frameon=False, fontsize=8)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ML_ROC_Curve.png'), dpi=300, bbox_inches='tight')
plt.close()

# RF importance
fig, ax = plt.subplots(figsize=(9, 7))
top30 = rf_importances.head(30)
colors = ['#E64B35' if g in tf_set else '#4DBBD5' for g in top30['gene']]
ax.barh(range(len(top30)), top30['importance'].values[::-1], color=colors[::-1])
ax.set_yticks(range(len(top30)))
ax.set_yticklabels(top30['gene'].values[::-1], fontsize=9)
ax.set_xlabel('Mean Decrease in Impurity', fontsize=11)
ax.set_title('Random Forest — Top 30 Discriminative Features', fontsize=13, fontweight='bold')
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], color='#E64B35', lw=4, label='Transcription Factor'),
    Line2D([0],[0], color='#4DBBD5', lw=4, label='Other Gene'),
], loc='lower right', frameon=False)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ML_RF_Importance.png'), dpi=300, bbox_inches='tight')
plt.close()

# Method comparison
fig, ax = plt.subplots(figsize=(8, 5))
method_names = ['LASSO', 'Random\nForest', 'SVM-RFE', 'Consensus\n(>=2)', 'Consensus\n(all 3)']
method_counts = [len(lasso_set), len(rf_set), len(svm_set), len(consensus_2of3), len(consensus_3of3)]
bar_colors = ['#4DBBD5', '#00A087', '#E64B35', '#F39B7F', '#B09C85']
ax.bar(method_names, method_counts, color=bar_colors, edgecolor='white', linewidth=0.8)
for i, v in enumerate(method_counts):
    ax.text(i, v + 0.5, str(v), ha='center', fontweight='bold', fontsize=11)
ax.set_ylabel('Number of Genes', fontsize=11)
ax.set_title('ML Biomarker Selection — Method Comparison', fontsize=13, fontweight='bold')
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ML_Method_Comparison.png'), dpi=300, bbox_inches='tight')
plt.close()

# Expression heatmap for consensus genes
if len(consensus_list) >= 3:
    consensus_present = [g for g in consensus_list if g in raw_ml.var_names][:20]
    raw_hypo = raw[raw.obs['brain_region'] == 'Hypothalamus'].copy()
    sc.pp.normalize_total(raw_hypo, target_sum=1e4)
    sc.pp.log1p(raw_hypo)

    fig, ax = plt.subplots(figsize=(7, max(4, len(consensus_present)*0.35)))
    heat_data = np.zeros((len(consensus_present), 3))
    for gi, gene in enumerate(consensus_present):
        for ci, code in enumerate(['A1', 'A2', 'A3']):
            m = raw_hypo.obs['condition_code'] == code
            heat_data[gi, ci] = raw_hypo[m, gene].X.toarray().mean()
    heat_z = (heat_data - heat_data.mean(axis=1, keepdims=True)) / (heat_data.std(axis=1, keepdims=True) + 1e-8)

    from matplotlib.colors import TwoSlopeNorm
    norm = TwoSlopeNorm(vcenter=0)
    im = ax.imshow(heat_z, aspect='auto', cmap='RdBu_r', norm=norm)
    ax.set_xticks([0, 1, 2])
    ax.set_xticklabels(['Normal (A1)', 'SD (A2)', 'Recovery (A3)'], fontsize=10)
    ax.set_yticks(range(len(consensus_present)))
    ax.set_yticklabels(consensus_present, fontsize=9)
    ax.set_title('Consensus Biomarkers Across Conditions\n(Z-score normalized)', fontsize=12, fontweight='bold')
    plt.colorbar(im, ax=ax, shrink=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'ML_Consensus_Heatmap.png'), dpi=300, bbox_inches='tight')
    plt.close()

print("  All ML figures saved\n")

# ============= 8. Summary =============
print("=" * 60)
print("ML BIOMARKER DISCOVERY — RESULTS")
print("=" * 60)
print(f"Input: {X.shape[0]} pseudo-bulk samples, {X.shape[1]} features")
print(f"LASSO selected: {len(lasso_set)} genes (AUC={cv_lasso.best_score_:.3f})")
print(f"RF top-30 AUC: {rf_cv_score.mean():.3f}")
print(f"Consensus: {len(consensus)} genes")
print(f"Panel AUC: {auc_full:.3f}")
print(f"\nKey Biomarkers (consensus):")
for g in sorted(consensus):
    tags = ' [TF]' if g in tf_set else ''
    g_auc = gene_aucs.get(g, 0)
    print(f"  {g:15s}  AUC={g_auc:.3f}{tags}")
print(f"\nSaved: {TBL}/ML_consensus_biomarkers.csv")
