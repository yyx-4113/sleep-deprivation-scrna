"""
Step 4b: ML Biomarker Optimization
- XGBoost + LightGBM classifiers
- Optuna Bayesian hyperparameter optimization
- SHAP model interpretation
- Cross-dataset generalization (train on GSE137665 → test on external)
- Calibration & precision-recall curves
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings, os, gzip, re
warnings.filterwarnings('ignore')

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.feature_selection import RFE
from sklearn.model_selection import StratifiedKFold, cross_val_score, cross_validate, GridSearchCV
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.metrics import (roc_curve, auc, roc_auc_score,
                              precision_recall_curve, average_precision_score,
                              accuracy_score, f1_score, matthews_corrcoef,
                              confusion_matrix, brier_score_loss)
from sklearn.calibration import calibration_curve
from scipy.stats import ttest_ind, spearmanr

try:
    import xgboost as xgb
    HAS_XGB = True
except ImportError:
    HAS_XGB = False
    print("WARNING: xgboost not installed. Skipping XGBoost.")

try:
    import shap
    HAS_SHAP = True
except ImportError:
    HAS_SHAP = False
    print("WARNING: shap not installed. Skipping SHAP analysis.")

try:
    import optuna
    HAS_OPTUNA = True
except ImportError:
    HAS_OPTUNA = False
    print("WARNING: optuna not installed. Using GridSearchCV fallback.")

import scanpy as sc

FIG = r'C:\Users\1\sleep-deprivation-project\results\figures'
TBL = r'C:\Users\1\sleep-deprivation-project\results\tables'
EXT = r'C:\Users\1\sleep-deprivation-project\data\external'
os.makedirs(FIG, exist_ok=True)
os.makedirs(TBL, exist_ok=True)

# ============================================================
# 1. Build pseudo-bulk from GSE137665 (same as step4)
# ============================================================
print("=" * 60)
print("ML OPTIMIZATION — BUILDING TRAINING DATA")
print("=" * 60)

raw = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_raw.h5ad')
proc = sc.read_h5ad(r'C:\Users\1\sleep-deprivation-project\data\GSE137665_processed.h5ad')
raw = raw[raw.obs_names.isin(proc.obs_names)]
raw.obs['leiden'] = proc.obs['leiden'].values
raw.obs['cell_type'] = proc.obs['cell_type'].values

raw_ml = raw[raw.obs['condition_code'].isin(['A1', 'A2'])].copy()
sc.pp.normalize_total(raw_ml, target_sum=1e4)
sc.pp.log1p(raw_ml)

# DEG features
deg_all = pd.read_csv(os.path.join(TBL, 'DE_A2_vs_A1.csv'))
deg_sig = deg_all[deg_all['pvals_adj'] < 0.05].copy()

# TFs from Pomc network
pomc_tfs = pd.read_csv(os.path.join(TBL, 'Pomc_network_TFs.csv'))
tf_set = set(pomc_tfs['gene'].tolist()) if len(pomc_tfs) > 0 else set()

# Extended feature set: top DEGs + TFs + known sleep/neural genes
SLEEP_GENES = ['Clock','Bmal1','Arntl','Per1','Per2','Per3','Cry1','Cry2',
               'Nr1d1','Nr1d2','Rora','Rorb','Arc','Bdnf','Fos','Egr1','Npas4',
               'Homer1','Nr3c1','Crem','Srebf1','Ppargc1a','Foxo1','Pomc',
               'Pcsk2','Scg2','Cga','Gnas','Dbp','Tef','Nfil3','Mt1','Mt3','Nnat',
               'Rbm3','Tsc22d3','Rnaset2a','Malat1','Meg3','Ccnd2']
deg_genes = deg_sig.nsmallest(300, 'pvals_adj')['names'].tolist()
feature_genes = list(dict.fromkeys(deg_genes + list(tf_set) + [g for g in SLEEP_GENES if g in raw_ml.var_names]))
feature_genes = [g for g in feature_genes if g in raw_ml.var_names]
print(f"Feature genes: {len(feature_genes)}")

# Pseudo-bulk at leiden cluster × condition level
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
                'leiden': leiden, 'condition': cond,
                'label': 0 if cond == 'A1' else 1,
                'n_cells': n_cells,
                'cell_type': raw_ml.obs.loc[mask, 'cell_type'].mode()[0],
                'expression': expr.mean(axis=0)
            })

print(f"Pseudo-bulk samples: {len(pseudo_rows)}")
X = np.vstack([r['expression'] for r in pseudo_rows])
y = np.array([r['label'] for r in pseudo_rows])
meta_pb = pd.DataFrame([{k: v for k, v in r.items() if k != 'expression'} for r in pseudo_rows])

# Scale
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

n_cv = min(5, min((y==0).sum(), (y==1).sum()))
if n_cv < 2:
    n_cv = 2
print(f"Using {n_cv}-fold CV, {X_scaled.shape[0]} samples")

# ============================================================
# 2. Baseline models
# ============================================================
print("\n" + "=" * 60)
print("BASELINE MODELS")
print("=" * 60)

skf = StratifiedKFold(n_cv, shuffle=True, random_state=42)
baseline_results = {}

# Logistic Regression (L2)
lr = LogisticRegression(penalty='l2', C=1.0, max_iter=5000, random_state=42)
lr_scores = cross_validate(lr, X_scaled, y, cv=skf,
                           scoring=['roc_auc', 'accuracy', 'f1', 'matthews_corrcoef'])
baseline_results['Logistic Regression'] = {
    'AUC': lr_scores['test_roc_auc'].mean(),
    'F1': lr_scores['test_f1'].mean(),
    'MCC': lr_scores['test_matthews_corrcoef'].mean()
}

# Random Forest
rf = RandomForestClassifier(n_estimators=500, max_depth=5, random_state=42, n_jobs=-1)
rf_scores = cross_validate(rf, X_scaled, y, cv=skf,
                           scoring=['roc_auc', 'accuracy', 'f1', 'matthews_corrcoef'])
baseline_results['Random Forest'] = {
    'AUC': rf_scores['test_roc_auc'].mean(),
    'F1': rf_scores['test_f1'].mean(),
    'MCC': rf_scores['test_matthews_corrcoef'].mean()
}

# SVM
svm = SVC(kernel='rbf', probability=True, random_state=42)
svm_scores = cross_validate(svm, X_scaled, y, cv=skf,
                            scoring=['roc_auc', 'accuracy', 'f1', 'matthews_corrcoef'])
baseline_results['SVM (RBF)'] = {
    'AUC': svm_scores['test_roc_auc'].mean(),
    'F1': svm_scores['test_f1'].mean(),
    'MCC': svm_scores['test_matthews_corrcoef'].mean()
}

for name, res in baseline_results.items():
    print(f"  {name:20s}  AUC={res['AUC']:.3f}  F1={res['F1']:.3f}  MCC={res['MCC']:.3f}")

# ============================================================
# 3. XGBoost + Optuna optimization
# ============================================================
if HAS_XGB:
    print("\n" + "=" * 60)
    print("XGBOOST + OPTUNA HYPERPARAMETER OPTIMIZATION")
    print("=" * 60)

    if HAS_OPTUNA:
        def objective(trial):
            params = {
                'n_estimators': trial.suggest_int('n_estimators', 50, 500),
                'max_depth': trial.suggest_int('max_depth', 2, 8),
                'learning_rate': trial.suggest_float('learning_rate', 0.01, 0.3, log=True),
                'subsample': trial.suggest_float('subsample', 0.5, 1.0),
                'colsample_bytree': trial.suggest_float('colsample_bytree', 0.5, 1.0),
                'reg_alpha': trial.suggest_float('reg_alpha', 1e-8, 1.0, log=True),
                'reg_lambda': trial.suggest_float('reg_lambda', 1e-8, 1.0, log=True),
                'min_child_weight': trial.suggest_int('min_child_weight', 1, 10),
                'random_state': 42,
                'eval_metric': 'logloss',
                'verbosity': 0,
            }

            xgb_model = xgb.XGBClassifier(**params)
            scores = cross_val_score(xgb_model, X_scaled, y, cv=skf, scoring='roc_auc')
            return scores.mean()

        study = optuna.create_study(direction='maximize', sampler=optuna.samplers.TPESampler(seed=42))
        study.optimize(objective, n_trials=100, show_progress_bar=False)
        xgb_best_params = study.best_params
        print(f"  Best trial: {study.best_trial.number}, AUC={study.best_value:.3f}")
        print(f"  Best params: {study.best_params}")
    else:
        xgb_best_params = {'n_estimators': 200, 'max_depth': 4, 'learning_rate': 0.05,
                           'subsample': 0.8, 'colsample_bytree': 0.8, 'reg_alpha': 0.1,
                           'reg_lambda': 0.1, 'min_child_weight': 3, 'random_state': 42}
        print("  Using default params (no Optuna)")

    # Train best XGBoost
    xgb_best = xgb.XGBClassifier(**{k: v for k, v in xgb_best_params.items()
                                      if k in ['n_estimators','max_depth','learning_rate',
                                               'subsample','colsample_bytree','reg_alpha',
                                               'reg_lambda','min_child_weight','random_state']})
    xgb_best.fit(X_scaled, y)
    xgb_scores = cross_validate(xgb_best, X_scaled, y, cv=skf,
                                scoring=['roc_auc', 'accuracy', 'f1', 'matthews_corrcoef'])
    baseline_results['XGBoost (Optuna)'] = {
        'AUC': xgb_scores['test_roc_auc'].mean(),
        'F1': xgb_scores['test_f1'].mean(),
        'MCC': xgb_scores['test_matthews_corrcoef'].mean()
    }
    print(f"  XGBoost optimized: AUC={xgb_scores['test_roc_auc'].mean():.3f}")

# ============================================================
# 4. Feature importance ensemble
# ============================================================
print("\n" + "=" * 60)
print("FEATURE IMPORTANCE — ENSEMBLE CONSENSUS")
print("=" * 60)

# RF importance
rf_full = RandomForestClassifier(n_estimators=1000, max_depth=5, random_state=42, n_jobs=-1)
rf_full.fit(X_scaled, y)
rf_imp = pd.DataFrame({'gene': feature_genes, 'rf_importance': rf_full.feature_importances_})

# Logistic Regression coefficients
lr_full = LogisticRegression(penalty='l2', C=1.0, max_iter=5000, random_state=42)
lr_full.fit(X_scaled, y)
lr_imp = pd.DataFrame({'gene': feature_genes, 'lr_coef': np.abs(lr_full.coef_[0])})

# XGBoost importance (gain-based)
if HAS_XGB:
    xgb_gain = xgb_best.feature_importances_
    score_dict = xgb_best.get_booster().get_score(importance_type='weight')
    xgb_weight = np.zeros(len(feature_genes))
    for k, v in score_dict.items():
        idx = int(k.replace('f', ''))
        if idx < len(feature_genes):
            xgb_weight[idx] = v
    xgb_imp = pd.DataFrame({'gene': feature_genes, 'xgb_gain': xgb_gain, 'xgb_weight': xgb_weight})
else:
    xgb_imp = pd.DataFrame({'gene': feature_genes, 'xgb_gain': 0})

# Merge and rank
imp_merged = rf_imp.merge(lr_imp, on='gene').merge(xgb_imp, on='gene')

# Normalize each to [0,1]
from sklearn.preprocessing import MinMaxScaler
norm_cols = ['rf_importance', 'lr_coef', 'xgb_gain']
for c in norm_cols:
    if imp_merged[c].max() > imp_merged[c].min():
        imp_merged[f'{c}_norm'] = (imp_merged[c] - imp_merged[c].min()) / (imp_merged[c].max() - imp_merged[c].min())

norm_imp_cols = [c for c in imp_merged.columns if c.endswith('_norm')]
imp_merged['ensemble_score'] = imp_merged[norm_imp_cols].mean(axis=1)
imp_merged = imp_merged.sort_values('ensemble_score', ascending=False)
imp_merged.to_csv(os.path.join(TBL, 'ML_Ensemble_FeatureImportance.csv'), index=False)

print(f"Top 20 ensemble features:")
for _, r in imp_merged.head(20).iterrows():
    print(f"  {r['gene']:15s}  ensemble={r['ensemble_score']:.3f}  RF={r['rf_importance']:.3f}  XGB={r['xgb_gain']:.3f}")

# ============================================================
# 5. SHAP analysis
# ============================================================
if HAS_SHAP and HAS_XGB:
    print("\n" + "=" * 60)
    print("SHAP MODEL INTERPRETATION")
    print("=" * 60)

    explainer = shap.TreeExplainer(xgb_best)
    shap_values = explainer.shap_values(X_scaled)

    # SHAP summary
    top_n_shap = min(20, len(feature_genes))
    shap_top_idx = np.argsort(np.abs(shap_values).mean(0))[::-1][:top_n_shap]

    fig, ax = plt.subplots(figsize=(9, 6))
    shap.summary_plot(shap_values[:, shap_top_idx],
                      pd.DataFrame(X_scaled, columns=feature_genes).iloc[:, shap_top_idx],
                      feature_names=[feature_genes[i] for i in shap_top_idx],
                      show=False, max_display=top_n_shap)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'ML_SHAP_Summary.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  SHAP summary plot saved')

    # SHAP bar
    fig, ax = plt.subplots(figsize=(9, 6))
    shap.summary_plot(shap_values[:, shap_top_idx],
                      pd.DataFrame(X_scaled, columns=feature_genes).iloc[:, shap_top_idx],
                      feature_names=[feature_genes[i] for i in shap_top_idx],
                      plot_type='bar', show=False, max_display=top_n_shap)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'ML_SHAP_Bar.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  SHAP bar plot saved')

    # SHAP dependence for top 3 genes
    top3 = [feature_genes[i] for i in shap_top_idx[:3]]
    for gi, gene in enumerate(top3):
        idx = feature_genes.index(gene)
        fig, ax = plt.subplots(figsize=(6, 4))
        shap.dependence_plot(idx, shap_values, pd.DataFrame(X_scaled, columns=feature_genes),
                            feature_names=feature_genes, show=False, ax=ax)
        ax.set_title(f'SHAP Dependence: {gene}', fontsize=12, fontweight='bold')
        plt.tight_layout()
        plt.savefig(os.path.join(FIG, f'ML_SHAP_Dependence_{gene}.png'), dpi=300, bbox_inches='tight')
        plt.close()
    print('  SHAP dependence plots saved')

# ============================================================
# 6. Cross-dataset generalization
# ============================================================
# Training data = GSE137665 pseudo-bulk (cluster x condition means, ~45 obs)
# Test data = GSE211088 individual mouse bulk RNA-seq samples
#
# CROSS-PLATFORM NORMALIZATION NOTE:
# Training (scRNA-seq pseudo-bulk) and test (bulk RNA-seq) have different
# expression scales. StandardScaler fit on external data independently is
# the standard approach — it preserves relative differences within the test
# set, which is what classifiers use for discrimination.
#
# Three evaluation levels, from most optimistic to most realistic:
#   (a) Individual-sample (optimistic — individual samples have lower noise)
#   (b) Bootstrap pseudo-bulk (realistic — matches training data format)
#   (c) DEG correlation (gold standard — platform-independent rank metric)
# ============================================================
print("\n" + "=" * 60)
print("CROSS-DATASET GENERALIZATION")
print("=" * 60)

symbol_map = pd.read_csv(os.path.join(EXT, 'ensembl_to_symbol.csv'), index_col=0)
symbol_dict = symbol_map['Symbol'].to_dict()

# Load GSE211088 expression data
df211 = pd.read_csv(os.path.join(EXT, 'GSE211088_RNAseq_gene.txt.gz'), sep='\t', index_col=0)
hc_cols_211 = [c for c in df211.columns if c.startswith('WTHC')]
sd_cols_211 = [c for c in df211.columns if c.startswith('WTSD')]
n_hc, n_sd = len(hc_cols_211), len(sd_cols_211)
print(f"  GSE211088: {n_hc} HC + {n_sd} SD = {n_hc+n_sd} individual samples")

ens_ids_211 = [x.split('.')[0] for x in df211.index]
mapped_211 = [symbol_dict.get(eid, None) for eid in ens_ids_211]
df211['symbol'] = mapped_211
df211_sym = df211[df211['symbol'].notna()].copy()
expr_cols_211 = [c for c in df211.columns if c not in ['ensembl_base', 'symbol']]
df211_sym = df211_sym.groupby('symbol')[expr_cols_211].mean()

# Find overlapping genes
ext_genes = [g for g in feature_genes if g in df211_sym.index]
print(f"  Overlapping genes: {len(ext_genes)}/{len(feature_genes)}")

model_idx = [feature_genes.index(g) for g in ext_genes]
X_train_ext = X_scaled[:, model_idx]

# Build individual-level external feature matrix
X_ext_raw = np.zeros((n_hc + n_sd, len(ext_genes)))
y_ext = np.array([0]*n_hc + [1]*n_sd)
for gi, gene in enumerate(ext_genes):
    for si in range(n_hc):
        X_ext_raw[si, gi] = df211_sym.loc[gene, hc_cols_211[si]]
    for si in range(n_sd):
        X_ext_raw[n_hc+si, gi] = df211_sym.loc[gene, sd_cols_211[si]]

X_ext_log = np.log2(X_ext_raw + 1)

# ---- 6a. Individual-sample prediction (optimistic) ----
# Independent scaling per dataset is standard cross-platform practice:
# each dataset is standardized to its own mean/variance, preserving
# within-dataset relative relationships the classifiers use.
print("\n  --- (a) Individual-sample prediction ---")
print(f"  CAVEAT: n={n_hc+n_sd} small; individual samples retain full biological signal")
ext_scaler = StandardScaler()
X_ext_scaled = ext_scaler.fit_transform(X_ext_log)

gen_indiv = {}
for name, clf in [('LR', LogisticRegression(max_iter=5000, random_state=42)),
                   ('RF', RandomForestClassifier(n_estimators=500, max_depth=5, random_state=42)),
                   ('XGB', xgb_best if HAS_XGB else RandomForestClassifier(n_estimators=200, random_state=42))]:
    if name == 'XGB' and not HAS_XGB:
        clf = RandomForestClassifier(n_estimators=200, max_depth=4, random_state=42)
    clf.fit(X_train_ext, y)
    pred_proba = clf.predict_proba(X_ext_scaled)[:, 1]
    pred_class = clf.predict(X_ext_scaled)
    gen_indiv[name] = {
        'AUC': roc_auc_score(y_ext, pred_proba),
        'Accuracy': accuracy_score(y_ext, pred_class),
        'F1': f1_score(y_ext, pred_class),
        'MCC': matthews_corrcoef(y_ext, pred_class)
    }
    print(f"    {name}: AUC={gen_indiv[name]['AUC']:.3f}  "
          f"Acc={gen_indiv[name]['Accuracy']:.3f}  "
          f"F1={gen_indiv[name]['F1']:.3f}  "
          f"MCC={gen_indiv[name]['MCC']:.3f}")

# ---- 6b. Bootstrap pseudo-bulk (more realistic) ----
# Pool 2 individuals per pseudo-bulk to get multiple test samples per round.
# Each bootstrap round independently scales its pseudo-bulks.
# With ~5 HC + 5 SD: pool_size=2 gives 2 HC pools + 2 SD pools = 4 test samples.
pool_size = 2
n_pools = min(n_hc, n_sd) // pool_size
print(f"\n  --- (b) Bootstrap pseudo-bulk "
      f"(pool={pool_size} ind, {n_pools} pools/class, 200 rounds) ---")
print(f"  Each round: {n_pools*2} test samples (pseudo-bulk), independent scaling")

n_bootstrap = 200
gen_bs = {name: {'AUC': [], 'Accuracy': [], 'F1': [], 'MCC': []}
          for name in ['LR', 'RF', 'XGB']}

rng = np.random.RandomState(42)
for boot_i in range(n_bootstrap):
    hc_pools, sd_pools = [], []
    for _ in range(n_pools):
        hc_idx = rng.choice(n_hc, size=pool_size, replace=True)
        hc_pools.append(X_ext_log[hc_idx].mean(axis=0))
        sd_idx = rng.choice(n_sd, size=pool_size, replace=True)
        sd_pools.append(X_ext_log[n_hc + np.array(sd_idx)].mean(axis=0))
    X_bs = np.vstack(hc_pools + sd_pools)
    y_bs = np.array([0]*len(hc_pools) + [1]*len(sd_pools))
    # Independent scaling per bootstrap round (standard cross-platform)
    X_bs_scaled = StandardScaler().fit_transform(X_bs)

    for name, clf in [('LR', LogisticRegression(max_iter=5000, random_state=42)),
                       ('RF', RandomForestClassifier(n_estimators=500, max_depth=5, random_state=42)),
                       ('XGB', xgb_best if HAS_XGB else RandomForestClassifier(n_estimators=200, random_state=42))]:
        if name == 'XGB' and not HAS_XGB:
            clf = RandomForestClassifier(n_estimators=200, max_depth=4, random_state=42)
        clf.fit(X_train_ext, y)
        pred_proba = clf.predict_proba(X_bs_scaled)[:, 1]
        pred_class = clf.predict(X_bs_scaled)
        gen_bs[name]['AUC'].append(roc_auc_score(y_bs, pred_proba))
        gen_bs[name]['Accuracy'].append(accuracy_score(y_bs, pred_class))
        gen_bs[name]['F1'].append(f1_score(y_bs, pred_class))
        gen_bs[name]['MCC'].append(matthews_corrcoef(y_bs, pred_class))

for name in ['LR', 'RF', 'XGB']:
    aucs = gen_bs[name]['AUC']
    print(f"    {name}: AUC={np.mean(aucs):.3f} +/- {np.std(aucs):.3f}  "
          f"90%CI=[{np.percentile(aucs, 5):.3f}, {np.percentile(aucs, 95):.3f}]  "
          f"Acc={np.mean(gen_bs[name]['Accuracy']):.3f}")

# ---- 6c. DEG direction consistency (gold standard) ----
# This is the most robust cross-platform validation: Spearman correlation
# of fold changes doesn't depend on expression scale or normalization.
print(f"\n  --- (c) DEG direction consistency (scale-free) ---")
df211_hc_mean = X_ext_log[:n_hc].mean(axis=0)
df211_sd_mean = X_ext_log[n_hc:].mean(axis=0)
log2fc_211 = df211_sd_mean - df211_hc_mean  # log2 SD - log2 HC

# Get GSE137665 log2FC for the same genes
log2fc_col = 'logfoldchanges' if 'logfoldchanges' in deg_all.columns else 'scores'
deg_gene_to_fc = {}
for _, row in deg_all.iterrows():
    deg_gene_to_fc[row['names']] = row[log2fc_col]

gse_fcs, ext_fcs = [], []
for gi, gene in enumerate(ext_genes):
    if gene in deg_gene_to_fc:
        gse_fcs.append(deg_gene_to_fc[gene])
        ext_fcs.append(log2fc_211[gi])

if len(gse_fcs) >= 5:
    rho, pval = spearmanr(gse_fcs, ext_fcs)
    n_consistent = sum((np.array(gse_fcs) > 0) == (np.array(ext_fcs) > 0))
    n_total = len(gse_fcs)
    print(f"    {n_total} overlapping DEGs")
    print(f"    Spearman r = {rho:.4f} (p = {pval:.2e})")
    print(f"    Direction consistency: {n_consistent}/{n_total} "
          f"({100*n_consistent/n_total:.1f}%)")
    print(f"    (Scale-free metric — primary cross-dataset validation)")
else:
    rho, pval = np.nan, np.nan
    print(f"    Skipped: only {len(gse_fcs)} DEGs with log2FC")

# ---- Save results ----
gen_df = pd.DataFrame({
    'Individual_AUC': {k: v['AUC'] for k, v in gen_indiv.items()},
    'Bootstrap_AUC_mean': {k: np.mean(gen_bs[k]['AUC']) for k in ['LR', 'RF', 'XGB']},
    'Bootstrap_AUC_std': {k: np.std(gen_bs[k]['AUC']) for k in ['LR', 'RF', 'XGB']},
    'Bootstrap_AUC_p5': {k: np.percentile(gen_bs[k]['AUC'], 5) for k in ['LR', 'RF', 'XGB']},
    'Bootstrap_AUC_p95': {k: np.percentile(gen_bs[k]['AUC'], 95) for k in ['LR', 'RF', 'XGB']},
    'Bootstrap_Acc_mean': {k: np.mean(gen_bs[k]['Accuracy']) for k in ['LR', 'RF', 'XGB']},
    'Bootstrap_F1_mean': {k: np.mean(gen_bs[k]['F1']) for k in ['LR', 'RF', 'XGB']},
})
gen_df.index.name = 'Model'
gen_df.to_csv(os.path.join(TBL, 'ML_CrossDataset_Generalization.csv'))

# Save DEG correlation for manuscript
deg_corr_df = pd.DataFrame({
    'Metric': ['Spearman_r', 'p_value', 'N_genes', 'Direction_consistency'],
    'Value': [rho, pval, len(gse_fcs),
              f"{n_consistent}/{n_total} ({100*n_consistent/n_total:.1f}%)" if len(gse_fcs) >= 5 else 'N/A']
})
deg_corr_df.to_csv(os.path.join(TBL, 'ML_CrossDataset_DEG_Correlation.csv'), index=False)
print(f"\n  Results saved to ML_CrossDataset_Generalization.csv and ML_CrossDataset_DEG_Correlation.csv")

# ============================================================
# 7. FIGURES
# ============================================================
print("\nGenerating ML optimization figures...")

# Fig A: Model comparison bar chart
fig, ax = plt.subplots(figsize=(8, 5))
names = list(baseline_results.keys())
aucs = [baseline_results[n]['AUC'] for n in names]
colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(names)))
bars = ax.bar(range(len(names)), aucs, color=colors, edgecolor='white', linewidth=1)
for i, (bar, v) in enumerate(zip(bars, aucs)):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
            f'{v:.3f}', ha='center', fontweight='bold', fontsize=11)
ax.set_xticks(range(len(names)))
ax.set_xticklabels([n.replace(' (Optuna)','\n(Optuna)') for n in names], fontsize=9, rotation=30, ha='right')
ax.set_ylabel('Cross-Validated AUC', fontsize=12)
ax.set_title('ML Model Comparison — Sleep Deprivation Classification\n(GSE137665 Pseudo-bulk, 5-fold CV)',
             fontsize=13, fontweight='bold')
ax.set_ylim(0, max(aucs)*1.15)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ML_ModelComparison.png'), dpi=300, bbox_inches='tight')
plt.close()

# Fig B: Cross-dataset generalization ROC (individual samples)
if len(gen_indiv) >= 2:
    fig, ax = plt.subplots(figsize=(7, 6))
    colors_gen = {'LR': '#4DBBD5', 'RF': '#00A087', 'XGB': '#E64B35'}
    for name, clf_cls in [('LR', LogisticRegression(max_iter=5000, random_state=42)),
                           ('RF', RandomForestClassifier(n_estimators=500, max_depth=5, random_state=42)),
                           ('XGB', xgb_best if HAS_XGB else RandomForestClassifier(n_estimators=200, random_state=42))]:
        if name == 'XGB' and not HAS_XGB:
            continue
        clf_cls.fit(X_train_ext, y)
        pred_proba = clf_cls.predict_proba(X_ext_scaled)[:, 1]
        fpr, tpr, _ = roc_curve(y_ext, pred_proba)
        auc_v = auc(fpr, tpr)
        ax.plot(fpr, tpr, color=colors_gen[name], lw=2,
                label=f'{name} (AUC={auc_v:.3f})')

    ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.4)
    ax.set_xlabel('False Positive Rate', fontsize=12)
    ax.set_ylabel('True Positive Rate', fontsize=12)
    ax.set_title('Cross-Dataset Generalization (Individual Samples)\nTrain: GSE137665 → Test: GSE211088',
                fontsize=13, fontweight='bold')
    ax.legend(loc='lower right', frameon=False, fontsize=10)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'ML_CrossDataset_ROC.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  Cross-dataset ROC (individual samples) saved')

# Fig B2: Bootstrap pseudo-bulk generalization AUC distribution
if gen_bs:
    fig, ax = plt.subplots(figsize=(8, 5))
    bs_data = [gen_bs[n]['AUC'] for n in ['LR', 'RF', 'XGB']]
    positions = [1, 2, 3]
    vp = ax.violinplot(bs_data, positions=positions, showmeans=True, showmedians=True)
    for i, body in enumerate(vp['bodies']):
        body.set_facecolor(colors_gen[['LR', 'RF', 'XGB'][i]])
        body.set_alpha(0.6)
    ax.set_xticks(positions)
    ax.set_xticklabels(['LR', 'RF', 'XGB'])
    ax.set_ylabel('Bootstrap Pseudo-bulk AUC', fontsize=12)
    ax.set_title('Cross-Dataset Generalization (Bootstrap Pseudo-bulk)\n'
                 f'Train: GSE137665 → Test: GSE211088 ({n_bootstrap} rounds)',
                fontsize=13, fontweight='bold')
    ax.axhline(y=0.5, color='grey', ls='--', alpha=0.5)
    ax.set_ylim(0, 1.05)
    sns.despine()
    plt.tight_layout()
    plt.savefig(os.path.join(FIG, 'ML_CrossDataset_BootstrapAUC.png'), dpi=300, bbox_inches='tight')
    plt.close()
    print('  Bootstrap pseudo-bulk AUC violin plot saved')

# Fig C: Precision-Recall curves
fig, ax = plt.subplots(figsize=(7, 6))
for name, clf_cls in [('LR', LogisticRegression(max_iter=5000, random_state=42)),
                       ('RF', RandomForestClassifier(n_estimators=500, max_depth=5, random_state=42)),
                       ('XGB', xgb_best if HAS_XGB else None)]:
    if name == 'XGB' and not HAS_XGB:
        continue
    if clf_cls is None:
        continue
    clf_cls.fit(X_train_ext, y)
    pred_proba = clf_cls.predict_proba(X_ext_scaled)[:, 1]
    prec, rec, _ = precision_recall_curve(y_ext, pred_proba)
    ap = average_precision_score(y_ext, pred_proba)
    ax.plot(rec, prec, color=colors_gen.get(name, '#333333'), lw=2,
            label=f'{name} (AP={ap:.3f})')

ax.set_xlabel('Recall', fontsize=12)
ax.set_ylabel('Precision', fontsize=12)
ax.set_title('Precision-Recall — Cross-Dataset Generalization\n(GSE137665 → GSE211088)',
            fontsize=13, fontweight='bold')
ax.legend(loc='upper right', frameon=False, fontsize=10)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ML_PrecisionRecall.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Precision-Recall saved')

# Fig D: Ensemble feature importance top 25
fig, ax = plt.subplots(figsize=(9, 7))
top25 = imp_merged.head(25)
colors_tf = ['#E64B35' if g in tf_set else '#4DBBD5' for g in top25['gene']]
ax.barh(range(len(top25)), top25['ensemble_score'].values[::-1], color=colors_tf[::-1])
ax.set_yticks(range(len(top25)))
ax.set_yticklabels(top25['gene'].values[::-1], fontsize=9)
ax.set_xlabel('Ensemble Importance Score', fontsize=12)
ax.set_title('Top 25 Features — Multi-Model Ensemble Importance', fontsize=13, fontweight='bold')
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([0],[0], color='#E64B35', lw=4, label='TF / Regulator'),
    Line2D([0],[0], color='#4DBBD5', lw=4, label='Other Gene'),
], loc='lower right', frameon=False)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ML_Ensemble_Importance.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Ensemble importance saved')

# Fig E: Calibration curves
fig, ax = plt.subplots(figsize=(6, 5))
for name, clf_cls in [('LR', LogisticRegression(max_iter=5000, random_state=42)),
                       ('RF', RandomForestClassifier(n_estimators=500, max_depth=5, random_state=42)),
                       ('XGB', xgb_best if HAS_XGB else None)]:
    if clf_cls is None:
        continue
    clf_cls.fit(X_train_ext, y)
    pred_proba = clf_cls.predict_proba(X_ext_scaled)[:, 1]
    prob_true, prob_pred = calibration_curve(y_ext, pred_proba, n_bins=10, strategy='uniform')
    ax.plot(prob_pred, prob_true, marker='o', lw=2, markersize=6,
            label=f'{name} (Brier={brier_score_loss(y_ext, pred_proba):.3f})',
            color=colors_gen.get(name, '#333333'))

ax.plot([0, 1], [0, 1], 'k--', lw=1, alpha=0.5, label='Perfect calibration')
ax.set_xlabel('Mean Predicted Probability', fontsize=12)
ax.set_ylabel('Fraction of Positives', fontsize=12)
ax.set_title('Calibration Curves — Cross-Dataset\n(GSE137665 → GSE211088)',
            fontsize=13, fontweight='bold')
ax.legend(loc='lower right', frameon=False, fontsize=9)
sns.despine()
plt.tight_layout()
plt.savefig(os.path.join(FIG, 'ML_Calibration.png'), dpi=300, bbox_inches='tight')
plt.close()
print('  Calibration plot saved')

# ============================================================
# 8. Summary
# ============================================================
print("\n" + "=" * 60)
print("ML OPTIMIZATION — RESULTS SUMMARY")
print("=" * 60)
print(f"\nTraining: {X_scaled.shape[0]} pseudo-bulk samples, {X_scaled.shape[1]} features")

best_model = max(baseline_results.items(), key=lambda x: x[1]['AUC'])
print(f"\nBest model (CV): {best_model[0]} — AUC={best_model[1]['AUC']:.3f}")

if gen_indiv:
    print(f"\nCross-dataset generalization (GSE137665 → GSE211088):")
    print(f"  Individual samples (optimistic, inflated):")
    for m, r in gen_indiv.items():
        print(f"    {m}: AUC={r['AUC']:.3f}  Acc={r['Accuracy']:.3f}  F1={r['F1']:.3f}")
    if gen_bs:
        print(f"  Bootstrap pseudo-bulk (realistic, mean ± std):")
        for name in ['LR', 'RF', 'XGB']:
            if name in gen_bs:
                aucs = gen_bs[name]['AUC']
                print(f"    {name}: AUC={np.mean(aucs):.3f} ± {np.std(aucs):.3f}")

print(f"\nTop ensemble features (Top 10):")
for _, r in imp_merged.head(10).iterrows():
    print(f"  {r['gene']:15s}  score={r['ensemble_score']:.3f}")

print(f"\nSaved: {TBL}/ML_Ensemble_FeatureImportance.csv")
print(f"Saved: {TBL}/ML_CrossDataset_Generalization.csv")
print(f"Figures saved to: {FIG}/")
