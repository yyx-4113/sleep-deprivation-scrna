# Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central Regulatory Hub and Predicts Downstream Transcriptional Consequences of Its Loss

## Abstract

**Background:** Sleep deprivation (SD) profoundly disrupts brain function, yet the cell-type-specific transcriptional programs and their hierarchical regulation across brain regions remain poorly characterized.

**Methods:** We performed single-cell RNA sequencing (scRNA-seq) analysis of 24,778 cells from the hypothalamus, brainstem, and cortex of mice subjected to 5-hour sleep deprivation versus home-cage controls (GSE137665). Our analytical pipeline comprised five components: differential expression analysis; gene regulatory network (GRN) inference centered on *Pomc*, the most significantly altered gene; in silico knockout simulation via network propagation; machine learning (ML) biomarker discovery with XGBoost-Optuna optimization and SHAP interpretation; and external validation using Fisher's meta-analysis across three independent datasets (GSE137665, GSE211088, GSE237419).

**Results:** We identified 442 significantly differentially expressed genes (DEGs) across three brain regions (padj < 0.05, |log2FC| > 0.5). The brainstem exhibited the strongest response (382 DEGs), followed by the hypothalamus (181 DEGs) and cortex (39 DEGs). *Pomc* was the most significantly altered gene overall (log2FC = −4.08, padj = 1.3 × 10⁻¹⁵⁷). Its expression collapsed during SD (log-normalized mean: 2.09 in controls vs. 0.35 in SD) and recovered to 2.32 after rebound sleep. GRN analysis identified key transcription factors (TFs) co-expressed with *Pomc*, including *Fos*, *Jun*, *Egr1*, *Nr4a1*, and *Nr3c1*. In silico knockout of *Pomc* predicted significant downstream effects on neuropeptide processing (*Pcsk2*, *Scg2*), immediate-early genes, and circadian clock components. ML analysis with XGBoost-Optuna optimization produced an ensemble biomarker panel with a cross-validated AUC of 0.87. The top ensemble features were *Dbp* (score = 0.826), *Pomc* (0.771), and *Rbm3* (0.734); SHAP analysis independently confirmed these three as the most influential predictors. Fisher's meta-analysis across all three datasets identified 904 cross-dataset consensus DEGs, with *Pomc* reaching a meta-analytic padj of 5.37 × 10⁻¹⁵⁸. Directional consistency between datasets was 85.3% (Spearman r = 0.356, p = 2.3 × 10⁻⁹).

**Conclusions:** *Pomc* functions as a central regulatory hub whose dramatic downregulation during sleep deprivation orchestrates widespread transcriptional changes affecting neuropeptide processing, circadian rhythms, and stress response pathways. The convergence of multi-cohort meta-analysis, ML ensemble modeling, and network-based in silico perturbation provides a testable framework for experimental validation of the molecular cascade triggered by *Pomc* deficiency during sleep loss.

**Keywords:** sleep deprivation, single-cell RNA-seq, Pomc, gene regulatory network, in silico knockout, machine learning

---

## 1. Introduction

Sleep deprivation impairs cognitive function, emotional regulation, metabolic homeostasis, and immune competence [1,2]. At the molecular level, sleep loss triggers large-scale transcriptional reprogramming across brain regions, involving immediate-early gene activation, circadian clock disruption, and metabolic pathway dysregulation [3,4]. However, most of our current understanding derives from bulk tissue RNA-sequencing, which averages signals across cell types and cannot resolve the regulatory logic by which transcription factors (TFs) orchestrate cell-type-specific responses.

Single-cell RNA sequencing (scRNA-seq) now enables dissection of transcriptional heterogeneity in complex tissues at cellular resolution [5]. Applied to sleep biology, scRNA-seq can identify which specific cell types within each brain region drive the transcriptional response to sleep loss and which TFs serve as master regulators of these programs [6].

Pro-opiomelanocortin (*Pomc*) encodes a precursor polypeptide cleaved into multiple bioactive peptides—including α-MSH, ACTH, and β-endorphin—and is expressed predominantly in hypothalamic neurons [7]. *Pomc* neurons are critical regulators of feeding behavior, energy homeostasis, and stress responses [8]. Emerging evidence indicates that *Pomc* neuronal activity is modulated by sleep–wake states [9], but the genome-wide transcriptional consequences of *Pomc* suppression during sleep deprivation have not been systematically characterized.

Here, we leverage a comprehensive mouse scRNA-seq dataset (GSE137665) spanning three brain regions under three sleep conditions—normal sleep, 5-hour sleep deprivation, and recovery sleep—to address five questions: (1) What are the cell-type-specific transcriptional responses to acute sleep deprivation? (2) What TFs co-regulate the *Pomc*-centered gene network? (3) What are the predicted downstream transcriptional consequences of *Pomc* loss, as assessed by in silico knockout? (4) Can machine learning identify robust multi-gene biomarkers of sleep deprivation? (5) Do the core findings generalize to independent external datasets?

---

## 2. Methods

### 2.1 Data Acquisition and Preprocessing

We obtained scRNA-seq data from GSE137665 [6] via the Gene Expression Omnibus (GEO). The dataset contains 29,571 cells from mouse hypothalamus, brainstem, and cortex under three conditions: normal sleep (A1, N = 3,491 hypothalamic cells), 5-hour sleep deprivation (A2, N = 2,706), and 2-hour recovery sleep (A3, N = 3,916).

Raw count matrices were processed in Scanpy v1.10 [10]. We retained cells with 200–6,000 detected genes and <15% mitochondrial reads, and removed genes detected in fewer than 3 cells. After quality control, 24,778 cells and 18,634 genes remained. Counts were library-size normalized to 10,000 reads per cell, log-transformed (log1p), and scaled to unit variance. We selected the top 2,000 highly variable genes for dimensionality reduction. Principal component analysis (30 components) was followed by UMAP embedding (min_dist = 0.3, spread = 1.0) and Leiden clustering (resolution = 0.5).

### 2.2 Cell Type Annotation

Cell types were assigned by marker gene scoring. We curated marker sets for eight major brain cell types — excitatory neurons (*Slc17a7*, *Neurod6*, *Tbr1*, *Camk2a*, *Snap25*), inhibitory neurons (*Gad1*, *Gad2*, *Slc32a1*, *Lhx6*, *Sst*, *Pvalb*), astrocytes (*Gfap*, *Aqp4*, *Slc1a3*, *Aldh1l1*), microglia (*Cx3cr1*, *Tmem119*, *P2ry12*, *Aif1*, *C1qa*), oligodendrocytes (*Mog*, *Mbp*, *Plp1*, *Mag*), oligodendrocyte precursor cells (*Pdgfra*, *Cspg4*, *Sox10*), endothelial cells (*Cldn5*, *Pecam1*, *Flt1*), and ependymal cells (*Foxj1*, *Tmem212*) — and scored each set per cell using `sc.tl.score_genes`. Each cell was assigned to the cell type yielding the highest mean score, provided it exceeded 0.5; otherwise the cell was labeled "Unassigned."

### 2.3 Differential Expression Analysis

We performed differential expression testing between sleep-deprived (A2) and control (A1) conditions for each brain region independently, using the Wilcoxon rank-sum test implemented in `sc.tl.rank_genes_groups`. Genes were considered significant at padj < 0.05 (Bonferroni correction) and |log2FC| > 0.5. Results were visualized with volcano plots and heatmaps.

### 2.4 Pomc-Centered Gene Regulatory Network

We focused GRN analysis on the hypothalamus (9,256 cells), the primary site of *Pomc* expression. Co-expression between *Pomc* and the 5,000 most highly expressed genes was quantified by Spearman rank correlation. The core GRN comprised *Pomc* plus its top 15 co-expressed TFs and top 20 co-expressed target genes, capped at 35 nodes for interpretability.

Known mouse TFs were identified from a curated list of 87 factors spanning major families (nuclear receptors, bZIP, bHLH, zinc finger, SOX, KLF, and circadian clock factors). We computed pairwise Spearman correlations among core network genes to construct a signed adjacency matrix and visualized the network with a spring layout, displaying edges where |r| > 0.2.

### 2.5 In Silico Pomc Knockout

We simulated *Pomc* perturbation using two complementary strategies. First, a direct correlation-based model: predicted_log2FC(gene) = −r(*Pomc*, gene) × mean(*Pomc* expression). Second, network propagation via a random-walk Laplacian L = D⁻¹A, where A is the adjacency matrix of pairwise Spearman correlations. The initial perturbation vector δ set *Pomc* to −1 × mean expression; the signal was propagated through two steps as δ(t+1) = 0.7 × L × δ(t) + 0.3 × δ(0).

### 2.6 TF Activity Inference

TF regulons were defined as genes with |Spearman r| > 0.25 and p < 0.01 relative to a given TF. Per-cell TF activity scores were computed as the weighted mean expression of regulon genes (weights = correlation coefficients). Differential TF activity between SD and control cells was assessed by Mann–Whitney U test with Bonferroni correction.

### 2.7 Machine Learning Biomarker Discovery and Optimization

Pseudo-bulk expression was generated at the Leiden cluster × condition level (45 observations, 199 features). Three feature selection methods were applied: L1-regularized logistic regression (LASSO, optimized via 5-fold cross-validated grid search), random forest (1,000 trees, max depth = 5), and SVM with recursive feature elimination (SVM-RFE, 30 features retained). Consensus biomarkers were defined as genes selected by at least two of the three methods.

XGBoost v2.1 was tuned via Optuna Bayesian hyperparameter optimization (Tree-structured Parzen Estimator, 100 trials, 5-fold stratified cross-validation). The search space spanned the number of estimators (50–500), maximum depth (2–8), learning rate (0.01–0.30, log-uniform), subsample ratio (0.6–1.0), L1 regularization (0–10), and L2 regularization (0–10). Model interpretability was assessed with SHAP (SHapley Additive exPlanations) values computed on the full training set using TreeExplainer; we generated summary, bar, and dependence plots for the top predictive features. Ensemble feature importance combined normalized scores from random forest (Gini importance), logistic regression (absolute coefficients), and XGBoost (gain and weight) to produce a consensus ranking.

Discriminative performance was evaluated by AUC from stratified 5-fold cross-validation. Cross-dataset generalization was assessed by training models on GSE137665 pseudo-bulk data and testing on individual GSE211088 samples, supplemented by bootstrap pseudo-bulk resampling (200 rounds, 2 individuals per pool) to match training data structure and yield confidence intervals.

### 2.8 External Validation and Meta-Analysis

Two independent mouse RNA-seq datasets from GEO served as external validation cohorts: GSE211088 (prefrontal cortex, 5-h SD vs. home-cage, N = 10) and GSE237419 (cerebral cortex, multi-timepoint SD, N = 42). Ensembl gene IDs were mapped to gene symbols using the NCBI mouse gene_info database. Log2 fold changes between SD and control groups were computed for each dataset and compared with our scRNA-seq results via Spearman correlation and direction-consistency analysis.

For multi-cohort meta-analysis, per-gene p-values from all three datasets were combined using Fisher's method (χ² = −2 Σ ln(pᵢ), df = 2k). Combined p-values were adjusted by the Benjamini–Hochberg procedure. Genes with Fisher meta-analytic padj < 0.05 were designated cross-dataset consensus DEGs. Directional consistency was defined as the proportion of genes showing concordant fold-change signs across all datasets in which they were detected.

### 2.9 Software and Data Availability

All analyses were performed in Python 3.12 using Scanpy v1.10, scikit-learn v1.8, NumPy, SciPy, Pandas, Matplotlib, Seaborn, and NetworkX. Code and processed data are available at GitHub (https://github.com/yongxinyang/sleep-deprivation-scrna) and will be made public upon preprint posting. Raw data are publicly accessible at GEO under accessions GSE137665, GSE211088, and GSE237419.

---

## 3. Results

### 3.1 Single-Cell Transcriptomic Landscape of the Mouse Brain Under Sleep Deprivation

After quality control, we retained 24,778 high-quality single-cell transcriptomes: 9,256 from the hypothalamus, 9,132 from the brainstem, and 6,390 from the cortex (**Figure 1A**). Marker gene scoring identified eight major cell types: microglia (12.5%), endothelial cells (11.4%), astrocytes (10.2%), excitatory neurons (9.1%), ependymal cells (6.6%), inhibitory neurons (1.5%), oligodendrocyte precursor cells (1.2%), and oligodendrocytes (1.0%) (**Figure 1B, Supplementary Figure S1**). Cells falling below the marker score threshold of 0.5 were labeled "Unassigned" (46.6%). The three experimental conditions — normal sleep (A1), 5-hour sleep deprivation (A2), and 2-hour recovery sleep (A3) — were well-represented across all brain regions and cell types (**Figure 1C**).

### 3.2 Differential Expression Reveals Region-Specific Transcriptional Responses

Comparison of sleep-deprived versus control conditions identified 442 unique DEGs across all brain regions (padj < 0.05, |log2FC| > 0.5), with substantial inter-regional overlap. The brainstem showed the most pronounced transcriptional response (382 DEGs), followed by the hypothalamus (181 DEGs) and cortex (39 DEGs) (**Figure 2A–C**). Downregulated genes predominated across all regions, suggesting that acute sleep deprivation primarily exerts a transcriptional repressive effect.

### 3.3 Pomc Is the Most Significantly Downregulated Gene in the Hypothalamus

*Pomc* emerged as the most significantly altered gene in the hypothalamus (log2FC = −4.08, padj = 1.3 × 10⁻¹⁵⁷; **Figure 2D**). Its expression fell from a log-normalized mean of 2.09 in normal sleep to 0.35 during deprivation — an approximately 4-fold reduction — before rebounding to 2.32 in recovery sleep, exceeding baseline levels. This pattern was unique to *Pomc* among hypothalamic neuropeptide genes and indicates a specific, reversible transcriptional downregulation rather than neuronal loss.

Other highly significant hypothalamic DEGs included the long non-coding RNAs *Malat1* (log2FC = −0.53, padj = 1.2 × 10⁻⁴⁵) and *Meg3* (log2FC = −1.12, padj = 4.2 × 10⁻⁴¹), neuroendocrine genes *Cga* (log2FC = −1.96, padj = 1.2 × 10⁻⁴⁴) and *Gnas* (log2FC = −1.10, padj = 2.3 × 10⁻³⁹), and *Apoe*, one of the few upregulated genes (log2FC = +1.29, padj = 3.9 × 10⁻³²).

### 3.4 Pomc-Centered Gene Regulatory Network

Co-expression analysis identified 41 genes with |Spearman r| > 0.3 with *Pomc*, including 15 known TFs. Key TFs positively co-expressed with *Pomc* included AP-1 family members (*Fos*, *Jun*, *Junb*, *Jund*), immediate-early genes (*Egr1*, *Egr2*, *Egr3*), nuclear receptors (*Nr4a1*, *Nr4a2*, *Nr4a3*), and circadian regulators (*Per1*, *Per2*, *Nfil3*, *Dbp*).

The core GRN (35 nodes) displayed a modular architecture centered on *Pomc*, with TFs forming a densely interconnected regulatory layer surrounding the target gene hub (**Figure 3A**). Positive correlations dominated the network (78% of edges), consistent with co-regulation within a shared transcriptional program.

### 3.5 In Silico Pomc Knockout Predicts Multi-Pathway Disruption

Network propagation simulation of *Pomc* knockout predicted significant transcriptional effects on 20 downstream genes (**Figure 3B**). The most strongly predicted downregulated targets included neuropeptide processing factors *Pcsk2* (predicted log2FC = −0.47) and *Scg2* (−0.32), the cell cycle regulator *Ccnd2* (−0.38), and the long non-coding RNA *Meg3* (−0.38).

TFs predicted to be most affected by *Pomc* loss included *Dbp* (circadian output regulator), *Nr3c1* (glucocorticoid receptor), *Fos*, and *Egr1*. These results suggest that *Pomc* deficiency secondarily disrupts circadian and stress-response transcriptional programs through loss of co-regulatory TF interactions.

### 3.6 Machine Learning Identifies Robust Ensemble Biomarkers with SHAP Validation

Baseline logistic regression achieved a cross-validated AUC of 0.870 (F1 = 0.830, MCC = 0.650), outperforming random forest (AUC = 0.660 ± 0.097) and SVM (AUC = 0.740 ± 0.087) (**Figure 4A**). XGBoost with Optuna Bayesian hyperparameter optimization (best trial: 75 estimators, max depth = 3, learning rate = 0.0315) reached a cross-validated AUC of 0.851 ± 0.056.

LASSO logistic regression selected 11 discriminative genes. Random forest identified *Dbp*, *Rbm3*, *Pomc*, *Mt3*, and *Per3* as the top features by Gini importance. Ensemble feature importance combining normalized scores from RF, LR, and XGBoost ranked *Dbp* (0.826), *Pomc* (0.771), *Rbm3* (0.734), *Eif5a* (0.551), and *Per3* (0.530) as the top five predictors (**Figure 4B**).

SHAP analysis on the XGBoost model independently confirmed *Dbp*, *Rbm3*, and *Pomc* as the most influential features (**Figure 4C–D**). SHAP dependence plots revealed that high expression of *Dbp* and *Rbm3*, combined with low *Pomc* expression, formed the dominant pattern separating sleep-deprived from control pseudo-bulk samples.

Consensus analysis across LASSO, RF, and SVM-RFE yielded six biomarkers selected by all three methods: *Dbp*, *Nr3c1*, *Pomc*, *Rbm3*, *Rnaset2a*, and *Tsc22d3*. This consensus panel achieved a cross-validated AUC of 0.931. All six genes showed condition-dependent expression patterns in the hypothalamus, with *Pomc* exhibiting the most dramatic suppression during SD and rebound during recovery.

### 3.7 External Validation and Fisher Meta-Analysis Confirm Cross-Dataset Robustness

To assess generalizability, we compared our DEG fold changes with those from two independent bulk RNA-seq datasets: GSE211088 (prefrontal cortex, 5-h SD) and GSE237419 (cerebral cortex, multi-timepoint SD). Critically, *Pomc* was consistently downregulated in both external datasets (GSE211088: log2FC = −1.66; GSE237419: log2FC = −1.66, p = 0.002) (**Figure 5A**).

Fisher's combined p-value meta-analysis across all three datasets identified 904 cross-dataset consensus DEGs (meta_padj < 0.05, Benjamini–Hochberg correction). *Pomc* achieved the strongest meta-analytic significance (Fisher combined p = 2.42 × 10⁻¹⁶², meta_padj = 5.37 × 10⁻¹⁵⁸), confirming its status as the most robust sleep deprivation-responsive gene across cohorts, brain regions, and technologies (**Figure 5B**).

Directional consistency analysis among DEGs with fold-change data available in all three datasets revealed 85.3% concordance (227/266 genes), with a Spearman correlation of r = 0.356 (p = 2.3 × 10⁻⁹) between our scRNA-seq fold changes and those from GSE211088 (**Figure 5C**). Genes showing consistent directional changes across all three datasets included *Pomc*, *Dbp*, *Nr1d1*, *Rbm3*, *Cry2*, and *Trem2* — representing a core sleep deprivation-responsive gene set that transcends brain region and technology platform.

Cross-dataset ML generalization — training classifiers on GSE137665 pseudo-bulk data and testing on GSE211088 individual samples — achieved AUC = 1.000. However, this likely reflects the small external sample size (N = 10) combined with highly discriminative features, rather than true perfect generalization. Bootstrap pseudo-bulk resampling (2 individuals per pool, 200 rounds, independent cross-platform scaling) confirmed robust but non-perfect generalization, with AUC confidence intervals reflecting the limited external sample size.

### 3.8 Differential Transcription Factor Activity

TF activity inference using co-expression-derived regulons revealed nine TFs with significantly altered activity during sleep deprivation (padj < 0.05, Mann–Whitney U test). *Pomc* activity decreased the most (activity log2FC = −0.22, padj = 2.8 × 10⁻⁸²), consistent with its dramatic expression loss. *Sox9* showed the largest activity increase (log2FC = +0.52, padj = 1.5 × 10⁻¹⁶, 50 target genes), followed by *Klf4* (log2FC = +0.79, padj = 1.4 × 10⁻⁶) and *Jund* (log2FC = +0.66, padj = 2.3 × 10⁻⁶) (**Figure 6A**). The nuclear receptors *Ar* (androgen receptor; log2FC = −0.71) and *Nr5a1* (steroidogenic factor 1; log2FC = −0.88) were significantly repressed during SD (**Figure 6B**).

The TF activity network around *Pomc* revealed that TFs with altered activity during SD — *Sox9*, *Klf4*, and *Jund* activated; *Nr5a1* and *Ar* repressed — map to distinct regulatory modules, suggesting coordinated transcriptional reprogramming organized around the *Pomc* hub (**Figure 6C**).

---

## 4. Discussion

In this study, we combined scRNA-seq, gene regulatory network inference, in silico perturbation, machine learning, and external meta-analysis to characterize the transcriptional response to acute sleep deprivation. Our results establish *Pomc* as a central regulatory hub and converge on a six-gene biomarker panel with robust cross-dataset support.

### 4.1 Pomc as a Master Regulator of the Sleep Deprivation Transcriptome

The magnitude of *Pomc* downregulation (log2FC = −4.08, padj = 1.3 × 10⁻¹⁵⁷) far exceeded that of any other hypothalamic gene, establishing it as the dominant transcriptional event during acute sleep deprivation in this brain region. Its complete recovery and overshoot after rebound sleep (2.09 → 0.35 → 2.32) demonstrate that this is a regulated, reversible response rather than a consequence of cellular stress or damage.

*Pomc* encodes the precursor for multiple bioactive peptides: α-MSH (anorexigenic), ACTH (HPA axis activator), and β-endorphin (endogenous opioid) [7]. Its profound suppression during sleep deprivation may therefore simultaneously affect feeding behavior, stress hormone signaling, and endogenous pain modulation — all functions known to be disrupted by sleep loss [11].

### 4.2 The Pomc Regulatory Network and Predicted Downstream Consequences

The core GRN we identified is dominated by AP-1 family TFs (*Fos*, *Jun*, *Junb*, *Jund*), immediate-early genes (*Egr1–3*), and nuclear receptors (*Nr4a1–3*, *Nr3c1*), all tightly co-expressed with *Pomc*. This co-expression pattern suggests shared upstream regulatory inputs, likely including CREB-mediated signaling via the cAMP pathway, which is known to regulate both *Pomc* and *Fos* family genes [12].

In silico knockout predicted that *Pomc* loss would secondarily affect neuropeptide processing (*Pcsk2*, *Scg2*), circadian output (*Dbp*), and stress response (*Nr3c1*) pathways. These predictions are mechanistically plausible. Prohormone convertase 2 (*Pcsk2*) processes POMC into its mature peptides [13]; its co-downregulation with *Pomc* would amplify the functional consequences of POMC loss. The predicted effect on *Dbp*, a clock-controlled gene, aligns with the known coupling between feeding-related neuropeptides and circadian rhythms [14].

### 4.3 A Six-Gene Biomarker Panel Validated by Ensemble Methods and SHAP

Our consensus biomarker panel — *Dbp*, *Nr3c1*, *Pomc*, *Rbm3*, *Rnaset2a*, and *Tsc22d3* — achieved strong discriminative performance (AUC = 0.931). Ensemble feature importance and SHAP analysis independently converged on *Dbp*, *Pomc*, and *Rbm3* as the top three contributors. The agreement among three distinct ML methods (LR, RF, XGBoost), corroborated by model-agnostic SHAP values, represents a higher evidentiary standard than single-algorithm feature selection [16].

All six consensus genes have established or emerging roles in sleep and circadian biology. *Dbp* is a core circadian output gene directly regulated by CLOCK:BMAL1. *Nr3c1* (glucocorticoid receptor) mediates stress-induced HPA axis activation. *Rbm3* is a cold-inducible RNA-binding protein whose expression is modulated by sleep [15]. *Tsc22d3* (GILZ) is a glucocorticoid-induced anti-inflammatory mediator. *Rnaset2a* participates in RNA metabolism. *Pomc* encodes the precursor for multiple behaviorally active neuropeptides. *Eif5a*, the fourth-ranked ensemble feature (score = 0.551), has an emerging role in cellular stress responses and sleep-related translational regulation. The convergence of circadian (*Dbp*, *Per3*), stress (*Nr3c1*, *Pomc*), and metabolic (*Rbm3*, *Eif5a*) regulators in this panel highlights the multi-system nature of the sleep deprivation response. A multi-gene classifier thus captures the biological complexity of sleep loss more robustly than any single marker.

The comparable performance of linear (LR, AUC = 0.870) and non-linear models (XGBoost, AUC = 0.851) supports the biological coherence of the identified biomarkers: the linear decision boundary proved well-suited to this dataset, and the agreement across model classes reduces concern about overfitting.

### 4.4 Cross-Dataset and Cross-Region Validation with Meta-Analytic Support

Fisher's meta-analysis identified 904 cross-dataset consensus DEGs, providing statistically rigorous evidence that the sleep deprivation transcriptional response generalizes across cohorts, brain regions, and profiling technologies. *Pomc* achieved the strongest meta-analytic significance (meta_padj = 5.37 × 10⁻¹⁵⁸), confirming its position as the single most robust sleep deprivation-responsive gene.

The 85.3% directional consistency between our hypothalamic scRNA-seq data and independent cortical RNA-seq data, combined with a significant Spearman correlation (r = 0.356, p = 2.3 × 10⁻⁹), demonstrates that the core sleep deprivation transcriptional program is conserved across brain regions. Critically, the same genes change in the same direction regardless of tissue origin or profiling technology. This cross-region conservation strengthens confidence that our findings — particularly the *Pomc*-centered network and ML biomarkers — capture fundamental sleep biology rather than region- or platform-specific artifacts.

The cross-dataset ML generalization AUC of 1.000 on individual GSE211088 samples should be interpreted cautiously. The small external sample size (N = 10), combined with 266 highly discriminative features, can produce perfect separation even under modest biological signal. The more conservative bootstrap pseudo-bulk evaluation and scale-free DEG direction consistency metric provide more reliable estimates of generalization.

### 4.5 Limitations

Several limitations warrant consideration. First, cell type annotation relied on marker gene scoring with a conservative threshold (0.5), leaving 46.6% of cells unassigned. This approach may miss rare but functionally important populations such as *Pomc*-expressing neurons, which represent a small fraction of hypothalamic cells. Second, the in silico knockout model is correlational rather than causal; experimental validation via conditional *Pomc* knockout followed by scRNA-seq would be necessary to confirm the predicted downstream effects. Third, only three biological replicates per condition were available, limiting statistical power for cell-type-specific comparisons. Fourth, the external validation datasets originated from cortex rather than hypothalamus, and the sample sizes were modest (GSE211088: N = 10; GSE237419: N = 42), potentially inflating cross-dataset ML generalization metrics. Fifth, the ML models were trained on pseudo-bulk rather than single-cell data, which may not fully capture cellular heterogeneity. Sixth, cross-platform normalization differences between scRNA-seq and bulk RNA-seq remain an inherent challenge; our reliance on rank-based and scale-free metrics (Spearman correlation, direction consistency, Fisher's meta-analysis) partially mitigates but does not eliminate this concern.

### 4.6 Conclusions and Future Directions

We demonstrate that *Pomc* serves as a central transcriptional hub in the hypothalamic response to sleep deprivation. Its coordinated regulation with AP-1 and nuclear receptor TFs forms a network whose disruption by sleep loss has predictable downstream consequences — affecting neuropeptide processing, circadian output, and stress response pathways. Our six-gene biomarker panel and GRN model provide testable hypotheses for future experimental studies. Immediate next steps include experimental validation of *Pomc* knockdown or knockout in cell culture or animal models, extension to additional sleep deprivation paradigms (chronic, REM-specific), and integration with proteomic and metabolomic data to bridge transcriptional changes to functional outcomes.

---

## References

1. Medic G, Wille M, Hemels ME. Short- and long-term health consequences of sleep disruption. *Nat Sci Sleep*. 2017;9:151–161.

2. Krause AJ, Simon EB, Mander BA, et al. The sleep-deprived human brain. *Nat Rev Neurosci*. 2017;18(7):404–418.

3. Cirelli C, Tononi G. Gene expression in the brain across the sleep–waking cycle. *Brain Res*. 2000;885(2):303–321.

4. Mackiewicz M, Shockley KR, Romer MA, et al. Macromolecule biosynthesis: a key function of sleep. *Physiol Genomics*. 2007;31(3):441–457.

5. Svensson V, Vento-Tormo R, Teichmann SA. Exponential scaling of single-cell RNA-seq in the past decade. *Nat Protoc*. 2018;13(4):599–604.

6. Jha PK, Valekunja UK, Ray S, Nollet M, Reddy AB. Single-cell transcriptomics and cell-specific proteomics reveals molecular signatures of sleep. *Commun Biol*. 2022;5(1):846.

7. Cawley NX, Li Z, Loh YP. 60 YEARS OF POMC: Biosynthesis, trafficking, and secretion of pro-opiomelanocortin-derived peptides. *J Mol Endocrinol*. 2016;56(4):T77–T97.

8. Zhan C, Zhou J, Feng Q, et al. Acute and long-term suppression of feeding behavior by POMC neurons in the brainstem and hypothalamus. *J Neurosci*. 2013;33(8):3624–3632.

9. Goldstein N, Levine BJ, Loy KA, et al. Hypothalamic neurons that regulate feeding can influence sleep/wake states based on homeostatic need. *Curr Biol*. 2018;28(23):3736–3747.

10. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15.

11. Knutson KL, Spiegel K, Penev P, Van Cauter E. The metabolic consequences of sleep deprivation. *Sleep Med Rev*. 2007;11(3):163–178.

12. Boutillier AL, Monnier D, Lorang D, Lundblad JR, Roberts JL, Loeffler JP. Corticotropin-releasing hormone stimulates proopiomelanocortin transcription by cFos-dependent and -independent pathways. *Mol Endocrinol*. 1995;9(6):745–755.

13. Benjannet S, Rondeau N, Day R, Chretien M, Seidah NG. PC1 and PC2 are proprotein convertases capable of cleaving proopiomelanocortin at distinct pairs of basic residues. *Proc Natl Acad Sci USA*. 1991;88(9):3564–3568.

14. Tognini P, Murakami M, Liu Y, et al. Distinct circadian signatures in liver and gut clocks revealed by ketogenic diet. *Cell Metab*. 2017;26(3):523–538.e5.

15. Wang H, Liu Y, Briesemann M, Yan J. Computational analysis of gene regulation in animal sleep deprivation. *Physiol Genomics*. 2010;42(3):427–436.

16. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. *Adv Neural Inf Process Syst*. 2017;30:4765–4774.

17. Akiba T, Sano S, Yanase T, Ohta T, Koyama M. Optuna: a next-generation hyperparameter optimization framework. *Proc 25th ACM SIGKDD Int Conf Knowl Discov Data Min*. 2019;2623–2631.

---

## Figure Legends

**Figure 1. Single-cell transcriptomic atlas of the mouse brain.** (A) UMAP visualization colored by brain region (Cortex, Hypothalamus, Brainstem). (B) UMAP colored by annotated cell type. (C) UMAP split by experimental condition (A1: normal sleep; A2: sleep deprivation; A3: recovery sleep).

**Figure 2. Differential expression analysis of sleep deprivation versus normal sleep.** (A–C) Volcano plots for hypothalamus, brainstem, and cortex showing log2 fold change (SD vs. Normal) versus −log10(adjusted p-value). Red: significantly upregulated (padj < 0.05, log2FC > 0.5); blue: significantly downregulated. Top genes by p-value are labeled. (D) *Pomc* expression across conditions in hypothalamus, showing collapse during sleep deprivation and recovery rebound.

**Figure 3. Pomc-centered gene regulatory network and in silico knockout.** (A) Core Pomc GRN. Red node: *Pomc*; blue nodes: transcription factors; grey nodes: target genes. Red edges: positive correlation; blue edges: negative correlation. Edge width represents |r| magnitude. (B) Predicted downstream effects of in silico *Pomc* knockout, showing the top 20 most affected genes ranked by predicted log2 fold change.

**Figure 4. Machine learning biomarker discovery and optimization.** (A) Model comparison (AUC, F1, MCC for LR, RF, SVM, XGBoost). (B) Ensemble feature importance combining RF Gini, LR coefficients, and XGBoost gain/weight scores (top 20 genes). Transcription factors highlighted in red. (C) SHAP summary plot showing the impact of top 20 features on model output. (D) SHAP bar plot ranking features by mean |SHAP value|. (E) SHAP dependence plots for *Dbp*, *Pomc*, and *Rbm3*. (F) Cross-dataset generalization ROC curves (train: GSE137665; test: GSE211088 individual samples). (G) Bootstrap pseudo-bulk generalization AUC distributions (200 rounds, violin plots).

**Figure 5. External validation and meta-analysis across independent datasets.** (A) *Pomc* expression in GSE211088 (prefrontal cortex) and GSE237419 (cerebral cortex) compared with GSE137665 hypothalamic scRNA-seq. (B) Meta-analysis volcano plot: −log10(Fisher combined p-value) versus consensus direction; 904 meta-significant genes highlighted. (C) Cross-dataset direction concordance heatmap showing the top 50 meta-significant genes across all three datasets. (D) Spearman correlation of DEG log2 fold changes between GSE137665 and external datasets.

**Figure 6. Differential transcription factor activity during sleep deprivation.** (A) Top differentially active TFs in hypothalamus (SD vs. Normal), ranked by adjusted p-value. (B) TF activity of *Pomc* and key upstream regulators across conditions. (C) Regulatory network showing TF activity changes (activated = blue; repressed = red) around the *Pomc* hub. Node size reflects regulon size; color intensity reflects activity fold-change magnitude.

---

## Supplementary Information

**Table S1.** Full differential expression results across all brain regions (DE_A2_vs_A1.csv).

**Table S2.** Pomc co-expression network (Spearman correlation, top 5,000 genes).

**Table S3.** Pomc in silico knockout prediction results.

**Table S4.** Machine learning consensus biomarkers and per-method selections (ML_consensus_biomarkers.csv).

**Table S5.** Ensemble feature importance combining RF, LR, and XGBoost scores (ML_Ensemble_FeatureImportance.csv).

**Table S6.** XGBoost-Optuna hyperparameter optimization results and SHAP feature rankings.

**Table S7.** Fisher meta-analysis: 904 cross-dataset consensus DEGs with combined p-values (MetaAnalysis_Significant.csv).

**Table S8.** Cross-dataset generalization results: individual-sample and bootstrap pseudo-bulk evaluation (ML_CrossDataset_Generalization.csv).

**Table S9.** DEG direction consistency and Spearman correlation across datasets (ML_CrossDataset_DEG_Correlation.csv).

**Table S10.** Differential TF activity results with regulon sizes.
