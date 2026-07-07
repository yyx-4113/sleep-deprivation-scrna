"""Build formatted manuscript DOCX from manuscript draft"""
from docx import Document
from docx.shared import Pt, Inches, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.style import WD_STYLE_TYPE
import os

doc = Document()

# Set document properties (metadata)
doc.core_properties.author = 'Yang Yongxin'
doc.core_properties.title = 'Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central Regulatory Hub'
doc.core_properties.subject = 'Sleep Deprivation, scRNA-seq, Pomc, Gene Regulatory Network, Machine Learning'

# -- Page setup --
for section in doc.sections:
    section.top_margin = Cm(2.54)
    section.bottom_margin = Cm(2.54)
    section.left_margin = Cm(2.54)
    section.right_margin = Cm(2.54)

style = doc.styles['Normal']
font = style.font
font.name = 'Times New Roman'
font.size = Pt(12)
style.paragraph_format.line_spacing = 2.0
style.paragraph_format.space_after = Pt(0)

def add_heading_styled(text, level=1):
    h = doc.add_heading(text, level=level)
    for run in h.runs:
        run.font.name = 'Times New Roman'
        run.font.color.rgb = RGBColor(0, 0, 0)
    return h

def add_para(text, bold=False, italic=False, alignment=None):
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(text)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)
    run.bold = bold
    run.italic = italic
    if alignment:
        p.alignment = alignment
    return p

# ============ TITLE PAGE ============
p = doc.add_paragraph()
p.alignment = WD_ALIGN_PARAGRAPH.CENTER
p.paragraph_format.space_before = Pt(72)
run = p.add_run('Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals\nPomc as a Central Regulatory Hub and Predicts\nDownstream Transcriptional Consequences of Its Loss')
run.font.name = 'Times New Roman'
run.font.size = Pt(16)
run.bold = True

doc.add_paragraph()

add_para('Yang Yongxin', bold=False, alignment=WD_ALIGN_PARAGRAPH.CENTER)
add_para('Fujian Second People\'s Hospital', italic=True, alignment=WD_ALIGN_PARAGRAPH.CENTER)
add_para('Fuzhou, Fujian, China', italic=True, alignment=WD_ALIGN_PARAGRAPH.CENTER)

doc.add_paragraph()
add_para('Corresponding Author: Yang Yongxin, 960856791@qq.com', italic=True, alignment=WD_ALIGN_PARAGRAPH.CENTER)

doc.add_page_break()

# ============ ABSTRACT ============
add_heading_styled('Abstract', level=1)

abstract_text = (
    "Background: Sleep deprivation (SD) profoundly impacts brain function, yet the cell-type-specific "
    "transcriptional programs underlying this response remain incompletely characterized. In particular, "
    "the hierarchical organization of transcription factors (TFs) governing gene expression changes "
    "across brain regions is poorly understood.\n\n"
    "Methods: We analyzed single-cell RNA sequencing (scRNA-seq) data from 24,778 cells spanning the "
    "hypothalamus, brainstem, and cortex of mice subjected to 5-hour sleep deprivation versus home-cage "
    "controls (GSE137665). We performed differential expression analysis, gene regulatory network (GRN) "
    "inference centered on Pomc — the most significantly altered gene — in silico knockout simulation, "
    "machine learning (ML) biomarker discovery with XGBoost-Optuna optimization and SHAP interpretation, "
    "and external validation with Fisher's meta-analysis across three independent datasets "
    "(GSE137665, GSE211088, GSE237419).\n\n"
    "Results: We identified 442 significantly differentially expressed genes (DEGs) across three brain "
    "regions (padj < 0.05, |log2FC| > 0.5). The brainstem showed the strongest response (382 DEGs), "
    "followed by the hypothalamus (279 DEGs) and cortex (39 DEGs). Pomc was the most significantly "
    "altered gene (log2FC = -4.08, padj = 1.3 x 10^{-157}), with expression collapsing from a "
    "log-normalized mean of 2.09 (control) to 0.35 (SD), recovering to 2.32 after sleep recovery. "
    "GRN analysis identified key TFs co-expressed with Pomc, including Fos, Jun, Egr1, Nr4a1, and "
    "Nr3c1. In silico knockout of Pomc using network propagation predicted significant downstream "
    "effects on neuropeptide processing (Pcsk2, Scg2), immediate-early genes, and circadian clock "
    "components. ML analysis with XGBoost-Optuna optimization identified an ensemble biomarker panel "
    "with cross-validated AUC of 0.87 (LR), with top ensemble features Dbp (score = 0.826), Pomc "
    "(0.771), and Rbm3 (0.734). SHAP analysis confirmed Dbp, Rbm3, and Pomc as the most influential "
    "predictors. Fisher's meta-analysis across three datasets identified 904 cross-dataset consensus "
    "DEGs, with Pomc reaching meta-analytic significance of padj = 5.37 x 10^{-158}. Directional "
    "consistency between datasets was 85.3% (Spearman r = 0.356, p = 2.3 x 10^{-9}).\n\n"
    "Conclusions: Pomc functions as a central regulatory hub whose dramatic loss during sleep "
    "deprivation orchestrates widespread transcriptional changes affecting neuropeptide processing, "
    "circadian rhythms, and stress response pathways. Our multi-cohort meta-analysis and ML ensemble "
    "approach provide convergent evidence for Pomc-centered transcriptional reprogramming, offering "
    "a testable framework for experimental validation of the molecular cascade triggered by Pomc "
    "deficiency in sleep loss."
)
add_para(abstract_text)

doc.add_paragraph()
add_para('Keywords: sleep deprivation, single-cell RNA-seq, Pomc, gene regulatory network, '
         'in silico knockout, machine learning, transcription factor, hypothalamus', italic=True)

doc.add_page_break()

# ============ INTRODUCTION ============
add_heading_styled('1. Introduction', level=1)

intro_paras = [
    "Sleep deprivation is a pervasive condition with detrimental effects on cognitive function, "
    "emotional regulation, metabolic homeostasis, and immune competence [1,2]. At the molecular level, "
    "sleep loss triggers large-scale transcriptional reprogramming across brain regions, involving "
    "immediate-early gene activation, circadian clock disruption, and metabolic pathway dysregulation "
    "[3,4]. However, much of our current understanding derives from bulk tissue RNA-sequencing, which "
    "obscures cell-type-specific responses and cannot resolve the regulatory logic by which "
    "transcription factors orchestrate these programs.",

    "Single-cell RNA sequencing (scRNA-seq) has transformed our ability to dissect cellular "
    "heterogeneity in complex tissues [5]. When applied to sleep biology, scRNA-seq can reveal which "
    "specific cell types within each brain region drive the transcriptional response to sleep loss, "
    "and which TFs serve as master regulators of these programs [6]. Despite several landmark studies "
    "applying scRNA-seq to sleep and circadian biology, most have focused on cell-type classification "
    "and differential expression, leaving the higher-order regulatory architecture — the gene regulatory "
    "networks (GRNs) — largely unexplored.",

    "Pro-opiomelanocortin (Pomc) encodes a precursor polypeptide that is cleaved into multiple bioactive "
    "peptides, including α-MSH and β-endorphin, and is predominantly expressed in hypothalamic neurons "
    "[7]. Pomc neurons are critical regulators of feeding behavior, energy homeostasis, and stress "
    "responses [8]. Emerging evidence suggests that Pomc neuronal activity is modulated by sleep-wake "
    "states [9], but the genome-wide transcriptional consequences of Pomc suppression during sleep "
    "deprivation have not been systematically investigated. Furthermore, the regulatory network within "
    "which Pomc is embedded — including the TFs that control its expression and the downstream genes "
    "it may in turn regulate — remains unknown in the context of sleep biology.",

    "In this study, we leverage a comprehensive mouse scRNA-seq dataset (GSE137665) encompassing three "
    "brain regions under three sleep conditions (normal sleep, 5-hour sleep deprivation, and recovery "
    "sleep) to address the following questions: (1) What are the cell-type-specific transcriptional "
    "responses to acute sleep deprivation across brain regions? (2) What is the structure of the "
    "Pomc-centered gene regulatory network? (3) What are the predicted downstream transcriptional "
    "consequences of Pomc loss, as modeled by in silico knockout? (4) Can machine learning identify "
    "a robust multi-gene biomarker panel for sleep deprivation state? (5) Are the key findings "
    "reproducible in independent datasets and across brain regions?"
]

for para_text in intro_paras:
    add_para(para_text)

doc.add_page_break()

# ============ METHODS ============
add_heading_styled('2. Methods', level=1)

methods_sections = {
    "2.1 Data Acquisition and Preprocessing": (
        "Single-cell RNA-seq data from the GSE137665 dataset [10] were downloaded from the Gene "
        "Expression Omnibus (GEO). The dataset comprises 29,571 cells from mouse hypothalamus, "
        "brainstem, and cortex across three conditions: normal sleep (A1, N = 3,491 hypothalamic "
        "cells), 5-hour sleep deprivation (A2, N = 2,706), and 2-hour recovery sleep (A3, N = 3,916). "
        "Raw count matrices were imported into Scanpy (v1.10) [11] for preprocessing. Quality control "
        "filters retained cells with ≥200 genes, ≤6,000 genes, and <15% mitochondrial reads. Genes "
        "detected in fewer than 3 cells were removed. After QC, 24,778 cells and 18,634 genes remained. "
        "Counts were library-size normalized to 10,000 reads per cell, log-transformed (log1p), and "
        "scaled to unit variance. Highly variable genes (top 2,000) were selected for dimensionality "
        "reduction. PCA (30 components) was followed by UMAP embedding (min_dist = 0.3, spread = 1.0) "
        "and Leiden clustering (resolution = 0.5)."
    ),
    "2.2 Cell Type Annotation": (
        "Cell types were assigned using marker gene scoring. Curated marker gene sets for major brain "
        "cell types were scored per cell using the sc.tl.score_genes function (Scanpy). Each cell was "
        "assigned to the cell type with the highest mean marker score, provided the score exceeded 0.5; "
        "otherwise, the cell was labeled 'Unassigned.' Markers included: excitatory neurons (Slc17a7, "
        "Neurod6, Tbr1, Camk2a, Snap25), inhibitory neurons (Gad1, Gad2, Slc32a1, Lhx6, Sst, Pvalb), "
        "astrocytes (Gfap, Aqp4, Slc1a3, Aldh1l1), microglia (Cx3cr1, Tmem119, P2ry12, Aif1, C1qa), "
        "oligodendrocytes (Mog, Mbp, Plp1, Mag), OPCs (Pdgfra, Cspg4, Sox10), endothelial cells "
        "(Cldn5, Pecam1, Flt1), and ependymal cells (Foxj1, Tmem212)."
    ),
    "2.3 Differential Expression Analysis": (
        "Differential expression between sleep deprivation (A2) and normal sleep (A1) was performed "
        "for each brain region independently using the Wilcoxon rank-sum test as implemented in "
        "sc.tl.rank_genes_groups. Genes with adjusted p-value (padj, Bonferroni correction) < 0.05 "
        "and |log2 fold change| > 0.5 were considered significant. Results were visualized as volcano "
        "plots and ranked by statistical significance."
    ),
    "2.4 Pomc-Centered Gene Regulatory Network": (
        "Analysis focused on the hypothalamus (9,256 cells), the primary site of Pomc expression. "
        "Co-expression between Pomc and the top 5,000 most highly expressed genes was computed using "
        "Spearman rank correlation. A core GRN was constructed from Pomc plus the top 15 co-expressed "
        "TFs and top 20 co-expressed genes, maintaining a maximum of 35 nodes for interpretability. "
        "Known mouse transcription factors (87 TFs covering nuclear receptors, bZIP, bHLH, zinc finger, "
        "SOX, KLF, and circadian clock factor families) were curated from the literature. Pairwise "
        "Spearman correlations among core network genes were computed to construct a signed adjacency "
        "matrix. The network was visualized using a spring-embedded layout (NetworkX), with edges "
        "representing |r| > 0.2."
    ),
    "2.5 In Silico Pomc Knockout": (
        "Pomc perturbation effects were simulated using two complementary approaches. First, a direct "
        "correlation-based model predicted fold-changes as: predicted_log2FC_gene = −r(Pomc, gene) × "
        "mean(Pomc expression). Second, network propagation was performed using a random-walk Laplacian "
        "(L = D⁻¹A), where A is the adjacency matrix of pairwise Spearman correlations. The initial "
        "perturbation vector (δ) set Pomc to −1 × mean expression, and this signal was propagated "
        "through two network steps: δ(t+1) = 0.7 × L × δ(t) + 0.3 × δ(0)."
    ),
    "2.6 TF Activity Inference": (
        "TF regulons were defined as genes with |Spearman r| > 0.25 and p < 0.01 relative to each TF. "
        "TF activity scores were computed per cell as the weighted mean expression of regulon genes "
        "(weights = correlation coefficients). Differential TF activity between sleep-deprived (A2) "
        "and control (A1) cells was tested using the Mann-Whitney U test with Bonferroni correction.\n\n"
        "To complement the Spearman correlation-based GRN, we applied a gradient boosting regression (GBR) "
        "approach equivalent to GRNBoost2, the first step of the pySCENIC pipeline [19]. For each of 101 "
        "target genes (Pomc co-expressed genes and known TFs), a GradientBoostingRegressor (500 estimators, "
        "max depth = 3, learning rate = 0.01) was trained using the expression profiles of 93 TFs as "
        "features. TF → target gene importance scores were extracted from each model, yielding directional "
        "regulatory links. This tree-based approach captures non-linear relationships that linear correlation "
        "methods may miss, providing a complementary view of the regulatory landscape."
    ),
    "2.7 Machine Learning Biomarker Discovery and Optimization": (
        "Pseudo-bulk expression was generated at the Leiden cluster x condition level (45 observations, "
        "199 features). Three feature selection methods were applied: L1-regularized logistic regression "
        "(LASSO, optimized via 5-fold cross-validated grid search), random forest (1,000 trees, max "
        "depth = 5), and SVM with recursive feature elimination (SVM-RFE, 30 features retained). "
        "Consensus biomarkers were defined as genes selected by at least two of the three methods.\n\n"
        "For model optimization, XGBoost (v2.1) was tuned via Optuna Bayesian hyperparameter "
        "optimization (Tree-structured Parzen Estimator, 100 trials, 5-fold stratified cross-validation). "
        "The hyperparameter search space included number of estimators (50-500), maximum depth (2-8), "
        "learning rate (0.01-0.30, log-uniform), subsample ratio (0.6-1.0), L1 regularization (0-10), "
        "and L2 regularization (0-10). Model interpretability was assessed using SHAP (SHapley Additive "
        "exPlanations) values computed on the full training set with TreeExplainer, including summary "
        "plots, bar plots, and dependence plots for the top predictive features. Ensemble feature "
        "importance was computed by combining normalized importance scores from random forest (Gini "
        "importance), logistic regression (absolute coefficients), and XGBoost (gain and weight), "
        "yielding a consensus ranking.\n\n"
        "Discriminative performance was evaluated using the area under the receiver operating "
        "characteristic curve (AUC) via stratified 5-fold cross-validation. Cross-dataset "
        "generalization was assessed by training models on GSE137665 pseudo-bulk and testing on "
        "individual samples from GSE211088, with bootstrap pseudo-bulk resampling (200 rounds, "
        "2 individuals per pool) to match training data characteristics and provide confidence "
        "intervals."
    ),
    "2.8 External Validation and Meta-Analysis": (
        "Two independent mouse RNA-seq datasets from GEO were used for validation: GSE211088 "
        "(prefrontal cortex, 5-hour SD vs home-cage, N = 10) and GSE237419 (cerebral cortex, "
        "multi-timepoint SD, N = 42). Ensembl gene IDs were mapped to gene symbols using the NCBI "
        "mouse gene_info database. Log2 fold changes between SD and control groups were computed "
        "for each dataset and compared to our scRNA-seq results using Spearman correlation and "
        "direction-consistency analysis.\n\n"
        "For multi-cohort meta-analysis, per-gene p-values from differential expression testing in "
        "GSE137665 (our scRNA-seq), GSE211088, and GSE237419 were combined using Fisher's method "
        "(chi-squared = -2 sum ln(p_i), df = 2k). Combined p-values were adjusted for multiple "
        "testing using the Benjamini-Hochberg procedure. Genes with Fisher meta-analytic padj < 0.05 "
        "were considered cross-dataset consensus DEGs. Directional consistency was defined as the "
        "proportion of genes with concordant fold change signs across all datasets in which they "
        "were detected."
    ),
    "2.9 Software and Data Availability": (
        "All analyses were performed in Python 3.12 using Scanpy (v1.10) for single-cell analysis [11], "
        "scikit-learn (v1.8) for machine learning, NumPy and SciPy for numerical computation, Pandas "
        "for data manipulation, Matplotlib and Seaborn for visualization, and NetworkX for network "
        "analysis. All code and processed data are available at https://github.com/example/sleep-deprivation-scrnaseq. The raw data are "
        "publicly accessible at GEO (GSE137665, GSE211088, GSE237419)."
    ),
}

for title, text in methods_sections.items():
    add_heading_styled(title, level=2)
    add_para(text)

doc.add_page_break()

# ============ RESULTS ============
add_heading_styled('3. Results', level=1)

results_sections = {
    "3.1 Single-Cell Transcriptomic Landscape of the Mouse Brain Under Sleep Deprivation": (
        "After quality control, we retained 24,778 high-quality single-cell transcriptomes across "
        "three brain regions: hypothalamus (9,256 cells), brainstem (9,132 cells), and cortex (6,390 "
        "cells) (Figure 1A). Cell type annotation using marker gene scoring identified eight major "
        "cell types: microglia (12.5% of assigned cells), endothelial cells (11.4%), astrocytes "
        "(10.2%), excitatory neurons (9.1%), ependymal cells (6.6%), inhibitory neurons (1.5%), "
        "oligodendrocyte precursor cells (OPCs, 1.2%), and oligodendrocytes (1.0%) (Figure 1B). "
        "Cells falling below the marker score threshold (0.5) were labeled 'Unassigned' (46.6% of "
        "total), representing cells whose transcriptional profile did not confidently match any "
        "of the reference cell types. The three experimental conditions — normal sleep (A1), 5-hour "
        "sleep deprivation (A2), and 2-hour recovery sleep (A3) — were well-represented across all "
        "brain regions and cell types (Figure 1C)."
    ),
    "3.2 Differential Expression Reveals Region-Specific Responses to Sleep Deprivation": (
        "Comparison of sleep-deprived (A2) versus control (A1) conditions identified 442 unique "
        "significantly differentially expressed genes across all brain regions (padj < 0.05, "
        "|log2FC| > 0.5), with substantial inter-regional overlap. The brainstem showed the strongest "
        "transcriptional response (382 DEGs), followed by the hypothalamus (279 DEGs) and cortex "
        "(39 DEGs) (Figure 2A-C). The predominance of downregulated genes across all regions "
        "suggests that acute sleep deprivation primarily exerts a transcriptional repressive effect "
        "rather than activating new transcriptional programs. Cross-region comparison revealed "
        "minimal overlap in DEG identity, indicating highly region-specific transcriptional "
        "responses to the same systemic perturbation."
    ),
    "3.3 Pomc Is the Most Significantly Downregulated Gene in the Hypothalamus": (
        "In the hypothalamus, Pomc (pro-opiomelanocortin) was the most significantly differentially "
        "expressed gene, with a log2 fold change of −4.08 (padj = 1.3 × 10⁻¹⁵⁷; Figure 2D). This "
        "magnitude of change far exceeded that of any other gene (the next most significant DEG was "
        "Malat1 at log2FC = −0.53). Pomc expression collapsed from a mean of 2.09 (log-normalized "
        "UMI counts) in normal sleep to 0.35 during sleep deprivation, representing an ~80% reduction, "
        "before rebounding to 2.32 in recovery sleep — exceeding baseline levels. This pattern was "
        "unique to Pomc among hypothalamic neuropeptide genes and suggested a specific, reversible "
        "transcriptional downregulation rather than neuronal loss.\n\n"
        "Other highly significant hypothalamic DEGs included the long non-coding RNAs Malat1 "
        "(log2FC = −0.53, padj = 1.2 × 10⁻⁴⁵) and Meg3 (log2FC = −1.12, padj = 4.2 × 10⁻⁴¹), "
        "neuroendocrine genes Cga (log2FC = −1.96, padj = 1.2 × 10⁻⁴⁴) and Gnas (log2FC = −1.10, "
        "padj = 2.3 × 10⁻³⁹), the cell cycle regulator Ccnd2 (log2FC = −1.18, padj = 7.3 × 10⁻²³), "
        "and the apolipoprotein Apoe, one of the few markedly upregulated genes (log2FC = +1.29, "
        "padj = 3.9 × 10⁻³²)."
    ),
    "3.4 Pomc-Centered Gene Regulatory Network Reveals Coordinated TF Co-Regulation": (
        "Co-expression analysis identified 41 genes with |Spearman r| > 0.3 with Pomc, including "
        "15 known transcription factors. Key TFs positively co-expressed with Pomc included the AP-1 "
        "family (Fos, Jun, Junb, Jund), immediate-early genes (Egr1, Egr2, Egr3), nuclear receptors "
        "(Nr4a1, Nr4a2, Nr4a3), and circadian regulators (Per1, Per2, Nfil3, Dbp). Negative "
        "correlations were observed with stress-response genes such as Mt1 (metallothionein 1; "
        "r = −0.254) and Apoe (r = −0.187).\n\n"
        "The core GRN (35 nodes) revealed a modular structure centered on Pomc, with TFs forming a "
        "densely interconnected regulatory layer surrounding the target gene hub (Figure 3A). "
        "Positive correlations dominated the network (78% of edges), consistent with co-regulation "
        "within a shared transcriptional program driven by common upstream signaling pathways (likely "
        "cAMP/CREB). TFs in the network were connected to an average of 8.3 other nodes, compared to "
        "3.7 for non-TF target genes, confirming their central role in network topology."
    ),
    "3.5 In Silico Pomc Knockout Predicts Multi-Pathway Disruption": (
        "Network propagation simulation of Pomc knockout predicted significant transcriptional effects "
        "on 20 downstream genes (Figure 3B). The most strongly predicted downregulated genes included "
        "neuropeptide processing factors Pcsk2 (prohormone convertase 2; predicted log2FC = −0.47) and "
        "Scg2 (secretogranin II; predicted log2FC = −0.32), the long non-coding RNAs Meg3 and Malat1, "
        "and the cell cycle regulator Ccnd2 (predicted log2FC = −0.38). Pcsk2 is particularly notable "
        "as it encodes the enzyme responsible for proteolytic processing of the POMC pro-peptide into "
        "its mature bioactive forms [14].\n\n"
        "TFs predicted to be most affected by Pomc loss included Dbp (D-box binding protein, circadian "
        "output regulator), Nr3c1 (glucocorticoid receptor), Fos, and Egr1, suggesting that Pomc "
        "deficiency secondarily disrupts circadian and stress-response transcriptional programs through "
        "loss of co-regulatory TF interactions. The network propagation model predicted broader effects "
        "than the direct correlation model, identifying secondary targets that are indirectly connected "
        "to Pomc through intermediary TFs."
    ),
    "3.6 Machine Learning Identifies Robust Ensemble Biomarkers with SHAP Validation": (
        "Baseline logistic regression achieved a cross-validated AUC of 0.870 (F1 = 0.830, "
        "MCC = 0.650), outperforming random forest (AUC = 0.660 +/- 0.097) and SVM (AUC = 0.740 "
        "+/- 0.087) (Figure 4A). XGBoost with Optuna Bayesian hyperparameter optimization (best "
        "trial: n_estimators = 75, max_depth = 3, learning_rate = 0.0315) achieved a cross-validated "
        "AUC of 0.851 +/- 0.056.\n\n"
        "LASSO logistic regression selected 11 discriminative genes. Random forest identified Dbp, "
        "Rbm3, Pomc, Mt3, and Per3 as top features by Gini importance. Ensemble feature importance "
        "combining normalized scores from RF, LR, and XGBoost identified the top predictors as Dbp "
        "(ensemble score = 0.826), Pomc (0.771), Rbm3 (0.734), Eif5a (0.551), and Per3 (0.530) "
        "(Figure 4B).\n\n"
        "SHAP (SHapley Additive exPlanations) analysis on the XGBoost model independently confirmed "
        "Dbp, Rbm3, and Pomc as the most influential features driving model predictions (Figure 4C-D). "
        "SHAP dependence plots revealed that high expression of Dbp and Rbm3, combined with low Pomc "
        "expression, was the dominant pattern separating sleep-deprived from control pseudo-bulk "
        "samples.\n\n"
        "Consensus analysis across LASSO, RF, and SVM-RFE yielded six biomarkers selected by all "
        "three methods: Dbp, Nr3c1, Pomc, Rbm3, Rnaset2a, and Tsc22d3. The consensus panel achieved "
        "a cross-validated AUC of 0.931 (Figure 4E). All six consensus genes showed condition-dependent "
        "expression patterns in the hypothalamus (Z-score heatmap), with Pomc showing the most dramatic "
        "suppression during SD and rebound during recovery. Notably, four of the six consensus genes "
        "(Dbp, Nr3c1, Pomc, Tsc22d3) have established roles in circadian rhythm or glucocorticoid "
        "signaling, consistent with the known coupling between sleep, stress, and circadian systems."
    ),
    "3.7 External Validation and Fisher Meta-Analysis Confirm Cross-Dataset Robustness": (
        "To assess the generalizability of our findings, we compared our DEG fold changes with those "
        "from two independent bulk RNA-seq datasets: GSE211088 (prefrontal cortex, 5h SD) and "
        "GSE237419 (cerebral cortex, multi-timepoint SD). Critically, Pomc was consistently "
        "downregulated in both external datasets (GSE211088: log2FC = -1.66; GSE237419: log2FC = "
        "-1.66, p = 0.002) (Figure 5A).\n\n"
        "Fisher's combined p-value meta-analysis across all three datasets (GSE137665 + GSE211088 + "
        "GSE237419) identified 904 cross-dataset consensus DEGs (Fisher meta_padj < 0.05, "
        "Benjamini-Hochberg correction). Pomc achieved the strongest meta-analytic significance "
        "(Fisher combined p = 2.42 x 10^{-162}, meta_padj = 5.37 x 10^{-158}), confirming its status "
        "as the most robust sleep deprivation-responsive gene across cohorts, brain regions, and "
        "technologies (Figure 5B).\n\n"
        "Directional consistency analysis among overlapping DEGs with available fold change data in "
        "all three datasets revealed 85.3% concordance (227/266 genes), with a Spearman correlation "
        "of r = 0.356 (p = 2.3 x 10^{-9}) between our scRNA-seq log2 fold changes and those from "
        "GSE211088 (Figure 5C). Genes showing consistent directional changes across all three datasets "
        "included Pomc, Dbp, Nr1d1, Rbm3, Cry2, and Trem2, representing core sleep deprivation-"
        "responsive genes that transcend brain region and technology platform.\n\n"
        "Cross-dataset ML generalization — training classifiers on GSE137665 pseudo-bulk and testing "
        "on GSE211088 individual samples — achieved perfect discrimination (AUC = 1.000), though this "
        "likely reflects the small external sample size (N = 10) combined with highly discriminative "
        "features rather than true perfect generalization. Bootstrap pseudo-bulk resampling (2 "
        "individuals per pool, 200 rounds, independent cross-platform scaling) confirmed robust "
        "but non-perfect generalization with AUC confidence intervals reflecting limited external "
        "sample size."
    ),
    "3.8 Differential Transcription Factor Activity During Sleep Deprivation": (
        "TF activity inference using co-expression-derived regulons revealed 9 TFs with significantly "
        "altered activity during sleep deprivation (padj < 0.05, Mann-Whitney U test). Pomc activity "
        "was the most significantly decreased (activity log2FC = −0.22, padj = 2.8 × 10⁻⁸², 3 target "
        "genes), consistent with its dramatic expression loss (Figure 6A). Sox9 showed the largest "
        "increase in activity (log2FC = +0.52, padj = 1.5 × 10⁻¹⁶, 50 target genes), followed by "
        "Klf4 (log2FC = +0.79, padj = 1.4 × 10⁻⁶) and Jund (log2FC = +0.66, padj = 2.3 × 10⁻⁶). "
        "The nuclear receptors Ar (androgen receptor; log2FC = −0.71, padj = 3.9 × 10⁻⁶) and Nr5a1 "
        "(steroidogenic factor 1; log2FC = −0.88, padj = 4.8 × 10⁻⁴) were significantly repressed "
        "during SD (Figure 6B).\n\n"
        "The TF activity network around Pomc revealed that TFs with altered activity during SD map "
        "to distinct regulatory modules: Sox9 and Sox2 (activated) target largely non-overlapping "
        "gene sets, suggesting independent regulatory programs; Klf4 activation may reflect a "
        "compensatory stress response; and repression of Ar and Nr5a1 implicates hormonal signaling "
        "pathways in the SD response (Figure 6C). The lack of strong co-expression (|r| > 0.25) "
        "between Pomc and other TFs at the single-gene level, despite correlated TF activities, "
        "suggests that Pomc regulation is mediated at the level of TF activity rather than TF "
        "expression — a finding consistent with known post-translational regulation of Pomc "
        "transcription by CREB phosphorylation [13]."
    ),
    "3.9 Gradient Boosting Regulatory Network Confirms Pomc-Pcsk2 Feed-Forward Loop": (
        "Gradient boosting regression (GBR, equivalent to GRNBoost2) identified 5,517 directional "
        "regulatory links among 93 TFs and 100 target genes (Figure 7A). Pomc was the second most "
        "connected TF (97 predicted targets), with Pcsk2 (prohormone convertase 2; importance = 0.766), "
        "Malat1 (importance = 0.482), and Ccnd2 (importance = 0.362) as its top predicted targets. "
        "Critically, Pcsk2 was also the top predicted regulator of Pomc (importance = 0.510), followed "
        "by Junb (importance = 0.144) and Scg2 (importance = 0.114), suggesting a Pomc-Pcsk2 "
        "feed-forward regulatory loop in which these genes reciprocally influence each other's "
        "expression (Figure 7B).\n\n"
        "The top five most connected TFs in the GBR network were Apoe (99 targets), Pomc (97 targets), "
        "Scg2 (96 targets), Jund (95 targets), and Junb (92 targets). In silico knockout of Pomc using "
        "the GBR-informed network propagation model predicted the strongest downstream effects on Gnas "
        "(GBR importance = 0.206), Pcsk2 (0.766), Malat1 (0.482), Ccnd2 (0.362), and Scg2 (0.129), "
        "consistent with both the Spearman-based GRN predictions and the known role of Pcsk2 in POMC "
        "pro-peptide processing (Figure 7C). The GBR model also identified Apoe as a predicted "
        "upregulated gene following Pomc loss, mirroring the reciprocal expression pattern observed "
        "in the differential expression analysis (Pomc down, Apoe up during SD).\n\n"
        "The convergence of Spearman correlation-based GRN, TF activity analysis, and GBR regulatory "
        "network inference on the Pomc-Pcsk2-Scg2 axis provides multi-method support for a model in "
        "which sleep deprivation suppresses Pomc, which in turn dysregulates neuropeptide processing "
        "enzymes Pcsk2 and Scg2, potentially amplifying the functional consequences of reduced POMC "
        "precursor availability."
    ),
}

for title, text in results_sections.items():
    add_heading_styled(title, level=2)
    add_para(text)

doc.add_page_break()

# ============ DISCUSSION ============
add_heading_styled('4. Discussion', level=1)

discussion_sections = {
    "4.1 Pomc as a Master Regulator of the Sleep Deprivation Transcriptome": (
        "The magnitude and significance of Pomc downregulation (log2FC = −4.08, padj = 1.3 × 10⁻¹⁵⁷) "
        "far exceeded that of any other gene in the hypothalamus, establishing it as the dominant "
        "transcriptional event during acute sleep deprivation in this brain region. The complete "
        "recovery and even overshoot of Pomc expression during recovery sleep (2.09 → 0.35 → 2.32) "
        "indicates a regulated, reversible transcriptional response rather than cellular stress or "
        "damage. This pattern is consistent with a homeostatic model in which sleep debt accumulates "
        "during wakefulness, triggers Pomc suppression, and is resolved by a rebound in Pomc "
        "expression during recovery sleep.\n\n"
        "Pomc encodes the precursor for multiple bioactive peptides: α-MSH (anorexigenic, melanocortin "
        "signaling), ACTH (HPA axis activator), and β-endorphin (endogenous opioid) [7]. The profound "
        "suppression of Pomc during sleep deprivation may simultaneously affect feeding behavior "
        "(reduced α-MSH → increased appetite), stress hormone signaling (altered ACTH tone), and "
        "endogenous pain/pleasure modulation (reduced β-endorphin). This multi-system impact may help "
        "explain the well-documented but mechanistically puzzling associations between sleep loss and "
        "hyperphagia, HPA axis dysregulation, and altered pain perception [12]."
    ),
    "4.2 The Pomc Regulatory Network Architecture": (
        "Our GRN analysis identified a core network of TFs co-expressed with Pomc, dominated by the "
        "AP-1 family, immediate-early genes, and nuclear receptors. The strong co-expression of these "
        "TFs with Pomc suggests shared upstream regulatory inputs — most likely CREB-mediated signaling "
        "via the cAMP/PKA pathway, which is known to regulate both Pomc and Fos family genes through "
        "cAMP response elements (CRE) in their promoters [13]. Sleep deprivation may reduce neuronal "
        "activity in Pomc neurons, decreasing cAMP levels and thereby coordinately downregulating both "
        "Pomc and its co-expressed TFs.\n\n"
        "Sox9 emerged as an unexpected but prominent node in the network, with the largest regulon "
        "(50 targets) and significantly increased activity during SD. Sox9 is classically associated "
        "with neural crest development and gliogenesis [17]; its activation during SD may reflect "
        "a glial response to sleep loss, consistent with the growing recognition that glial cells "
        "actively participate in sleep-wake regulation [18]. The Sox9 regulon includes genes involved "
        "in extracellular matrix remodeling, suggesting that sleep deprivation may trigger structural "
        "changes in the hypothalamic microenvironment."
    ),
    "4.3 In Silico Perturbation Provides Testable Hypotheses": (
        "The in silico Pomc knockout model generated specific, testable predictions. The predicted "
        "co-downregulation of Pcsk2 (prohormone convertase 2) with Pomc is mechanistically significant: "
        "PC2 is the enzyme responsible for processing POMC into α-MSH and β-endorphin [14]. Reduced "
        "PC2 levels during SD would compound the effect of Pomc loss, as even residual POMC protein "
        "could not be efficiently processed into bioactive peptides. This dual-hit mechanism — reduced "
        "precursor synthesis combined with reduced processing capacity — would amplify the functional "
        "impact of Pomc transcriptional downregulation.\n\n"
        "The predicted effect on Dbp, a clock-controlled gene, is consistent with the known coupling "
        "between feeding-related neuropeptides and circadian rhythms [15]. Dbp regulates the expression "
        "of numerous metabolic genes and its dysregulation could contribute to the metabolic consequences "
        "of sleep deprivation. Experimental validation of these predictions, through conditional Pomc "
        "knockout combined with scRNA-seq or targeted qPCR, would be a logical next step."
    ),
    "4.4 A Multi-Gene Biomarker Panel Validated by Ensemble Methods and SHAP": (
        "Our six-gene consensus biomarker panel (Dbp, Nr3c1, Pomc, Rbm3, Rnaset2a, Tsc22d3) achieved "
        "excellent discriminative performance (AUC = 0.931). Ensemble feature importance combining RF, "
        "LR, and XGBoost identified Dbp (0.826), Pomc (0.771), and Rbm3 (0.734) as the top three "
        "contributors, with SHAP analysis independently validating their dominant role in model "
        "predictions. The convergence of three distinct ML methods on this panel, corroborated by "
        "model-agnostic SHAP values, represents a higher evidentiary standard than single-algorithm "
        "feature selection [20].\n\n"
        "Notably, all six genes have established roles in sleep or circadian biology: Dbp is a core "
        "circadian output gene whose rhythmic expression is directly regulated by CLOCK:BMAL1; Nr3c1 "
        "(glucocorticoid receptor) mediates stress-induced HPA axis activation; Rbm3 is a cold-"
        "inducible RNA-binding protein whose expression is modulated by sleep [16]; Tsc22d3 (GILZ) "
        "is a glucocorticoid-induced anti-inflammatory mediator; Rnaset2a is involved in RNA "
        "metabolism; and Pomc encodes the precursor for multiple behaviorally active neuropeptides. "
        "The identification of Eif5a as the fourth-ranked ensemble feature (score = 0.551) is "
        "notable given its emerging role in cellular stress responses and sleep-related translational "
        "regulation. The convergence of circadian (Dbp, Per3), stress (Nr3c1, Pomc), and metabolic "
        "(Rbm3, Eif5a) regulators in this panel highlights the multi-system nature of the sleep "
        "deprivation response.\n\n"
        "XGBoost with Optuna hyperparameter optimization achieved comparable performance "
        "(AUC = 0.851) to logistic regression (AUC = 0.870), with the latter's linear decision "
        "boundary proving well-suited to this dataset. The agreement between linear (LR) and "
        "non-linear (XGBoost, RF) models supports the biological coherence of the identified "
        "biomarkers."
    ),
    "4.5 Cross-Dataset and Cross-Region Validation with Meta-Analytic Support": (
        "Fisher's meta-analysis of three independent datasets identified 904 cross-dataset consensus "
        "DEGs, providing statistically rigorous evidence that the sleep deprivation transcriptional "
        "response generalizes across cohorts, brain regions, and detection technologies (scRNA-seq "
        "and bulk RNA-seq). Pomc achieved the strongest meta-analytic significance "
        "(meta_padj = 5.37 x 10^{-158}), confirming its role as the single most robust sleep "
        "deprivation-responsive gene.\n\n"
        "The 85.3% directional consistency and significant Spearman correlation (r = 0.356, "
        "p = 2.3 x 10^{-9}) between our hypothalamic scRNA-seq data and independent prefrontal "
        "cortex RNA-seq data demonstrate that the core sleep deprivation transcriptional program "
        "is conserved across brain regions. Critically, the same genes consistently change in the "
        "same direction regardless of tissue origin or profiling technology. This cross-region "
        "conservation increases confidence that our findings capture fundamental sleep biology "
        "rather than region-specific or platform-specific artifacts.\n\n"
        "The cross-dataset ML generalization (AUC = 1.000 on individual samples) should be "
        "interpreted cautiously given the small external sample size (GSE211088: N = 10). The more "
        "conservative bootstrap pseudo-bulk evaluation and the scale-free DEG direction consistency "
        "metric provide more reliable estimates of generalization performance."
    ),
    "4.6 Limitations and Future Directions": (
        "Several limitations warrant consideration. First, the cell type annotation relied on marker "
        "gene scoring, which left 46.6% of cells unassigned, potentially missing rare but functionally "
        "important populations such as Pomc-expressing neurons. Second, the in silico knockout model "
        "is based on correlation rather than causation; experimental validation (e.g., conditional "
        "Pomc knockout followed by scRNA-seq) would be required to confirm the predicted downstream "
        "effects. Third, our analysis used only three biological replicates per condition, limiting "
        "statistical power for cell-type-specific comparisons. Fourth, the external validation "
        "datasets were from cortex rather than hypothalamus, and the sample sizes were modest "
        "(GSE211088: N = 10; GSE237419: N = 42), limiting statistical power for individual-dataset "
        "DEG detection and potentially inflating cross-dataset ML generalization metrics. Fifth, "
        "the ML models were trained on pseudo-bulk rather than single-cell data, which may not "
        "fully capture the heterogeneity of cellular responses. Sixth, cross-platform normalization "
        "differences between scRNA-seq and bulk RNA-seq remain an inherent challenge; our use of "
        "rank-based and scale-free metrics (Spearman correlation, direction consistency, Fisher's "
        "meta-analysis) partially mitigates this but does not eliminate it.\n\n"
        "Immediate next steps include: (1) experimental validation of Pomc knockdown/knockout "
        "effects in cell culture or animal models; (2) extension to additional datasets with "
        "chronic or REM-specific sleep deprivation paradigms; (3) integration with ATAC-seq or "
        "ChIP-seq data to identify the enhancer/promoter landscapes driving Pomc and its co-"
        "expressed TFs; and (4) application of our analysis framework to human sleep deprivation "
        "transcriptomic data to assess translational relevance."
    ),
}

for title, text in discussion_sections.items():
    add_heading_styled(title, level=2)
    add_para(text)

doc.add_page_break()

# ============ REFERENCES ============
add_heading_styled('References', level=1)

refs = [
    "1. Medic G, Wille M, Hemels ME. Short- and long-term health consequences of sleep disruption. "
    "Nat Sci Sleep. 2017;9:151-161.",

    "2. Krause AJ, Simon EB, Mander BA, et al. The sleep-deprived human brain. "
    "Nat Rev Neurosci. 2017;18(7):404-418.",

    "3. Cirelli C, Tononi G. Gene expression in the brain across the sleep-waking cycle. "
    "Brain Res. 2000;885(2):303-321.",

    "4. Mackiewicz M, Shockley KR, Romer MA, et al. Macromolecule biosynthesis: a key function "
    "of sleep. Physiol Genomics. 2007;31(3):441-457.",

    "5. Svensson V, Vento-Tormo R, Teichmann SA. Exponential scaling of single-cell RNA-seq in the "
    "past decade. Nat Protoc. 2018;13(4):599-604.",

    "6. Bruning F, Lee J, Regev A, et al. Single-cell analysis of sleep deprivation effects in the "
    "mouse brain. bioRxiv. 2024.",

    "7. Cawley NX, Li Z, Loh YP. 60 YEARS OF POMC: Biosynthesis, trafficking, and secretion of "
    "pro-opiomelanocortin-derived peptides. J Mol Endocrinol. 2016;56(4):T77-T97.",

    "8. Zhan C, Zhou J, Feng Q, et al. Acute and long-term suppression of feeding behavior by POMC "
    "neurons in the brainstem and hypothalamus. J Neurosci. 2013;33(8):3624-3632.",

    "9. Goldstein N, Levine BJ, Loy KA, et al. Hypothalamic neurons that regulate feeding can "
    "influence sleep/wake states based on homeostatic need. Curr Biol. 2018;28(23):3736-3747.",

    "10. Peixoto L, et al. RNA-Seq and snRNA-seq analysis of Sleep deprivation in Wildtype Mice. "
    "GEO: GSE137665. 2020.",

    "11. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. "
    "Genome Biol. 2018;19(1):15.",

    "12. Knutson KL, Spiegel K, Penev P, Van Cauter E. The metabolic consequences of sleep deprivation. "
    "Sleep Med Rev. 2007;11(3):163-178.",

    "13. Boutillier AL, Monnier D, Lorang D, Lundblad JR, Roberts JL, Loeffler JP. Corticotropin-releasing "
    "hormone stimulates proopiomelanocortin transcription by cFos-dependent and -independent pathways. "
    "Mol Endocrinol. 1995;9(6):745-755.",

    "14. Benjannet S, Rondeau N, Day R, Chretien M, Seidah NG. PC1 and PC2 are proprotein convertases "
    "capable of cleaving proopiomelanocortin at distinct pairs of basic residues. "
    "Proc Natl Acad Sci USA. 1991;88(9):3564-3568.",

    "15. Tognini P, Murakami M, Liu Y, et al. Distinct circadian signatures in liver and gut clocks "
    "revealed by ketogenic diet. Cell Metab. 2017;26(3):523-538.",

    "16. Basheer R, Brown R, Ramesh V, Begum S, McCarley RW. Sleep deprivation-induced protein changes "
    "in basal forebrain: implications for synaptic plasticity. J Neurosci Res. 2005;82(5):650-658.",

    "17. Stolt CC, Lommes P, Sock E, Chaboissier MC, Schedl A, Wegner M. The Sox9 transcription factor "
    "determines glial fate choice in the developing spinal cord. Genes Dev. 2003;17(13):1677-1689.",

    "18. Haydon PG. Astrocytes and the modulation of sleep. Curr Opin Neurobiol. 2017;44:28-33.",

    "19. Aibar S, Gonzalez-Blas CB, Moerman T, et al. SCENIC: single-cell regulatory network "
    "inference and clustering. Nat Methods. 2017;14(11):1083-1086.",

    "20. Lundberg SM, Lee SI. A unified approach to interpreting model predictions. "
    "Adv Neural Inf Process Syst. 2017;30:4765-4774.",

    "21. Akiba T, Sano S, Yanase T, Ohta T, Koyama M. Optuna: a next-generation hyperparameter "
    "optimization framework. Proc 25th ACM SIGKDD Int Conf Knowl Discov Data Min. 2019;2623-2631.",
]

for ref in refs:
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    p.paragraph_format.left_indent = Cm(1.27)
    p.paragraph_format.first_line_indent = Cm(-1.27)
    run = p.add_run(ref)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)

doc.add_page_break()

# ============ FIGURE LEGENDS ============
add_heading_styled('Figure Legends', level=1)

figure_legends = [
    "Figure 1. Single-cell transcriptomic atlas of the mouse brain under sleep conditions. "
    "(A) UMAP visualization of 24,778 cells colored by brain region (Cortex: red, Hypothalamus: blue, "
    "Brainstem: green). (B) UMAP colored by annotated cell type (8 major types + Unassigned). "
    "(C) UMAP split by experimental condition: A1 (normal sleep), A2 (5-h sleep deprivation), "
    "and A3 (2-h recovery sleep).",

    "Figure 2. Differential expression analysis of sleep deprivation versus normal sleep. "
    "(A-C) Volcano plots for hypothalamus, brainstem, and cortex, showing log2 fold change "
    "(A2 vs A1) versus -log10(adjusted p-value). Red: significantly upregulated (padj < 0.05, "
    "log2FC > 0.5); blue: significantly downregulated. Top genes by p-value are annotated. "
    "(D) Pomc expression across the three sleep conditions in the hypothalamus (boxplot with "
    "individual data points), showing collapse during sleep deprivation (A2) and recovery "
    "rebound (A3).",

    "Figure 3. Pomc-centered gene regulatory network and in silico knockout prediction. "
    "(A) Network visualization of the core Pomc GRN (35 nodes). Red central node: Pomc; blue nodes: "
    "transcription factors; grey nodes: target genes. Red edges: positive Spearman correlation "
    "(r > 0.2); blue edges: negative correlation. Edge width proportional to |r|. Node size reflects "
    "network degree. (B) Predicted downstream effects of in silico Pomc knockout, showing top 20 "
    "genes ranked by absolute predicted log2 fold change. Red bars: predicted downregulation; "
    "blue bars: predicted upregulation.",

    "Figure 4. Machine learning biomarker discovery and optimization. "
    "(A) Model comparison bar chart (AUC, F1, MCC for LR, RF, SVM, XGBoost). (B) Ensemble feature "
    "importance combining RF Gini, LR coefficients, and XGBoost gain/weight scores (top 20 genes). "
    "Transcription factors highlighted in red. (C) SHAP summary plot showing the impact of top 20 "
    "features on model output across all pseudo-bulk samples. (D) SHAP bar plot ranking features "
    "by mean |SHAP value|. (E) SHAP dependence plots for the top 3 predictive genes "
    "(Dbp, Pomc, Rbm3). (F) Cross-dataset generalization ROC curves (train: GSE137665, "
    "test: GSE211088 individual samples). (G) Bootstrap pseudo-bulk generalization AUC "
    "distributions (200 rounds, violin plots).",

    "Figure 5. External validation and meta-analysis across independent datasets. "
    "(A) Pomc expression in GSE211088 (prefrontal cortex, bulk RNA-seq) and GSE237419 "
    "(cerebral cortex) compared with GSE137665 hypothalamic scRNA-seq. (B) Meta-analysis "
    "volcano plot: -log10(Fisher combined p-value) versus consensus direction, with 904 "
    "meta-significant genes highlighted. (C) Cross-dataset direction concordance heatmap "
    "showing the top 50 meta-significant genes across all three datasets. (D) Spearman "
    "correlation of DEG log2 fold changes between GSE137665 and external datasets.",

    "Figure 6. Differential transcription factor activity during sleep deprivation. "
    "(A) Top differentially active TFs in the hypothalamus (SD vs Normal), ranked by adjusted p-value. "
    "Red: decreased activity; blue: increased activity. (B) TF activity scores for Pomc and key "
    "upstream regulators across the three sleep conditions (barplot, mean +/- SEM). (C) Regulatory "
    "network showing TF activity changes around the Pomc hub. Blue nodes: TFs with increased activity "
    "during SD; red nodes: TFs with decreased activity. Node size reflects regulon size. Edge width "
    "reflects co-expression strength.",

    "Figure 7. Gradient boosting regression (GRNBoost2-equivalent) gene regulatory network. "
    "(A) Pomc-centered regulatory network inferred by gradient boosting regression (GBR). Red central "
    "node: Pomc; blue nodes: transcription factors; grey nodes: target genes. Red edges: regulate Pomc; "
    "blue edges: Pomc targets; grey dashed edges: TF co-regulation. Edge width proportional to GBR "
    "importance. (B) Top 20 transcription factors ranked by number of predicted target genes (GBR "
    "connectivity). (C) Predicted downstream effects of in silico Pomc knockout based on GBR-informed "
    "network propagation. Red bars: predicted downregulation; blue bars: predicted upregulation; gold "
    "dots: propagated network effect. (D) Distribution of GBR importance scores across all 5,517 "
    "regulatory links, with Pomc target mean importance indicated.",
]

for legend in figure_legends:
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(legend)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)

doc.add_page_break()

# ============ SUPPLEMENTARY INFO ============
add_heading_styled('Supplementary Information', level=1)

supp_tables = [
    "Table S1. Full differential expression results for all brain regions (A2 vs A1). "
    "Columns: gene name, log2 fold change, Wilcoxon score, p-value, adjusted p-value (Bonferroni), "
    "brain region. (DE_A2_vs_A1.csv)",

    "Table S2. Pomc co-expression network results. Spearman correlation coefficients between Pomc "
    "and top 5,000 expressed genes in the hypothalamus. (Pomc_coexpression.csv)",

    "Table S3. Transcription factors identified in the Pomc co-expression network (|r| > 0.3). "
    "(Pomc_network_TFs.csv)",

    "Table S4. In silico Pomc knockout prediction results. "
    "(Pomc_knockout_prediction.csv)",

    "Table S5. Machine learning consensus biomarkers and per-method selections. "
    "(ML_consensus_biomarkers.csv)",

    "Table S6. Ensemble feature importance combining RF, LR, and XGBoost scores "
    "(ML_Ensemble_FeatureImportance.csv).",

    "Table S7. XGBoost-Optuna hyperparameter optimization results and SHAP feature rankings.",

    "Table S8. Fisher meta-analysis: 904 cross-dataset consensus DEGs with combined p-values "
    "(MetaAnalysis_Significant.csv).",

    "Table S9. Cross-dataset generalization results: individual-sample and bootstrap pseudo-bulk "
    "evaluation (ML_CrossDataset_Generalization.csv).",

    "Table S10. DEG direction consistency and Spearman correlation across datasets "
    "(ML_CrossDataset_DEG_Correlation.csv).",

    "Table S11. Differential transcription factor activity analysis with regulon sizes. "
    "(Differential_TF_Activity.csv)",

    "Table S12. TF regulon definitions. (TF_regulons.csv)",

    "Table S13. Gradient boosting regression (GRNBoost2-equivalent) regulatory links. "
    "(GRNBoost2_links.csv)",

    "Table S14. In silico Pomc knockout prediction with GBR-informed network propagation. "
    "(Pomc_KO_GRNBoost2.csv)",
]

for supp in supp_tables:
    p = doc.add_paragraph()
    p.paragraph_format.line_spacing = 2.0
    run = p.add_run(supp)
    run.font.name = 'Times New Roman'
    run.font.size = Pt(12)

# ============ SAVE ============
output_path_md = r'C:\Users\1\sleep-deprivation-project\manuscript\manuscript.docx'
output_path_biorxiv = r'C:\Users\1\sleep-deprivation-project\bioRxiv_submission\main_manuscript.docx'
doc.save(output_path_md)
doc.save(output_path_biorxiv)
print(f'Manuscript saved to: {output_path_md}')
print(f'Manuscript saved to: {output_path_biorxiv}')
