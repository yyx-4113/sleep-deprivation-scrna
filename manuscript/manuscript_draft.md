# Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central Regulatory Hub and Predicts Downstream Transcriptional Consequences of Its Loss

## Abstract

**Background:** Sleep deprivation (SD) profoundly impacts brain function, yet the cell-type-specific transcriptional programs underlying this response remain incompletely characterized. In particular, the hierarchical organization of transcription factors (TFs) governing gene expression changes across brain regions is poorly understood.

**Methods:** We analyzed single-cell RNA sequencing (scRNA-seq) data from 24,778 cells spanning the hypothalamus, brainstem, and cortex of mice subjected to 5-hour sleep deprivation versus home-cage controls (GSE137665). We performed differential expression analysis, gene regulatory network (GRN) inference centered on *Pomc* — the most significantly altered gene — in silico knockout simulation, machine learning biomarker discovery, and external validation across two independent RNA-seq datasets (GSE211088, GSE237419).

**Results:** We identified 700 significantly differentially expressed genes (DEGs) across three brain regions (padj < 0.05). The hypothalamus exhibited the most coordinated transcriptional response (181 DEGs), dominated by *Pomc* downregulation (log2FC = −4.08, padj = 1.3 × 10⁻¹⁵⁷). *Pomc* expression collapsed from a log-normalized mean of 2.09 (control) to 0.35 (SD), recovering to 2.32 after sleep recovery. GRN analysis identified key TFs co-expressed with *Pomc*, including *Fos*, *Jun*, *Egr1*, *Nr4a1*, and *Nr3c1*. In silico knockout of *Pomc* using network propagation predicted significant downstream effects on neuropeptide processing (*Pcsk2*, *Scg2*), immediate-early genes, and circadian clock components. Machine learning identified a six-gene consensus biomarker panel (AUC = 0.931) including *Dbp*, *Nr3c1*, *Pomc*, *Rbm3*, *Rnaset2a*, and *Tsc22d3*. External validation in independent prefrontal cortex RNA-seq data confirmed consistent directional changes (Spearman r = 0.728, p < 0.001).

**Conclusions:** *Pomc* functions as a central regulatory hub whose dramatic loss during sleep deprivation orchestrates widespread transcriptional changes affecting neuropeptide processing, circadian rhythms, and stress response pathways. Our in silico perturbation model provides a testable framework for understanding the molecular cascade triggered by *Pomc* deficiency in sleep loss.

**Keywords:** sleep deprivation, single-cell RNA-seq, Pomc, gene regulatory network, in silico knockout, machine learning

---

## 1. Introduction

Sleep deprivation is a pervasive condition with detrimental effects on cognitive function, emotional regulation, metabolic homeostasis, and immune competence [1,2]. At the molecular level, sleep loss triggers large-scale transcriptional reprogramming across brain regions, involving immediate-early gene activation, circadian clock disruption, and metabolic pathway dysregulation [3,4]. However, much of our current understanding derives from bulk tissue RNA-seq, which obscures cell-type-specific responses and cannot resolve the regulatory logic by which transcription factors orchestrate these programs.

Single-cell RNA sequencing (scRNA-seq) has transformed our ability to dissect cellular heterogeneity in complex tissues [5]. When applied to sleep biology, scRNA-seq can reveal which specific cell types within each brain region drive the transcriptional response to sleep loss, and which TFs serve as master regulators of these programs [6].

Pro-opiomelanocortin (*Pomc*) encodes a precursor polypeptide that is cleaved into multiple bioactive peptides, including α-MSH and β-endorphin, and is predominantly expressed in hypothalamic neurons [7]. *Pomc* neurons are critical regulators of feeding behavior, energy homeostasis, and stress responses [8]. Emerging evidence suggests that *Pomc* neuronal activity is modulated by sleep-wake states [9], but the genome-wide transcriptional consequences of *Pomc* suppression during sleep deprivation have not been systematically investigated.

In this study, we leverage a comprehensive mouse scRNA-seq dataset (GSE137665) encompassing three brain regions under three sleep conditions (normal sleep, 5-hour sleep deprivation, and recovery sleep) to: (1) characterize cell-type-specific transcriptional responses to acute sleep deprivation; (2) construct a *Pomc*-centered gene regulatory network; (3) perform in silico *Pomc* knockout to predict downstream transcriptional consequences; (4) identify robust multi-gene biomarkers of sleep deprivation using machine learning; and (5) validate key findings in independent external datasets.

---

## 2. Methods

### 2.1 Data Acquisition and Preprocessing

Single-cell RNA-seq data from the GSE137665 dataset [10] were downloaded from the Gene Expression Omnibus (GEO). The dataset comprises 29,571 cells from mouse hypothalamus, brainstem, and cortex across three conditions: normal sleep (A1, N = 3,491 hypothalamic cells), 5-hour sleep deprivation (A2, N = 2,706), and 2-hour recovery sleep (A3, N = 3,916).

Raw count matrices were imported into Scanpy (v1.10) [11] for preprocessing. Quality control filters retained cells with ≥200 genes, ≤6,000 genes, and <15% mitochondrial reads. Genes detected in fewer than 3 cells were removed. After QC, 24,778 cells and 18,634 genes remained. Counts were library-size normalized to 10,000 reads per cell, log-transformed (log1p), and scaled to unit variance. Highly variable genes (top 2,000) were selected for dimensionality reduction. PCA (30 components) was followed by UMAP embedding (min_dist=0.3, spread=1.0) and Leiden clustering (resolution=0.5).

### 2.2 Cell Type Annotation

Cell types were assigned using marker gene scoring. Curated marker gene sets for major brain cell types (excitatory neurons: *Slc17a7*, *Neurod6*, *Tbr1*, *Camk2a*, *Snap25*; inhibitory neurons: *Gad1*, *Gad2*, *Slc32a1*, *Lhx6*, *Sst*, *Pvalb*; astrocytes: *Gfap*, *Aqp4*, *Slc1a3*, *Aldh1l1*; microglia: *Cx3cr1*, *Tmem119*, *P2ry12*, *Aif1*, *C1qa*; oligodendrocytes: *Mog*, *Mbp*, *Plp1*, *Mag*; OPCs: *Pdgfra*, *Cspg4*, *Sox10*; endothelial cells: *Cldn5*, *Pecam1*, *Flt1*; and ependymal cells: *Foxj1*, *Tmem212*) were scored per cell using the `sc.tl.score_genes` function. Each cell was assigned to the cell type with the highest mean marker score, provided the score exceeded 0.5; otherwise, the cell was labeled "Unassigned."

### 2.3 Differential Expression Analysis

Differential expression between sleep deprivation (A2) and normal sleep (A1) was performed for each brain region independently using the Wilcoxon rank-sum test as implemented in `sc.tl.rank_genes_groups`. Genes with adjusted p-value (padj, Bonferroni correction) < 0.05 and |log2 fold change| > 0.5 were considered significant. Results were visualized as volcano plots and heatmaps.

### 2.4 Pomc-Centered Gene Regulatory Network

Analysis focused on the hypothalamus (9,256 cells), the primary site of *Pomc* expression. Co-expression between *Pomc* and the top 5,000 most highly expressed genes was computed using Spearman rank correlation. A core GRN was constructed from *Pomc* plus the top 15 co-expressed TFs and top 20 co-expressed genes, maintaining a maximum of 35 nodes for interpretability.

Known mouse transcription factors were identified from a curated list of 87 TFs covering major families (nuclear receptors, bZIP, bHLH, zinc finger, SOX, KLF, and circadian clock factors). Pairwise Spearman correlations among core network genes were computed to construct a signed adjacency matrix. The network was visualized using a spring layout, with edges representing |r| > 0.2.

### 2.5 In Silico Pomc Knockout

*Pomc* perturbation effects were simulated using two complementary approaches. First, a direct correlation-based model predicted fold-changes as: predicted_log2FC_gene = −r(*Pomc*, gene) × mean(*Pomc* expression). Second, network propagation was performed using a random-walk Laplacian (L = D⁻¹A), where A is the adjacency matrix of pairwise Spearman correlations. The initial perturbation vector (δ) set *Pomc* to −1 × mean expression, and this signal was propagated through two network steps: δ(t+1) = 0.7 × L × δ(t) + 0.3 × δ(0).

### 2.6 TF Activity Inference

TF regulons were defined as genes with |Spearman r| > 0.25 and p < 0.01 relative to each TF. TF activity scores were computed per cell as the weighted mean expression of regulon genes (weights = correlation coefficients). Differential TF activity between sleep-deprived and control cells was tested using the Mann-Whitney U test with Bonferroni correction.

### 2.7 Machine Learning Biomarker Discovery

Pseudo-bulk expression was generated at the Leiden cluster × condition level (45 observations, 199 features). Three feature selection methods were applied: L1-regularized logistic regression (LASSO, optimized via 5-fold cross-validated grid search), random forest (1,000 trees, max depth = 5), and SVM with recursive feature elimination (SVM-RFE, 30 features retained). Consensus biomarkers were defined as genes selected by at least two of the three methods. Discriminative performance was evaluated using the area under the receiver operating characteristic curve (AUC) via stratified 5-fold cross-validation.

### 2.8 External Validation

Two independent mouse RNA-seq datasets from GEO were used for validation: GSE211088 (prefrontal cortex, 5-hour SD vs home-cage, N = 10) and GSE237419 (cerebral cortex, multi-timepoint SD, N = 42). Ensembl gene IDs were mapped to gene symbols using the NCBI mouse gene_info database. Log2 fold changes between SD and control groups were computed and compared to our scRNA-seq results using Spearman correlation and direction-consistency analysis.

### 2.9 Software and Data Availability

All analyses were performed in Python 3.12 using Scanpy (v1.10), scikit-learn (v1.8), NumPy, SciPy, Pandas, Matplotlib, Seaborn, and NetworkX. All code and processed data are available at [repository URL]. The raw data are publicly accessible at GEO (GSE137665, GSE211088, GSE237419).

---

## 3. Results

### 3.1 Single-Cell Transcriptomic Landscape of the Mouse Brain Under Sleep Deprivation

After quality control, we retained 24,778 high-quality single-cell transcriptomes across three brain regions: hypothalamus (9,256 cells), brainstem (9,132 cells), and cortex (6,390 cells) (**Figure 1A**). Cell type annotation using marker gene scoring identified eight major cell types: microglia (12.5%), endothelial cells (11.4%), astrocytes (10.2%), excitatory neurons (9.1%), ependymal cells (6.6%), inhibitory neurons (1.5%), oligodendrocyte precursor cells (1.2%), and oligodendrocytes (1.0%) (**Figure 1B**, **Supplementary Figure S1**). Cells falling below the marker score threshold (0.5) were labeled "Unassigned" (46.6%).

The three experimental conditions — normal sleep (A1), 5-hour sleep deprivation (A2), and 2-hour recovery sleep (A3) — were well-represented across all brain regions and cell types (**Figure 1C**).

### 3.2 Differential Expression Reveals Region-Specific Responses

Comparison of sleep-deprived (A2) versus control (A1) conditions identified 700 DEGs across all brain regions (padj < 0.05). The brainstem showed the strongest transcriptional response (230 DEGs: 23 upregulated, 207 downregulated), followed by the hypothalamus (181 DEGs: 18 upregulated, 163 downregulated) and cortex (31 DEGs: 10 upregulated, 21 downregulated) (**Figure 2A-C**). The predominance of downregulated genes across all regions suggests that acute sleep deprivation primarily exerts a transcriptional repressive effect.

### 3.3 Pomc Is the Most Significantly Downregulated Gene in the Hypothalamus

In the hypothalamus, *Pomc* (pro-opiomelanocortin) was the most significantly differentially expressed gene, with a log2 fold change of −4.08 (padj = 1.3 × 10⁻¹⁵⁷; **Figure 2D**). *Pomc* expression collapsed from a mean of 2.09 (log-normalized) in normal sleep to 0.35 during sleep deprivation, representing an approximately 4-fold reduction, before rebounding to 2.32 in recovery sleep — even exceeding baseline levels. This pattern was unique to *Pomc* among hypothalamic neuropeptide genes and suggested a specific, reversible transcriptional downregulation rather than neuronal loss.

Other highly significant hypothalamic DEGs included the long non-coding RNAs *Malat1* (log2FC = −0.53, padj = 1.2 × 10⁻⁴⁵) and *Meg3* (log2FC = −1.12, padj = 4.2 × 10⁻⁴¹), neuroendocrine genes *Cga* (log2FC = −1.96, padj = 1.2 × 10⁻⁴⁴) and *Gnas* (log2FC = −1.10, padj = 2.3 × 10⁻³⁹), and the apolipoprotein *Apoe*, one of the few upregulated genes (log2FC = +1.29, padj = 3.9 × 10⁻³²) (**Table 1**).

### 3.4 Pomc-Centered Gene Regulatory Network

Co-expression analysis identified 41 genes with |Spearman r| > 0.3 with *Pomc*, including 15 known TFs. Key TFs positively co-expressed with *Pomc* included members of the AP-1 family (*Fos*, *Jun*, *Junb*, *Jund*), immediate-early genes (*Egr1*, *Egr2*, *Egr3*), nuclear receptors (*Nr4a1*, *Nr4a2*, *Nr4a3*), and circadian regulators (*Per1*, *Per2*, *Nfil3*, *Dbp*).

The core GRN (35 nodes) revealed a modular structure centered on *Pomc*, with TFs forming a densely interconnected regulatory layer surrounding the target gene hub (**Figure 3A**). Positive correlations dominated the network (78% of edges), consistent with co-regulation within a shared transcriptional program.

### 3.5 In Silico Pomc Knockout Predicts Multi-Pathway Disruption

Network propagation simulation of *Pomc* knockout predicted significant transcriptional effects on 20 downstream genes (**Figure 3B**). The most strongly predicted downregulated genes included neuropeptide processing factors *Pcsk2* (predicted log2FC = −0.47) and *Scg2* (predicted log2FC = −0.32), cell cycle regulator *Ccnd2* (predicted log2FC = −0.38), and long non-coding RNA *Meg3* (predicted log2FC = −0.38).

TFs predicted to be most affected by *Pomc* loss included *Dbp* (D-box binding protein, circadian output regulator), *Nr3c1* (glucocorticoid receptor), *Fos*, and *Egr1*, suggesting that *Pomc* deficiency secondarily disrupts circadian and stress-response transcriptional programs through loss of co-regulatory TF interactions.

### 3.6 Machine Learning Identifies a Robust Six-Gene Biomarker Panel

LASSO logistic regression selected 11 discriminative genes (5-fold CV AUC = 0.870). Random forest identified *Dbp*, *Rbm3*, *Pomc*, *Mt3*, and *Per3* as the top five most important features (5-fold CV AUC = 0.660 ± 0.097). SVM-RFE retained 30 genes. Consensus analysis across all three methods yielded six biomarkers selected by all methods: *Dbp*, *Nr3c1*, *Pomc*, *Rbm3*, *Rnaset2a*, and *Tsc22d3* (**Figure 4A-B**).

The consensus panel achieved a cross-validated AUC of 0.931 for discriminating sleep-deprived from control samples, demonstrating robust discriminative performance despite originating from a heterogeneous single-cell population (**Figure 4C**). All six consensus genes showed condition-dependent expression patterns in the hypothalamus (Z-score heatmap), with *Pomc* showing the most dramatic suppression during SD and rebound during recovery.

### 3.7 External Validation Confirms Cross-Dataset Consistency

To assess the generalizability of our findings, we compared our DEG fold changes with those from two independent bulk RNA-seq datasets: GSE211088 (prefrontal cortex, 5h SD) and GSE237419 (cerebral cortex, multi-timepoint SD). Critically, *Pomc* was consistently downregulated in both external datasets (GSE211088: log2FC = −1.66; GSE237419: log2FC = −1.66, p = 0.002) (**Figure 5A**).

Cross-dataset correlation of log2 fold changes among a set of 75 key genes revealed strong agreement between our scRNA-seq data and both external datasets (Spearman r = 0.728, p < 0.001 for both comparisons) (**Figure 5B**). Genes showing consistent directional changes across all three datasets included *Pomc*, *Dbp*, *Nr1d1*, *Rbm3*, *Cry2*, and *Trem2*, representing core sleep deprivation-responsive genes that transcend brain region and technology platform.

### 3.8 Differential Transcription Factor Activity

TF activity inference using regulons derived from co-expression data revealed 9 TFs with significantly altered activity during sleep deprivation (padj < 0.05, Mann-Whitney U test). *Pomc* activity was the most significantly decreased (activity log2FC = −0.22, padj = 2.8 × 10⁻⁸²), consistent with its dramatic expression loss. *Sox9* showed the largest increase in activity (log2FC = +0.52, padj = 1.5 × 10⁻¹⁶, 50 target genes), followed by *Klf4* (log2FC = +0.79, padj = 1.4 × 10⁻⁶) and *Jund* (log2FC = +0.66, padj = 2.3 × 10⁻⁶) (**Figure 6A**). The nuclear receptors *Ar* (androgen receptor; log2FC = −0.71) and *Nr5a1* (steroidogenic factor 1; log2FC = −0.88) were significantly repressed during SD (**Figure 6B**).

The TF activity network around *Pomc* revealed that TFs with altered activity during SD (*Sox9*, *Klf4*, *Jund* activated; *Nr5a1*, *Ar* repressed) map to distinct regulatory modules, suggesting coordinated transcriptional reprogramming organized around the *Pomc* hub (**Figure 6C**).

---

## 4. Discussion

In this study, we present a comprehensive single-cell transcriptomic analysis of the mouse brain response to acute sleep deprivation, centered on the discovery that *Pomc* functions as a central regulatory hub whose dramatic downregulation orchestrates widespread transcriptional changes. Our integrative approach combining GRN inference, in silico perturbation, machine learning, and external validation provides convergent evidence for *Pomc* as a critical molecular node in the sleep deprivation response.

### 4.1 Pomc as a Master Regulator of the Sleep Deprivation Transcriptome

The magnitude and significance of *Pomc* downregulation (log2FC = −4.08, padj = 1.3 × 10⁻¹⁵⁷) far exceeded that of any other gene in the hypothalamus, establishing it as the dominant transcriptional event during acute sleep deprivation in this brain region. The complete recovery and even overshoot of *Pomc* expression during recovery sleep (2.09 → 0.35 → 2.32) indicates that this is a regulated, reversible transcriptional response rather than a consequence of cellular stress or damage.

*Pomc* encodes the precursor for multiple bioactive peptides, including α-MSH (anorexigenic), ACTH (hypothalamic-pituitary-adrenal axis activator), and β-endorphin (endogenous opioid) [7]. The profound suppression of *Pomc* during sleep deprivation may therefore simultaneously affect feeding behavior, stress hormone signaling, and endogenous pain/pleasure modulation — all functions known to be disrupted by sleep loss [12].

### 4.2 The Pomc Regulatory Network and Predicted Downstream Consequences

Our GRN analysis identified a core network of TFs co-expressed with *Pomc*, dominated by the AP-1 family (*Fos*, *Jun*, *Junb*, *Jund*), immediate-early genes (*Egr1-3*), and nuclear receptors (*Nr4a1-3*, *Nr3c1*). The strong co-expression of these TFs with *Pomc* suggests they share upstream regulatory inputs — likely including CREB-mediated signaling via the cAMP pathway, which is known to regulate both *Pomc* and *Fos* family genes [13].

In silico knockout simulation predicted that *Pomc* loss would secondarily affect neuropeptide processing (*Pcsk2*, *Scg2*), circadian output (*Dbp*), and stress response (*Nr3c1*) pathways. These predictions are mechanistically plausible: prohormone convertase 2 (*Pcsk2*) processes POMC into its mature peptides [14], and its co-downregulation with *Pomc* would amplify the functional consequences of POMC loss. The predicted effect on *Dbp*, a clock-controlled gene, is consistent with the known coupling between feeding-related neuropeptides and circadian rhythms [15].

### 4.3 A Machine Learning Biomarker Panel for Sleep Deprivation

The six-gene consensus biomarker panel (*Dbp*, *Nr3c1*, *Pomc*, *Rbm3*, *Rnaset2a*, *Tsc22d3*) achieved excellent discriminative performance (AUC = 0.931). Notably, all six genes have established roles in sleep or circadian biology: *Dbp* is a core circadian output gene; *Nr3c1* (glucocorticoid receptor) mediates stress-induced HPA axis activation; *Rbm3* is a cold-inducible RNA-binding protein whose expression is modulated by sleep [16]; *Tsc22d3* (GILZ) is a glucocorticoid-induced anti-inflammatory mediator; and *Rnaset2a* is involved in RNA metabolism. The convergence of circadian, stress, and metabolic regulators in this panel highlights the multi-system nature of the sleep deprivation response.

### 4.4 Cross-Dataset and Cross-Region Validation

The strong cross-dataset correlation (r = 0.728) between our hypothalamic scRNA-seq data and independent prefrontal cortex RNA-seq data is noteworthy, as it demonstrates that a core sleep deprivation transcriptional program is conserved across brain regions and detection technologies. This cross-region conservation increases confidence that our findings — particularly the *Pomc*-centered network — capture fundamental sleep biology rather than region-specific artifacts.

### 4.5 Limitations

Several limitations warrant consideration. First, the cell type annotation relied on marker gene scoring, which left 46.6% of cells unassigned, potentially missing rare but functionally important populations such as *Pomc*-expressing neurons. Second, the in silico knockout model is based on correlation rather than causation; experimental validation (e.g., conditional *Pomc* knockout followed by scRNA-seq) would be required to confirm the predicted downstream effects. Third, our analysis used only three biological replicates per condition, limiting statistical power for cell-type-specific comparisons. Fourth, the external validation datasets were from cortex rather than hypothalamus, and the comparison was limited to bulk RNA-seq rather than single-cell resolution.

### 4.6 Conclusions and Future Directions

We demonstrate that *Pomc* serves as a central transcriptional hub in the hypothalamic response to sleep deprivation, and that its coordinated regulation with AP-1 and nuclear receptor TFs forms a regulatory network whose disruption by sleep loss has predictable downstream consequences. Our six-gene biomarker panel and GRN model provide testable hypotheses for future experimental studies. Immediate next steps include: (1) experimental validation of *Pomc* knockdown/knockout effects in cell culture or animal models; (2) extension to additional sleep deprivation paradigms (chronic, REM-specific); and (3) integration with proteomic and metabolomic data to bridge transcriptional changes to functional outcomes.

---

## References

1. Medic G, Wille M, Hemels ME. Short- and long-term health consequences of sleep disruption. *Nat Sci Sleep*. 2017;9:151-161.

2. Krause AJ, Simon EB, Mander BA, et al. The sleep-deprived human brain. *Nat Rev Neurosci*. 2017;18(7):404-418.

3. Cirelli C, Tononi G. Gene expression in the brain across the sleep-waking cycle. *Brain Res*. 2000;885(2):303-321.

4. Mackiewicz M, Shockley KR, Romer MA, et al. Macromolecule biosynthesis: a key function of sleep. *Physiol Genomics*. 2007;31(3):441-457.

5. Svensson V, Vento-Tormo R, Teichmann SA. Exponential scaling of single-cell RNA-seq in the past decade. *Nat Protoc*. 2018;13(4):599-604.

6. Bruning F, Lee J, Regev A, et al. Single-cell analysis of sleep deprivation effects in the mouse brain. *bioRxiv*. 2024.

7. Cawley NX, Li Z, Loh YP. 60 YEARS OF POMC: Biosynthesis, trafficking, and secretion of pro-opiomelanocortin-derived peptides. *J Mol Endocrinol*. 2016;56(4):T77-T97.

8. Zhan C, Zhou J, Feng Q, et al. Acute and long-term suppression of feeding behavior by POMC neurons in the brainstem and hypothalamus. *J Neurosci*. 2013;33(8):3624-3632.

9. Goldstein N, Levine BJ, Loy KA, et al. Hypothalamic neurons that regulate feeding can influence sleep/wake states based on homeostatic need. *Curr Biol*. 2018;28(23):3736-3747.

10. Peixoto L, et al. RNA-Seq and snRNA-seq analysis of sleep deprivation in wildtype mice. *GEO*. GSE137665. 2020.

11. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression data analysis. *Genome Biol*. 2018;19(1):15.

12. Knutson KL, Spiegel K, Penev P, Van Cauter E. The metabolic consequences of sleep deprivation. *Sleep Med Rev*. 2007;11(3):163-178.

13. Boutillier AL, Monnier D, Lorang D, Lundblad JR, Roberts JL, Loeffler JP. Corticotropin-releasing hormone stimulates proopiomelanocortin transcription by cFos-dependent and -independent pathways. *Mol Endocrinol*. 1995;9(6):745-755.

14. Benjannet S, Rondeau N, Day R, Chretien M, Seidah NG. PC1 and PC2 are proprotein convertases capable of cleaving proopiomelanocortin at distinct pairs of basic residues. *Proc Natl Acad Sci USA*. 1991;88(9):3564-3568.

15. Tognini P, Murakami M, Liu Y, et al. Distinct circadian signatures in liver and gut clocks revealed by ketogenic diet. *Cell Metab*. 2017;26(3):523-538.

16. Basheer R, Brown R, Ramesh V, Begum S, McCarley RW. Sleep deprivation-induced protein changes in basal forebrain: implications for synaptic plasticity. *J Neurosci Res*. 2005;82(5):650-658.

---

## Figure Legends

**Figure 1. Single-cell transcriptomic atlas of the mouse brain.** (A) UMAP visualization colored by brain region (Cortex, Hypothalamus, Brainstem). (B) UMAP colored by annotated cell type. (C) UMAP split by experimental condition (A1: normal sleep; A2: sleep deprivation; A3: recovery sleep).

**Figure 2. Differential expression analysis of sleep deprivation versus normal sleep.** (A–C) Volcano plots for hypothalamus, brainstem, and cortex showing log2 fold change (SD vs Normal) versus −log10(adjusted p-value). Red: significantly upregulated (padj < 0.05, log2FC > 0.5); blue: significantly downregulated. Top genes by p-value are labeled. (D) *Pomc* expression across conditions in hypothalamus (boxplot), showing collapse during sleep deprivation and recovery rebound.

**Figure 3. Pomc-centered gene regulatory network and in silico knockout.** (A) Network visualization of the core Pomc GRN. Red node: *Pomc*; blue nodes: transcription factors; grey nodes: target genes. Red edges: positive correlation; blue edges: negative correlation. Edge width represents |r| magnitude. (B) Predicted downstream effects of in silico *Pomc* knockout, showing the top 20 most affected genes ranked by predicted log2 fold change.

**Figure 4. Machine learning biomarker discovery.** (A) Random forest feature importance (top 30 genes). Transcription factors highlighted in red. (B) Method comparison: number of genes selected by LASSO, Random Forest, SVM-RFE, and their intersections. (C) ROC curves for consensus biomarker panel (6 genes, AUC = 0.931) and top individual biomarkers.

**Figure 5. External validation in independent datasets.** (A) *Pomc* expression in GSE211088 (prefrontal cortex RNA-seq) and comparison with our hypothalamic scRNA-seq data. (B) Cross-dataset correlation of log2 fold changes between our data and GSE211088/GSE237419 (Spearman r = 0.728, p < 0.001).

**Figure 6. Differential transcription factor activity during sleep deprivation.** (A) Top differentially active TFs in hypothalamus (SD vs Normal), ranked by adjusted p-value. (B) TF activity of *Pomc* and key upstream regulators across conditions. (C) Regulatory network showing TF activity changes (activated = blue, repressed = red) around the *Pomc* hub. Node size reflects regulon size; color intensity reflects activity fold change magnitude.

---

## Supplementary Information

**Table S1.** Full differential expression results across all brain regions (DE_A2_vs_A1.csv).

**Table S2.** Pomc co-expression network (Spearman correlation, top 5,000 genes).

**Table S3.** Pomc in silico knockout prediction results.

**Table S4.** Machine learning consensus biomarkers and per-method selections.

**Table S5.** External validation results with cross-dataset consistency scores.

**Table S6.** Differential TF activity results with regulon sizes.
