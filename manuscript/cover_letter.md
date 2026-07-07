# Cover Letter

**Date:** 2026-05-18

**To:** The Editor
Frontiers in Neuroscience
EPFL Innovation Park, Building I, 1015 Lausanne, Switzerland

**Re:** Submission of Original Research Article: "Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central Regulatory Hub and Predicts Downstream Transcriptional Consequences of Its Loss"

---

Dear Editor,

We are pleased to submit our manuscript entitled **"Single-Cell Transcriptomic Analysis of Sleep Deprivation Reveals Pomc as a Central Regulatory Hub and Predicts Downstream Transcriptional Consequences of Its Loss"** for consideration for publication in *Frontiers in Neuroscience*.

## Significance and Novelty

Sleep deprivation affects millions of people worldwide and is associated with cognitive impairment, metabolic dysfunction, and increased disease risk. Despite decades of research, the cell-type-specific transcriptional programs and gene regulatory hierarchies underlying the brain's response to sleep loss remain poorly characterized.

Our study makes the following **novel contributions**:

1. **Discovery of *Pomc* as a central transcriptional hub**: We identify *Pomc* as the most significantly altered gene in the hypothalamus during acute sleep deprivation (log2FC = −4.08, padj = 1.3 × 10⁻¹⁵⁷), with expression collapsing by ~80% and fully recovering upon sleep — a pattern that implicates it as a specific, regulated transcriptional switch.

2. **First in silico *Pomc* knockout in sleep biology**: Using network propagation through a *Pomc*-centered gene regulatory network, we predict downstream consequences of *Pomc* loss on neuropeptide processing, circadian output, and stress-response pathways. This provides a computational framework for hypothesis generation that can guide future experimental validation.

3. **Robust multi-gene biomarker panel**: Machine learning (LASSO, Random Forest, SVM-RFE consensus) identified a six-gene panel for classifying sleep-deprived versus control states with an AUC of 0.931, integrating circadian (*Dbp*), stress (*Nr3c1*, *Tsc22d3*), and metabolic (*Pomc*, *Rbm3*) regulators.

4. **Cross-dataset, cross-technology validation**: Independent validation in two external RNA-seq datasets (GSE211088, GSE237419) confirmed consistent directional changes for key sleep deprivation-responsive genes (Spearman r = 0.728, p < 0.001), demonstrating the generalizability of our findings beyond the primary dataset.

5. **Comprehensive analysis pipeline**: Our study integrates single-cell differential expression, co-expression network analysis, transcription factor activity inference, machine learning, and external validation — providing a reproducible computational framework applicable to other sleep and circadian transcriptomic studies.

## Data and Code Availability

All raw data analyzed in this study are publicly available from the Gene Expression Omnibus (GSE137665, GSE211088, GSE237419). All analysis code, processed data files, and figure-generation scripts are available at [repository URL to be provided upon acceptance]. The complete analysis pipeline is implemented in Python 3.12 using open-source libraries (Scanpy, scikit-learn, NetworkX, etc.), ensuring full reproducibility.

## Suggested Reviewers

1. Dr. Chiara Cirelli — University of Wisconsin-Madison (sleep and transcriptomics)
2. Dr. Sara J. Aton — University of Michigan (sleep and synaptic plasticity)
3. Dr. Lucia Peixoto — Washington State University (sleep deprivation genomics)
4. Dr. Lior Pachter — California Institute of Technology (single-cell computational methods)
5. Dr. Paul Franken — University of Lausanne (sleep genetics)

## Competing Interests

The authors declare no competing interests.

## Author Contributions

Yang Yongxin: Conceptualization, methodology, software, formal analysis, writing — original draft, writing — review and editing.

We confirm that this manuscript has not been published elsewhere and is not under consideration by another journal. All authors have approved the manuscript and agree with its submission to *Frontiers in Neuroscience*.

We thank you and the reviewers for your time and consideration.

Sincerely,

Yang Yongxin
Fujian Second People's Hospital
960856791@qq.com

---

### Suggested Target Journals (Ranked by Fit)

| Rank | Journal | IF (2024) | Rationale |
|------|---------|-----------|-----------|
| 1 | **Sleep** | ~5.0 | Top sleep journal; perfect thematic fit; single-cell + computational acceptable |
| 2 | **Frontiers in Neuroscience** | ~3.5 | Open access; accepts computational studies; relatively fast review |
| 3 | **Journal of Sleep Research** | ~3.5 | Strong sleep focus; European sleep research society |
| 4 | **Scientific Reports** | ~4.0 | Broad scope; accepts bioinformatics; high acceptance rate |
| 5 | **BMC Genomics** | ~3.5 | Good for computational genomics; open access |
| 6 | **Neurobiology of Sleep and Circadian Rhythms** | ~2.5 | Newer journal; very specific fit |

**Recommendation for a first submission**: Start with **Frontiers in Neuroscience** — good balance of impact factor, scope fit, review speed, and acceptance rate for computational work.
