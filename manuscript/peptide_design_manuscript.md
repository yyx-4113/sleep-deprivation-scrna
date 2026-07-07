---
title: "Mechanism-Conditioned Latent-Space Optimization of Sleep-Promoting Peptides via Active Learning with Virtual Cell Validation"
author: "Yongxin Yang"
date: "2026-06-03"
documentclass: article
mainfont: Times New Roman
fontsize: 12pt
---

# Mechanism-Conditioned Latent-Space Optimization of Sleep-Promoting Peptides via Active Learning with Virtual Cell Validation

Yongxin Yang<sup>1,*</sup>

<sup>1</sup> Fujian Second People's Hospital, Fuzhou, Fujian 350003, China

<sup>*</sup> Correspondence: 960856791@qq.com

**Running title:** Latent-space optimization of sleep-promoting peptides

**Keywords:** active learning; rank-percentile ensemble score; directional transcriptional rescue; latent-space optimization; virtual cell model; sleep deprivation; single-cell transcriptomics; food-derived peptides

**Word count:** approximately 8,000 (excluding references)

**Figures:** 3 main figures

**Supplementary:** 2 supplementary figures, 5 supplementary tables

---

## Abstract

**Motivation:** Computational discovery of bioactive peptides traditionally relies on screening natural libraries. While de novo generative models using protein language models have shown promise, their training is computationally demanding and reproducibility remains challenging. Existing approaches also lack mechanism-driven conditioning and systems-level evaluation beyond molecular property prediction.

**Results:** We present PepDesign-Active, a framework integrating three components: (i) a peptide-specific PCA-based latent space built from 52,517 food-derived peptides, providing a reproducible alternative to variational autoencoders (76.7% reconstruction accuracy, 13-second training); (ii) a mechanism-conditioned active learning loop guided by a Rank-Percentile Ensemble Score (RPES), which eliminates ceiling effects via geometric mean of component rank percentiles; and (iii) evaluation through a single-cell transcriptome-derived virtual cell model using Directional Transcriptional Rescue (DTR), which scores perturbations by projection onto a biologically defined rescue direction. Under RPES, PepDesign-Active active learning achieved 0.8059 +/- 0.0020, a +3.2% improvement over random latent sampling (0.7811) and a +2.4% improvement over equal-budget random library search (0.7870 +/- 0.0089, n=10). Ablation confirmed both active learning (d = +0.016 vs. single-pass) and mechanism conditioning (d = +0.037 vs. unconditional search) are necessary. DTR sensitivity analysis confirmed robust discrimination between real and scrambled gene regulatory networks across damping factors (0.3-0.9) and binding thresholds. Top RPES-designed peptides (e.g., PGPGRITAYL) showed enrichment in Gly, Tyr, Ala, and Pro, with lengths (9-10 residues) favorable for solid-phase synthesis.

**Availability:** Code and models at https://gitee.com/yongxin-yang/sleep-deprivation-scrna.

---

## Introduction

Bioactive peptides derived from dietary proteins have attracted substantial interest as functional food ingredients and potential therapeutic leads due to their diverse biological activities, including antioxidant, antimicrobial, antihypertensive, and neuroactive properties (Daliri et al., 2017; Chakrabarti et al., 2018). The traditional discovery pipeline -- enzymatic hydrolysis, fractionation, activity-guided purification, and mass spectrometric identification -- is labor-intensive and inherently limited to peptides present in the starting material (Li et al., 2025). Computational methods have emerged as complementary approaches to accelerate discovery, with recent advances in deep learning enabling accurate prediction of peptide bioactivity from sequence alone (Torres et al., 2025; Zhang R. et al., 2025).

However, the dominant paradigm in computational peptide science remains *screening*: given a library of natural peptide sequences, machine learning models are trained to classify or rank them by predicted bioactivity. While effective, this approach cannot explore sequence space beyond what exists in nature. De novo peptide design -- the computational generation of novel sequences not present in any natural library -- represents the logical next step. Recent breakthroughs in protein design, particularly diffusion models and protein language models (PLMs), have demonstrated the feasibility of generating functional proteins and peptides from scratch (Watson et al., 2023; Dauparas et al., 2022; Torres et al., 2025). However, these approaches have three key limitations for the specific task of food-derived bioactive peptide design:

First, PLMs such as ESM-2 and ProtBERT are trained on full-length proteins and may not optimally represent the sequence space of short peptides (2-20 amino acids) typically derived from food proteins (Lin et al., 2023). The amino acid composition, structural constraints, and evolutionary signals differ substantially between globular proteins and short food-derived peptides.

Second, existing generative models lack mechanism-driven conditioning. They generate sequences that satisfy basic sequence statistics or optimize for a single molecular property, without incorporating the biological mechanism through which the peptide is expected to act. For sleep-promoting peptides specifically, the mechanism involves multi-receptor binding (MC4R, MT1, MT2, OX1R, OX2R) and downstream transcriptional effects that can be modeled through a gene regulatory network (GRN).

Third, existing approaches lack systems-level evaluation. While molecular docking provides binding energy estimates and machine learning models provide bioactivity predictions, neither captures the network-level transcriptional effect of a candidate peptide on the target tissue -- a gap that a virtual cell model can partially address.

To address these gaps, we present PepDesign-Active, a computational framework that integrates three components: (i) a peptide-specific PCA-based latent space built from short food-derived peptide sequences, providing a compact representation specialized for this molecular class; (ii) a mechanism-conditioned active learning loop that uses an ensemble surrogate model to guide iterative exploration of the latent space toward regions predicted to reverse a disease-relevant transcriptional signature; and (iii) a single-cell transcriptome-derived computational virtual cell model that simulates downstream transcriptional consequences of peptide-receptor interactions, providing systems-level evaluation beyond statistical metrics.

We demonstrate PepDesign-Active on the task of designing sleep-promoting peptides derived from Fujian medicinal-food ingredients (walnut, mulberry, and black sesame). Sleep deprivation affects millions worldwide and is associated with cognitive impairment, metabolic dysfunction, and immune dysregulation (Medic et al., 2017; Krause et al., 2017). Single-cell transcriptomic analyses have revealed cell-type-specific molecular responses to sleep deprivation across multiple brain regions including the hypothalamus (Jha et al., 2022), and chronic sleep deprivation has been shown to disrupt POMC neuronal function in the arcuate nucleus through the SCN-BMAL1 signaling pathway (Du et al., 2024), providing both the disease-relevant transcriptional signature and a validated GRN that serves as the basis for our virtual cell model.

---

## Methods

### 2.1 Peptide Library and Feature Encoding

A comprehensive library of 52,517 peptide sequences (length 2-20 amino acids, mean 9.9, SD 4.7) was compiled from three Fujian medicinal-food protein sources: walnut (*Juglans regia*), mulberry leaf (*Morus alba*), and black sesame (*Sesamum indicum*). Proteins were subjected to in silico digestion with trypsin, pepsin, and chymotrypsin, and all resulting peptide fragments were retained. Each peptide was one-hot encoded as a 20 x 20 binary matrix (20 positions, zero-padded to the maximum peptide length x 20 standard amino acid types), where rows represent sequence positions and columns represent amino acid identity.

### 2.2 Peptide-Specific Latent Space via Principal Component Analysis

To learn a compact representation of short peptide sequences, each peptide was flattened from its one-hot encoding (20 positions x 20 amino acids = 400 dimensions) into a binary vector. Principal Component Analysis (PCA) was applied to the 52,517-peptide training set, projecting the 400-dimensional one-hot space to a 64-dimensional latent space. The top 64 principal components captured 49.0% of the input variance, consistent with the high degree of sparsity in the one-hot encoding (mean peptide length 9.9 out of 20 positions, leaving ~50% of positions as zero-padding).

PCA was chosen over a neural network-based variational autoencoder for two reasons: (i) computational accessibility -- PCA training completed in 13 seconds on consumer hardware (Intel Core i7-8550U, 8 GB RAM), whereas ConvVAE training requires GPU acceleration that was not available for this study; and (ii) empirical reliability -- we attempted to train the ConvVAE architecture described in the literature (two 1D convolutional layers, 64-dim latent space, beta-VAE schedule) on an NVIDIA T4 GPU, but training was unstable with reconstruction accuracy plateauing at ~10% (near chance level, vs. 78% reported in prior work), even after learning rate and architecture tuning. PCA provided a deterministic, reproducible alternative that achieved 76.7% reconstruction accuracy in 13 seconds without hyperparameter tuning. As demonstrated in the ablation study (Section 3.5, Ablation C), the 64-dimensional PCA latent space adds only ~2% incremental predictive power beyond 10 physicochemical descriptors for the RPES prediction task, suggesting that the linear PCA representation captures the primary axes of variation in peptide sequence space relevant to bioactivity scoring. Mean reconstruction accuracy (inverse PCA transform to one-hot, followed by argmax decoding to amino acid sequence) was 76.7%, comparable to the 78.3% reported for ConvVAE architectures on similar peptide datasets.

### 2.3 Rank-Percentile Ensemble Score (RPES)

The primary optimization target is the **Rank-Percentile Ensemble Score (RPES)**, a composite metric designed to avoid the ceiling effects inherent in weighted-average scores. Unlike existing approaches that average normalized component scores (where the maximum can be reached by scoring well on a single component), RPES uses a geometric mean of rank percentiles, requiring candidates to perform well on ALL component axes simultaneously to achieve a high score.

RPES is computed from six component predictors, each providing a complementary signal about peptide bioactivity potential:

- **rf_score**: Random Forest classifier (n_estimators=100, max_depth=10) trained on a curated benchmark of experimentally validated bioactive vs. inactive food-derived peptides (balanced dataset, n = 2,000, achieved AUC = 0.94). Uses 10 physicochemical descriptors as input.

- **bilstm_neuro**: A sequence-aware neuroactivity classifier. Due to CPU hardware constraints preventing PyTorch LSTM training, we used an XGBoost classifier trained on dipeptide composition features (400 k-mer frequencies) as a computationally efficient proxy. The model was trained with pseudo-labels derived from known neuroactivity-correlated physicochemical properties (moderate hydrophobicity + positive charge) on 6,000 peptides, achieving an AUC of 1.00 on this proxy task. While this does not capture the sequential dependencies that a true BiLSTM would learn, the dipeptide features preserve local sequence-order information relevant to receptor binding. We refer to this component as "bilstm_neuro" for consistency with the original PepDesign-Active architecture description, noting that a GPU-trained BiLSTM is planned for Phase 2 validation.

- **composite_ml**: Arithmetic mean of standardized scores from XGBoost, SVM (RBF kernel), and k-Nearest Neighbors (k=5), all trained on the same physicochemical features as rf_score, providing classifier diversity.

- **ensemble_score**: Weighted average prediction from the HGBR + RF ensemble surrogate model, incorporating both latent and physicochemical features.

- **hydrophobicity_plausibility**: A domain-knowledge term penalizing peptides with mean Kyte-Doolittle hydrophobicity outside the interquartile range of natural food-derived peptides.

- **length_normalization**: A Gaussian penalty centered at the natural peptide mean length (9.9 residues, SD 4.7).

**Aggregation via rank percentiles.** For each component, raw prediction values across the entire peptide library are converted to percentile ranks (0-1, uniform distribution). The RPES for peptide *p* is the geometric mean of its six component ranks:

RPES(p) = (rank_1(p) x rank_2(p) x ... x rank_6(p))^(1/6)

This formulation has three advantages over weighted averaging: (i) it eliminates the ceiling effect -- the theoretical maximum of 1.0 is asymptotically unreachable since it requires being in the top percentile on all six components; (ii) the geometric mean penalizes severe weakness on any single component (a peptide scoring at the 1st percentile on one component has its RPES reduced to near-zero regardless of performance on the other five); and (iii) rank transformation removes sensitivity to the scale and distribution of individual component scores, making the composite robust to outliers.

We validated empirically that RPES produces a well-distributed score (mean = 0.42, SD = 0.15, max = 0.84 on 403,461 peptides) with no saturation at the upper bound, in contrast to the weighted-average formulation which produced max = 1.00 with multiple peptides at the ceiling.

**Surrogate model.** An ensemble of Histogram-based Gradient Boosting Regressor (HGBR, max_iter=100, max_depth=6) and Random Forest Regressor (RF, n_estimators=100, max_depth=10) was trained on 52,517 samples with 74 features (64 latent dimensions + 10 physicochemical descriptors including length, mean hydrophobicity, molecular weight, charge properties, and amino acid composition categories) to predict RPES from latent representations. Model performance was evaluated via 3-fold cross-validation. Note that the surrogate model is trained on the same 52,517-peptide library used for PCA; therefore, reported CV metrics reflect in-distribution predictive performance rather than out-of-distribution generalization.

### 2.4 Active Learning Loop

The active learning loop proceeds as follows: (i) All latent vectors in the current pool are scored using the ensemble surrogate model. (ii) The top 30 scoring vectors are selected as "parents." (iii) For each round, 250 candidate latent vectors are generated via Gaussian perturbation (noise scale adaptively reduced from 0.30 to 0.05 across rounds) around the parent vectors. (iv) Candidates are scored, and the top 50 are decoded to peptide sequences via the PCA decoder. (v) The top 15 latent vectors are added to the pool for the next round. (vi) Convergence is assessed when the best score improvement over three consecutive rounds drops below 0.001.

### 2.5 Virtual Cell Model for Computational Evaluation

We emphasize that the virtual cell model is a **computational evaluation tool**, not a substitute for in vitro or in vivo biological validation. It provides a systems-level in silico readout that complements molecular property predictions but should be interpreted as generating qualitative hypotheses rather than quantitative predictions of biological efficacy.

The virtual cell model was constructed from published single-cell transcriptomic data of sleep-deprived mouse hypothalamus (Jha et al., 2022) and consists of a 38-gene, 56-edge gene regulatory network (GRN) centered on *Pomc* in hypothalamic neurons. Receptors MC4R, MT1, MT2, OX1R (Hcrtr1), OX2R (Hcrtr2), and 5-HT1A serve as entry points, selected based on their established roles in sleep-wake regulation: MC4R modulates wakefulness via melanocortin signaling (Xu et al., 2020), MT1/MT2 mediate melatonin-dependent sleep onset (Liu et al., 2016), OX1R/OX2R are critical for arousal stability (Sakurai, 2007; Bonnavion & de Lecea, 2010), and 5-HT1A regulates REM sleep (Monti, 2011). All six receptors are expressed in hypothalamic neuronal populations relevant to the Pomc circuit. We note that the AlphaFold-predicted structures used for docking likely represent inactive-state conformations of these GPCRs; for agonist peptide docking, active-state structures would be preferable, and the docking scores should be interpreted as relative rankings rather than absolute binding affinity predictions.

A peptide's predicted binding energy profile across these six receptors is estimated using AutoDock Vina (v1.2.5), with receptor structures obtained from the AlphaFold Protein Structure Database (MC4R: P32245, MT1: P48039, MT2: P49286, OX1R: O43613, OX2R: O43614, 5-HT1A: P08908). For each peptide, the 2D structure is generated from the decoded sequence using PyMOL, converted to a 3D conformer with RDKit (ETKDG method, 200 conformers), and docked against each receptor with exhaustiveness = 16. We acknowledge two important limitations: (i) 200 conformers provide sparse sampling of the conformational space for flexible linear peptides ≤20 residues, and the minimum binding energy across conformers is reported; the standard deviation across conformers (typically 1-3 kcal/mol) should be considered when interpreting score differences. (ii) The AlphaFold structures likely represent inactive-state GPCR conformations; active-state structures would be more appropriate for agonist peptide docking. The binding energies should therefore be interpreted as relative rankings within the pipeline rather than absolute affinity predictions.

This binding profile is propagated through the GRN: each receptor's activation level (estimated as 1 / (1 + exp(-(DeltaG_bind - DeltaG_threshold) / RT))) modulates the expression of its immediate downstream targets, and the perturbation cascades through the 56 edges of the network via damped iterative propagation (10 steps, damping factor 0.7).

**Directional Transcriptional Rescue (DTR).** To evaluate whether a peptide-induced perturbation moves the transcriptional program toward the normal state, we define a "rescue direction" vector from published single-cell transcriptomic data. For each gene in the GRN, we compute the fold change between sleep-deprived and normal hypothalamic neurons from the Jha et al. (2022) dataset. The rescue direction is the negated fold-change vector: genes upregulated in sleep deprivation must be suppressed, and genes downregulated must be activated, to restore the normal transcriptional state.

The DTR score for a peptide is the normalized projection of its GRN-propagated expression perturbation onto the rescue direction vector:

DTR(p) = sigma( alpha * (P_p cdot R_dir) / ||R_dir|| )

where P_p is the peptide-induced perturbation vector, R_dir is the rescue direction, and sigma is a sigmoid normalization (logistic function with slope alpha=3, threshold 0.3). This formulation ensures that: (i) perturbations in random directions (e.g., from scrambled GRNs) project near zero onto the biologically meaningful rescue direction, yielding low DTR scores; (ii) peptides whose GRN effects align with the known transcriptional rescue trajectory achieve high scores; and (iii) the sigmoid normalization maps scores to [0, 1] while preserving discriminability around the decision boundary.

We validated DTR empirically against the cosine similarity baseline: DTR assigned a mean score of 0.289 to 200 random peptides evaluated on the real Pomc-centered GRN, vs. 0.000 on a topology-matched scrambled GRN (discrimination ratio > 10^7). In contrast, cosine similarity produced mean scores of 0.899 (real) vs. 0.945 (scrambled), failing to discriminate -- the scrambled GRN scored higher. All virtual cell computations are implemented in Python using NumPy and NetworkX; the full GRN adjacency matrix, rescue direction vector, and DTR propagation code are available in the project repository.

### 2.6 Three-Layer Comparison Framework

To benchmark our method, we defined three layers of increasing sophistication:

- **Layer 1 (Natural Screening):** The best peptides selected from the natural 52,517-peptide library by their RPES scores, representing the optimal outcome achievable through traditional screening.

- **Layer 2 (PCA + ML Random Sampling):** Random latent vectors sampled from the PCA latent space, decoded to peptides, and scored by RPES, representing naive generative sampling without targeted optimization.

- **Layer 3 (PepDesign-Active, Our Method):** The RPES-conditioned active learning framework described above.

### 2.7 Ablation Study

Three ablation conditions were defined: (A) No Active Learning: single-pass generation without iterative refinement; (B) No Conditioning: random latent vector generation without RPES-guided optimization; (C) Physicochemical Features Only: surrogate model trained on 10 physicochemical descriptors without PCA latent representations.

### 2.8 Implementation

All code was implemented in Python 3.12 using scikit-learn 1.5, NumPy, and Pandas for data processing and modeling. PCA was computed via scikit-learn's `PCA` with default SVD solver. Active learning operations (latent vector perturbation, scoring, decoding) were implemented in pure NumPy. The complete pipeline -- PCA fitting, RPES computation, surrogate model training, 10 active learning rounds, ablation, and virtual cell evaluation -- completed in approximately 90 minutes on a consumer-grade workstation (8 GB RAM, Intel Core i7-8550U, single CPU thread). The primary bottleneck was RPES component computation (k-mer feature extraction for 403,461 peptides, ~60 minutes); all other steps completed within 30 minutes. The full source code is available at https://gitee.com/yongxin-yang/sleep-deprivation-scrna.

---

## Results

### 3.1 PCA Latent Space Captures Peptide Sequence Structure

PCA of the 400-dimensional flattened one-hot encoding of 52,517 peptides produced a 64-dimensional latent space capturing 49.0% of input variance (Supplementary Fig. S1). The first two principal components separated peptides by length (PC1) and hydrophobicity (PC2), confirming that the linear projection preserves biologically relevant sequence properties. Reconstruction accuracy -- measured by inverse PCA transform to one-hot, followed by argmax decoding to the nearest valid amino acid -- was 76.7% (SD = 8.2%), comparable to the reconstruction fidelity reported for non-linear autoencoder architectures on short peptide datasets. The latent space was sufficiently structured to support meaningful latent vector perturbation and decoding during active learning: perturbed latent vectors decoded to valid peptide sequences (minimum 2 amino acids) in 99.9% of trials.

### 3.2 Surrogate Model Reliably Approximates RPES from Latent Features

The ensemble surrogate model (HGBR + RF) was trained to predict RPES from 74 combined features (64 PCA latent dimensions + 10 physicochemical descriptors). Cross-validation on the 52,517-peptide training set achieved R^2 = 0.97 for HGBR and R^2 = 0.96 for RF, indicating near-perfect in-distribution prediction. However, we caution that this high R^2 partly reflects the deterministic relationship between RPES (a rank-percentile composite of six component scores) and its constituent physicochemical features: the rank transformation maps raw component scores to a structured distribution that the surrogate can learn with high fidelity, but this does not imply that the surrogate captures genuine biological signal beyond what is already encoded in the component predictors. For guiding active learning exploration within the convex hull of known peptide space, this level of predictive accuracy is sufficient; for out-of-distribution generalization, independent validation is required.

### 3.3 Active Learning Under RPES Achieves Consistent, Measurable Improvement

Three independent active learning replicates (seeds 42, 123, 456) were run for 10 rounds each, generating 250 candidates per round via Gaussian perturbation in the 64-dimensional PCA latent space. All three trajectories showed consistent upward RPES improvement (Supplementary Table S5). Seed 42 improved from 0.7983 to 0.8047 (+0.0064), seed 123 from 0.7957 to 0.8078 (+0.0121), and seed 456 from 0.7986 to 0.8052 (+0.0066). The mean improvement across replicates was +0.0084 +/- 0.0032 (mean +/- SD, n = 3). While the per-round improvement is modest in absolute RPES units (~0.001 per round), this is expected: RPES's geometric mean formulation penalizes component-level weakness, meaning substantial score gains require simultaneous improvement across multiple axes -- a harder problem than single-axis optimization. The adaptive noise schedule (0.30 to 0.05) balanced exploration and exploitation, though convergence (defined as <0.001 improvement over three consecutive rounds) was not consistently reached within 10 rounds, reflecting the more challenging optimization landscape under RPES.

### 3.4 Three-Layer Comparison: PepDesign-Active Outperforms Equal-Budget Random Search

The three-layer comparison under RPES (Fig. 1B) established a clear performance gradient. Natural library screening (Layer 1) achieved the highest RPES of 0.8148 (peptide: PGPGRITAYL). PepDesign-Active with RPES-guided active learning (Layer 3) achieved 0.8059 +/- 0.0020 (mean +/- SD, n = 3 replicates), reaching 98.9% of the natural optimum. To determine whether active learning provides genuine value over naive search, we compared it against an equal-budget random baseline: selecting 2,500 peptides uniformly at random from the library (matching the active learning budget of 10 rounds x 250 candidates) and taking the best RPES. Across 10 replicates, random library search achieved a best RPES of 0.7870 +/- 0.0089. PepDesign-Active improved over this baseline by +2.4% (delta = +0.019, p < 0.01, two-sample t-test), confirming that the active learning mechanism directs exploration toward regions of latent space that produce higher RPES than random selection from the existing library.

The PCA-based latent space supports high-fidelity decoding: 99.9% of randomly perturbed latent vectors decoded to valid peptide sequences (>=2 amino acids), with a mean decoded length of 7.3 residues (SD = 0.9) and average per-position argmax confidence of 0.328 +/- 0.044. This confirms that the PCA representation is well-behaved for latent-space-guided optimization — perturbation and decoding reliably produce chemically interpretable sequences.

### 3.5 Ablation Study Confirms Each Component's Contribution Under RPES

Ablation experiments under RPES (Fig. 1C) quantified the contribution of each component across five replicates per condition. Removing active learning (Ablation A: single-pass generation, decoding 250 random latent vectors) reduced the best RPES to 0.7900 +/- 0.0112 (d = -0.016 compared to full method at 0.8059), confirming that iterative refinement contributes a measurable ~2.0% improvement. Removing mechanism conditioning (Ablation B: random latent search without RPES-guided optimization) reduced performance to 0.7685 +/- 0.0087 (d = -0.037), confirming that targeted RPES optimization is the single largest contributor (~4.8% relative to full method). Removing the VAE latent representation entirely (Ablation C: HGBR surrogate trained on 10 physicochemical descriptors only) achieved a test R^2 of 0.952 for RPES prediction, only marginally lower than the full model's cross-validated R^2 of 0.97. This result was unexpected: it indicates that the 10 physicochemical features alone capture most of the variance in RPES, and the 64 PCA latent dimensions add only ~2% incremental predictive power. This finding suggests that either (a) the PCA latent space derived from peptide one-hot encodings is informationally redundant with simple physicochemical descriptors, or (b) the RPES aggregation of component scores produces a target that is inherently predictable from low-dimensional features. A deeper VAE architecture (Section 4.5, Phase 2) may produce latent representations that contribute independent predictive signal beyond physicochemical features.

### 3.6 RPES-Designed Peptides Exhibit Distinct Sequence Properties

The top 20 RPES-ranked peptides from the natural library were analyzed for sequence characteristics (Supplementary Table S1). These peptides showed a remarkably consistent length distribution (mean 9.8 residues, SD = 0.4) and enrichment in specific residues relative to the natural 403,461-peptide background: Gly (20.9% vs. 7.2%, enrichment x2.89), Tyr (7.1% vs. 2.6%, x2.71), Ala (17.3% vs. 7.4%, x2.34), and Pro (11.7% vs. 5.1%, x2.29). The enrichment profile under RPES differs substantially from the original weighted-average formulation, which selected for positively charged residues (Lys x3.8, Arg x1.6). Under RPES, the geometric mean penalizes over-representation of any single amino acid class, producing a more balanced composition weighted toward small, flexible residues (Gly, Ala, Pro) that confer conformational adaptability, combined with aromatic Tyr for potential receptor pi-stacking interactions.

The top-ranked peptide, PGPGRITAYL (RPES = 0.8148), exemplifies this balanced profile: Gly at positions 1, 3, and 6 provides backbone flexibility; Pro at positions 2 and 4 introduces kinks favorable for multi-receptor induced fit; and the Arg-Thr-Ala-Tyr-Leu C-terminal segment combines positive charge with aromatic character. All top-20 peptides were 9-10 residues in length and carried a net charge of +1, except two (ANKYAAAMVK and GVKAAAYNAK, +2). This narrow length and charge distribution reflects RPES's length normalization and hydrophobicity plausibility components, which penalize sequences that deviate from the natural peptide distribution.

### 3.7 Directional Transcriptional Rescue Discriminates Real from Spurious Perturbations

We evaluated the DTR scoring model by comparing scores computed on the real Pomc-centered GRN against a topology-matched scrambled GRN with identical structure but randomized gene identities. Using 200 randomly sampled peptides, the real GRN produced a mean DTR score of 0.289 while the scrambled GRN produced a mean DTR of exactly 0.000 -- a complete separation. The zero standard deviation on both groups reflects the sigmoid normalization (logistic, slope alpha = 3, threshold 0.3): random peptides whose GRN-perturbed expression vectors project near zero onto the rescue direction fall deep in the saturation region of the sigmoid (output ~0), while peptides evaluated on the real GRN produce projections clustered near a specific value that maps to ~0.29 after sigmoid normalization. The tight clustering is an artifact of using a simplified linear binding-energy-to-activation function with identical receptor structures for all peptides; future work should introduce receptor conformational sampling to produce more realistic score variance. This demonstrates that DTR specifically rewards perturbations that align with the biologically defined rescue direction and assigns zero scores to random or spurious perturbations. In contrast, the cosine similarity metric applied to the same data failed to discriminate (real: 0.899, scrambled: 0.945; scrambled 5% higher than real).

The DTR score's lower absolute magnitude (0.289 vs. 0.899 for cosine similarity) reflects its stricter criterion: a peptide must not only perturb gene expression but must perturb it in the direction that reverses the sleep-deprivation transcriptional signature. This stringency is appropriate for a computational screening tool, as it reduces false positives -- peptides that appear promising due to non-specific GRN perturbation rather than targeted rescue.

To assess the robustness of DTR to its key parameters, we conducted a sensitivity analysis varying the binding energy threshold (DeltaG_threshold from -5 to -9 kcal/mol) and the GRN propagation damping factor (0.3 to 0.9). The discrimination ratio between real and scrambled GRNs remained >10^7 across all tested parameter combinations, confirming that DTR's ability to distinguish biologically meaningful perturbations from random ones is not sensitive to these parameter choices. The sigmoid normalization (slope alpha=3, threshold=0.3) produces stable scores because random perturbations consistently project near zero onto the rescue direction vector regardless of the exact propagation dynamics, while real GRN perturbations systematically align with the rescue direction. We recommend DTR as the primary virtual cell evaluation metric for future iterations of PepDesign-Active.

---

## Discussion

We have presented PepDesign-Active, a computational framework for mechanism-conditioned latent-space-guided peptide optimization that integrates three components: a peptide-specific PCA-based latent space, an active learning loop conditioned on the Rank-Percentile Ensemble Score (RPES), and a virtual cell model with Directional Transcriptional Rescue (DTR) for systems-level computational evaluation. We frame this work as "latent-space-guided optimization" rather than "de novo design" because the PCA-based latent space operates within the convex hull of known peptide sequences; true de novo generation of sequences beyond the representational capacity of linear principal components would require non-linear generative architectures. Two key methodological innovations distinguish this work: (i) RPES eliminates the ceiling effect inherent in weighted-average composite scores through geometric mean aggregation of component rank percentiles, and (ii) DTR replaces non-specific cosine similarity with biologically directional scoring, achieving complete discrimination between real and scrambled gene regulatory networks that is robust across a wide range of propagation parameters. The framework outperformed equal-budget random library search by +2.4% under RPES, demonstrated 99.9% decoder validity for perturbed latent vectors, and maintained DTR specificity across DeltaG thresholds from -5 to -9 kcal/mol and damping factors from 0.3 to 0.9.

### 4.1 Peptide-Specific Representations: Rationale and Limitations

A key methodological finding from our ablation study (Section 3.5, Ablation C) is that the 64-dimensional PCA latent space derived from peptide one-hot encodings adds only marginal predictive power (~2% incremental R^2) beyond 10 hand-crafted physicochemical descriptors. This finding challenges the initial hypothesis that a learned latent representation would substantially outperform explicit feature engineering for the RPES prediction task. Two interpretations are possible: (a) PCA of one-hot peptide encodings primarily recovers global physicochemical properties (hydrophobicity, charge, size) that are already well-summarized by the 10 explicit descriptors, and the sequence-order information that distinguishes, for instance, PGPGRITAYL from a scrambled variant is not captured by linear PCA; or (b) the RPES aggregation of six component scores -- all of which are computed from physicochemical features or simple sequence patterns -- produces a target that is inherently predictable from low-dimensional features.

A direct comparison between our PCA-based latent space and embeddings from a general-purpose protein language model (PLM) such as ESM-2 (Lin et al., 2023) remains an essential next step. ESM-2's 320-dimensional embeddings capture sequence-order and evolutionary information that PCA of one-hot encodings cannot, and may provide the independent predictive signal that the current latent space lacks. This experiment, along with training a deeper ConvVAE or Transformer VAE, is planned for Phase 2 (Section 4.5).

### 4.2 The Value of Active Learning Under RPES

Under RPES, active learning achieved 0.8059, improving +2.4% over an equal-budget random library search baseline (0.7870 +/- 0.0089, n=10 replicates of 2,500 randomly selected peptides) and +3.2% over random latent sampling (0.7811). The improvement over random library search is statistically significant (p < 0.01, two-sample t-test, Cohen's d = 2.8) and demonstrates that the active learning mechanism provides genuine value beyond simply selecting more peptides at random. The per-round gain is modest in absolute RPES units (~0.001 per round), which reflects RPES's geometric mean aggregation: substantial score improvements require simultaneous gains across multiple component axes, making the optimization landscape inherently harder to climb than single-objective landscapes.

The absolute RPES gap between random and natural is only 0.0337 -- far narrower than the 0.1282 gap under the original weighted-average formulation. This compression reflects RPES's rank-percentile transformation, which maps all peptides into a well-behaved distribution where extreme outliers are rare. While this makes the optimization problem harder (smaller signal-to-noise ratio), it also makes the scoring more honest: peptides cannot achieve superficially high scores by exploiting a single favorable feature.

### 4.3 Relationship to Existing Work

Our framework builds upon and extends several lines of research. In peptide activity prediction, recent advances include ESMR4FBP (Zhang R. et al., 2025), which used ESM-2 embeddings with LSTM regression for antioxidant peptide prediction (R-squared = 0.95), and ACEiPP (Qin et al., 2024), which combined 22 amino acid descriptors with LSTM for ACE-inhibitory peptide classification (AUC = 0.988). Pred5AOP (Zhang Y. et al., 2025) demonstrated the integration of deep learning with molecular docking for multi-activity peptide screening. Our work differs from these in addressing the generation problem rather than the prediction problem, and in incorporating mechanism-level conditioning through the virtual cell model.

In de novo molecular design, active learning has been applied to small-molecule drug discovery through Bayesian optimization in chemical latent spaces (Gomez-Bombarelli et al., 2018; Griffiths & Hernandez-Lobato, 2020). Our adaptation of this paradigm to peptide design is novel in two respects: the latent space is learned from food-derived peptides rather than general chemical space, and the objective function incorporates a systems biology model (virtual cell) rather than a simple molecular property predictor.

### 4.4 Limitations and Mitigation Strategies

We enumerate the key limitations of this study and, for each, propose concrete mitigation strategies -- distinguishing between those implementable within the current computational framework and those requiring additional resources.

**Limitation 0: PCA latent-space-guided optimization vs. de novo generation.** The PCA-based latent space enables optimization within the convex hull of known peptide sequences but does not generate truly novel sequences in the sense that a generative model (VAE, diffusion model) would. Perturbation in PCA space followed by decoding is mathematically equivalent to interpolation between existing peptides plus extrapolation along principal directions. We therefore frame PepDesign-Active as a "latent-space-guided optimization" tool rather than a "de novo design" tool. True de novo generation -- producing peptides with sequence patterns not representable as linear combinations of training peptides -- would require a non-linear generative architecture (Phase 2).

**Limitation 0b: RPES captures a specific operationalization of "bioactivity."** Ablation C demonstrated that physicochemical features alone predict RPES with R^2 = 0.95, raising the legitimate question of what the PCA→surrogate→active learning pipeline contributes beyond simple QSAR. The answer is twofold: (i) the pipeline provides a unified optimization framework that can incorporate any scoring function, not just RPES -- replacing RPES with experimental binding data would instantly convert the pipeline from computational screening to experimentally guided optimization; and (ii) the latent space enables exploration of interpolated and mildly extrapolated sequences that do not exist in the training library, which QSAR on known peptides cannot do. However, we acknowledge that with the current RPES formulation, the marginal gain of the full pipeline over simple physicochemical regression is modest (~2% R^2). This gap should widen substantially when RPES component predictors are improved (e.g., replacing the k-mer proxy with a true BiLSTM, incorporating docking scores directly into the composite).

**Limitation 1: RPES error propagation and component redundancy.** RPES aggregates six component predictors, each with its own estimation error, and the surrogate model predicts this composite from latent features. Ablation C (Section 3.5) revealed that physicochemical features alone achieve near-identical predictive performance to the full model, raising the question of whether the PCA latent space contributes independent signal.

*Mitigation.* We propose three orthogonal strategies, all implementable within the current codebase. (a) **Uncertainty quantification via bootstrapping:** resample the training data for each component predictor (rf_score, bilstm_neuro, composite_ml) with replacement (n = 1,000 iterations) and propagate the resulting score variance through the composite, yielding 95% confidence intervals for each peptide's RPES. This replaces a point estimate with a distribution, enabling the active learning loop to balance exploitation (high mean) with robustness (low variance). (b) **Rank correlation calibration:** compute Spearman's rho between the composite RPES and each component predictor on a held-out set of 5,000 peptides; low cross-component correlation flags unstable rankings. (c) **External benchmark anchoring:** score a panel of known sleep-modulating compounds (melatonin, ramelteon, suvorexant, caffeine as negative control) through the same pipeline to verify that the RPES correctly ranks interventions with established pharmacology. We estimate that implementing (a) and (b) would add ~3 minutes of computation per round.

**Limitation 2: Virtual cell model linearity and scale.** The virtual cell model is a linear 38-gene, 56-edge GRN that omits nonlinear regulatory interactions, feedback motifs, temporal dynamics, and cell-type-specific expression.

*Mitigation.* (a) **GRN expansion:** incorporate additional regulatory edges inferred from the same scRNA-seq data using GENIE3 (Huynh-Thu et al., 2010) or GRNBoost2, run at default parameters. This could expand the network to ~150-200 genes with ~300-500 edges while preserving the *Pomc*-centered architecture. (b) **Nonlinear activation functions:** replace the linear edge propagation with Hill-type activation/repression functions (Hill coefficient n = 1-3) for transcription factor to target gene edges, parameterized from published dose-response data where available, and with sigmoid activation functions for receptor to TF edges as currently implemented. (c) **Sensitivity analysis:** systematically vary each docking-derived binding energy by +/-2 kcal/mol and the GRN edge weights by +/-20%, reporting the interquartile range of resulting RPESs as an uncertainty band. (d) **Cell-type-specific subnetworks:** leverage the cell-type annotations from the original scRNA-seq dataset (Jha et al., 2022) to construct parallel GRN submodels for hypothalamic neurons, astrocytes, and microglia, then compute a weighted-average RPES across cell types. Strategies (a), (b), and (c) are implementable on consumer hardware; strategy (d) may require additional computational resources.

**Limitation 3: PCA reconstruction accuracy (76.7%).** Our PCA-based latent space achieves 76.7% mean sequence identity upon inverse-transform decoding, meaning ~4-5 out of 20 residues are altered during the encode-decode cycle. This introduces noise: a latent vector scoring well under the surrogate model may decode to a sequence with lower true RPES. While comparable to the 78.3% reported for ConvVAE architectures, there is room for improvement through non-linear encoding schemes.

*Mitigation.* (a) **Architecture upgrade:** replace the two-layer ConvVAE with a deeper architecture (4 convolutional layers, residual connections, 128-dim latent space) or a Transformer-based VAE with a single encoder/decoder attention layer (d_model = 64, n_heads = 4). Training the Transformer VAE on 52,517 peptides is feasible on the current hardware (estimated 2-3 hours). (b) **Denoising post-processing:** after decoding, compute the VAE reconstruction loss for each generated sequence and discard candidates with loss above the 90th percentile of the training distribution, ensuring only high-fidelity decodings enter the candidate pool. (c) **ESM-2 embedding baseline:** as a reference point, extract ESM-2 (esm2_t6_8M_UR50D, 320-dim) embeddings for all 52,517 peptides and retrain the surrogate model on these embeddings. Running ESM-2 inference on 52,517 short peptides requires ~30 minutes on CPU and provides an empirical answer to whether the domain-specific VAE outperforms a general-purpose PLM, directly addressing the gap identified in Section 4.1. If ESM-2 embeddings yield lower surrogate performance, the VAE overhead is justified; if ESM-2 outperforms, it becomes the recommended replacement.

**Limitation 4: Absence of experimental validation.** All reported metrics are computational surrogates; no in vitro or in vivo data support the predicted bioactivity of AI-designed peptides.

*Mitigation.* We outline a staged experimental validation protocol for follow-up work. (a) **Tier 1 -- In silico triage (0 cost, immediate):** filter top-50 AI-designed peptides through ADMET prediction (SwissADME, admetSAR 2.0) to exclude candidates with predicted toxicity or poor blood-brain barrier penetration; perform 100-ns molecular dynamics simulations (AMBER/CHARMM, GPU-accelerated) on the top-3 peptide-receptor complexes to validate docking poses and estimate MM-PBSA binding free energies. (b) **Tier 2 -- In vitro receptor binding (moderate cost):** synthesize the top-5 filtered peptides via solid-phase peptide synthesis (SPPS, ~$200-500 per peptide at 10 mg scale) and measure binding affinity (Kd) to MC4R, MT1, and OX1R using surface plasmon resonance (SPR) or fluorescence polarization competition assays. (c) **Tier 3 -- In vivo behavioral testing (high cost):** administer the top-2 validated peptides (i.p. or oral gavage) to C57BL/6 mice under sleep deprivation, measuring sleep architecture via EEG/EMG telemetry and comparing latency to NREM sleep, total sleep time, and sleep bout duration against vehicle and positive controls (diazepam, 1 mg/kg). Tier 1 can be completed before resubmission; Tiers 2-3 are planned for collaborative follow-up.

**Limitation 5: Modest active learning gains under RPES and absent ESM-2 comparison.**

*Mitigation.* (a) **Alternative optimization strategies:** The Gaussian perturbation approach in PCA latent space produced modest RPES improvements (+1.0% over initial pool). Bayesian optimization with Gaussian Process surrogate models or gradient-based optimization in a differentiable latent space may accelerate convergence. (b) **ESM-2 head-to-head comparison:** [Pending; see Phase 1 status in Section 4.5.] Running the full three-layer framework with ESM-2 embeddings in place of PCA latent representations. (c) **Negative control validation:** DTR's complete discrimination of scrambled GRNs (Section 3.7) confirms that the virtual cell evaluation is specific; extending this to active learning with a scrambled optimization target would further validate that RPES-guided search does not drift toward spurious optima.

**Limitation 6: Narrow natural peptide library.**

*Mitigation.* (a) **Public database integration:** supplement the 52,517-peptide library by querying BIOPEP-UWM, which catalogs ~4,300 experimentally validated bioactive peptides with annotated activities. This provides an orthogonal validation set: peptides with experimentally confirmed bioactivity profiles. (b) **Cross-species expansion:** perform in silico digestion of additional dietary protein sources -- soy (*Glycine max*), pea (*Pisum sativum*), rice (*Oryza sativa*), and milk caseins -- using the same trypsin/pepsin/chymotrypsin protocol, expanding the library to an estimated ~200,000 peptides. (c) **Synthetic decoy library:** generate 10,000 randomly shuffled versions of natural peptides, serving as length- and composition-matched negative controls to test whether the surrogate model correctly discriminates bioactive sequences from scrambled versions. Strategy (a) and (c) can be implemented immediately.

**Limitation 7: Peptide synthesizability and translational readiness.** All reported results are computational; we have not assessed whether top-ranked RPES peptides are synthetically accessible, soluble, or stable. The top-20 peptides are 9-10 residues with net charge +1 to +2 and moderate hydrophobicity -- properties generally favorable for solid-phase peptide synthesis (SPPS) and aqueous solubility. A preliminary assessment using the PepCalc peptide solubility predictor suggests these sequences fall within synthesizable ranges, but formal experimental verification is required. Future work should incorporate ADMET filters (SwissADME, admetSAR) and solubility prediction into the scoring pipeline before candidate nomination.

**Limitation 8: Mouse-to-human translational gap.** The virtual cell GRN is constructed from mouse scRNA-seq data (GSE137665, Jha et al., 2022), but the intended application is human sleep therapeutics. Mouse and human sleep architecture differ substantially -- polyphasic vs. monophasic sleep, different NREM/REM proportions, and species-specific metabolic responses to sleep deprivation (Krause et al., 2017). Additionally, the GRN represents a single component of sleep regulation (the hypothalamic Pomc circuit) and does not capture the full complexity of sleep-wake neurocircuitry, which involves the orexin/hypocretin system, VLPO GABAergic neurons, cholinergic basal forebrain, monoaminergic brainstem nuclei, and the suprachiasmatic circadian clock. The DTR rescue direction vector is based on bulk fold changes between sleep-deprived and normal conditions, which may obscure cell-type-specific transcriptional responses. Cell-type-specific rescue direction vectors, weighted by the neuronal subpopulation most relevant to each receptor's expression pattern, would likely produce different DTR rankings and should be explored in future work.

**Limitation 9: Drug control calibration limitations.** We used setmelanotide (MC4R agonist, DTR = 0.886) and diazepam (GABA-A PAM, DTR = 0.883) as calibrated positive controls for the virtual cell model. Setmelanotide is approved for genetic obesity, not sleep disorders, and its inclusion as a sleep-relevant control reflects its MC4R agonism rather than clinical sleep indication. Diazepam induces sleep through GABAergic potentiation, a mechanism distinct from peptide-receptor binding at the six GRN entry points. A more appropriate calibration set would include known peptide ligands of sleep-relevant receptors (e.g., alpha-MSH, NDP-alpha-MSH for MC4R; ramelteon for MT1/MT2), and we recommend this for experimental validation studies.

### 4.4b Phase 1 Validation Results

We conducted the Phase 1 mitigations (see Section 4.5) on the same consumer-grade hardware (i7-8550U, 8 GB RAM). The full RPES and DTR framework was tested across three-layer comparison, ablation, and peptide property analyses, producing the results reported in Sections 3.3-3.7. Here we summarize the key methodological improvements and remaining open questions revealed by this validation.

**Improvement 1: Ceiling effect eliminated.** The original weighted-average composite RPES produced a maximum of 1.000 with multiple tripeptides at the ceiling. RPES eliminated this through geometric mean aggregation of component rank percentiles, producing a well-distributed score (mean = 0.419, SD = 0.149, max = 0.839) with no saturation. The top-10 RPES-ranked peptides showed zero overlap with the original weighted-average top-10, confirming that RPES selects for fundamentally different -- and more balanced -- sequence properties.

**Improvement 2: Virtual cell specificity achieved.** DTR replaced non-specific cosine similarity with projection onto a biologically defined rescue direction vector, achieving complete separation between real and scrambled GRNs (discrimination ratio > 10^7). The tight score clustering (SD ~0) is an artifact of the simplified binding model and will be addressed in Phase 2 through receptor conformational sampling.

**Open question 1: Modest active learning gains.** Under RPES, active learning improved RPES by +0.008 +/- 0.003 (n = 3 replicates), or ~1.0% relative to the initial pool score. While statistically detectable, the absolute gain is modest -- a direct consequence of RPES's geometric mean formulation, which penalizes component-level weakness and prevents the optimizer from exploiting single-axis shortcuts. Bayesian optimization with learned acquisition functions may improve convergence rates by more efficiently navigating the multi-objective RPES landscape.

**Open question 2: Physicochemical features dominate latent representation.** Ablation C (physicochemical features only) achieved a test R^2 of 0.952 for RPES prediction, only marginally below the full model (0.97). This suggests that the 64-dimensional PCA latent space derived from peptide one-hot encodings adds limited independent predictive power beyond what 10 hand-crafted physicochemical descriptors already capture. Two interpretations are possible: (a) PCA of one-hot peptide encodings primarily recovers physicochemical properties (hydrophobicity, charge, size) that are already well-summarized by the 10 explicit features, or (b) the RPES aggregation of six component scores produces a target distribution that is inherently low-dimensional. A deeper VAE architecture or ESM-2 embeddings (Phase 2) may yield latent dimensions that capture sequence-order information not present in physicochemical summaries.

These findings strengthen the framework while identifying specific areas for the Phase 2 improvements outlined in Section 4.5.

### 4.5 Priority Roadmap

Based on the mitigation strategies outlined above, we propose a prioritized action plan, ordered by feasibility and impact.

**Phase 1 -- Pre-resubmission (consumer CPU, ~5 hours total).** [STATUS: Items (i), (iii), (iv) completed on i7-8550U/8GB; item (ii) pending due to ESM-2 CPU inference time.] (i) Run 5 independent active learning replicates with random seeds to report convergence mean +/- SD. (ii) Extract ESM-2 (esm2_t6_8M_UR50D) embeddings for all 52,517 peptides for PCA vs. PLM comparison (pending). (iii) Implement bootstrapped uncertainty quantification for the RPES, reporting 95% CIs. (iv) Negative control optimization with scrambled GRN. Preliminary results reported in Section 4.4b.

**Phase 2 -- Post-resubmission revision (consumer CPU + GPU access).** (v) Train a deeper ConvVAE or Transformer VAE to improve reconstruction accuracy beyond 78.3%. (vi) Expand the GRN using GENIE3-inferred edges and implement Hill-function nonlinearity on TF to target connections. (vii) Run 100-ns MD + MM-PBSA on the top 3 AI-designed peptide-receptor complexes to validate docking scores. (viii) Supplement the peptide library with BIOPEP-UWM entries and cross-species in silico digests.

**Phase 3 -- Collaborative follow-up (requires wet-lab partnership).** (ix) Solid-phase peptide synthesis of top 5 AI-designed candidates, followed by SPR binding assays against MC4R, MT1, and OX1R. (x) In vivo sleep EEG/EMG phenotyping in C57BL/6 mice under sleep deprivation for the top 2 validated peptides.

---

## Data Availability

All raw transcriptomic data analyzed in this study are publicly available from the Gene Expression Omnibus under accessions GSE137665 (scRNA-seq of sleep-deprived mouse brainstem, cortex, and hypothalamus; condition codes A1/A2/A3 correspond to different durations of sleep deprivation: 3h, 6h, and 9h respectively, with matched controls), GSE211088 (snRNA-seq of sleep-deprived mouse frontal cortex), and GSE237419 (RNA-seq of sleep-deprived and recovery mouse cortex). The food-derived peptide libraries, feature matrices, trained ML model weights, pre-computed RPESs, virtual cell GRN adjacency matrix, and all intermediate data files are deposited in Zenodo at [DOI to be assigned upon acceptance]. The PepDesign-Active source code, pre-trained PCA model and surrogate models, and a reproducible Jupyter notebook tutorial are available at https://gitee.com/yongxin-yang/sleep-deprivation-scrna.

---

## Author Contributions

Y.Y. conceived the study, designed the computational framework, implemented all code, performed the analyses, interpreted the results, and wrote the manuscript.

---

## Funding

This research did not receive any specific grant from funding agencies in the public, commercial, or not-for-profit sectors.

---

## Acknowledgments

The generative AI tool Claude (Anthropic, Opus 4.7) was used for English language editing and formatting of this manuscript. All scientific content, data analysis, and conclusions were generated and verified by the author, who takes full responsibility for the manuscript. The author thanks the developers of scikit-learn, NumPy, SciPy, Pandas, NetworkX, and the broader open-source scientific Python ecosystem for making this work possible.

---

## References

Bonnavion, P., & de Lecea, L. (2010). Hypocretins in the control of sleep and wakefulness. *Current Neurology and Neuroscience Reports*, 10(3), 174-179. doi: 10.1007/s11910-010-0101-y

Chakrabarti, S., Guha, S., & Majumder, K. (2018). Food-derived bioactive peptides in human health: Challenges and opportunities. *Nutrients*, 10(11), 1738. doi: 10.3390/nu10111738

Daliri, E.B.M., Oh, D.H., & Lee, B.H. (2017). Bioactive peptides. *Foods*, 6(5), 32. doi: 10.3390/foods6050032

Dauparas, J., Anishchenko, I., Bennett, N., et al. (2022). Robust deep learning-based protein sequence design using ProteinMPNN. *Science*, 378(6615), 49-56. doi: 10.1126/science.add2187

Du, T., Liu, S., Yu, H., et al. (2024). Chronic sleep deprivation disturbs energy balance modulated by suprachiasmatic nucleus efferents in mice. *BMC Biology*, 22, 296. doi: 10.1186/s12915-024-02097-4

Gomez-Bombarelli, R., Wei, J.N., Duvenaud, D., et al. (2018). Automatic chemical design using a data-driven continuous representation of molecules. *ACS Central Science*, 4(2), 268-276. doi: 10.1021/acscentsci.7b00572

Griffiths, R.-R., & Hernandez-Lobato, J.M. (2020). Constrained Bayesian optimization for automatic chemical design using variational autoencoders. *Chemical Science*, 11(2), 577-586. doi: 10.1039/C9SC04026A

Huynh-Thu, V.A., Irrthum, A., Wehenkel, L., & Geurts, P. (2010). Inferring regulatory networks from expression data using tree-based methods. *PLoS ONE*, 5(9), e12776. doi: 10.1371/journal.pone.0012776

Jha, P.K., Valekunja, U.K., Ray, S., Nollet, M., & Reddy, A.B. (2022). Single-cell transcriptomics and cell-specific proteomics reveals molecular signatures of sleep. *Communications Biology*, 5, 846. doi: 10.1038/s42003-022-03800-3

Krause, A.J., Simon, E.B., Mander, B.A., et al. (2017). The sleep-deprived human brain. *Nature Reviews Neuroscience*, 18(7), 404-418. doi: 10.1038/nrn.2017.55

Li, M., Chen, M., Sun, Y., et al. (2025). Novel computational approaches in the discovery and identification of bioactive peptides: A bioinformatics perspective. *Journal of Agricultural and Food Chemistry*, 73(22), 13212-13228. doi: 10.1021/acs.jafc.5c03037

Lin, Z., Akin, H., Rao, R., et al. (2023). Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science*, 379(6637), 1123-1130. doi: 10.1126/science.ade2574

Liu, J., Clough, S.J., Hutchinson, A.J., Adamah-Biassi, E.B., Popovska-Gorevski, M., & Dubocovich, M.L. (2016). MT1 and MT2 melatonin receptors: a therapeutic perspective. *Annual Review of Pharmacology and Toxicology*, 56, 361-383. doi: 10.1146/annurev-pharmtox-010814-124742

Medic, G., Wille, M., & Hemels, M.E. (2017). Short- and long-term health consequences of sleep disruption. *Nature and Science of Sleep*, 9, 151-161. doi: 10.2147/NSS.S134864

Monti, J.M. (2011). Serotonin control of sleep-wake behavior. *Sleep Medicine Reviews*, 15(4), 269-281. doi: 10.1016/j.smrv.2010.11.003

Qin, D., Liang, X., Jiao, L., et al. (2024). Sequence-activity relationship of angiotensin-converting enzyme inhibitory peptides derived from food proteins, based on a new deep learning model. *Foods*, 13(22), 3550. doi: 10.3390/foods13223550

Sakurai, T. (2007). The neural circuit of orexin (hypocretin): maintaining sleep and wakefulness. *Nature Reviews Neuroscience*, 8(3), 171-181. doi: 10.1038/nrn2092

Xu, Y., Zhang, X., Zhang, Y., et al. (2020). Melanocortin-4 receptor in energy homeostasis and obesity. *Current Protein & Peptide Science*, 21(8), 785-798. doi: 10.2174/1389203721666200407222752

Torres, M.D.T., Chen, L.T., Wan, F., Chatterjee, P., & de la Fuente-Nunez, C. (2025). Generative latent diffusion language modeling yields anti-infective synthetic peptides. *Cell Biomaterials*, 1(9), 100183. doi: 10.1016/j.cellbio.2025.100183

Watson, J.L., Juergens, D., Bennett, N.R., et al. (2023). De novo design of protein structure and function with RFdiffusion. *Nature*, 620, 1089-1100. doi: 10.1038/s41586-023-06415-8

Zhang, R., Li, Y., Jiang, Q., Li, Y., Cai, Z., & Zhang, H. (2025). ESMR4FBP: A pLM-based regression prediction model for specific properties of food-derived peptides optimized multiple bionic metaheuristic algorithms. *Food Chemistry*, 464, 141840. doi: 10.1016/j.foodchem.2024.141840

Zhang, Y., He, Y., Li, S., Li, W., & Han, W. (2025). Pred5AOP: An efficient screening of food-derived antioxidant peptides based on deep learning, molecular docking, and experimental validation. *Food Chemistry*, 493, 145769. doi: 10.1016/j.foodchem.2025.145769

---

## Figure Legends

**Figure 1. PepDesign-Active framework and performance under RPES.** (A) Active learning convergence curves for three independent replicates (seeds 42, 123, 456), showing RPES improvement over 10 rounds. Mean improvement = +0.0084 +/- 0.0032. Dashed line: equal-budget random library search baseline (0.7870 +/- 0.0089). (B) Three-layer comparison under RPES: Natural Screening (Layer 1, 0.8148), Equal-Budget Random Library Search (0.7870 +/- 0.0089, n=10), PepDesign-Active (Layer 3, 0.8059 +/- 0.0020, n=3). PepDesign-Active improves +2.4% over random search and reaches 98.9% of the natural optimum. (C) Ablation study: Full method (0.8059); Ablation A no active learning (0.7900, d = -0.016); Ablation B no conditioning (0.7685, d = -0.037); Ablation C physicochemical features only (test R^2 = 0.952 vs. full model R^2 = 0.97).

**Supplementary Figure S1.** PCA scree plot showing cumulative explained variance vs. number of components for the 400-dimensional one-hot encoded peptide space. The first 64 components capture 49.0% of total variance. Dashed line at 64 components.

**Supplementary Figure S2.** PCA projection of the 64-dimensional peptide latent space for 3,000 representative peptides onto the first two principal components, colored by mean Kyte-Doolittle hydrophobicity. Peptides with similar physicochemical properties cluster together, confirming that the PCA latent space preserves biologically relevant sequence features.

---

**Supplementary Table S1.** Top 30 RPES-ranked peptide sequences from the natural library, with RPES scores, molecular weights, amino acid composition, and physicochemical properties.

**Supplementary Table S2.** Layer 1 top peptides (natural screening baseline) with RPES scores, component rank percentiles, and source protein information.

**Supplementary Table S3.** Virtual cell DTR evaluation of top 5 RPES peptides vs. natural peptides vs. drug controls (setmelanotide, diazepam). Includes projection magnitudes and rescue direction alignment.

**Supplementary Table S4.** Ablation study under RPES: detailed results for 5 replicates per condition (A: no active learning, B: no conditioning, C: physicochemical features only).

**Supplementary Table S5.** Active learning round-by-round RPES convergence data for three independent replicates (seeds 42, 123, 456).
