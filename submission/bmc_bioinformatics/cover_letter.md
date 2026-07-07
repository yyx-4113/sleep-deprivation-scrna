Dear Editor,

I am pleased to submit our manuscript entitled **"Mechanism-Conditioned Latent-Space Optimization of Sleep-Promoting Peptides via Active Learning with Virtual Cell Validation"** for consideration for publication in *BMC Bioinformatics*.

**Why BMC Bioinformatics?** This manuscript reports a computational framework that integrates dimensionality reduction, ensemble machine learning, active learning, and single-cell transcriptomics-derived gene regulatory network modeling — spanning multiple subfields of computational biology that fall squarely within the journal's scope. The methods are fully reproducible on consumer hardware (8 GB RAM, single CPU thread), aligning with BMC Bioinformatics's emphasis on accessible computational methodology.

**What is novel?** We present two methodological innovations: (i) a Rank-Percentile Ensemble Score (RPES) that eliminates the ceiling effect inherent in conventional weighted-average composite scoring through geometric mean aggregation of component rank percentiles; and (ii) a Directional Transcriptional Rescue (DTR) metric that replaces non-specific cosine similarity with biologically informed projection onto a disease-relevant transcriptomic rescue direction. We also document, transparently, our failed attempts to reproduce a published ConvVAE architecture and our decision to use PCA as a deterministic, computationally accessible alternative — a reproducibility note that we believe strengthens the paper's scientific rigor.

**Key findings:**
- Active learning with RPES outperforms equal-budget random library search by +2.4% (p < 0.01, Cohen's d = 2.8)
- DTR achieves complete discrimination between real and scrambled gene regulatory networks (>10^7 ratio), robust across binding thresholds and damping factors
- PCA-based latent space achieves 76.7% reconstruction accuracy with 99.9% decoder validity
- All experiments completed on consumer hardware in ~90 minutes; full source code publicly available

**Declarations:** This manuscript has not been published elsewhere and is not under consideration by another journal. All authors have approved the manuscript and agree with its submission to BMC Bioinformatics. The author declares no competing interests. This research received no specific grant funding. The generative AI tool Claude (Anthropic, Opus 4.7) was used for English language editing; all scientific content was generated and verified by the author.

**Suggested reviewers:**
1. Dr. Pawan K. Jha — University of California, San Diego (sleep scRNA-seq expertise; GSE137665 dataset author)
2. Dr. Marcus D. T. Torres — University of Pennsylvania (de novo peptide design; generative peptide models)
3. Dr. Jose Manuel Hernandez-Lobato — Microsoft Research / University of Cambridge (Bayesian optimization; active learning in chemical space)

Thank you for your time and consideration.

Sincerely,
Yongxin Yang
Fujian Second People's Hospital
No. 282 Wusi Road, Gulou District
Fuzhou, Fujian 350003, China
E-mail: 960856791@qq.com
