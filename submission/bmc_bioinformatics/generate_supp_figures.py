"""generate_supp_figures.py — Generate Supplementary Figures S1 and S2."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os, random, sys

sys.path.insert(0, r'D:\projects\sleep-deprivation-project\pepdesign')
from phase1_pipeline import (
    MODEL_DIR, OUTPUT_DIR, PCALatentSpace, load_peptide_library, flatten_onehot, log
)
from sklearn.decomposition import PCA as PCA_sk

OUT = r'D:\projects\sleep-deprivation-project\submission\bmc_bioinformatics\figures'
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 10,
    'axes.labelsize': 11, 'axes.titlesize': 12,
    'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

# Load data
peptides_df = load_peptide_library()
unique_peptides = peptides_df['peptide'].tolist()
vae_peptides = random.sample(unique_peptides, 52517)

# Fit PCA
oh_flat = flatten_onehot(vae_peptides)
pca_full = PCA_sk(n_components=64, random_state=42)
pca_full.fit(oh_flat)

# ---- Fig S1: PCA Scree Plot ----
fig, ax = plt.subplots(figsize=(7, 4.5))
cumvar = np.cumsum(pca_full.explained_variance_ratio_)
ax.plot(range(1, 65), cumvar, 'b-', linewidth=1.5)
ax.axvline(x=64, color='gray', linestyle='--', alpha=0.5, label=f'64 components ({cumvar[-1]*100:.0f}% variance)')
ax.axhline(y=0.49, color='red', linestyle=':', alpha=0.5)
ax.set_xlabel('Number of Principal Components')
ax.set_ylabel('Cumulative Explained Variance')
ax.set_title('Supplementary Figure S1: PCA Scree Plot', fontweight='bold')
ax.legend(fontsize=9)
ax.set_xlim(0, 65)
ax.set_ylim(0, 1.05)
plt.tight_layout()
plt.savefig(f'{OUT}/FigS1_pca_scree.png')
plt.close()
print('Fig S1 saved.')

# ---- Fig S2: PCA Latent Space Projection ----
pca_vis = PCA_sk(n_components=2)
n_vis = min(3000, len(vae_peptides))
vis_idx = np.random.choice(len(vae_peptides), n_vis, replace=False)
oh_vis = oh_flat[vis_idx]
lat_2d = pca_vis.fit_transform(oh_vis)

KD = {'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8,
      'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'Y':-1.3}
colors = [np.mean([KD.get(aa, 0.0) for aa in vae_peptides[i]]) for i in vis_idx]

fig, ax = plt.subplots(figsize=(8, 6))
sc = ax.scatter(lat_2d[:, 0], lat_2d[:, 1], c=colors, cmap='RdYlBu', s=4, alpha=0.6, edgecolors='none')
cbar = plt.colorbar(sc, ax=ax)
cbar.set_label('Mean Kyte-Doolittle Hydrophobicity', fontsize=9)
ax.set_xlabel(f'PC1 ({pca_vis.explained_variance_ratio_[0]*100:.1f}%)')
ax.set_ylabel(f'PC2 ({pca_vis.explained_variance_ratio_[1]*100:.1f}%)')
ax.set_title('Supplementary Figure S2: PCA Latent Space of Peptide Sequences', fontweight='bold')
plt.tight_layout()
plt.savefig(f'{OUT}/FigS2_latent_pca.png')
plt.close()
print('Fig S2 saved.')
print(f'\nAll supplementary figures saved to {OUT}/')
