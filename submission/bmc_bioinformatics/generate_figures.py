"""generate_figures.py — Generate publication-quality figures for BMC Bioinformatics submission."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import json, os

OUT = r'D:\projects\sleep-deprivation-project\submission\bmc_bioinformatics\figures'
os.makedirs(OUT, exist_ok=True)

# BMC style: 300 dpi, Arial/sans-serif, no gridlines, clean
plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 10,
    'axes.labelsize': 11, 'axes.titlesize': 12,
    'figure.dpi': 300, 'savefig.dpi': 300,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.1,
})

# Colors
C1, C2, C3 = '#2196F3', '#FF9800', '#4CAF50'
CGRAY = '#9E9E9E'
CRED = '#F44336'

# ==========================================
# Fig 1A: Active Learning Convergence
# ==========================================
# Three replicates from RPES experiment
seeds = ['Seed 42', 'Seed 123', 'Seed 456']
histories = [
    (42,  [0.7983, 0.8002, 0.8011, 0.8018, 0.8025, 0.8031, 0.8036, 0.8040, 0.8044, 0.8047]),
    (123, [0.7957, 0.8015, 0.8037, 0.8048, 0.8057, 0.8064, 0.8070, 0.8074, 0.8077, 0.8078]),
    (456, [0.7986, 0.7998, 0.8023, 0.8032, 0.8038, 0.8043, 0.8047, 0.8050, 0.8051, 0.8052]),
]

fig, ax = plt.subplots(figsize=(7, 4.5))
for (seed, hist), c, ls in zip(histories, [C1, C2, C3], ['-', '--', '-.']):
    ax.plot(range(1, len(hist)+1), hist, color=c, linestyle=ls, linewidth=1.5,
            marker='o', markersize=4, label=f'Seed {seed}')

# Random baseline
rand_mean, rand_std = 0.7870, 0.0089
ax.axhline(y=rand_mean, color=CGRAY, linestyle=':', linewidth=1.2, alpha=0.8)
ax.fill_between([0.5, 10.5], rand_mean-rand_std, rand_mean+rand_std,
                color=CGRAY, alpha=0.12)
ax.text(5.5, rand_mean-0.004, 'Random library baseline', color=CGRAY, fontsize=9, ha='center')

ax.set_xlabel('Active Learning Round')
ax.set_ylabel('Best RPES')
ax.set_xlim(0.5, 10.5)
ax.legend(fontsize=8, loc='lower right')
ax.set_title('A', fontweight='bold', loc='left')
plt.tight_layout()
plt.savefig(f'{OUT}/Fig1A_al_convergence.png')
plt.close()
print('Fig 1A saved.')

# ==========================================
# Fig 1B: Three-Layer Comparison
# ==========================================
fig, ax = plt.subplots(figsize=(6, 4.5))
layers = ['Natural\nScreening\n(Layer 1)', 'Random Library\nSearch\n(Equal Budget)',
          'PepDesign-\nActive\n(Layer 3)']
values = [0.8148, 0.7870, 0.8059]
errors = [0, 0.0089, 0.0020]
colors = [C1, CGRAY, C3]

bars = ax.bar(layers, values, color=colors, edgecolor='white', linewidth=0.8, width=0.5)
ax.errorbar(layers, values, yerr=errors, fmt='none', ecolor='black', capsize=5, linewidth=1.2)

# Annotate values
for bar, val in zip(bars, values):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.008, f'{val:.4f}',
            ha='center', fontsize=10, fontweight='bold')

# Improvement arrow
ax.annotate('+2.4%', xy=(2, 0.8059), xytext=(1.5, 0.815),
            fontsize=10, fontweight='bold', color=C3,
            arrowprops=dict(arrowstyle='->', color=C3, lw=1.5))

ax.set_ylabel('Best RPES')
ax.set_ylim(0.74, 0.84)
ax.set_title('B', fontweight='bold', loc='left')
ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
plt.tight_layout()
plt.savefig(f'{OUT}/Fig1B_three_layer.png')
plt.close()
print('Fig 1B saved.')

# ==========================================
# Fig 1C: Ablation Study
# ==========================================
fig, ax = plt.subplots(figsize=(6, 4.5))
conditions = ['Full\nMethod', 'A: No Active\nLearning', 'B: No\nConditioning',
              'C: Physchem\nOnly']
values = [0.8059, 0.7900, 0.7685, 0.8059]
errors_c = [0.0020, 0.0112, 0.0087, 0.0050]
colors_c = [C3, C1, C2, CGRAY]

# Ablation C uses R² instead of RPES, show on secondary axis conceptually
bars = ax.bar(conditions[:3], values[:3], color=colors_c[:3], edgecolor='white', linewidth=0.8, width=0.45)
ax.errorbar(conditions[:3], values[:3], yerr=errors_c[:3], fmt='none', ecolor='black', capsize=5)

for bar, val in zip(bars, values[:3]):
    ax.text(bar.get_x() + bar.get_width()/2, val + 0.012, f'{val:.4f}',
            ha='center', fontsize=9, fontweight='bold')

# Annotate deltas
ax.annotate('', xy=(0, 0.8059), xytext=(1, 0.7900),
            arrowprops=dict(arrowstyle='<->', color=CRED, lw=1.2))
ax.text(0.5, 0.796, f'd = -0.016', ha='center', fontsize=8, color=CRED, fontweight='bold')

ax.annotate('', xy=(0, 0.8059), xytext=(2, 0.7685),
            arrowprops=dict(arrowstyle='<->', color=CRED, lw=1.2))
ax.text(1.5, 0.781, f'd = -0.037', ha='center', fontsize=8, color=CRED, fontweight='bold')

# Ablation C as text
ax.text(3, 0.8059, f'Test R^2\n= 0.952\n(vs 0.97 full)',
        ha='center', fontsize=9, fontweight='bold', color=CGRAY,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=CGRAY, alpha=0.5))

ax.set_ylabel('Best RPES')
ax.set_ylim(0.72, 0.84)
ax.set_title('C', fontweight='bold', loc='left')
plt.tight_layout()
plt.savefig(f'{OUT}/Fig1C_ablation.png')
plt.close()
print('Fig 1C saved.')

print(f'\nAll figures saved to {OUT}/')
