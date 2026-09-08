#!/usr/bin/env python
"""Figures for the NBSCR data commentary (manuscript_commentary.md).

Only uses values already computed in results/tables (no new analysis).
Outputs PNGs to figures/.
"""
import csv
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 600,  # Elsevier line-art requirement (>=600 dpi at final size)
})

NS_C, SD_C, RS_C = "#3C6E9F", "#C0504D", "#6E9F5C"
GREY = "#8C8C8C"


def fig1_design():
    fig, ax = plt.subplots(figsize=(7.0, 3.2))
    regions = ["Hypothalamus", "Brainstem", "Cortex"]
    conds = [("NS", NS_C), ("SD", SD_C), ("RS", RS_C)]
    for ri, r in enumerate(regions):
        for ci, (c, col) in enumerate(conds):
            x = ci * 3 + ri
            ax.add_patch(plt.Rectangle((x, 0), 0.86, 1.0, facecolor=col,
                                       edgecolor="white", linewidth=2))
            ax.text(x + 0.43, 0.5, f"{c}\n1 library\n(3 mice pooled)",
                    ha="center", va="center", color="white", fontsize=7.5)
    ax.set_xlim(-0.3, 9.3)
    ax.set_ylim(-0.35, 1.5)
    ax.set_yticks([])
    ax.set_xticks([1.0, 4.0, 7.0])
    ax.set_xticklabels(regions, fontsize=9)
    ax.set_title("GSE137665: 9 libraries = 3 regions x 3 conditions, 0 replicate libraries",
                 fontsize=10, loc="left")
    ax.text(0, -0.28,
            "Terminal tissue collection -> NS, SD, RS are different animals -> all contrasts are "
            "between pooled groups.\nSample-level n = 1 per condition: no between-library variance "
            "is estimable.",
            fontsize=8, color="#333333", va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig1_design.png"))
    plt.close(fig)


def fig2_pseudoreplication():
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.9))
    ax = axes[0]
    ax.bar(["Cell-level\nWilcoxon", "Pseudo-bulk\nCPM ratio"], [-4.08, -3.10],
           color=[GREY, NS_C], width=0.6)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("log2 fold change (SD vs NS)")
    ax.set_title("A  Two effect sizes, one contrast", fontsize=9.5, loc="left")
    ax.text(0, -3.6, "cell-level p\n= read sampling", fontsize=7, ha="center")
    ax.text(1, -2.7, "n = 1 library\nper condition", fontsize=7, ha="center")
    ax.set_ylim(-5, 0.6)

    ax = axes[1]
    libs = ["NS (A1)", "SD (A2)", "RS (A3)"]
    cpm = [22369, 2605, 21251]
    cols = [NS_C, SD_C, RS_C]
    b = ax.bar(libs, cpm, color=cols, width=0.6)
    for r, v in zip(b, cpm):
        ax.text(r.get_x() + r.get_width() / 2, v + 600, f"{v:,}", ha="center", fontsize=8)
    ax.set_ylabel("*Pomc* CPM (hypothalamus)")
    ax.set_title("B  Three point observations, no variance", fontsize=9.5, loc="left")
    ax.text(0.02, 0.62,
            "log2 SD/NS = -3.10 (8.6x)\nRS = 0.95x NS\nno CI, no test available",
            transform=ax.transAxes, fontsize=8, color="#333333", va="top")
    ax.set_ylim(0, 27000)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_pseudoreplication.png"))
    plt.close(fig)


def fig3_composition():
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.7))
    ax = axes[0]
    ax.bar(["any\ndetected", ">=2 aux\nmarkers", ">=3 aux\nmarkers"], [45.0, 30.6, 7.5],
           color=[GREY, NS_C, "#2E5E88"], width=0.65)
    ax.set_ylabel("% hypothalamic cells *Pomc*+")
    ax.set_title("A  Definition matters", fontsize=9, loc="left")
    ax.set_ylim(0, 52)

    ax = axes[1]
    ax.bar(["NS", "SD", "RS"], [35.1, 5.1, 44.4], color=[NS_C, SD_C, RS_C], width=0.65)
    ax.set_ylabel("% *Pomc*+ (multi-marker)")
    ax.set_title("B  Composition tracks condition", fontsize=9, loc="left")
    ax.set_ylim(0, 52)

    ax = axes[2]
    ax.bar(["Neuronal", "Non-neuronal\n(ambient)"], [63.5, 36.5],
           color=[NS_C, "#B07AA1"], width=0.65)
    ax.set_ylabel("% cells detecting *Pomc*")
    ax.set_title("C  Ambient RNA", fontsize=9, loc="left")
    ax.set_ylim(0, 78)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig3_composition.png"))
    plt.close(fig)


def fig4_leakage():
    fig, ax = plt.subplots(figsize=(4.2, 2.9))
    ax.bar(["Cell-level\nCV (leaky)", "Library-holdout\nCV", "Permuted\nlabels"],
           [0.93, 0.50, 0.50], color=["#B07AA1", NS_C, GREY], width=0.6)
    ax.axhline(0.5, color="black", linestyle="--", linewidth=0.9)
    ax.text(2.0, 0.53, "chance", fontsize=7.5, ha="center")
    ax.set_ylabel("AUC")
    ax.set_ylim(0, 1.1)
    ax.set_title("Six-gene classifier: leakage vs correct holdout", fontsize=9.5, loc="left")
    ax.text(0.02, 0.72, "library-holdout vs\npermutation: P = 0.57",
            transform=ax.transAxes, fontsize=8, color="#333333")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_leakage.png"))
    plt.close(fig)


def fig5_external():
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    ax.bar(["GSE137665\n(unreplicated)", "GSE243489\n(5 libraries/condition)"],
           [3.03, -1.24], color=[RS_C, "#B07AA1"], width=0.55)
    ax.axhline(0, color="black", linewidth=0.9)
    ax.set_ylabel("log2 FC *Pomc* (recovery vs sleep deprivation / comparator)")
    ax.text(0, 3.15, "rebound", ha="center", fontsize=8)
    ax.text(1, -1.45, "no rebound", ha="center", fontsize=8)
    ax.set_ylim(-2.2, 4.0)
    ax.set_title("Directionally discordant for *Pomc*", fontsize=9.5, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_external.png"))
    plt.close(fig)


def fig6_inventory():
    # curated inventory (Table 1)
    data = [
        ("GSE137665", 1, "Commun Biol"),
        ("GSE214337", 1, "Nature"),
        ("GSE213496", 1, "Cell"),
        ("GSE218624", 1, "Nat Neurosci"),
        ("GSE280145", 1, "unpublished"),
        ("GSE289089", 2, "unpublished"),
        ("GSE211088", 3, "iScience"),
        ("GSE245537", 3, "Science"),
        ("GSE256140", 4, "eLife"),
        ("GSE243489", 5, "Cell Rep"),
    ]
    fig, ax = plt.subplots(figsize=(6.6, 3.1))
    names = [d[0] for d in data]
    vals = [d[1] for d in data]
    cols = [SD_C if v == 1 else NS_C for v in vals]
    b = ax.barh(names, vals, color=cols)
    for r, d in zip(b, data):
        ax.text(r.get_width() + 0.08, r.get_y() + r.get_height() / 2,
                f"n={d[1]}  ({d[2]})", va="center", fontsize=7.5)
    ax.set_xlabel("Independent libraries per condition")
    ax.set_xlim(0, 7.2)
    ax.invert_yaxis()
    ax.set_title("5 of 10 sleep single-cell datasets have one library per condition",
                 fontsize=9.5, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig6_inventory.png"))
    plt.close(fig)


def fig7_checklist():
    items = [
        "1. How many independent libraries per condition?",
        "2. What is the unit of inference?",
        "3. Is pooling disclosed?",
        "4. Were conditions measured in the same animals?",
        "5. Is the effect size computed at library level?",
        "6. Has composition been separated from expression?",
        "7. Is ambient RNA quantified?",
        "8. Are cell identities defined by multiple markers?",
        "9. Was machine learning split by library?",
        "10. Is there independent replication for the specific claim?",
    ]
    fig, ax = plt.subplots(figsize=(6.6, 3.4))
    ax.axis("off")
    ax.set_title("Box 1 — Ten questions to ask of a single-cell dataset with a condition contrast",
                 fontsize=10, loc="left", pad=10)
    for i, t in enumerate(items):
        ax.text(0.01, 0.90 - i * 0.093, t, fontsize=9, va="top", color="#222222")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig7_checklist.png"))
    plt.close(fig)


if __name__ == "__main__":
    fig1_design()
    fig2_pseudoreplication()
    fig3_composition()
    fig4_leakage()
    fig5_external()
    fig6_inventory()
    fig7_checklist()
    print("figures written to", FIG)
    for f in sorted(os.listdir(FIG)):
        print("  ", f)
