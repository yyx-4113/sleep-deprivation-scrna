"""
Generate NSFC-quality mechanism diagrams (v3 -- matching template style).
Template: modern, clean, white bg, colored accent boxes, light tints.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import os

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['SimHei', 'Microsoft YaHei', 'Arial'],
    'axes.unicode_minus': False,
})
DPI = 300

# Color palette -- extracted from template
C = {
    'white':      '#FCFEFE',
    'warm_bg':    '#FCF4E8',
    'blue_bg':    '#E7EBFE',
    'blue':       '#454FFE',
    'blue_mid':   '#94A6FD',
    'coral':      '#DD4B5C',
    'coral_lt':   '#E9C5B8',
    'orange':     '#FEAF63',
    'mint':       '#4FFCDA',
    'mint_lt':    '#B8FEEB',
    'green':      '#4CE53E',
    'green_lt':   '#A4FCA0',
    'purple':     '#827199',
    'yellow':     '#E8FF4B',
    'dark':       '#323039',
    'gray':       '#827199',
    'lgray':      '#C8C8D0',
}

def box_solid(ax, x, y, w, h, color, text='', fontsize=9, bold=True, z=4):
    """Solid colored box with white text."""
    b = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.05',
                       facecolor=color, edgecolor='none',
                       linewidth=0, zorder=z)
    ax.add_patch(b)
    if text:
        ax.text(x+w/2, y+h/2, text, ha='center', va='center',
                fontsize=fontsize, fontweight='bold' if bold else 'normal',
                color='white', zorder=z+1)

def box_tinted(ax, x, y, w, h, bg_color, border_color, text='',
               fontsize=9, bold=True, text_color=None, z=4):
    """Light-tinted box with colored border."""
    b = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.05',
                       facecolor=bg_color, edgecolor=border_color,
                       linewidth=1.5, zorder=z)
    ax.add_patch(b)
    if text:
        ax.text(x+w/2, y+h/2, text, ha='center', va='center',
                fontsize=fontsize, fontweight='bold' if bold else 'normal',
                color=text_color or border_color, zorder=z+1)

def box_outline(ax, x, y, w, h, color, text='', fontsize=9, bold=True, z=4):
    """White box with colored border."""
    b = FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.05',
                       facecolor='white', edgecolor=color,
                       linewidth=2, zorder=z)
    ax.add_patch(b)
    if text:
        ax.text(x+w/2, y+h/2, text, ha='center', va='center',
                fontsize=fontsize, fontweight='bold' if bold else 'normal',
                color=color, zorder=z+1)

def arrow(ax, x1, y1, x2, y2, color=C['gray'], lw=1.8, z=2, ls='-'):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw, ls=ls),
                zorder=z)

def arrow_d(ax, x, y1, y2, color=C['gray'], lw=1.8, z=2):
    arrow(ax, x, y1, x, y2, color, lw, z)

def arrow_r(ax, x1, x2, y, color=C['gray'], lw=1.8, z=2):
    arrow(ax, x1, y, x2, y, color, lw, z)


# ===================================================================
# FIG 1: NSFC -- Pomc-Pcsk2 mechanism (template style)
# ===================================================================
def draw_nsfc(output_path):
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 11)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # Title (clean, modern)
    ax.text(8, 10.55, 'Pomc-Pcsk2 正反馈环路在睡眠剥夺中的断裂机制',
            ha='center', fontsize=15, fontweight='bold', color=C['dark'])
    ax.text(8, 10.2, 'Sleep Deprivation Disrupts the Pomc-Pcsk2 Positive Feedback Loop',
            ha='center', fontsize=7.5, color=C['gray'], style='italic')

    # Top trigger bar
    box_solid(ax, 5.5, 9.35, 5.0, 0.5, C['coral'], 'Sleep Deprivation  睡眠剥夺', fontsize=10)
    arrow_d(ax, 8, 9.35, 8.95, C['coral'], lw=2.5)

    # Main cell panel (light warm bg)
    cell = FancyBboxPatch((1.0, 2.3), 14.0, 6.4, boxstyle='round,pad=0.25',
                          facecolor=C['warm_bg'], edgecolor='none', zorder=0)
    ax.add_patch(cell)
    ax.text(1.5, 8.4, '下丘脑弓状核  Pomc+ 神经元', fontsize=9,
            fontweight='bold', color=C['gray'])

    # Nucleus compartment (light blue bg)
    nuc = FancyBboxPatch((1.4, 3.0), 5.8, 5.0, boxstyle='round,pad=0.15',
                         facecolor=C['blue_bg'], edgecolor=C['blue_mid'],
                         linewidth=1.5, zorder=1)
    ax.add_patch(nuc)
    ax.text(4.3, 7.78, '细胞核  Nucleus', fontsize=8, color=C['blue'],
            fontweight='bold', ha='center', zorder=2)

    # Nuclear elements
    box_outline(ax, 2.0, 6.9, 2.0, 0.45, C['blue'], 'CREB / STAT3', fontsize=8)
    arrow_d(ax, 3.0, 6.9, 6.5, C['blue'], lw=1.5)
    box_outline(ax, 2.0, 6.0, 2.0, 0.45, C['blue'], 'Pomc 基因', fontsize=8)
    arrow_d(ax, 3.0, 6.0, 5.6, C['blue'], lw=1.5)

    # mRNA
    box_tinted(ax, 1.6, 5.05, 2.8, 0.48, C['coral_lt'], C['coral'],
               'Pomc mRNA  DOWN', fontsize=8.5, text_color=C['coral'])
    ax.annotate('log2FC = -4.08\npadj = 1.3e-157', xy=(4.4, 5.29),
                fontsize=6.5, color=C['coral'], ha='left', va='center', style='italic')

    # Cytoplasm label
    ax.text(11.5, 7.78, '细胞质  Cytoplasm', fontsize=8, color=C['orange'],
            fontweight='bold', ha='center')

    # mRNA -> POMC
    arrow_r(ax, 4.4, 8.0, 5.29, C['gray'], lw=1.8)
    box_tinted(ax, 8.0, 6.1, 2.8, 0.48, C['coral_lt'], C['coral'],
               'POMC 前体蛋白  DOWN', fontsize=8.5, text_color=C['coral'])

    # POMC -> Pcsk2
    arrow_d(ax, 9.4, 6.1, 5.68, C['gray'], lw=1.8)
    # GRNBoost2 importance label on Pomc -> Pcsk2 edge
    ax.text(10.0, 5.63, 'GBR imp=0.766\n(Top 1 target)', fontsize=5.5,
            color=C['coral'], ha='center', style='italic', zorder=6)
    box_tinted(ax, 8.0, 5.15, 2.8, 0.48, C['coral_lt'], C['coral'],
               'Pcsk2 蛋白酶  DOWN', fontsize=8.5, text_color=C['coral'])

    # Pcsk2 promoter TF binding predictions (new!)
    ax.text(8.0, 4.35, 'Pcsk2 Promoter: CREB1 (-291) Fos/JunB (-607 TRE)',
            fontsize=5.5, color=C['blue'], ha='left', style='italic', zorder=6)
    # Virtual KO predicted downstream targets
    ax.text(8.0, 4.15, 'Predicted targets: Scg2(-0.32) Ccnd2(-0.38) Dbp(↓)',
            fontsize=5.2, color=C['purple'], ha='left', style='italic', zorder=6)

    # Peptides
    arrow(ax, 9.4, 5.15, 8.3, 4.5, C['gray'], lw=1.5)
    arrow(ax, 9.4, 5.15, 10.5, 4.5, C['gray'], lw=1.5)
    box_tinted(ax, 7.3, 4.05, 2.0, 0.42, C['coral_lt'], C['coral'],
               'a-MSH  DOWN', fontsize=8, text_color=C['coral'])
    box_tinted(ax, 9.6, 4.05, 2.2, 0.42, C['coral_lt'], C['coral'],
               'b-endorphin  DOWN', fontsize=8, text_color=C['coral'])

    # Feedback loop
    ax.annotate('', xy=(3.0, 6.2), xytext=(8.8, 5.15),
                arrowprops=dict(arrowstyle='->', color=C['orange'], lw=4,
                                connectionstyle='arc3,rad=-0.55'),
                zorder=5)
    # GBR importance on feedback edge
    ax.text(4.5, 3.65, 'GBR imp=0.510 (Top 1 TF)\nPcsk2 -> Pomc',
            fontsize=5.5, color=C['orange'], ha='center', style='italic', zorder=6)
    ax.text(5.8, 3.3, '正反馈环路  Positive Feedback', fontsize=8.5,
            fontweight='bold', color=C['orange'], ha='center',
            bbox=dict(boxstyle='round,pad=0.25', facecolor='white',
                      edgecolor=C['orange'], alpha=0.9))
    ax.text(5.8, 2.95, '[X] 环路断裂  (双打击放大器)', fontsize=8.5,
            fontweight='bold', color=C['coral'], ha='center',
            bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                      edgecolor=C['coral'], lw=1.8))

    # Bottom: Outcomes
    outcomes = [
        (2.5, '睡眠-觉醒\n节律紊乱', C['coral']),
        (5.5, '氧化应激\n增加', C['orange']),
        (8.5, '神经炎症\nTNF-a/IL-1b', C['purple']),
        (11.5, 'HPA轴激素\n节律异常', C['blue']),
        (14.5, '下游靶基因\n表达紊乱', C['green']),
    ]
    for ox, otext, oc in outcomes:
        box_tinted(ax, ox-0.9, 1.3, 1.8, 0.6, C['white'], oc, otext, fontsize=7,
                   text_color=oc)

    arrow_d(ax, 8, 2.3, 2.0, C['dark'], lw=2)

    # Right sidebar
    sx = 14.5
    for i, (mnum, mdesc, mc) in enumerate([
        ('I', '多层面验证\nmRNA/蛋白/肽', C['blue']),
        ('II', '因果验证\nsiRNA+Rescue', C['orange']),
        ('III', '机制解析\nChIP-qPCR', C['purple']),
    ]):
        y = 7.8 - i * 0.8
        box_solid(ax, sx-0.35, y, 0.7, 0.6, mc, mnum, fontsize=7)
        ax.text(sx+0.5, y+0.3, mdesc, fontsize=6, color=C['dark'], va='center')

    plt.tight_layout(pad=0.3)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=DPI, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)
    print(f'[OK] NSFC: {output_path}')


# ===================================================================
# FIG 2: TCM -- AI peptide pipeline (template style)
# ===================================================================
def draw_tcm(output_path):
    fig, ax = plt.subplots(1, 1, figsize=(13, 10))
    ax.set_xlim(0, 17)
    ax.set_ylim(0, 12.5)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor('white')

    # Title
    ax.text(8.5, 12.1, 'AI 驱动的闽产药食同源食材助眠活性肽高效发现管道',
            ha='center', fontsize=14, fontweight='bold', color=C['dark'])

    # Three ingredients (top)
    foods = [
        (2.8, '建瓯核桃', 'Juglans regia\n褪黑素+Trp  SSA-LSTM验证', C['blue']),
        (8.5, '武夷桑葚', 'Morus alba\n助眠肽空白  创新点最大', C['purple']),
        (14.2, '闽北黑芝麻', 'Sesamum indicum\n蛋白22%  ACE肽已验证', C['green']),
    ]
    for fx, fname, fdesc, fc in foods:
        box_solid(ax, fx-1.5, 11.0, 3.0, 0.58, fc, fname, fontsize=10)
        ax.text(fx, 10.68, fdesc, fontsize=6.2, color=C['gray'], ha='center')
        arrow_d(ax, fx, 11.0, 10.78, C['lgray'], lw=1.2)

    # Merge arrows
    arrow(ax, 2.8, 10.78, 6.5, 10.35, C['lgray'], lw=1, ls='dotted')
    arrow(ax, 14.2, 10.78, 6.5, 10.35, C['lgray'], lw=1, ls='dotted')
    arrow_d(ax, 8.5, 10.6, 10.35, C['dark'], lw=2.5)

    # Modules
    MW, MH = 6.5, 0.72
    MX_L, MX_R = 1.5, 9.0
    R1, R2, R3 = 9.5, 8.25, 7.0

    mods = [
        (MX_L, R1, MW, MH, '一  虚拟胃肠酶切', 'In Silico Digestion',
         '胃蛋白酶/胰蛋白酶/糜蛋白酶/木瓜蛋白酶  |  <=2漏切, 2-20AA -> 去冗余',
         '->  5,000+ 肽段库', C['blue'], '第1-3月'),
        (MX_R, R1, MW, MH, '二  AI 深度学习筛选', 'BiLSTM + Attention + RF',
         'AAC+PseAAC+二肽组成->440维  |  多标签:神经保护/抗氧化/ACE/抗炎',
         '->  60 条候选肽段', C['green'], '第3-6月'),
        (1.5, R2, 14.0, 0.8, '三  多靶点分子对接 + MD 模拟验证',
         'AutoDock Vina + GROMACS 50ns',
         'GABA_AR(6D6T) | OX2R(7TAC) | MC4R(6W25) | 5-HT1A(7E2Z) | Keap1(2FLU)  |  60肽x5靶点=300对接  |  <-8kcal/mol -> TOP10 -> MD验证TOP5',
         '->  10 条高结合候选肽', C['orange'], '第6-9月'),
        (MX_L, R3, MW, MH, '四  虚拟细胞模型验证', 'Pomc Neuron Virtual Cell',
         'Pomc+神经元GRN(正常/SD/恢复)  |  肽段->受体->TF->网络传播->基因表达',
         '->  逆转分数排序', C['purple'], '第9-12月'),
        (MX_R, R3, MW, MH, '五  体外实验初步验证', 'HT22 Oxidative Stress Model',
         'H2O2诱导(200uM)+肽段梯度(10/50/100uM)  |  CCK-8活力+ROS探针',
         '->  3 条验证候选肽', C['coral'], '第12-18月'),
    ]

    for mx, my, mw, mh, mtitle, msub, mitems, mout, mc, mtime in mods:
        b = FancyBboxPatch((mx, my), mw, mh, boxstyle='round,pad=0.1',
                           facecolor='white', edgecolor=mc, linewidth=2, zorder=3)
        ax.add_patch(b)
        # Left accent
        accent = Rectangle((mx, my), 0.06, mh, facecolor=mc, zorder=4)
        ax.add_patch(accent)

        ax.text(mx+0.2, my+mh-0.12, mtitle, fontsize=9, fontweight='bold',
                color=mc, va='top', zorder=5)
        ax.text(mx+0.2, my+mh-0.33, msub, fontsize=6.2, color=C['gray'],
                va='top', style='italic', zorder=5)
        ax.text(mx+0.25, my+mh-0.47, mitems, fontsize=6.5, color=C['dark'],
                va='top', zorder=5)
        ax.text(mx+mw-0.15, my+0.06, mout, fontsize=7.5, fontweight='bold',
                color=mc, va='bottom', ha='right', zorder=5)
        ax.text(mx+mw-0.15, my+mh-0.12, mtime, fontsize=6, color=C['lgray'],
                va='top', ha='right', zorder=5)

    # Connectors
    arrow_d(ax, 4.75, R1, R2+0.78, C['lgray'], lw=2)
    arrow_d(ax, 12.25, R1, R2+0.78, C['lgray'], lw=2)
    arrow_d(ax, 8.5, R2, R3+0.7, C['lgray'], lw=2)
    arrow_d(ax, 4.75, R2, R3+0.7, C['lgray'], lw=2)
    arrow_d(ax, 12.25, R2, R3+0.7, C['lgray'], lw=2)

    # Target receptor panel
    py = 4.15
    panel = FancyBboxPatch((1.0, py), 15.0, 2.4, boxstyle='round,pad=0.15',
                           facecolor=C['blue_bg'], edgecolor='none', zorder=0)
    ax.add_patch(panel)
    ax.text(8.5, py+2.1, '五靶点协同信号传导网络  |  Multi-Target Signaling Network',
            fontsize=9, fontweight='bold', color=C['dark'], ha='center')

    targets = [
        ('GABA_AR', 'Ca2+ down\nPKA down', 'CREB down', '镇静\nPomc down', C['blue']),
        ('OX2R', 'Ca2+ up\nMAPK up', 'CREB up', '觉醒 down\n(拮抗)', C['green']),
        ('MC4R', 'cAMP up\nPKA up', 'CREB/STAT3 up', 'Pomc up\nPcsk2 up', C['orange']),
        ('5-HT1A', 'cAMP down\nPKA down', 'CREB down', '节律调节', C['purple']),
        ('Keap1', 'Nrf2转位 up\nROS down', 'Nrf2/ARE up', 'HO-1 up\nSOD up', C['coral']),
    ]
    for i, (tn, tp, ttf, te, tc) in enumerate(targets):
        tx = 1.3 + i * 2.95
        box_solid(ax, tx, py+1.2, 2.2, 0.55, tc, tn, fontsize=8)
        ax.text(tx+1.1, py+0.9, tp, fontsize=5.8, color=C['gray'], ha='center')
        box_outline(ax, tx, py+0.15, 2.2, 0.45, tc, ttf, fontsize=7)
        ax.text(tx+1.1, py-0.05, te, fontsize=6, color=tc, ha='center', fontweight='bold')

    # Connectors
    arrow_d(ax, 4.75, 5.6, py+2.4, C['gray'], lw=2)
    arrow_d(ax, 12.25, 5.6, py+2.4, C['gray'], lw=2)

    # Bottom output
    box_solid(ax, 3.5, 1.2, 10.0, 0.55, C['dark'],
              '产出: SCI 1-2篇  |  助眠肽 3条  |  AI筛选代码 1套  |  肽段库 5000+条', fontsize=8)
    arrow_d(ax, 8.5, 3.9, 1.85, C['dark'], lw=2.5)

    # Innovation badge
    ax.text(15.8, 10.2,
            '创新点\n======\n1 首次AI筛选\n  助眠活性肽\n2 虚拟细胞\n  验证药效\n3 "逆转分数"\n  新概念\n4 聚焦闽产\n  地域特色',
            fontsize=5.5, color=C['gray'], va='top', ha='left', linespacing=1.3,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=C['lgray']))

    plt.tight_layout(pad=0.3)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=DPI, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)
    print(f'[OK] TCM: {output_path}')


# ===================================================================
# FIG 3: NSFC Research Roadmap (技术路线图)
# ===================================================================
def draw_nsfc_roadmap(output_path):
    fig, ax = plt.subplots(1, 1, figsize=(12, 14.5))
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 16)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor('white')

    MX = 8.0  # center x
    MW = 13.5  # module width
    ML = MX - MW/2  # module left

    # ---- Title ----
    ax.text(MX, 15.5, '研究技术路线',
            ha='center', fontsize=16, fontweight='bold', color=C['dark'])
    ax.text(MX, 15.05, 'Roadmap: From Computational Discovery to Mechanistic Validation',
            ha='center', fontsize=7.5, color=C['gray'], style='italic')

    # ---- Step 0: Pre-existing findings ----
    y0 = 14.0
    h0 = 0.55
    box_solid(ax, ML, y0, MW, h0, C['blue'],
              '前期计算发现  —  bioRxiv 2026 预印本', fontsize=10)
    ax.text(MX, y0-0.2,
            'Pomc 为睡眠剥夺下丘脑转录组核心枢纽基因 (log2FC = -4.08, padj = 1.3e-157)  |  GRN 35节点  |  虚拟敲除算法',
            fontsize=7, color=C['gray'], ha='center')
    arrow_d(ax, MX, y0, y0-0.45, C['dark'], lw=2.2)

    # ---- Helper: draw a module box ----
    def draw_module(my, mh, mcolor, mnum, mtitle, mtime, items_left, items_right):
        """Draw a module panel with two-column content."""
        b = FancyBboxPatch((ML, my), MW, mh, boxstyle='round,pad=0.12',
                           facecolor='white', edgecolor=mcolor, linewidth=2, zorder=3)
        ax.add_patch(b)
        # Left color accent bar
        accent = Rectangle((ML+0.05, my+0.05), 0.08, mh-0.1,
                           facecolor=mcolor, zorder=4)
        ax.add_patch(accent)
        # Module number badge
        box_solid(ax, ML+0.25, my+mh-0.5, 0.55, 0.38, mcolor, mnum, fontsize=9)
        # Title
        ax.text(ML+0.95, my+mh-0.22, mtitle, fontsize=10.5, fontweight='bold',
                color=C['dark'], va='center')
        # Timeline
        ax.text(MX+MW/2-0.15, my+mh-0.22, mtime, fontsize=7.5, color=C['lgray'],
                va='center', ha='right')
        # Left items
        for j, item in enumerate(items_left):
            ax.text(ML+0.35, my+mh-0.6-j*0.28, item, fontsize=7.2, color=C['dark'], va='center')
        # Right items
        for j, item in enumerate(items_right):
            ax.text(MX+0.15, my+mh-0.6-j*0.28, item, fontsize=7.2, color=C['dark'], va='center')
        # Divider line
        ax.plot([MX, MX], [my+0.1, my+mh-0.65], color=C['lgray'], lw=0.8, ls='--', zorder=5)

    # ---- Module I ----
    y1 = 11.2
    h1 = 2.5
    draw_module(y1, h1, C['coral'], 'I', '多层面表达验证  (mRNA → 蛋白 → 功能肽)',
                '第 1-9 月',
                ['qPCR 验证 Pomc / Pcsk2 及 GRN 中',
                 '  TOP 10 连接基因的 mRNA 表达',
                 '免疫荧光 (IF) 定位 Pomc 在',
                 '  弓状核神经元中的表达变化'],
                ['Western Blot 检测 POMC 前体蛋白',
                 '  及 Pcsk2 蛋白酶的表达水平',
                 'ELISA 定量下游功能肽 (a-MSH、',
                 '  b-内啡肽) 脑脊液/血清浓度'])
    arrow_d(ax, MX, y1, y1+h1+0.18, C['coral'], lw=2.2)

    # ---- Module II ----
    y2 = 7.9
    h2 = 2.5
    draw_module(y2, h2, C['orange'], 'II', '因果验证  (siRNA 敲低 → RNA-seq → Rescue)',
                '第 9-24 月',
                ['小鼠弓状核立体定位注射 Pomc siRNA',
                 '  实现 Pomc 体内敲低',
                 'shRNA 慢病毒转导 HT22 神经元',
                 '  构建 Pomc 稳定敲低细胞系'],
                ['RNA-seq 检测 Pomc KD 后全转录组变化',
                 '  对比虚拟 KO 预测结果',
                 'Rescue: Pomc KD 后回补 a-MSH /',
                 '  POMC 蛋白 → 验证靶基因恢复'])
    arrow_d(ax, MX, y2, y2+h2+0.18, C['orange'], lw=2.2)

    # ---- Module III ----
    y3 = 4.6
    h3 = 2.5
    draw_module(y3, h3, C['purple'], 'III', '机制初探  (Pomc–Pcsk2 正反馈环路解析)',
                '第 12-30 月',
                ['生信分析 Pomc / Pcsk2 启动子区',
                 '  转录因子结合位点',
                 'Pcsk2 siRNA KD → qPCR 检测',
                 '  Pomc 变化，验证环路反向调控'],
                ['ChIP-qPCR 验证关键转录因子 (CREB/',
                 '  STAT3) 对 Pomc/Pcsk2 启动子结合',
                 '整合计算预测 + 湿实验数据 →',
                 '  构建 Pomc–Pcsk2 环路调控模型'])
    arrow_d(ax, MX, y3, y3+h3+0.18, C['purple'], lw=2.2)

    # ---- Output ----
    y_out = 2.8
    h_out = 1.1
    b_out = FancyBboxPatch((ML+1.5, y_out), MW-3, h_out, boxstyle='round,pad=0.15',
                           facecolor=C['dark'], edgecolor='none', zorder=3)
    ax.add_patch(b_out)
    ax.text(MX, y_out+0.65, '预期成果', fontsize=11, fontweight='bold', color='white', ha='center')
    ax.text(MX, y_out+0.2,
            'SCI 论文 2-3 篇  |  阐明 Pomc-Pcsk2 环路在 SD 中的作用机制  |  RNA-seq+ChIP 数据集  |  后续国自然申报基础',
            fontsize=7, color=C['lgray'], ha='center')
    arrow_d(ax, MX, y_out+h_out+0.08, y_out+h_out+0.5, color=C['dark'], lw=2.5)

    # ---- Side annotation: innovation points ----
    ax.text(15.5, 9.0,
            '创新\n====\n1 首次验证\n  Pomc->Pcsk2\n  因果关系\n2 首次构建\n  正反馈环路\n  调控模型\n3 虚拟KO\n  预测+湿\n  实验互验',
            fontsize=5.2, color=C['gray'], va='center', ha='left', linespacing=1.2,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=C['lgray']))

    # ---- Method labels on right ----
    methods = [
        (14.8, 12.4, 'qPCR', C['blue']),
        (14.8, 11.8, 'Western\nBlot', C['coral']),
        (14.8, 11.2, 'IF\nELISA', C['green']),
        (14.8, 9.1, 'siRNA\nKD', C['orange']),
        (14.8, 8.5, 'RNA-seq', C['purple']),
        (14.8, 7.9, 'Rescue', C['coral']),
        (14.8, 5.8, 'ChIP-\nqPCR', C['blue']),
        (14.8, 5.2, 'Pcsk2\nKD', C['green']),
    ]
    for mx, my, mt, mc in methods:
        box_outline(ax, mx-0.5, my-0.2, 1.0, 0.4, mc, mt, fontsize=5.2)

    plt.tight_layout(pad=0.3)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=DPI, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)
    print(f'[OK] NSFC Roadmap: {output_path}')


# ===================================================================
# FIG 4: TCM Research Roadmap (技术路线图)
# ===================================================================
def draw_tcm_roadmap(output_path):
    fig, ax = plt.subplots(1, 1, figsize=(14, 17))
    ax.set_xlim(0, 18)
    ax.set_ylim(0, 18)
    ax.set_aspect('equal')
    ax.axis('off')
    fig.patch.set_facecolor('white')

    MX = 9.0  # center x
    FW = 7.2  # full module width
    FL = MX - FW/2
    HW = 7.2  # half-module width (for two-module rows)
    GAP = 0.6  # gap between two modules

    # ---- Title ----
    ax.text(MX, 17.5, 'AI 驱动的闽产药食同源食材助眠活性肽高效发现管道',
            ha='center', fontsize=15, fontweight='bold', color=C['dark'])
    ax.text(MX, 17.02, 'Technical Roadmap: Integrated Computational-Experimental Pipeline for Sleep-Enhancing Peptide Discovery',
            ha='center', fontsize=7, color=C['gray'], style='italic')

    # ---- Three ingredients (top row) ----
    foods = [
        (4.0, '建瓯核桃', 'Juglans regia', C['blue']),
        (9.0, '武夷桑葚', 'Morus alba', C['purple']),
        (14.0, '闽北黑芝麻', 'Sesamum indicum', C['green']),
    ]
    for fx, fname, flat, fc in foods:
        box_solid(ax, fx-1.4, 16.02, 2.8, 0.52, fc, fname, fontsize=9.5)
        ax.text(fx, 15.8, flat, fontsize=6.2, color=C['gray'], ha='center', style='italic')
        arrow_d(ax, fx, 16.02, 15.82, C['lgray'], lw=1.2)

    # Converge arrows
    arrow(ax, 4.0, 15.82, MX, 15.35, C['lgray'], lw=1.5, ls='dotted')
    arrow(ax, 14.0, 15.82, MX, 15.35, C['lgray'], lw=1.5, ls='dotted')
    arrow_d(ax, MX, 15.35, 14.95, C['dark'], lw=2.5)

    # ---- Helper: two-module row ----
    def draw_module_row(my, mh, left_num, left_title, left_sub, left_items, left_out,
                        left_color, left_time,
                        right_num, right_title, right_sub, right_items, right_out,
                        right_color, right_time):
        """Two modules side by side in one row."""
        for is_left, (lx, mnum, mtitle, msub, mitems, mout, mc, mtime) in enumerate([
            (MX-HW-GAP/2, left_num, left_title, left_sub, left_items, left_out,
             left_color, left_time),
            (MX+GAP/2, right_num, right_title, right_sub, right_items, right_out,
             right_color, right_time),
        ]):
            b = FancyBboxPatch((lx, my), HW, mh, boxstyle='round,pad=0.1',
                               facecolor='white', edgecolor=mc, linewidth=2, zorder=3)
            ax.add_patch(b)
            accent = Rectangle((lx+0.05, my+0.05), 0.07, mh-0.1,
                               facecolor=mc, zorder=4)
            ax.add_patch(accent)
            box_solid(ax, lx+0.22, my+mh-0.44, 0.48, 0.33, mc, mnum, fontsize=8)
            ax.text(lx+0.82, my+mh-0.22, mtitle, fontsize=9.5, fontweight='bold',
                    color=C['dark'], va='center')
            ax.text(lx+0.82, my+mh-0.42, msub, fontsize=6.2, color=C['gray'],
                    va='center', style='italic')
            for j, item in enumerate(mitems):
                ax.text(lx+0.3, my+mh-0.62-j*0.25, item, fontsize=6.8, color=C['dark'],
                        va='center')
            ax.text(lx+HW-0.1, my+0.06, mout, fontsize=7.2, fontweight='bold', color=mc,
                    va='bottom', ha='right')
            ax.text(lx+HW-0.1, my+mh-0.22, mtime, fontsize=6.2, color=C['lgray'],
                    va='center', ha='right')

    # ---- Row 2: Module 1 & 2 ----
    y_row2 = 12.7
    h_row2 = 2.0
    draw_module_row(y_row2, h_row2,
                    '一', '虚拟胃肠酶切', 'In Silico Digestion',
                    ['- UniProt 获取三种食材全部蛋白序列',
                     '- 虚拟酶切: 胃蛋白酶/胰蛋白酶/',
                     '  糜蛋白酶/木瓜蛋白酶',
                     '- <=2漏切位点, 肽段2-20AA',
                     '- 去冗余 -> 理论肽段库'],
                    '-> 5,000+ 肽段库', C['blue'], '第 1-3 月',
                    '二', 'AI 深度学习筛选', 'BiLSTM + Attention + RF',
                    ['- 特征工程: AAC+PseAAC+二肽组成',
                     '  -> ~440 维特征向量',
                     '- 随机森林 (参考Qian 2024 QSAR)',
                     '- BiLSTM+Attention (参考SSA-LSTM-VD)',
                     '- 多标签预测: 神经保护/抗氧化/',
                     '  ACE抑制/抗炎'],
                    '-> 60 条候选肽', C['green'], '第 3-6 月')

    arrow_d(ax, MX-HW-GAP/2, y_row2, y_row2+h_row2+0.15, C['lgray'], lw=2)
    arrow_d(ax, MX+HW+GAP/2, y_row2, y_row2+h_row2+0.15, C['lgray'], lw=2)

    # ---- Row 3: Module 3 (full width) ----
    y_row3 = 9.6
    h_row3 = 2.5
    mcolor3 = C['orange']
    b3 = FancyBboxPatch((1.2, y_row3), 15.6, h_row3, boxstyle='round,pad=0.12',
                        facecolor='white', edgecolor=mcolor3, linewidth=2, zorder=3)
    ax.add_patch(b3)
    accent3 = Rectangle((1.25, y_row3+0.05), 0.08, h_row3-0.1, facecolor=mcolor3, zorder=4)
    ax.add_patch(accent3)
    box_solid(ax, 1.55, y_row3+h_row3-0.52, 0.55, 0.38, mcolor3, '三', fontsize=9)
    ax.text(2.25, y_row3+h_row3-0.28, '多靶点分子对接 + MD 模拟验证',
            fontsize=10.5, fontweight='bold', color=C['dark'], va='center')
    ax.text(2.25, y_row3+h_row3-0.55, 'AutoDock Vina + GROMACS 50ns',
            fontsize=7, color=C['gray'], va='center', style='italic')
    ax.text(15.5, y_row3+h_row3-0.28, '第 6-9 月', fontsize=7, color=C['lgray'], va='center', ha='right')

    # Receptor targets
    tars = ['GABA_AR', 'OX2R', 'MC4R', '5-HT1A', 'Keap1']
    tar_colors = [C['blue'], C['green'], C['orange'], C['purple'], C['coral']]
    for ti, (tn, tc) in enumerate(zip(tars, tar_colors)):
        tx = 1.6 + ti * 3.05
        box_solid(ax, tx, y_row3+1.2, 2.5, 0.5, tc, tn, fontsize=8)
        ax.text(tx+1.25, y_row3+0.9, '60肽 x 5靶点', fontsize=6, color=C['gray'], ha='center')
        ax.text(tx+1.25, y_row3+0.65, '= 300次对接', fontsize=6, color=C['gray'], ha='center')

    # Sub-results
    ax.text(MX, y_row3+0.25, '<-8 kcal/mol -> TOP 10 -> GROMACS 50ns MD -> TOP 5 复合物稳定性验证 (RMSD/RMSF/MM-PBSA)',
            fontsize=7, color=C['dark'], ha='center')
    ax.text(16.5, y_row3+0.06, '-> 10 条\n高结合肽', fontsize=7, fontweight='bold', color=mcolor3,
            va='bottom', ha='right')

    arrow_d(ax, 4.5, 9.6, 9.1, C['gray'], lw=2)
    arrow_d(ax, 13.5, 9.6, 9.1, C['gray'], lw=2)
    arrow_d(ax, MX, 9.1, 8.7, C['dark'], lw=2.5)

    # ---- Row 4: Module 4 & 5 ----
    y_row4 = 6.0
    h_row4 = 2.4
    draw_module_row(y_row4, h_row4,
                    '四', '虚拟细胞模型验证', 'Pomc Neuron Virtual Cell',
                    ['- Pomc+神经元GRN (正常/SD/恢复)',
                     '- 肽段->受体->TF活性改变->GRN',
                     '  传播->全基因表达变化',
                     '- 计算逆转分数(Rescue Score)',
                     '  = -cos(peptide, SD)'],
                    '-> 逆转分数排序', C['purple'], '第 9-12 月',
                    '五', '体外实验初步验证', 'HT22 Oxidative Stress Model',
                    ['- H2O2诱导 (200uM, 24h)',
                     '- 肽段梯度处理 (10/50/100 uM)',
                     '- CCK-8 检测细胞活力',
                     '- ROS 荧光探针检测活性氧水平',
                     '- 活力>=20% + ROS<=30% -> 候选'],
                    '-> 3 条验证肽', C['coral'], '第 12-18 月')

    arrow_d(ax, MX-HW-GAP/2, y_row4, y_row4+h_row4+0.15, C['gray'], lw=2)
    arrow_d(ax, MX+HW+GAP/2, y_row4, y_row4+h_row4+0.15, C['gray'], lw=2)
    arrow_d(ax, MX, 5.65, 5.25, C['dark'], lw=2.5)

    # ---- Target signaling panel ----
    py = 3.5
    panel = FancyBboxPatch((0.8, py), 16.4, 1.6, boxstyle='round,pad=0.15',
                           facecolor=C['blue_bg'], edgecolor=C['blue_mid'], linewidth=1.5, zorder=0)
    ax.add_patch(panel)
    ax.text(MX, py+1.32, '五靶点协同信号传导网络验证  |  Multi-Target Signaling Network',
            fontsize=9.5, fontweight='bold', color=C['dark'], ha='center')

    sig_targets = [
        ('GABA_AR', 'Ca2+ down\nPKA down', 'CREB down', '镇静\nPomc down', C['blue']),
        ('OX2R', 'Ca2+ up\nMAPK up', 'CREB up', '觉醒 down\n(拮抗)', C['green']),
        ('MC4R', 'cAMP up\nPKA up', 'CREB/STAT3 up', 'Pomc up\nPcsk2 up', C['orange']),
        ('5-HT1A', 'cAMP down\nPKA down', 'CREB down', '节律调节', C['purple']),
        ('Keap1', 'Nrf2转位 up\nROS down', 'Nrf2/ARE up', 'HO-1 up\nSOD up', C['coral']),
    ]
    for i, (tn, tp, ttf, te, tc) in enumerate(sig_targets):
        tx = 1.1 + i * 3.3
        box_solid(ax, tx, py+0.55, 2.4, 0.5, tc, tn, fontsize=8)
        ax.text(tx+1.2, py+0.3, tp, fontsize=5.8, color=C['gray'], ha='center')
        box_outline(ax, tx, py-0.15, 2.4, 0.4, tc, ttf, fontsize=7)
        ax.text(tx+1.2, py-0.35, te, fontsize=5.8, color=tc, ha='center', fontweight='bold')

    arrow_d(ax, MX, 2.7, 2.3, C['dark'], lw=2.5)

    # ---- Output bar ----
    box_solid(ax, 2.5, 1.5, 13.0, 0.6, C['dark'],
              '产出: SCI 1-2 篇  |  助眠活性肽 3 条  |  AI筛选管道代码 1 套  |  肽段库 5,000+ 条',
              fontsize=9)

    # ---- Innovation badge ----
    ax.text(17.2, 15.5,
            '创新点\n======\n1 首次AI筛选\n  助眠活性肽\n2 虚拟细胞验\n  证药效\n3 "逆转分数"\n  新概念\n4 聚焦闽产\n  地域特色',
            fontsize=5.3, color=C['gray'], va='top', ha='left', linespacing=1.3,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=C['lgray']))

    plt.tight_layout(pad=0.3)
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.savefig(output_path, dpi=DPI, facecolor='white', edgecolor='none',
                bbox_inches='tight', pad_inches=0.15)
    plt.close(fig)
    print(f'[OK] TCM Roadmap: {output_path}')


# ===================================================================
if __name__ == '__main__':
    draw_nsfc('D:/sleep-deprivation-project/grant/mechanism_diagram_nsfc.png')
    draw_tcm('D:/药膳/mechanism_diagram_tcm.png')
    draw_nsfc_roadmap('D:/sleep-deprivation-project/grant/roadmap_nsfc.png')
    draw_tcm_roadmap('D:/药膳/roadmap_tcm.png')
    print('\nDone.')
