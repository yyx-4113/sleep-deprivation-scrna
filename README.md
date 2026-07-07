# Sleep Deprivation scRNA-seq Analysis

## 项目概述
利用公共数据库的单细胞转录组数据，分析睡眠剥夺对脑组织的转录影响。通过 GRN 推断识别关键转录因子，虚拟敲除预测其下游效应，多数据集验证结果。

## 数据来源

| Accession | Tissue | Method | Year |
|-----------|--------|--------|------|
| GSE137665 | Brainstem, Cortex, Hypothalamus | scRNA-seq (10x) | 2022 |
| GSE211088 | Frontal Cortex | snRNA-seq + bulk | 2023 |
| GSE214337 | Cerebral Cortex | snRNA-seq | 2022 |
| GSE243489 | Preoptic Hypothalamus | snRNA-seq (10x) | 2024 |

## 环境限制与策略

| 限制 | 原因 | 解决方案 |
|------|------|----------|
| 本机 8GB 内存 | R 加载多个大包会崩溃 (exit 139) | 重型计算全走 Colab |
| NCBI 无法访问 | SSL CRYPT_E_REVOCATION_OFFLINE | Colab 云端下载 |
| R 库路径混用 | 系统库与个人库 DLL 冲突 | `.libPaths` 只用个人库 |

## 工作流（Colab-first）

### 第一步：Google Colab（核心分析）
上传到 https://colab.research.google.com 按顺序运行：

1. **`Colab_Full_Pipeline.ipynb`** — 数据下载 → 预处理 → 聚类 → 差异表达 → 富集 → 出图
2. **`03_pySCENIC_CellOracle.ipynb`** — GRN 推断 → 虚拟敲除

### 第二步：本地 RStudio（精调出图，可选）
打开 `C:\Program Files\RStudio\rstudio.exe`：
- `00_install_packages.R` — 装剩余 Bioc 包（部分大包可能需 RStudio GUI）
- `01_seurat_analysis.R` — 备选方案（小数据集可在 RStudio 中跑）
- `02_diff_expression.R` — 本地差异分析（从 Colab 导出的表达矩阵）
- `04_external_validation.R` — 独立数据集验证
- `05_machine_learning.R` — ML 筛选标志基因

### 本地库路径注意
所有 R 脚本需在开头加：`.libPaths("C:/Users/1/R/library")`

## 目录结构

```
sleep-deprivation-project/
├── data/
├── scripts/
│   ├── 00_install_packages.R
│   ├── 00_Colab_Download_and_Preprocess.ipynb
│   ├── Colab_Full_Pipeline.ipynb        ★ 主文件
│   ├── 01_seurat_analysis.R
│   ├── 02_diff_expression.R
│   ├── 03_pySCENIC_CellOracle.ipynb     ★ 虚拟敲除
│   ├── 04_external_validation.R
│   └── 05_machine_learning.R
├── results/
│   └── figures/
└── README.md
```

## 睡眠剥夺可关注的转录因子

| TF | 睡眠剥夺角色 |
|----|------------|
| Nr3c1 (GR) | 应激核心，HPA 轴 |
| Crem | cAMP 响应，节律调控 |
| Per1/Per2 | 昼夜节律核心 |
| Fos/Jun | 神经元活动标志 |
| Npas4 | 突触可塑性 |
| Srebf1 | 脂代谢，SD 后上调 |
| Nr4a1/Nr4a2 | 孤核受体，应激 |
