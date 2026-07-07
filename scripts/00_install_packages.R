# ============================================================
# 一键安装所有需要的 R 包
# 运行方式：在 RStudio 中打开此文件，Ctrl+A 全选，Ctrl+Enter 运行
# 重要：只用个人库路径，不要混用系统库（混用会导致 DLL 冲突崩溃）
# ============================================================

# 设置个人库路径
.libPaths("C:/Users/1/R/library")

# --- 先装 BiocManager（Bioconductor 的包管理器）---
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", lib = "C:/Users/1/R/library")

# --- CRAN 包 ---
cran_pkgs <- c(
  "tidyverse",      # 数据处理 + 绑图（ggplot2/dplyr/tidyr）
  "Seurat",         # 单细胞核心包
  "patchwork",      # 拼图
  "ggsci",          # SCI 期刊配色
  "ggpubr",         # 发表级统计图
  "ggplot2",        # 绑图核心
  "pheatmap",       # 热图
  "RColorBrewer",   # 配色
  "viridis",        # 色盲友好配色
  "ggrepel",        # 标签防重叠
  "survival",       # 生存分析
  "survminer",      # 生存曲线绑图
  "glmnet",         # LASSO 回归
  "pROC",           # ROC 曲线
  "rms",            # Nomogram
  "forestplot",     # 森林图
  "corrplot",       # 相关性图
  "circlize",       # 圈图
  "enrichplot",     # 富集图
  "DOSE",           # 疾病本体富集
  "msigdbr",        # MSigDB 基因集
  "WGCNA",          # 加权基因共表达网络
  "ggalluvial"      # 桑基图/冲积图
)

# 安装 CRAN 包
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# --- Bioconductor 包 ---
bioc_pkgs <- c(
  "GEOquery",         # 下载 GEO 数据
  "limma",            # 差异表达分析
  "DESeq2",           # 差异表达分析（count数据）
  "edgeR",            # 差异表达分析
  "clusterProfiler",  # GO/KEGG 富集分析
  "org.Mm.eg.db",     # 小鼠基因注释
  "org.Hs.eg.db",     # 人类基因注释
  "ComplexHeatmap",   # 高级热图
  "SingleR",          # 单细胞自动注释
  "celldex",          # 参考表达矩阵
  "AUCell",           # 基因集活性评分
  "GSVA",             # 基因集变异分析
  "sva",              # 批次校正
  "EnhancedVolcano",  # 火山图
  "fgsea"             # 快速 GSEA
)

BiocManager::install(bioc_pkgs, ask = FALSE)

cat("\n========================\n所有包安装完成！\n========================\n")
