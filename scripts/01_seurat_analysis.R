# ============================================================
# 01 - 单细胞 Seurat 标准分析流程
# 数据：GSE137665 - 睡眠剥夺小鼠多脑区 scRNA-seq
# ============================================================

library(Seurat)
library(tidyverse)
library(patchwork)
library(ggsci)

# -------------------- 0. 参数设置 --------------------
dir_data   <- "../data"
dir_result <- "../results"
dir_fig    <- "../results/figures"

# -------------------- 1. 下载数据 --------------------
# GSE137665 包含 brainstem, cortex, hypothalamus 三个区域
# 每个区域有 Control 和 Sleep Deprivation (SD) 两组

# 从 GEO 获取元数据
# 实际数据（count矩阵）需要从 GEO 的 supplementary files 下载
# 或从 SRA 下载 fastq 后重新比对

# 假设我们已下载 10x 格式的三个矩阵：
# data/GSE137665_brainstem/
# data/GSE137665_cortex/
# data/GSE137665_hypothalamus/

# -------------------- 2. 读取数据 --------------------
# 这里以脑干区域为例，其他区域同理
read_10x_data <- function(sample_dir, sample_name) {
  counts <- Read10X(data.dir = sample_dir)
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  return(obj)
}

# 示例（实际路径需根据下载后调整）
# brainstem_ctrl <- read_10x_data("data/brainstem_ctrl", "Brainstem_Control")
# brainstem_sd   <- read_10x_data("data/brainstem_sd",   "Brainstem_SD")

# 如果你用的是从 GEO 下载的 h5 文件：
# brainstem_ctrl <- Read10X_h5("data/brainstem_ctrl.h5")

# -------------------- 3. 合并样本 --------------------
# merged <- merge(brainstem_ctrl, y = brainstem_sd,
#                 add.cell.ids = c("ctrl", "sd"))

# -------------------- 4. 质控 (QC) --------------------
qc_and_filter <- function(obj) {
  # 计算线粒体基因比例
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  # 计算核糖体基因比例
  obj[["percent.ribo"]] <- PercentageFeatureSet(obj, pattern = "^Rp[sl]")

  # QC 过滤（参数需根据数据调整）
  obj <- subset(obj,
    nFeature_RNA > 200  & nFeature_RNA < 6000 &
    nCount_RNA   > 500  & nCount_RNA   < 30000 &
    percent.mt   < 20
  )
  return(obj)
}

# merged <- qc_and_filter(merged)

# -------------------- 5. 标准化与降维 --------------------
process_seurat <- function(obj) {
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = 2000)
  obj <- ScaleData(obj, vars.to.regress = c("percent.mt"))
  obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 30)
  obj <- RunUMAP(obj, dims = 1:20)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.5)
  return(obj)
}

# merged <- process_seurat(merged)

# -------------------- 6. 细胞注释 --------------------
library(SingleR)
library(celldex)

annotate_cells <- function(obj) {
  # 使用小鼠大脑参考数据集
  ref <- MouseRNAseqData()  # Immunogen 参考数据
  # 或者用 Allen Brain Atlas 参考

  pred <- SingleR(
    test  = GetAssayData(obj, layer = "data"),
    ref   = ref,
    labels = ref$label.main,
    clusters = obj$seurat_clusters
  )

  # 把注释映射回 Seurat 对象
  idx <- match(obj$seurat_clusters, rownames(pred))
  obj$cell_type <- pred$labels[idx]

  return(obj)
}

# merged <- annotate_cells(merged)

# -------------------- 7. 找出各群 marker 基因 --------------------
# markers <- FindAllMarkers(merged, only.pos = TRUE,
#                           logfc.threshold = 0.25,
#                           min.pct = 0.25)
# write.csv(markers, file.path(dir_result, "all_markers.csv"))

# -------------------- 8. 绑图 --------------------
# UMAP 聚类图
plot_umap_clusters <- function(obj) {
  DimPlot(obj, group.by = "seurat_clusters", label = TRUE,
          repel = TRUE, pt.size = 0.5) +
    scale_color_igv() +
    ggtitle("UMAP by Clusters") +
    theme_minimal()
}

# UMAP 分组图（对照 vs 睡眠剥夺）
plot_umap_groups <- function(obj) {
  DimPlot(obj, group.by = "orig.ident", pt.size = 0.5) +
    scale_color_npg() +
    ggtitle("UMAP by Condition") +
    theme_minimal()
}

# UMAP 细胞类型图
plot_umap_celltypes <- function(obj) {
  DimPlot(obj, group.by = "cell_type", pt.size = 0.5,
          repel = TRUE, label = FALSE) +
    scale_color_igv() +
    ggtitle("UMAP by Cell Type") +
    theme_minimal()
}

# 顶级 marker 基因热图
plot_marker_heatmap <- function(obj, markers, top_n = 10) {
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(n = top_n, order_by = avg_log2FC)

  DoHeatmap(obj, features = unique(top_markers$gene)) +
    scale_fill_gradientn(colors = c("navy", "white", "firebrick3"))
}

# -------------------- 9. 保存对象 --------------------
# saveRDS(merged, file.path(dir_data, "seurat_merged.rds"))

cat("Seurat 分析流程模板加载完毕。\n")
cat("请根据实际下载的数据修改文件路径和参数。\n")
