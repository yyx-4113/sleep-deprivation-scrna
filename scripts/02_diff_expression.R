# ============================================================
# 02 - 差异表达 + GO/KEGG/GSEA 富集分析
# 比较：睡眠剥夺 vs 对照，各脑区 × 各细胞类型
# ============================================================

library(Seurat)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(limma)
library(EnhancedVolcano)
library(enrichplot)
library(msigdbr)
library(fgsea)
library(ggpubr)

dir_result <- "../results"
dir_fig    <- "../results/figures"

# -------------------- 1. 加载数据 --------------------
# obj <- readRDS("../data/seurat_merged.rds")

# -------------------- 2. 比较设置 --------------------
# 按 [脑区 × 细胞类型] 做差异表达

run_de_analysis <- function(obj, cell_type, condition_col = "orig.ident",
                            ident.1 = "SD", ident.2 = "Control") {

  # 取该细胞类型的子集
  sub <- subset(obj, subset = cell_type == !!cell_type)
  Idents(sub) <- condition_col

  # limma 差异表达
  de <- FindMarkers(sub, ident.1 = ident.1, ident.2 = ident.2,
                    test.use = "MAST",    # 单细胞推荐用 MAST 或 wilcox
                    logfc.threshold = 0,
                    min.pct = 0.1)
  de$gene <- rownames(de)
  return(de)
}

# -------------------- 3. 批量差异表达 --------------------
# cell_types <- unique(obj$cell_type)
# de_list <- map(cell_types, ~ run_de_analysis(obj, .x))
# names(de_list) <- cell_types

# -------------------- 4. GO 富集分析 --------------------
run_go_enrichment <- function(de_genes, direction = "up",
                              logfc_cut = 0.25, pval_cut = 0.05) {
  if (direction == "up") {
    genes <- de_genes %>% filter(avg_log2FC > logfc_cut, p_val_adj < pval_cut)
  } else {
    genes <- de_genes %>% filter(avg_log2FC < -logfc_cut, p_val_adj < pval_cut)
  }
  genes <- genes$gene

  ego <- enrichGO(
    gene          = genes,
    OrgDb         = org.Mm.eg.db,
    keyType       = "SYMBOL",
    ont           = "BP",           # Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  return(simplify(ego, cutoff = 0.7))  # 去掉冗余 GO term
}

# -------------------- 5. KEGG 富集分析 --------------------
run_kegg_enrichment <- function(de_genes, logfc_cut = 0.25, pval_cut = 0.05) {
  # 先把 SYMBOL 转 ENTREZID
  genes <- de_genes %>%
    filter(abs(avg_log2FC) > logfc_cut, p_val_adj < pval_cut) %>%
    pull(gene)

  entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)

  ekegg <- enrichKEGG(
    gene         = entrez$ENTREZID,
    organism     = "mmu",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05
  )
  return(ekegg)
}

# -------------------- 6. GSEA (不截断阈值，用全部基因) --------------------
run_gsea <- function(de_genes, msigdb_sets = "H") {
  # 准备 ranked gene list
  de_genes <- de_genes %>% filter(!is.na(avg_log2FC))
  ranks <- de_genes$avg_log2FC
  names(ranks) <- de_genes$gene
  ranks <- sort(ranks, decreasing = TRUE)

  # 取 MSigDB Hallmark 基因集
  m_df <- msigdbr(species = "Mus musculus", category = msigdb_sets) %>%
    dplyr::select(gs_name, gene_symbol)

  fgsea_res <- fgsea(
    pathways = split(m_df$gene_symbol, m_df$gs_name),
    stats    = ranks,
    minSize  = 15,
    maxSize  = 500
  )
  return(fgsea_res %>% arrange(desc(NES)))
}

# -------------------- 7. 绑图函数 --------------------

# 火山图
plot_volcano <- function(de_genes, cell_type) {
  de_genes %>%
    mutate(sig = case_when(
      avg_log2FC > 0.25  & p_val_adj < 0.05 ~ "Up",
      avg_log2FC < -0.25 & p_val_adj < 0.05 ~ "Down",
      TRUE ~ "NS"
    )) %>%
    ggplot(aes(avg_log2FC, -log10(p_val_adj), color = sig)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey70")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    theme_minimal(base_size = 12) +
    labs(title = paste0(cell_type, ": SD vs Control"),
         x = "log2 Fold Change", y = "-log10(padj)") +
    theme(legend.position = "top")
}

# 气泡图（GO 富集）
plot_go_bubble <- function(ego, n = 20) {
  ego@result %>%
    slice_min(p.adjust, n = n) %>%
    mutate(gene_ratio = parse_ratio(GeneRatio)) %>%
    ggplot(aes(gene_ratio, fct_reorder(Description, gene_ratio),
               size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "#E64B35", high = "grey70") +
    scale_size(range = c(2, 8)) +
    theme_minimal(base_size = 10) +
    labs(x = "Gene Ratio", y = "", color = "p.adjust", size = "Count")
}

# GSEA 山脊图
plot_gsea_ridge <- function(fgsea_res, n = 15) {
  fgsea_res %>%
    filter(padj < 0.05) %>%
    slice_max(abs(NES), n = n) %>%
    ggplot(aes(NES, fct_reorder(pathway, NES), fill = -log10(padj))) +
    geom_col() +
    scale_fill_gradient(low = "grey80", high = "#E64B35") +
    theme_minimal(base_size = 10) +
    labs(y = "", fill = "-log10(padj)")
}

# -------------------- 8. 跨细胞类型比较：扰动响应得分 --------------------
# 计算每个细胞类型对睡眠剥夺的响应程度

calc_response_score <- function(de_list) {
  map_dfr(names(de_list), function(ct) {
    de <- de_list[[ct]]
    tibble(
      cell_type = ct,
      n_up   = sum(de$avg_log2FC > 0.25 & de$p_val_adj < 0.05),
      n_down = sum(de$avg_log2FC < -0.25 & de$p_val_adj < 0.05),
      total_de = n_up + n_down,
      response_score = n_up - n_down  # 净上调信号
    )
  })
}

# plot_response_scores <- function(score_df) {
#   score_df %>%
#     pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "n") %>%
#     ggplot(aes(fct_reorder(cell_type, total_de), n, fill = direction)) +
#     geom_col() +
#     coord_flip() +
#     scale_fill_npg() +
#     theme_minimal()
# }

cat("差异表达 + 富集分析脚本加载完毕。\n")
