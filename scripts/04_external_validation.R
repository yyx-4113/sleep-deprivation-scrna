# ============================================================
# 04 - 外部验证：独立 GEO 数据集 + 人类相关性
# ============================================================
# 策略：
#   1. GEO 独立数据集验证（GSE211088, GSE214337）
#   2. 人类同源基因比对（小鼠 TF → 人类直系同源基因）
#   3. 文献挖掘已知睡眠基因集做交集
#   4. 可用 GTEx 人脑表达做跨物种验证

library(GEOquery)
library(tidyverse)
library(limma)
library(ggpubr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(org.Hs.eg.db)

dir_data   <- "../data"
dir_result <- "../results"

# -------------------- 1. 下载独立验证数据集 --------------------
# GSE211088: snRNA-seq frontal cortex, 睡眠剥夺
# GSE214337: snRNA-seq cerebral cortex, 睡眠需求

download_geo <- function(gse_id) {
  gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
  return(gse)
}

# gse211088 <- download_geo("GSE211088")
# gse214337 <- download_geo("GSE214337")

# -------------------- 2. 提取表达矩阵 --------------------
extract_expr_and_pheno <- function(gse_obj) {
  expr <- exprs(gse_obj[[1]])
  pheno <- pData(gse_obj[[1]])
  return(list(expr = expr, pheno = pheno))
}

# -------------------- 3. 验证我们找到的关键基因 --------------------
# 从 Seurat 和 pySCENIC 结果中提取关键 TF 列表
validate_genes <- function(key_genes, validation_expr, pheno, group_col) {

  # 取交集（有的基因可能在验证集中没有）
  common_genes <- intersect(key_genes, rownames(validation_expr))

  results <- map_dfr(common_genes, function(gene) {
    # 对每个基因做 t 检验（SD vs Control）
    vals <- as.numeric(validation_expr[gene, ])
    groups <- pheno[[group_col]]

    test <- t.test(vals ~ groups)
    fc <- log2(mean(vals[groups == unique(groups)[1]]) /
               mean(vals[groups == unique(groups)[2]]))

    tibble(gene = gene, log2FC = fc, pval = test$p.value)
  })

  results$padj <- p.adjust(results$pval, method = "BH")
  return(results %>% arrange(padj))
}

# -------------------- 4. 小鼠基因 → 人类同源基因 --------------------
mouse_to_human <- function(mouse_genes) {
  homologs <- bitr(mouse_genes,
    fromType = "SYMBOL", toType = c("SYMBOL", "HGNC_SYMBOL"),
    OrgDb = org.Mm.eg.db
  )
  # 更精确的同源映射（1:1 ortholog）
  # 可用 biomaRt 做更准确的映射
  return(homologs)
}

# 如果要做 GTEx 人脑表达验证：
# library(biomaRt)
# human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# mouse  <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
# orthologs <- getLDS(
#   attributes = c("mgi_symbol"),
#   filters = "mgi_symbol", values = mouse_genes,
#   mart = mouse,
#   attributesL = c("hgnc_symbol"),
#   martL = human
# )

# -------------------- 5. 文献已知睡眠基因集 --------------------
sleep_genes <- c(
  # 昼夜节律
  "Clock", "Bmal1", "Per1", "Per2", "Per3", "Cry1", "Cry2",
  "Nr1d1", "Nr1d2", "Rora", "Rorb",
  # 突触可塑性
  "Arc", "Bdnf", "Fos", "Egr1", "Npas4", "Homer1", "Dlg4",
  # 应激/代谢
  "Nr3c1", "Crem", "Srebf1", "Ppargc1a", "Foxo1",
  # 免疫/胶质
  "Tnf", "Il1b", "Aif1", "Gfap", "Trem2",
  # Homer1 信号
  "Homer2", "Homer3", "Shank1", "Shank3"
)

# 把我们从 scRNA-seq 中找到的差异基因跟已知睡眠基因做交集
# sleep_overlap <- intersect(de_sig_genes, sleep_genes)

# -------------------- 6. 可视化验证结果 --------------------
plot_validation <- function(val_results) {
  val_results %>%
    filter(padj < 0.05) %>%
    mutate(sig = ifelse(log2FC > 0, "Up in SD", "Down in SD")) %>%
    ggplot(aes(log2FC, -log10(padj), color = sig, label = gene)) +
    geom_point(aes(size = abs(log2FC)), alpha = 0.7) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    scale_color_manual(values = c("Up in SD" = "#E64B35",
                                  "Down in SD" = "#4DBBD5")) +
    scale_size(range = c(1, 6)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal(base_size = 12) +
    labs(title = "Independent Dataset Validation", x = "log2FC", y = "-log10(adj.p)")
}

# -------------------- 7. 多数据集一致性评分 --------------------
# 计算我们发现的基因在各独立数据集中的效应方向是否一致
calc_consistency <- function(results_list) {
  # results_list: 命名列表，每个元素是一个数据集的验证结果
  bind_rows(results_list, .id = "dataset") %>%
    filter(padj < 0.05) %>%
    group_by(gene) %>%
    summarise(
      n_datasets = n(),
      sign_consistent = n_distinct(sign(log2FC)) == 1,
      mean_log2FC = mean(log2FC)
    ) %>%
    filter(n_datasets >= 2, sign_consistent) %>%
    arrange(desc(n_datasets), desc(abs(mean_log2FC)))
}

cat("外部验证脚本加载完毕。\n")
cat("关键点：至少用 2 个独立数据集验证，方向一致才有说服力。\n")
