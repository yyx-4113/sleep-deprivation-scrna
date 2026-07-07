# ============================================================
# 05 - 机器学习筛选睡眠剥夺关键标志基因
# 方法：LASSO + SVM-RFE + 随机森林
# 输出：一组可发表的高可信度基因 / TF
# ============================================================

library(tidyverse)
library(glmnet)      # LASSO
library(randomForest) # 随机森林
library(pROC)         # ROC 评估
library(caret)        # 交叉验证

dir_result <- "../results"
dir_fig    <- "../results/figures"

# -------------------- 1. 准备输入 --------------------
# 从 Seurat 伪批量表达矩阵开始
# 或从差异表达结果中取 top 差异基因作为特征

# 伪批量：按样本聚合成 pseudo-bulk
# pseudo <- AggregateExpression(obj, group.by = c("orig.ident", "cell_type"),
#                                assays = "RNA", slot = "data")
# expr_mat <- as.data.frame(pseudo$RNA)

# 标准化后用于机器学习
# 格式：行 = 样本（n个对照 + n个SD），列 = 特征基因

# -------------------- 2. LASSO 回归 --------------------
run_lasso <- function(expr_mat, group_labels) {
  set.seed(42)
  x <- t(expr_mat)  # glmnet 要求 样本×特征
  y <- ifelse(group_labels == "SD", 1, 0)

  cv <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 10)

  # 最优 lambda (1se)
  coefs <- coef(cv, s = "lambda.1se")
  selected <- rownames(coefs)[which(coefs[,1] != 0)][-1]  # 去掉 intercept

  return(list(
    model = cv,
    genes = selected,
    lambda_min = cv$lambda.min,
    lambda_1se  = cv$lambda.1se
  ))
}

# -------------------- 3. SVM-RFE (递归特征消除) --------------------
run_svm_rfe <- function(expr_mat, group_labels, size = 30) {
  x <- as.data.frame(t(expr_mat))
  y <- factor(group_labels)

  # 使用 caret 的 RFE
  ctrl <- rfeControl(
    functions = caretFuncs,
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    verbose = FALSE
  )

  subsets <- c(5, 10, 15, 20, 30, 50)

  rfe_res <- rfe(x, y, sizes = subsets,
                 rfeControl = ctrl,
                 method = "svmLinear")

  return(rfe_res)
}

# -------------------- 4. 随机森林 --------------------
run_rf <- function(expr_mat, group_labels, ntree = 1000) {
  x <- as.data.frame(t(expr_mat))
  y <- factor(group_labels)

  rf <- randomForest(x, y, ntree = ntree, importance = TRUE)
  imp <- importance(rf, type = 1)  # MeanDecreaseAccuracy
  imp <- imp[order(imp[,1], decreasing = TRUE), , drop = FALSE]

  return(list(model = rf, importance = imp))
}

# -------------------- 5. 多方法交集 --------------------
# 取 LASSO + RF + SVM-RFE 的交集
get_consensus_genes <- function(lasso_genes, rf_top_genes, svm_genes,
                                rf_top_n = 30) {
  rf_top <- names(rf_top_genes)[1:rf_top_n]
  intersect(intersect(lasso_genes, rf_top), svm_genes)
}

# -------------------- 6. ROC 评估 --------------------
evaluate_auc <- function(expr_mat, group_labels, selected_genes) {
  # 用选出的基因建一个多基因评分
  expr_sel <- expr_mat[selected_genes, , drop = FALSE]
  scores <- colMeans(expr_sel)  # 简单均值

  roc_res <- roc(group_labels, scores)
  return(roc_res)
}

# -------------------- 7. 可视化 --------------------

# LASSO 路径图
plot_lasso_cv <- function(cv_model) {
  plot(cv_model)
}

# ROC 曲线
plot_roc <- function(roc_obj) {
  ggroc(roc_obj, color = "#E64B35", size = 1.2) +
    geom_abline(linetype = "dashed", alpha = 0.3) +
    theme_minimal(base_size = 12) +
    annotate("text", x = 0.3, y = 0.2,
             label = paste0("AUC = ", round(auc(roc_obj), 3))) +
    labs(title = "ROC Curve", x = "1 - Specificity", y = "Sensitivity")
}

# 特征重要性热图（RF）
plot_rf_importance <- function(rf_imp, top_n = 30) {
  df <- data.frame(gene = rownames(rf_imp), importance = rf_imp[,1])
  df %>%
    slice_max(importance, n = top_n) %>%
    ggplot(aes(importance, fct_reorder(gene, importance))) +
    geom_col(fill = "#4DBBD5") +
    theme_minimal(base_size = 10) +
    labs(x = "Mean Decrease Accuracy", y = "", title = "Random Forest Importance")
}

cat("机器学习筛选标志基因脚本加载完毕。\n")
cat("建议：至少用 2-3 个独立数据集的验证 ROC > 0.7 才有说服力。\n")
