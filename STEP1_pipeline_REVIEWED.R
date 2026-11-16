# =============================================================================
# STEP 1 (REVIEWED) — Cuproptosis-Related Genes (CRGs) panel in TCGA-LUAD
# Mentor revision: 23/05/2026
# Project:  Cuproptosis Prognostic Signature in LUAD
# Author:   <student name>     |  Mentor: Thao
#
# Mục tiêu của STEP 1:
#   1. Tải TCGA-LUAD RNA-seq + clinical data
#   2. Normalize bằng TMM + log2-CPM (theo edgeR best practice)
#   3. Trích 19 CRGs (Tsvetkov et al., Science 2022)
#   4. Sinh 4 figure cho Figure 1: 1A violin DE | 1B correlation | 1C UMAP | 1D heatmap+OS
#
# Quy ước:
#   - Output figure  -> outputs/Figures/
#   - Output table   -> outputs/Tables/
#   - Object trung gian -> outputs/Rdata/
# =============================================================================

# ---- 0. Thiết lập môi trường --------------------------------------------------
# R ≥ 4.1 yêu cầu cho native pipe |>
stopifnot(getRversion() >= "4.1.0")

suppressPackageStartupMessages({
  library(TCGAbiolinks)          # tải data từ GDC
  library(SummarizedExperiment)  # cấu trúc dữ liệu của TCGA
  library(edgeR)                 # TMM, filterByExpr, cpm
  library(limma)                 # chuẩn bị cho Step 2 (volcano)
  library(ggplot2)               # plotting
  library(pheatmap)              # heatmap
  library(corrplot)              # correlation matrix
  library(uwot)                  # UMAP
  library(survival)              # Surv, survdiff
  library(survminer)             # ggsurvplot (KM curve có pval, risk table)
  library(dplyr)                 # data manipulation
  library(tibble)                # rownames_to_column
  library(tidyr)                 # pivot_longer
  library(RColorBrewer)
})

# Tạo cây thư mục output ngay từ đầu, tránh xả file linh tinh
dir.create("outputs/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Rdata",   recursive = TRUE, showWarnings = FALSE)

set.seed(123)   # seed cho mọi randomness (kmeans, umap)

# ---- 1. Tải TCGA-LUAD --------------------------------------------------------
# STAR - Counts là dạng raw count, dùng cho edgeR/DESeq2
query_rna <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
# files.per.chunk thấp giúp ổn định khi mạng yếu; method="api" mặc định
GDCdownload(query_rna, method = "api", files.per.chunk = 20)
rna_data <- GDCprepare(query_rna)

clin <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")

# Lưu raw object để nếu chạy lại không phải download
saveRDS(rna_data, "outputs/Rdata/rna_data_raw.rds")
saveRDS(clin,     "outputs/Rdata/clin_raw.rds")

# ---- 2. Pre-processing & Normalization ---------------------------------------
# 2.1. Ma trận count thô (unstranded là default cho STAR)
counts_raw <- assay(rna_data, "unstranded")

# 2.2. Lọc gene & chuyển Ensembl ID -> gene symbol AN TOÀN
#     Vấn đề của make.unique(): nếu 2 Ensembl ID khác nhau cùng map về 1 symbol
#     (vd CDKN2A có nhiều transcript), gene thứ 2 sẽ thành "CDKN2A.1" và bị mất
#     khi match với cup_genes.
#     Giải pháp: aggregate counts theo gene_symbol (cộng tổng counts).

gene_symbols <- rowData(rna_data)$gene_name
gene_types   <- rowData(rna_data)$gene_type

# (Tuỳ chọn) chỉ giữ protein-coding để giảm noise (~20k gene)
keep_pc <- gene_types == "protein_coding" & !is.na(gene_symbols) & gene_symbols != ""
counts_pc    <- counts_raw[keep_pc, ]
symbols_pc   <- gene_symbols[keep_pc]

# Aggregate: nếu nhiều Ensembl ID cùng symbol thì cộng counts
counts <- rowsum(counts_pc, group = symbols_pc, na.rm = TRUE)
cat("Sau aggregate: ", nrow(counts), " unique gene symbols\n")

# 2.3. Xác định sample type (theo TCGA barcode position 14-15)
#     "01" = primary tumour | "11" = solid tissue normal | "02" = recurrent (sẽ loại)
sample_codes  <- substr(colnames(counts), 14, 15)
n_other       <- sum(!sample_codes %in% c("01", "11"))
cat("Primary tumour (01):", sum(sample_codes == "01"), "\n",
    "Normal       (11):", sum(sample_codes == "11"), "\n",
    "Other        :", n_other, " (recurrent/metastatic — loại bỏ)\n", sep = "")

keep_sample   <- sample_codes %in% c("01", "11")
counts        <- counts[, keep_sample]
sample_codes  <- sample_codes[keep_sample]

# 2.4. Lọc gene biểu hiện thấp TRƯỚC khi TMM (best practice edgeR)
group <- factor(ifelse(sample_codes == "01", "Tumor", "Normal"),
                levels = c("Normal", "Tumor"))
dge   <- DGEList(counts = counts, group = group)
keep  <- filterByExpr(dge, group = group)
dge   <- dge[keep, , keep.lib.sizes = FALSE]    # phải set FALSE để recompute lib size
cat("Sau filterByExpr: ", nrow(dge), " genes giữ lại / ", nrow(counts), "\n")

# 2.5. TMM normalization + log2-CPM
#     prior.count = 2 (mặc định, an toàn cho gene biểu hiện thấp)
dge      <- calcNormFactors(dge, method = "TMM")
expr_mat <- cpm(dge, log = TRUE, prior.count = 2)

# 2.6. Tách tumour-only và tumour+normal
tumor_samples  <- colnames(expr_mat)[sample_codes == "01"]
normal_samples <- colnames(expr_mat)[sample_codes == "11"]
expr_tumor     <- expr_mat[, tumor_samples]
expr_all       <- expr_mat[, c(tumor_samples, normal_samples)]

# Save normalized matrices
saveRDS(expr_mat,   "outputs/Rdata/expr_mat_logCPM.rds")
saveRDS(expr_tumor, "outputs/Rdata/expr_tumor_logCPM.rds")

# ---- 3. 19 CRG panel ---------------------------------------------------------
cup_genes <- c("FDX1","LIAS","LIPT1","LIPT2","DLD","DLAT","PDHA1","PDHB",
               "MTF1","GLS","CDKN2A","SLC31A1","ATP7A","ATP7B","ATOX1",
               "CCS","SCO1","SCO2","COX17")
activator_genes  <- c("FDX1","LIAS","LIPT1","LIPT2","DLD","DLAT",
                      "PDHA1","PDHB","MTF1")    # 9 activator chính
suppressor_genes <- setdiff(cup_genes, activator_genes)

# Kiểm tra coverage — báo ngay nếu thiếu gene nào
missing_genes <- setdiff(cup_genes, rownames(expr_mat))
if (length(missing_genes) > 0) {
  warning("Các CRG sau KHÔNG có trong expression matrix: ",
          paste(missing_genes, collapse = ", "))
}
cup_genes_present <- intersect(cup_genes, rownames(expr_mat))
cat("CRGs có trong data: ", length(cup_genes_present), " / 19\n")

expr_cup <- expr_mat[cup_genes_present, , drop = FALSE]

# ---- 4. Figure 1A — Tumour vs Normal (violin + boxplot + sig stars) ----------
plot_df <- as.data.frame(expr_cup) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "sample", values_to = "expr") |>
  mutate(type = ifelse(substr(sample, 14, 15) == "01", "Tumor", "Normal"))

# Wilcoxon + effect size (log2FC = mean_T - mean_N vì expr đã ở log2 scale)
pvals <- plot_df |>
  group_by(gene) |>
  summarise(
    p       = wilcox.test(expr ~ type, exact = FALSE)$p.value,
    mean_T  = mean(expr[type == "Tumor"]),
    mean_N  = mean(expr[type == "Normal"]),
    log2FC  = mean_T - mean_N,
    .groups = "drop"
  ) |>
  mutate(
    padj = p.adjust(p, method = "BH"),                   # BH (FDR) theo convention
    sig  = case_when(
      padj < 1e-4 ~ "****",
      padj < 1e-3 ~ "***",
      padj < 1e-2 ~ "**",
      padj < 5e-2 ~ "*",
      TRUE        ~ "ns"
    )
  )

write.csv(pvals, "outputs/Tables/Table_S1_DE_Tumor_vs_Normal.csv", row.names = FALSE)

# Vị trí label sig: tính per-gene để label nằm gọn trên đỉnh từng violin
label_pos <- plot_df |>
  group_by(gene) |>
  summarise(y = max(expr) + 0.5, .groups = "drop")
pvals_lab <- left_join(pvals, label_pos, by = "gene")

p1 <- ggplot(plot_df, aes(gene, expr, fill = type)) +
  geom_violin(trim = FALSE, alpha = 0.85, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white",
               alpha = 0.6, linewidth = 0.3,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(Tumor = "#e74c3c", Normal = "#3498db")) +
  geom_text(data = pvals_lab,
            aes(x = gene, y = y, label = sig),
            inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  scale_x_discrete(limits = cup_genes_present) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  labs(title = "Figure 1A — Differential expression of 19 CRGs in TCGA-LUAD",
       subtitle = sprintf("n_tumor = %d, n_normal = %d  |  Wilcoxon + BH-FDR",
                          length(tumor_samples), length(normal_samples)),
       y = expression(Expression~(log[2]~CPM)), x = NULL, fill = NULL)

ggsave("outputs/Figures/Fig1A_CRG_violin_TumourVsNormal.png",
       p1, width = 13, height = 5.5, dpi = 300)

# ---- 5. Figure 1B — Correlation heatmap (tumour only) ------------------------
tumor_mat <- t(expr_cup[, tumor_samples])   # samples × genes
cor_mat   <- cor(tumor_mat, method = "spearman")

png("outputs/Figures/Fig1B_CRG_correlation_heatmap.png",
    width = 11, height = 9, units = "in", res = 300)
corrplot(cor_mat, method = "color", type = "upper",
         tl.cex = 0.9, tl.col = "black",
         col = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(200),
         addCoef.col = "black", number.cex = 0.55,
         title = "Figure 1B — Spearman correlation among CRGs (tumour samples)",
         mar = c(0, 0, 2, 0))
dev.off()

# ---- 6. Figure 1C — UMAP + K-means clustering --------------------------------
#  - Cluster TRÊN ma trận scaled X (không trên UMAP embedding) — đúng practice
#  - Set seed riêng cho UMAP (uwot dùng RNG khác)
X  <- scale(tumor_mat)

# Đánh giá k tối ưu bằng silhouette (chỉ informational, k=2 đã định trước)
# library(factoextra); fviz_nbclust(X, kmeans, method = "silhouette", k.max = 6)

km <- kmeans(X, centers = 2, nstart = 25)

set.seed(123)
emb <- umap(X, n_neighbors = 20, min_dist = 0.25,
            metric = "euclidean", n_threads = 4)

# Gán label "high/low" cuproptosis dựa trên TB của activator genes
activator_present <- intersect(activator_genes, colnames(tumor_mat))
cluster_means     <- tapply(rowMeans(tumor_mat[, activator_present]), km$cluster, mean)
high_id           <- as.integer(names(which.max(cluster_means)))
labels            <- ifelse(km$cluster == high_id,
                            "C1 (Cuproptosis-high)",
                            "C2 (Cuproptosis-low)")

umap_df <- data.frame(UMAP1 = emb[, 1], UMAP2 = emb[, 2], cluster = labels)
p3 <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(size = 1.6, alpha = 0.78) +
  scale_color_manual(values = c("C1 (Cuproptosis-high)" = "#d7191c",
                                "C2 (Cuproptosis-low)"  = "#2c7bb6")) +
  theme_bw(base_size = 11) +
  labs(title = "Figure 1C — UMAP: two cuproptosis subtypes (k-means, k=2)",
       color = NULL)
ggsave("outputs/Figures/Fig1C_UMAP_subtypes.png",
       p3, width = 8, height = 7, dpi = 300)

# ---- 7. Survival: KM curve + stratified heatmap (Figure 1D) ------------------
clin_subset <- clin |>
  dplyr::select(submitter_id, days_to_last_follow_up, days_to_death, vital_status) |>
  mutate(
    OS_time  = ifelse(vital_status == "Dead",
                      as.numeric(days_to_death),
                      as.numeric(days_to_last_follow_up)),
    OS_event = ifelse(vital_status == "Dead", 1, 0)
  )

# Patient barcode (12 ký tự) để join với clinical
tumor_ids12 <- substr(tumor_samples, 1, 12)

surv_df <- data.frame(sample  = tumor_samples,
                      barcode = tumor_ids12,
                      cluster = labels) |>
  left_join(clin_subset, by = c("barcode" = "submitter_id"))

# Báo cáo missing clinical
n_missing_clin <- sum(is.na(surv_df$OS_time))
cat("Sample tumour không có OS_time:", n_missing_clin, "\n")

# 7a. KM curve (Figure 1D-left) — KHÔNG dùng median split kiểu "Long/Short OS"
surv_df_complete <- surv_df |> filter(!is.na(OS_time) & OS_time > 0)
fit_km <- survfit(Surv(OS_time, OS_event) ~ cluster, data = surv_df_complete)

p_km <- ggsurvplot(
  fit_km, data = surv_df_complete,
  pval = TRUE, conf.int = TRUE, risk.table = TRUE,
  palette = c("#d7191c", "#2c7bb6"),
  xlab = "Time (days)", ylab = "Overall survival",
  title = "Figure 1D (left) — KM survival by cuproptosis subtype",
  legend.title = "Subtype", risk.table.height = 0.28
)
ggsave("outputs/Figures/Fig1D_KM_curve.png",
       plot = p_km$plot, width = 7, height = 5.5, dpi = 300)
ggsave("outputs/Figures/Fig1D_KM_risk_table.png",
       plot = p_km$table, width = 7, height = 2, dpi = 300)

# Tính p log-rank để báo vào title heatmap
fit_lr <- survdiff(Surv(OS_time, OS_event) ~ cluster, data = surv_df_complete)
p_lr   <- 1 - pchisq(fit_lr$chisq, df = length(fit_lr$n) - 1)

# 7b. Stratified heatmap (Figure 1D-right)
expr_z <- t(scale(t(expr_cup[, tumor_samples])))

# Annotation: chỉ dùng cluster + vital status (KHÔNG dùng Long/Short OS giả tạo)
ann <- data.frame(
  Cluster = labels,
  Vital   = ifelse(surv_df$OS_event[match(tumor_samples, surv_df$sample)] == 1,
                   "Deceased", "Alive/Censored"),
  row.names = tumor_samples
)
ann$Vital[is.na(ann$Vital)] <- "Unknown"

ann_colors <- list(
  Cluster = c("C1 (Cuproptosis-high)" = "#d7191c",
              "C2 (Cuproptosis-low)"  = "#2c7bb6"),
  Vital   = c("Deceased" = "#34495e",
              "Alive/Censored" = "#2ecc71",
              "Unknown" = "grey80")
)

png("outputs/Figures/Fig1D_stratified_heatmap.png",
    width = 15, height = 7, units = "in", res = 300)
pheatmap(
  expr_z,
  annotation_col   = ann,
  annotation_colors = ann_colors,
  show_colnames    = FALSE,
  cluster_cols     = TRUE,
  cluster_rows     = TRUE,
  color            = colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(100),
  breaks           = seq(-2.5, 2.5, length.out = 101),
  main = sprintf("Figure 1D (right) — CRG expression by subtype  (log-rank p = %.3g)",
                 p_lr)
)
dev.off()

# ---- 8. Lưu objects cuối cùng + sessionInfo ---------------------------------
saveRDS(
  list(
    expr_mat   = expr_mat,
    expr_cup   = expr_cup,
    cup_genes_present = cup_genes_present,
    cluster    = labels,
    km         = km,
    umap_emb   = emb,
    pvals_DE   = pvals,
    surv_df    = surv_df,
    log_rank_p = p_lr
  ),
  "outputs/Rdata/STEP1_objects.rds"
)

writeLines(capture.output(sessionInfo()),
           "outputs/sessionInfo_STEP1.txt")

cat("\n========================================\n")
cat(" STEP 1 completed.\n")
cat(" Figures saved to outputs/Figures/\n")
cat(" Tables  saved to outputs/Tables/\n")
cat(" Log-rank p (C1 vs C2):", signif(p_lr, 4), "\n")
cat("========================================\n")
