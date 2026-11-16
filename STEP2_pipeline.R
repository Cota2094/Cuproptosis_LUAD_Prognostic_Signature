# =============================================================================
# STEP 2 (REVIEWED) — Cuproptosis prognostic SIGNATURE in TCGA-LUAD
# Mentor revision: 23/06/2026
# Project:  Cuproptosis Prognostic Signature in LUAD
# Author:   <student name>     |  Mentor: Thao
#
# Mục tiêu của STEP 2 (toàn bộ phần xây dựng signature):
#   1. Univariate Cox screening 19 CRGs  -> Figure 2A forest plot + Table_S2
#   2. LASSO-Cox để chọn gene & tránh overfit -> Figure 2B (path) + 2C (CV)
#   3. Multivariate Cox của gene được chọn -> Table_S3
#   4. Risk score, chia high/low, KM + time-dependent ROC -> Fig 2D, 2E
#   5. Lưu objects + sessionInfo
#
# Đầu vào: outputs/Rdata/STEP1_objects.rds  (sinh ra ở STEP 1)
#          -> chứa: expr_mat, expr_cup, cup_genes_present, surv_df, ...
#
# Quy ước (giữ nguyên từ STEP 1):
#   - Output figure  -> outputs/Figures/
#   - Output table   -> outputs/Tables/
#   - Object trung gian -> outputs/Rdata/
#   - Multiple testing: BH-FDR (consistent với STEP 1)
#   - Survival endpoint: Overall Survival (OS)
# =============================================================================

# ---- 0. Thiết lập môi trường --------------------------------------------------
stopifnot(getRversion() >= "4.1.0")

suppressPackageStartupMessages({
  library(survival)     # coxph, Surv, cox.zph
  library(survminer)    # ggsurvplot, surv_cutpoint (nếu cần)
  library(glmnet)       # LASSO-Cox (family = "cox")
  library(timeROC)      # time-dependent ROC (AUC 1/3/5 năm)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(RColorBrewer)
})

dir.create("outputs/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Rdata",   recursive = TRUE, showWarnings = FALSE)

set.seed(123)   # seed toàn cục; LASSO-CV sẽ set seed riêng trước cv.glmnet

# Tham số có thể chỉnh — khai báo 1 chỗ, không "magic number" rải rác
ALPHA_SIG    <- 0.05   # ngưỡng "prognostic" sau FDR để báo cáo
CANDIDATE_P  <- 0.10   # ngưỡng univariate (raw p) để đưa gene vào LASSO
LAMBDA_RULE  <- "lambda.min"   # hoặc "lambda.1se" (parsimonious hơn)
ROC_TIMES    <- c(365, 1095, 1825)   # 1, 3, 5 năm (đơn vị: ngày)

# ---- 1. Nạp objects của STEP 1 ----------------------------------------------
step1 <- readRDS("outputs/Rdata/STEP1_objects.rds")

expr_cup <- step1$expr_cup               # 19 CRG × tất cả sample (log2-CPM)
surv_df  <- step1$surv_df                # đã có sample, barcode, OS_time, OS_event
cup_genes_present <- step1$cup_genes_present

# Chỉ giữ tumour (primary, "01") — surv_df ở STEP1 đã là tumour-only
tumor_samples <- surv_df$sample
expr_tumor    <- expr_cup[, tumor_samples, drop = FALSE]

# ---- 2. Làm sạch dữ liệu sống còn -------------------------------------------
# 2.1. Loại sample thiếu OS hoặc OS_time <= 0 (không hợp lệ cho Cox)
surv_clean <- surv_df |>
  filter(!is.na(OS_time), !is.na(OS_event), OS_time > 0)

# 2.2. QUAN TRỌNG: 1 bệnh nhân / 1 mẫu.
#      TCGA-LUAD có vài bệnh nhân được sequencing >1 lần (vial A/B) -> barcode
#      bệnh nhân (12 ký tự) trùng -> nếu không gộp sẽ vi phạm giả định độc lập
#      của Cox. Giữ mẫu đầu tiên theo thứ tự barcode.
dup_n <- sum(duplicated(surv_clean$barcode))
if (dup_n > 0) {
  cat("Phát hiện", dup_n, "mẫu trùng bệnh nhân -> giữ 1 mẫu/bệnh nhân\n")
}
surv_clean <- surv_clean |>
  arrange(barcode, sample) |>
  distinct(barcode, .keep_all = TRUE)

# 2.3. Đồng bộ ma trận biểu hiện với danh sách mẫu cuối
keep_samp  <- surv_clean$sample
expr_use   <- expr_tumor[, keep_samp, drop = FALSE]

cat("Số bệnh nhân vào phân tích survival:", ncol(expr_use), "\n")

# 2.4. Z-score TỪNG GENE (mean 0, sd 1) trên các tumour
#      Lý do: Cox trên biến liên tục cho HR theo "1 đơn vị" của biến.
#      log2-CPM mỗi gene có dynamic range khác nhau -> HR KHÔNG so sánh được
#      giữa các gene. Z-score giúp HR = "thay đổi nguy cơ trên mỗi 1 SD tăng
#      biểu hiện" -> forest plot mới có ý nghĩa khi đặt cạnh nhau.
expr_z <- t(scale(t(expr_use)))          # gene × sample, mỗi gene chuẩn hoá
stopifnot(!any(is.na(expr_z)))           # nếu có gene sd=0 sẽ ra NaN -> dừng

# Bảng phân tích: hàng = bệnh nhân, cột = gene (đã z-score) + OS
cox_df <- as.data.frame(t(expr_z))
cox_df$OS_time  <- surv_clean$OS_time
cox_df$OS_event <- surv_clean$OS_event

# ---- 3. Univariate Cox regression cho 19 CRGs (Figure 2A) -------------------
uni_cox <- lapply(cup_genes_present, function(g) {
  # Dùng formula theo tên gene (backtick an toàn) để cox.zph tái lập đúng model frame
  form <- reformulate(sprintf("`%s`", g), response = "Surv(OS_time, OS_event)")
  fit  <- coxph(form, data = cox_df)
  sm   <- summary(fit)
  # Kiểm tra giả định proportional hazards (Schoenfeld residuals)
  zph_p <- tryCatch(cox.zph(fit)$table[1, "p"], error = function(e) NA_real_)
  data.frame(
    gene   = g,
    HR     = sm$coefficients[1, "exp(coef)"],
    HR_low = sm$conf.int[1, "lower .95"],
    HR_high= sm$conf.int[1, "upper .95"],
    z      = sm$coefficients[1, "z"],
    p      = sm$coefficients[1, "Pr(>|z|)"],
    ph_p   = zph_p,                       # ph_p < 0.05 => vi phạm PH, cần lưu ý
    stringsAsFactors = FALSE
  )
}) |> bind_rows()

# BH-FDR cho 19 test (consistent với STEP 1)
uni_cox <- uni_cox |>
  mutate(
    padj      = p.adjust(p, method = "BH"),
    direction = ifelse(HR > 1, "Risk (HR > 1)", "Protective (HR < 1)"),
    sig_raw   = p    < 0.05,
    sig_fdr   = padj < ALPHA_SIG
  ) |>
  arrange(desc(HR))

n_sig_raw <- sum(uni_cox$sig_raw)
n_sig_fdr <- sum(uni_cox$sig_fdr)

write.csv(uni_cox, "outputs/Tables/Table_S2_univariate_cox_results.csv",
          row.names = FALSE)

cat(sprintf("Univariate Cox: %d/19 có raw p<0.05 | %d/19 sống sót BH-FDR<%.2f\n",
            n_sig_raw, length(cup_genes_present), n_sig_fdr, ALPHA_SIG))

# ---- 4. Forest plot Figure 2A ------------------------------------------------
# Solid = sống sót FDR (báo cáo trung thực), nhạt = phần còn lại.
fp_df <- uni_cox |>
  mutate(
    gene  = factor(gene, levels = rev(gene)),     # HR cao trên cùng
    label = sprintf("%.2f (%.2f–%.2f)", HR, HR_low, HR_high),
    plab  = ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)),
    alpha = ifelse(sig_fdr, 1, 0.45)
  )

x_text <- max(fp_df$HR_high) * 1.35
p_fig2a <- ggplot(fp_df, aes(x = HR, y = gene, color = direction, alpha = alpha)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = HR_low, xmax = HR_high), height = 0.25, linewidth = 0.7) +
  geom_point(shape = 18, size = 3.4) +
  geom_text(aes(x = x_text, label = label), color = "grey30",
            hjust = 0, size = 3, show.legend = FALSE) +
  geom_text(aes(x = x_text * 1.45, label = plab),
            color = "grey45", hjust = 0, size = 3, show.legend = FALSE) +
  scale_x_log10(limits = c(min(fp_df$HR_low) * 0.9, x_text * 1.9),
                breaks = c(0.3, 0.5, 0.75, 1, 1.5, 2, 3)) +
  scale_color_manual(values = c("Protective (HR < 1)" = "#2c7bb6",
                                "Risk (HR > 1)"        = "#c0392b")) +
  scale_alpha_identity() +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position    = "top",
        axis.text.y        = element_text(face = "bold")) +
  labs(
    title = "Figure 2A — Univariate Cox Regression: Prognostic CRGs in TCGA-LUAD",
    subtitle = sprintf(
      "n = %d patients | per 1-SD increase | %d/%d raw p<0.05, %d/%d BH-FDR<%.2f",
      nrow(cox_df), n_sig_raw, length(cup_genes_present),
      n_sig_fdr, length(cup_genes_present), ALPHA_SIG),
    x = "Hazard Ratio (95% CI, log scale)", y = NULL, color = NULL)

ggsave("outputs/Figures/Fig2A_univariate_cox_forest.png",
       p_fig2a, width = 11, height = 9, dpi = 300)

# ---- 5. LASSO-Cox: chọn gene & tránh overfit (Figure 2B, 2C) -----------------
# Đưa gene CANDIDATE (raw p < CANDIDATE_P) vào LASSO. Nếu quá ít gene survive,
# nới ngưỡng để LASSO còn nguyên liệu — và GHI RÕ điều này trong methods.
candidate_genes <- uni_cox$gene[uni_cox$p < CANDIDATE_P]
if (length(candidate_genes) < 2) {
  warning("Quá ít candidate (", length(candidate_genes),
          ") — dùng toàn bộ 19 CRG cho LASSO và nêu rõ tính exploratory.")
  candidate_genes <- cup_genes_present
}
cat("Số gene đưa vào LASSO:", length(candidate_genes), "\n")

x_mat <- t(expr_z[candidate_genes, , drop = FALSE])     # sample × gene
y_surv <- Surv(cox_df$OS_time, cox_df$OS_event)

set.seed(123)   # seed cho cross-validation (cv.glmnet dùng RNG)
cvfit <- cv.glmnet(x_mat, y_surv, family = "cox",
                   alpha = 1, nfolds = 10, type.measure = "deviance")
fit_lasso <- cvfit$glmnet.fit
lambda_sel <- cvfit[[LAMBDA_RULE]]

# Fig 2B — coefficient path
png("outputs/Figures/Fig2B_LASSO_coefficient_path.png",
    width = 8, height = 6, units = "in", res = 300)
plot(fit_lasso, xvar = "lambda", label = TRUE)
abline(v = log(lambda_sel), lty = 2, col = "red")
title("Figure 2B — LASSO coefficient path", line = 2.5)
dev.off()

# Fig 2C — cross-validation curve
png("outputs/Figures/Fig2C_LASSO_cv_curve.png",
    width = 8, height = 6, units = "in", res = 300)
plot(cvfit)
title("Figure 2C — LASSO 10-fold CV (partial likelihood deviance)", line = 2.5)
dev.off()

# Hệ số khác 0 tại lambda được chọn
coef_lasso <- as.matrix(coef(fit_lasso, s = lambda_sel))   # ép về matrix thường
sel_idx    <- which(coef_lasso[, 1] != 0)
sel_genes  <- rownames(coef_lasso)[sel_idx]
cat("Gene được LASSO chọn:", paste(sel_genes, collapse = ", "), "\n")
if (length(sel_genes) == 0)
  stop("LASSO không chọn gene nào ở lambda này — thử LAMBDA_RULE='lambda.min'.")

# ---- 6. Multivariate Cox của gene được chọn (Table_S3) ----------------------
multi_form <- as.formula(
  paste("Surv(OS_time, OS_event) ~", paste(sprintf("`%s`", sel_genes), collapse = " + ")))
fit_multi  <- coxph(multi_form, data = cox_df)
sm_multi   <- summary(fit_multi)

multi_tab <- data.frame(
  gene   = sel_genes,
  coef   = sm_multi$coefficients[, "coef"],
  HR     = sm_multi$coefficients[, "exp(coef)"],
  HR_low = sm_multi$conf.int[, "lower .95"],
  HR_high= sm_multi$conf.int[, "upper .95"],
  p      = sm_multi$coefficients[, "Pr(>|z|)"],
  row.names = NULL
)
write.csv(multi_tab, "outputs/Tables/Table_S3_multivariate_cox_signature.csv",
          row.names = FALSE)

# Kiểm tra PH cho mô hình đa biến
zph_multi <- cox.zph(fit_multi)
capture.output(print(zph_multi)) |>
  writeLines("outputs/Tables/Table_S3b_PH_assumption_check.txt")

# ---- 7. Risk score + phân nhóm high/low (Figure 2D) -------------------------
# Risk score = linear predictor của mô hình Cox đa biến (β·z).
# Dùng hệ số LASSO sẽ "shrunk"; ở đây dùng coxph đa biến trên gene đã chọn —
# nêu rõ công thức trong methods để tái lập được.
risk_score <- predict(fit_multi, type = "lp")
cox_df$risk_score <- as.numeric(risk_score)

# Cắt ở MEDIAN (đơn giản, tái lập được). Có thể đổi sang surv_cutpoint nếu muốn,
# nhưng phải báo cáo cách chọn để tránh "optimal cutpoint" bias.
cut_val <- median(cox_df$risk_score)
cox_df$risk_group <- factor(ifelse(cox_df$risk_score > cut_val, "High", "Low"),
                            levels = c("Low", "High"))

fit_km <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = cox_df)
p_km <- ggsurvplot(
  fit_km, data = cox_df,
  pval = TRUE, conf.int = TRUE, risk.table = TRUE,
  palette = c("#2c7bb6", "#c0392b"),
  legend.title = "Risk", legend.labs = c("Low", "High"),
  xlab = "Time (days)", ylab = "Overall survival",
  title = "Figure 2D — KM survival by cuproptosis risk score",
  risk.table.height = 0.28
)
ggsave("outputs/Figures/Fig2D_KM_riskscore.png",
       p_km$plot, width = 7, height = 5.5, dpi = 300)
ggsave("outputs/Figures/Fig2D_KM_riskscore_table.png",
       p_km$table, width = 7, height = 2, dpi = 300)

# ---- 8. Time-dependent ROC (Figure 2E) --------------------------------------
roc <- timeROC(T = cox_df$OS_time, delta = cox_df$OS_event,
               marker = cox_df$risk_score, cause = 1,
               times = ROC_TIMES, iid = TRUE)

png("outputs/Figures/Fig2E_timeROC.png",
    width = 7, height = 7, units = "in", res = 300)
plot(roc, time = ROC_TIMES[1], col = "#2c7bb6", title = FALSE)
plot(roc, time = ROC_TIMES[2], col = "#f39c12", add = TRUE)
plot(roc, time = ROC_TIMES[3], col = "#c0392b", add = TRUE)
legend("bottomright", bty = "n",
       legend = sprintf("%d-yr AUC = %.3f", c(1, 3, 5), roc$AUC),
       col = c("#2c7bb6", "#f39c12", "#c0392b"), lwd = 2)
title("Figure 2E — Time-dependent ROC of risk score")
dev.off()

cat("Time-ROC AUC (1/3/5-yr):", paste(round(roc$AUC, 3), collapse = " / "), "\n")

# ---- 9. Lưu objects + sessionInfo -------------------------------------------
saveRDS(
  list(
    uni_cox      = uni_cox,
    candidate_genes = candidate_genes,
    cvfit        = cvfit,
    lambda_sel   = lambda_sel,
    sel_genes    = sel_genes,
    fit_multi    = fit_multi,
    multi_tab    = multi_tab,
    cox_df       = cox_df,           # gồm risk_score + risk_group
    cut_val      = cut_val,
    timeROC_AUC  = roc$AUC,
    params       = list(ALPHA_SIG = ALPHA_SIG, CANDIDATE_P = CANDIDATE_P,
                        LAMBDA_RULE = LAMBDA_RULE, ROC_TIMES = ROC_TIMES)
  ),
  "outputs/Rdata/STEP2_objects.rds"
)

writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_STEP2.txt")

cat("\n========================================\n")
cat(" STEP 2 completed.\n")
cat(" Signature genes:", paste(sel_genes, collapse = ", "), "\n")
cat(" 5-yr AUC:", round(roc$AUC[3], 3), "\n")
cat(" Figures -> outputs/Figures/  | Tables -> outputs/Tables/\n")
cat("========================================\n")
