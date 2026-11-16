# =============================================================================
# STEP 4 (FINAL) — Downstream Characterization của Cuproptosis Prognostic Signature
# Bản gốc: học sinh (Thanh) | Mentor hoàn thiện: 02/08/2026
#
# >>> ĐÂY LÀ BẢN CỦA HỌC SINH ĐÃ ĐƯỢC MENTOR SỬA <<<
#   Giữ nguyên phần em làm tốt: dd<<- global env, regplot nomogram, fun.at thủ
#   công, subtitles=FALSE cho calibration, gộp mức smoking hiếm, join clin_raw.rds
#   bằng submitter_id, lazy-load + SKIP có message.
#   Mọi chỗ mentor sửa được đánh dấu:  # [MENTOR FIX #n] ...
#
#   5 lỗi đã sửa (chi tiết trong STEP4_REVIEW_for_student.md):
#     #1 [CRITICAL] ESTIMATE: TumorPurity không có khi platform="illumina"
#     #2 [MAJOR]    ESTIMATE: read.table đổi '-'->'.' làm merge rớt hết mẫu
#     #3 [MAJOR]    native pipe '_' cần R>=4.2 -> viết lại không dùng '_'
#     #4 [MAJOR]    DCA logistic bỏ censoring -> chuyển DCA time-to-event
#     #5 [MAJOR]    GSVA/msigdbr đổi API -> ssgseaParam() + collection=
#   + vài MINOR: GSEA rank theo t-stat, enrichplot vào guard, calibration m,
#     path đưa lên khối CONFIG đầu file.
#
# Package "nặng" (rms/regplot/estimate/GSVA/clusterProfiler/msigdbr/limma/dcurves)
# lazy-load trong hàm — thiếu thì SKIP có message, không crash cả script.
# =============================================================================

# ---- 0. Môi trường -----------------------------------------------------------
stopifnot(getRversion() >= "4.2.0")   # [MENTOR FIX #3] khớp việc dùng native pipe

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

dir.create("outputs/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Rdata",   recursive = TRUE, showWarnings = FALSE)

set.seed(123)
FDR_ALPHA <- 0.05
PALETTE   <- c(Low = "#2c7bb6", High = "#c0392b")   # khớp STEP2/STEP3

# ---- CONFIG: mọi đường dẫn khai báo 1 chỗ (dễ sửa khi đổi máy) --------------
# [MENTOR FIX MINOR] gom absolute path lên đây thay vì rải rác giữa file.
STEP1_RDS <- "/Users/aithanh/Cuproptosis_LUAD/Step1/outputs_step1_new/Rdata/STEP1_objects.rds"
STEP2_RDS <- "/Users/aithanh/Cuproptosis_LUAD/Step_2/outputs/Rdata/STEP2_objects.rds"
CLIN_RDS  <- "/Users/aithanh/Cuproptosis_LUAD/Step1/outputs_step1_new/Rdata/clin_raw.rds"
GMT_IMMUNE <- "outputs/immune_28_celltypes.gmt"   # <-- 28 immune cell gene set (Charoentong)

# ---- 1. Nạp objects "đóng băng" từ STEP 1 + STEP 2 --------------------------
step1 <- readRDS(STEP1_RDS)
step2 <- readRDS(STEP2_RDS)

expr_mat  <- as.matrix(step1$expr_mat)   # gene x sample, log2-CPM
cox_df    <- step2$cox_df                # rownames = sample; có risk_score, risk_group
sel_genes <- step2$sel_genes

stopifnot(all(c("risk_score", "risk_group") %in% colnames(cox_df)))
cat("  Samples (từ STEP 2, đã 1 bệnh nhân/1 mẫu):", nrow(cox_df), "\n")
cat("  Risk group:", paste(names(table(cox_df$risk_group)),
                           table(cox_df$risk_group), sep = "=", collapse = " | "), "\n")

# Lookup dùng chung: sample (full barcode) -> risk_group
risk_lookup <- cox_df |> rownames_to_column("sample") |> select(sample, risk_group)

# ---- 2. Ghép clinical (barcode 12 ký tự, giống cách join ở STEP 1) ----------
if (!file.exists(CLIN_RDS))
  stop("Không tìm thấy clin_raw.rds tại '", CLIN_RDS, "' — kiểm tra CONFIG ở đầu file.")
clin_df <- readRDS(CLIN_RDS)   # output của GDCquery_clinic() ở STEP 1
cat("\n>> colnames(clin_df) (một phần):\n"); print(head(colnames(clin_df), 40))

COL_STAGE   <- "ajcc_pathologic_stage"
COL_AGE     <- "age_at_index"
COL_SEX     <- "gender"
COL_BARCODE <- "submitter_id"

# Cột smoking không đảm bảo luôn có -> tự dò, chỉ lấy cột PHÂN LOẠI (status/history),
# bỏ biến liên tục (pack_years) để không nổ dummy variable.
smoking_candidates <- grep("smok|tobacco", colnames(clin_df), ignore.case = TRUE, value = TRUE)
cat(">> Cột liên quan smoking:", paste(smoking_candidates, collapse = ", "), "\n")
status_pref <- grep("status|history", smoking_candidates, ignore.case = TRUE, value = TRUE)
COL_SMOKING <- if (length(status_pref) > 0) status_pref[1] else NA_character_
if (is.na(COL_SMOKING))
  message("Không có cột smoking dạng phân loại -> BỎ QUA smoking khỏi mô hình.")

stopifnot(all(c(COL_STAGE, COL_AGE, COL_SEX, COL_BARCODE) %in% colnames(clin_df)))

short_barcode <- function(x) substr(x, 1, 12)
cox_df$barcode12  <- short_barcode(rownames(cox_df))
clin_df$barcode12 <- short_barcode(clin_df[[COL_BARCODE]])

# [MENTOR FIX MINOR] đặt tên rõ cho phần tử đầu của select vector
select_cols <- c(barcode12 = "barcode12", stage = COL_STAGE, age = COL_AGE, sex = COL_SEX)
if (!is.na(COL_SMOKING)) select_cols <- c(select_cols, smoking = COL_SMOKING)

model_df <- cox_df |>
  rownames_to_column("sample") |>
  inner_join(
    clin_df |> select(all_of(select_cols)) |> distinct(barcode12, .keep_all = TRUE),
    by = "barcode12"
  )
cat("\n>> Ghép clinical: ", nrow(model_df), "/", nrow(cox_df),
    " mẫu có đủ lâm sàng (", nrow(cox_df) - nrow(model_df), " mẫu rớt khi merge)\n", sep = "")
if (nrow(model_df) < 0.5 * nrow(cox_df))
  warning("Hơn 50% mẫu rớt khi merge — kiểm COL_BARCODE / định dạng barcode.")

# [MENTOR FIX #3] dọn stage KHÔNG dùng placeholder '_' (tương thích R rộng hơn)
model_df$stage_group <- factor(
  gsub("[ABC]$", "", gsub("Stage ", "", model_df$stage)),
  levels = c("I", "II", "III", "IV"))
model_df$sex <- factor(model_df$sex)
model_df$age <- suppressWarnings(as.numeric(model_df$age))
if ("smoking" %in% colnames(model_df)) {
  model_df$smoking <- factor(model_df$smoking)
  lvl_n <- table(model_df$smoking)
  if (length(lvl_n) > 6) {
    rare <- names(lvl_n[lvl_n < 10])
    if (length(rare) > 0) {
      levels(model_df$smoking)[levels(model_df$smoking) %in% rare] <- "Other"
      message("Smoking >6 mức -> gộp mức <10 mẫu thành 'Other'.")
    }
  }
}
model_df <- model_df |> filter(!is.na(stage_group), !is.na(age), !is.na(sex))
cat(">> Sau loại NA covariate: ", nrow(model_df), " mẫu cho Phần A\n", sep = "")

# =============================================================================
# PHẦN A — Independent prognosis + Nomogram (Figure 4A-4D)
# =============================================================================

# ---- A1. Univariate + Multivariate Cox (Figure 4A, Table S6) ---------------
covariates <- c("risk_score", "stage_group", "age", "sex")
if ("smoking" %in% colnames(model_df)) covariates <- c(covariates, "smoking")

uni_list <- lapply(covariates, function(v) {
  fit <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~", v)), data = model_df)
  sm  <- summary(fit)
  data.frame(variable = rownames(sm$coefficients),
             HR = sm$coefficients[, "exp(coef)"],
             HR_low = sm$conf.int[, "lower .95"], HR_high = sm$conf.int[, "upper .95"],
             p = sm$coefficients[, "Pr(>|z|)"], model = "Univariate", row.names = NULL)
}) |> bind_rows()

fit_multi4 <- coxph(as.formula(paste("Surv(OS_time, OS_event) ~",
                                     paste(covariates, collapse = " + "))), data = model_df)
sm4 <- summary(fit_multi4)
multi_tab <- data.frame(variable = rownames(sm4$coefficients),
                        HR = sm4$coefficients[, "exp(coef)"],
                        HR_low = sm4$conf.int[, "lower .95"], HR_high = sm4$conf.int[, "upper .95"],
                        p = sm4$coefficients[, "Pr(>|z|)"], model = "Multivariate", row.names = NULL)

table_s6 <- bind_rows(uni_list, multi_tab)
write.csv(table_s6, "outputs/Tables/Table_S6_independent_prognosis.csv", row.names = FALSE)

risk_p_multi <- multi_tab$p[multi_tab$variable == "risk_score"]
cat(sprintf("\n>> risk_score multivariate p = %.3g -> %s\n", risk_p_multi,
            ifelse(risk_p_multi < 0.05,
                   "ĐỘC LẬP với stage/age/sex/smoking (p<0.05)",
                   "KHÔNG còn ý nghĩa sau hiệu chỉnh — báo cáo trung thực + bàn Discussion")))

clean_varname <- function(v) dplyr::case_when(
  v == "risk_score"        ~ "Risk score",
  v == "age"               ~ "Age",
  grepl("^stage_group", v) ~ paste("Stage", sub("^stage_group", "", v)),
  grepl("^sex", v)         ~ paste0("Sex: ", tools::toTitleCase(sub("^sex", "", v))),
  grepl("^smoking", v)     ~ paste0("Smoking: ", sub("^smoking", "", v)),
  TRUE                     ~ v)

fp4 <- multi_tab |>
  mutate(variable_clean = clean_varname(variable),
         variable_clean = factor(variable_clean, levels = rev(variable_clean)),
         label_full = sprintf("%.2f (%.2f-%.2f), p %s", HR, HR_low, HR_high,
                              ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))))
x_text4 <- max(fp4$HR_high, na.rm = TRUE) * 1.3
p_fig4a <- ggplot(fp4, aes(x = HR, y = variable_clean)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = HR_low, xmax = HR_high), height = 0.2) +
  geom_point(shape = 18, size = 3.2, color = "#c0392b") +
  geom_text(aes(x = x_text4, label = label_full), hjust = 0, size = 3.1) +
  scale_x_log10(limits = c(min(fp4$HR_low, na.rm = TRUE) * 0.9, x_text4 * 3.4)) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  labs(title = "Figure 4A — Multivariate Cox: risk score vs. clinical covariates",
       subtitle = sprintf("n = %d | risk_score p = %.3g", nrow(model_df), risk_p_multi),
       x = "Hazard Ratio (95% CI, log scale)", y = NULL)
fig4a_h <- max(4.5, 0.35 * nrow(fp4) + 1.5)
ggsave("outputs/Figures/Fig4A_independent_prognosis_forest.png",
       p_fig4a, width = 11, height = fig4a_h, dpi = 300, bg = "white")
cat(">> Figure 4A saved.\n")

# ---- A2. Nomogram (Fig4B) + Calibration (Fig4C) ----------------------------
run_nomogram <- function(df) {
  if (!requireNamespace("rms", quietly = TRUE)) {
    message("Chưa cài 'rms' -> BỎ QUA Fig4B/4C."); return(invisible(NULL)) }
  suppressPackageStartupMessages(library(rms))
  dd <<- datadist(df); options(datadist = "dd")   # dd ở GLOBAL env (bẫy của cph)
  df$OS_time_yr <- df$OS_time / 365.25

  f_cph <- cph(Surv(OS_time_yr, OS_event) ~ risk_score + stage_group + age + sex,
               data = df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 5)
  surv_fun <- Survival(f_cph)
  fun_at_list <- list(
    c(0.95,0.9,0.85,0.8,0.7,0.6,0.5,0.4,0.3,0.2),
    c(0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1),
    c(0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1))
  nom <- nomogram(f_cph,
    fun = list(function(x) surv_fun(1, x), function(x) surv_fun(3, x), function(x) surv_fun(5, x)),
    funlabel = c("1-year OS", "3-year OS", "5-year OS"),
    fun.at = fun_at_list, lp = FALSE)

  used_regplot <- FALSE
  if (requireNamespace("regplot", quietly = TRUE)) {
    used_regplot <- tryCatch({
      png("outputs/Figures/Fig4B_nomogram.png", width = 10, height = 7, units = "in", res = 300, bg = "white")
      regplot::regplot(f_cph, observation = df[1, ], points = TRUE, odds = FALSE,
                       droplines = TRUE, failtime = c(1, 3, 5), prfail = FALSE,
                       title = "Figure 4B — Nomogram (1/3/5-year OS)")
      dev.off(); TRUE
    }, error = function(e) { message("regplot lỗi (", conditionMessage(e), ") -> rms::nomogram.")
      if (dev.cur() > 1) dev.off(); FALSE })
  }
  if (!used_regplot) {
    png("outputs/Figures/Fig4B_nomogram.png", width = 11, height = 7.5, units = "in", res = 300, bg = "white")
    plot(nom, xfrac = 0.35, cex.axis = 0.85, cex.var = 1)
    title("Figure 4B — Nomogram (1/3/5-year OS)", line = 3); dev.off()
  }
  cat(">> Figure 4B saved (", ifelse(used_regplot, "regplot", "rms::nomogram"), ").\n")

  # Calibration 1/3/5-year (subtitles=FALSE để caption không đè trục)
  png("outputs/Figures/Fig4C_calibration.png", width = 11, height = 4, units = "in", res = 300, bg = "white")
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 2) + 0.1)
  cal_summary <- list()
  for (yr in c(1, 3, 5)) {
    f_u <- cph(Surv(OS_time_yr, OS_event) ~ risk_score + stage_group + age + sex,
               data = df, x = TRUE, y = TRUE, surv = TRUE, time.inc = yr, u = yr)
    cal <- tryCatch(
      calibrate(f_u, method = "boot", u = yr, m = max(floor(nrow(df) / 5), 20), B = 200),  # [MENTOR FIX MINOR] m=n/5
      error = function(e) { message("Calibration ", yr, "-yr lỗi: ", conditionMessage(e)); NULL })
    if (!is.null(cal)) {
      plot(cal, xlab = paste0("Predicted ", yr, "-year OS"),
           ylab = paste0("Observed ", yr, "-year OS"), main = paste0(yr, "-year"), subtitles = FALSE)
      cal_summary[[paste0(yr, "yr")]] <- cal
    }
  }
  dev.off(); cat(">> Figure 4C saved.\n")
  options(datadist = NULL); if (exists("dd", envir = .GlobalEnv)) rm(dd, envir = .GlobalEnv)
  list(cph_fit = f_cph, nomogram = nom, calibration = cal_summary)
}
nomogram_res <- run_nomogram(model_df)

# ---- A3. Decision Curve Analysis (Figure 4D) — time-to-event ---------------
run_dca <- function(df) {
  if (!requireNamespace("dcurves", quietly = TRUE)) {
    message("Chưa cài 'dcurves' -> BỎ QUA Fig4D (không bắt buộc)."); return(invisible(NULL)) }
  suppressPackageStartupMessages(library(dcurves))
  # [MENTOR FIX #4] DCA time-to-event tại mốc 5 năm (1825 ngày), KHÔNG logistic.
  # Xác suất TỬ VONG 5 năm = 1 - S(1825) từ mô hình Cox đã fit.
  pred <- tryCatch({
    sf <- survfit(fit_multi4, newdata = df)
    s5 <- summary(sf, times = 1825)$surv
    1 - as.numeric(s5)
  }, error = function(e) { message("Không tính được pred 5-yr: ", conditionMessage(e)); NULL })
  if (is.null(pred)) return(invisible(NULL))
  df$pred_5yr <- pred
  p_dca <- dca(Surv(OS_time, OS_event) ~ pred_5yr, data = df, time = 1825) |>
    plot(smooth = TRUE) + labs(title = "Figure 4D — Decision Curve Analysis (5-year)")
  ggsave("outputs/Figures/Fig4D_DCA.png", p_dca, width = 7, height = 5, dpi = 300, bg = "white")
  cat(">> Figure 4D saved.\n"); invisible(p_dca)
}
dca_res <- run_dca(model_df)

# =============================================================================
# PHẦN B — Immune Microenvironment (Figure 5A-5C)
# =============================================================================

# ---- B1. ESTIMATE scores (Figure 5A) ----------------------------------------
run_estimate <- function(expr_mat, risk_df) {
  if (!requireNamespace("estimate", quietly = TRUE)) {
    message("Chưa cài 'estimate' (R-Forge) -> BỎ QUA Fig5A."); return(invisible(NULL)) }
  suppressPackageStartupMessages(library(estimate))
  tmp_in <- tempfile(fileext = ".txt"); tmp_filt <- tempfile(fileext = ".gct"); tmp_out <- tempfile(fileext = ".gct")
  write.table(data.frame(GeneSymbol = rownames(expr_mat), expr_mat, check.names = FALSE),
              tmp_in, sep = "\t", quote = FALSE, row.names = FALSE)
  filterCommonGenes(input.f = tmp_in, output.f = tmp_filt, id = "GeneSymbol")
  estimateScore(tmp_filt, tmp_out, platform = "illumina")

  # [MENTOR FIX #2] check.names=FALSE để giữ dấu '-' trong tên mẫu (không -> '.')
  sc <- read.table(tmp_out, skip = 2, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  sc <- as.data.frame(t(sc[, -1])); sc$sample <- rownames(sc)

  # [MENTOR FIX #1] platform="illumina" KHÔNG có TumorPurity -> lấy score động
  score_cols <- intersect(c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"), colnames(sc))

  merged <- inner_join(risk_df, sc, by = "sample")
  cat(">> Merge ESTIMATE:", nrow(merged), "/", nrow(risk_df), "mẫu khớp\n")   # [MENTOR FIX #2] kiểm chứng
  if (nrow(merged) < 0.5 * nrow(risk_df))
    warning("ESTIMATE: <50% mẫu khớp — kiểm tên mẫu (dấu '-' vs '.').")

  long <- merged |> pivot_longer(all_of(score_cols), names_to = "score_type", values_to = "value")
  stat <- long |> group_by(score_type) |>
    summarise(p = wilcox.test(value ~ risk_group)$p.value, .groups = "drop") |>
    mutate(padj = p.adjust(p, method = "BH"))

  p_fig5a <- ggplot(long, aes(risk_group, value, fill = risk_group)) +
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~ score_type, scales = "free_y") +
    scale_fill_manual(values = PALETTE) + theme_bw(base_size = 12) +
    labs(title = "Figure 5A — ESTIMATE scores by risk group",
         subtitle = "Wilcoxon test, BH-FDR", x = NULL, y = "Score")
  ggsave("outputs/Figures/Fig5A_ESTIMATE_scores.png", p_fig5a, width = 8, height = 6, dpi = 300, bg = "white")
  write.csv(stat, "outputs/Tables/Table_S7a_ESTIMATE_stat.csv", row.names = FALSE)
  cat(">> Figure 5A saved.\n"); list(scores = merged, stat = stat)
}
estimate_res <- run_estimate(expr_mat, risk_lookup)

# ---- B2. ssGSEA 28 loại tế bào miễn dịch (Figure 5B) ------------------------
run_ssgsea_immune <- function(expr_mat, risk_df) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    message("Chưa cài 'GSVA' -> BỎ QUA Fig5B."); return(invisible(NULL)) }
  if (!file.exists(GMT_IMMUNE)) {
    message("Chưa có ", GMT_IMMUNE, " (28 immune cell gene sets) -> BỎ QUA Fig5B."); return(invisible(NULL)) }
  suppressPackageStartupMessages(library(GSVA))
  gene_sets <- GSEABase::getGmt(GMT_IMMUNE)

  # [MENTOR FIX #5] API mới: ssgseaParam() + gsva(param); KHÔNG gsva(mat, method=)
  ss_par <- GSVA::ssgseaParam(as.matrix(expr_mat), gene_sets, normalize = TRUE)
  ssgsea_scores <- GSVA::gsva(ss_par)
  ssgsea_df <- as.data.frame(t(ssgsea_scores)) |> rownames_to_column("sample")

  merged <- inner_join(risk_df, ssgsea_df, by = "sample")
  cat(">> Merge ssGSEA:", nrow(merged), "/", nrow(risk_df), "mẫu khớp\n")
  cell_cols <- setdiff(colnames(merged), c("sample", "risk_group"))

  stat <- lapply(cell_cols, function(cc)
    data.frame(cell_type = cc, p = wilcox.test(merged[[cc]] ~ merged$risk_group)$p.value)) |>
    bind_rows() |> mutate(padj = p.adjust(p, method = "BH")) |> arrange(padj)
  write.csv(stat, "outputs/Tables/Table_S7b_immune_cell_stat.csv", row.names = FALSE)

  long <- merged |> pivot_longer(all_of(cell_cols), names_to = "cell_type", values_to = "score")
  sig_cells <- head(stat$cell_type[stat$padj < FDR_ALPHA], 12)
  long_sig <- long |> filter(cell_type %in% sig_cells)
  p_fig5b <- ggplot(long_sig, aes(risk_group, score, fill = risk_group)) +
    geom_boxplot(outlier.size = 0.4) + facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = PALETTE) + theme_bw(base_size = 10) +
    labs(title = "Figure 5B — Immune cell infiltration (ssGSEA) by risk group",
         subtitle = "Top significant cell types, BH-FDR < 0.05", x = NULL, y = "ssGSEA score")
  ggsave("outputs/Figures/Fig5B_immune_infiltration.png", p_fig5b, width = 11, height = 8, dpi = 300, bg = "white")
  cat(">> Figure 5B saved.\n"); list(scores = merged, stat = stat)
}
ssgsea_res <- run_ssgsea_immune(expr_mat, risk_lookup)

# ---- B3. Immune checkpoint expression (Figure 5C) ---------------------------
checkpoint_genes <- c(PDCD1 = "PD-1", CD274 = "PD-L1", CTLA4 = "CTLA-4",
                      LAG3 = "LAG-3", HAVCR2 = "TIM-3", TIGIT = "TIGIT")
present_cp <- intersect(names(checkpoint_genes), rownames(expr_mat))
if (length(present_cp) == 0) {
  message("Không có gene checkpoint trong expr_mat -> BỎ QUA Fig5C.")
} else {
  cp_df <- as.data.frame(t(expr_mat[present_cp, , drop = FALSE])) |>
    rownames_to_column("sample") |> inner_join(risk_lookup, by = "sample") |>
    pivot_longer(all_of(present_cp), names_to = "gene", values_to = "expr") |>
    mutate(gene_label = checkpoint_genes[gene])
  cp_stat <- cp_df |> group_by(gene) |>
    summarise(p = wilcox.test(expr ~ risk_group)$p.value, .groups = "drop") |>
    mutate(padj = p.adjust(p, method = "BH"))
  write.csv(cp_stat, "outputs/Tables/Table_S7c_checkpoint_stat.csv", row.names = FALSE)
  p_fig5c <- ggplot(cp_df, aes(risk_group, expr, fill = risk_group)) +
    geom_boxplot(outlier.size = 0.5) + facet_wrap(~ gene_label, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = PALETTE) + theme_bw(base_size = 11) +
    labs(title = "Figure 5C — Immune checkpoint expression by risk group",
         subtitle = "Wilcoxon test, BH-FDR", x = NULL, y = "log2-CPM")
  ggsave("outputs/Figures/Fig5C_immune_checkpoints.png", p_fig5c,
         width = 2 + 2 * length(present_cp), height = 4.5, dpi = 300, bg = "white")
  cat(">> Figure 5C saved.\n")
}

# =============================================================================
# PHẦN C — Functional Enrichment (Figure 6A-6B)
# =============================================================================
run_gsea <- function(expr_mat, risk_df) {
  if (!requireNamespace("limma", quietly = TRUE) ||
      !requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("msigdbr", quietly = TRUE)) {
    message("Thiếu limma/clusterProfiler/msigdbr -> BỎ QUA Fig6."); return(invisible(NULL)) }
  suppressPackageStartupMessages({ library(limma); library(clusterProfiler); library(msigdbr) })

  common <- intersect(colnames(expr_mat), risk_df$sample)
  em <- expr_mat[, common, drop = FALSE]
  grp <- factor(risk_df$risk_group[match(common, risk_df$sample)], levels = c("Low", "High"))
  design <- model.matrix(~ grp)
  fit <- eBayes(lmFit(em, design))
  deg <- topTable(fit, coef = 2, number = Inf) |> rownames_to_column("gene")
  write.csv(deg, "outputs/Tables/Table_S8_DEG_high_vs_low.csv", row.names = FALSE)

  # [MENTOR FIX MINOR] rank theo moderated t-statistic (có dấu), ổn định hơn logFC
  ranked <- sort(setNames(deg$t, deg$gene), decreasing = TRUE)

  # [MENTOR FIX #5] msigdbr >=10 dùng collection=; <10 dùng category=
  hallmark <- if (utils::packageVersion("msigdbr") >= "10.0.0")
    msigdbr(species = "Homo sapiens", collection = "H")[, c("gs_name", "gene_symbol")]
  else
    msigdbr(species = "Homo sapiens", category = "H")[, c("gs_name", "gene_symbol")]

  gsea_res <- GSEA(ranked, TERM2GENE = hallmark, pvalueCutoff = 1, seed = TRUE)
  gsea_tab <- as.data.frame(gsea_res) |> mutate(padj = p.adjust(pvalue, method = "BH"))
  write.csv(gsea_tab, "outputs/Tables/Table_S8b_GSEA_hallmark.csv", row.names = FALSE)

  if (nrow(gsea_tab) > 0) {
    p_fig6a <- dotplot(gsea_res, showCategory = 15, split = ".sign") + facet_grid(. ~ .sign) +
      labs(title = "Figure 6A — GSEA (HALLMARK), High vs Low risk")
    ggsave("outputs/Figures/Fig6A_GSEA_hallmark_dotplot.png", p_fig6a, width = 10, height = 8, dpi = 300, bg = "white")
    # [MENTOR FIX MINOR] Fig6B chỉ vẽ khi có enrichplot (không làm hàm chết)
    if (requireNamespace("enrichplot", quietly = TRUE)) {
      p_fig6b <- enrichplot::gseaplot2(gsea_res, geneSetID = gsea_tab$ID[1],
                                       title = paste0("Figure 6B — ", gsea_tab$ID[1]))
      ggsave("outputs/Figures/Fig6B_GSEA_top_pathway.png", p_fig6b, width = 8, height = 6, dpi = 300, bg = "white")
    } else message("Chưa cài 'enrichplot' -> bỏ Fig6B (Fig6A vẫn có).")
    cat(">> Figure 6A(-6B) saved.\n")
  }
  list(deg = deg, gsea = gsea_tab)
}
gsea_res_out <- run_gsea(expr_mat, risk_lookup)

# =============================================================================
# 7. Lưu objects + sessionInfo
# =============================================================================
saveRDS(list(model_df = model_df, table_s6 = table_s6, fit_multi4 = fit_multi4,
             nomogram_res = nomogram_res, dca_res = if (exists("dca_res")) TRUE else NULL,
             estimate_res = estimate_res, ssgsea_res = ssgsea_res,
             gsea_res = gsea_res_out),
        "outputs/Rdata/STEP4_objects.rds")
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_STEP4.txt")

cat("\n========================================\n")
cat(" STEP 4 completed (phần thiếu package/file input đã SKIP có message).\n")
cat(" risk_score multivariate p =", signif(risk_p_multi, 4), "\n")
cat(" Figures -> outputs/Figures/  | Tables -> outputs/Tables/\n")
cat("========================================\n")
