# =============================================================================
# STEP 3 (REVIEWED) — External validation + Drug sensitivity of the
#                     Cuproptosis prognostic SIGNATURE (LUAD)
# Mentor revision: 05/07/2026
# Project:  Cuproptosis Prognostic Signature in LUAD
# Author:   <student name>     |  Mentor: Thao
#
# Mục tiêu STEP 3:
#   PHẦN A — VALIDATION (Figure 3A–3C)
#     1. Nạp signature (gene + hệ số β) đã "đóng băng" từ STEP 2.
#     2. Tải 2 cohort GEO độc lập: GSE31210 (n≈226) + GSE72094 (n≈398).
#     3. Map probe -> gene signature, áp ĐÚNG công thức risk score của STEP 2.
#     4. Chia high/low, vẽ KM + time-dependent ROC cho TỪNG cohort.
#   PHẦN B — DRUG SENSITIVITY (Figure 3D)
#     5. Dự đoán IC50 (oncoPredict, GDSC2) cho TCGA-LUAD.
#     6. So sánh IC50 giữa nhóm High vs Low risk (Wilcoxon + BH-FDR).
#   7. Lưu objects + sessionInfo.
#
# Đầu vào: outputs/Rdata/STEP2_objects.rds  (sinh ra ở STEP 2)
#          -> chứa: sel_genes, fit_multi, multi_tab, cox_df (risk_score),
#                   cut_val, params ...
#          outputs/Rdata/STEP1_objects.rds  (cần expr TCGA cho drug sensitivity)
#
# Quy ước (giữ nguyên từ STEP 1–2):
#   - Output figure  -> outputs/Figures/     table -> outputs/Tables/
#   - Object trung gian -> outputs/Rdata/     Multiple testing: BH-FDR
#   - Survival endpoint: Overall Survival (OS)
#   - Để CODE tự sinh mọi con số trên figure (KHÔNG hard-code).
#
# NGUYÊN TẮC VÀNG CỦA VALIDATION (đọc kỹ):
#   * "Đóng băng" (freeze) toàn bộ mô hình ở training: gene, hệ số β, và
#     cách tính risk score. Validation TUYỆT ĐỐI không được fit lại Cox để
#     chọn lại gene / hệ số -> nếu fit lại là data leakage, không còn là
#     "external validation" nữa.
#   * Chỉ có 1 thứ được phép ước lượng trong cohort validation: ngưỡng chia
#     high/low. Ở đây dùng MEDIAN nội cohort và KHAI BÁO rõ (xem mục A5).
# =============================================================================

# ---- 0. Thiết lập môi trường --------------------------------------------------
stopifnot(getRversion() >= "4.1.0")

suppressPackageStartupMessages({
  library(GEOquery)     # tải GSE + fData annotation
  library(survival)     # coxph, Surv
  library(survminer)    # ggsurvplot
  library(timeROC)      # time-dependent ROC
  library(oncoPredict)  # calcPhenotype -> dự đoán IC50 (Part B)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

dir.create("outputs/Figures", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Rdata",   recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/GEO_raw", recursive = TRUE, showWarnings = FALSE)

set.seed(123)

# Tham số khai báo 1 chỗ
ROC_TIMES  <- c(365, 1095, 1825)   # 1,3,5 năm (ngày)
FDR_ALPHA  <- 0.05
DRUG_TOP_N <- 12                   # số thuốc khác biệt nhất để vẽ boxplot

# ---- 1. Nạp signature ĐÃ ĐÓNG BĂNG từ STEP 2 --------------------------------
step2 <- readRDS("outputs/Rdata/STEP2_objects.rds")
step1 <- readRDS("outputs/Rdata/STEP1_objects.rds")

sel_genes <- step2$sel_genes          # gene trong signature
fit_multi <- step2$fit_multi          # mô hình Cox đa biến (nguồn hệ số β)

# Trích hệ số β của signature (thứ tự khớp tên gene, backtick-safe)
beta_vec <- coef(fit_multi)
names(beta_vec) <- gsub("`", "", names(beta_vec))
beta_vec <- beta_vec[sel_genes]       # đảm bảo đúng thứ tự
stopifnot(!any(is.na(beta_vec)))
cat("Signature:", paste(sprintf("%s(%.3f)", names(beta_vec), beta_vec),
                        collapse = ", "), "\n")

# Hàm tính risk score = β · z (z-score TỪNG GENE trong cohort đang xét).
# Đây CHÍNH là công thức STEP 2, chỉ khác dữ liệu đầu vào.
compute_risk <- function(expr_gene_sample, beta) {
  # expr_gene_sample: matrix gene(signature) × sample, thang log2 tương đương
  z <- t(scale(t(expr_gene_sample)))                 # z-score từng gene
  z[is.na(z)] <- 0                                   # gene sd=0 -> đóng góp 0
  as.numeric(beta %*% z[names(beta), , drop = FALSE]) # linear predictor
}

# =============================================================================
# PHẦN A — EXTERNAL VALIDATION
# =============================================================================

# ---- A1. Hàm tải & chuẩn hoá 1 cohort GEO -----------------------------------
# LƯU Ý: mỗi series GEO đặt tên cột lâm sàng KHÁC nhau trong pData. Hàm dưới
# nhận vào tên biến OS_time / OS_event tuỳ series -> em phải MỞ pData ra xem
# trước (str(pData(gset))) rồi truyền đúng tên. Đừng tin mù các mặc định.
load_geo_cohort <- function(gse_id, time_field, event_field,
                            time_unit = c("day", "month"),
                            probe_symbol_col = NULL) {
  time_unit <- match.arg(time_unit)
  gset <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = TRUE,
                 destdir = "outputs/GEO_raw")[[1]]

  expr <- exprs(gset)                       # probe × sample
  pdat <- pData(gset)
  fdat <- fData(gset)

  # 1) Map probe -> gene symbol. Tên cột symbol khác nhau giữa platform:
  #    GPL570 hay dùng "Gene Symbol"; GPL15048 (GSE72094) dùng cột khác.
  if (is.null(probe_symbol_col)) {
    cand <- grep("symbol", tolower(colnames(fdat)))
    if (length(cand) == 0)
      stop(gse_id, ": không tìm thấy cột gene symbol trong fData -> khai báo probe_symbol_col thủ công.")
    probe_symbol_col <- colnames(fdat)[cand[1]]
  }
  sym <- as.character(fdat[[probe_symbol_col]])
  sym <- sub(" ///.*$", "", sym)            # 1 probe nhiều symbol -> lấy cái đầu

  # 2) Collapse nhiều probe/1 gene: giữ probe có TRUNG BÌNH biểu hiện cao nhất
  #    (probe "sáng" nhất). Cách này tái lập được và hay dùng trong paper.
  keep <- !is.na(sym) & sym != ""
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  ord  <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[ord, ]; sym <- sym[ord]
  dup  <- duplicated(sym)
  expr <- expr[!dup, , drop = FALSE]
  rownames(expr) <- sym[!dup]

  # 3) log2 nếu dữ liệu ở thang tuyến tính (heuristic: max lớn -> chưa log)
  if (max(expr, na.rm = TRUE) > 100) {
    expr <- log2(expr + 1)
    message(gse_id, ": phát hiện thang tuyến tính -> đã log2(x+1).")
  }

  # 4) Survival: đọc từ pData theo tên cột em truyền vào
  os_time <- suppressWarnings(as.numeric(pdat[[time_field]]))
  os_evt  <- pdat[[event_field]]
  # Chuẩn hoá event về 0/1 (GEO hay ghi "alive"/"dead", "0"/"1", TRUE/FALSE)
  os_evt  <- ifelse(grepl("dead|1|yes|true", tolower(as.character(os_evt))), 1L,
              ifelse(grepl("alive|0|no|false", tolower(as.character(os_evt))), 0L, NA))
  if (time_unit == "month") os_time <- os_time * 30.44   # -> ngày

  df <- data.frame(sample = colnames(expr),
                   OS_time = os_time, OS_event = os_evt,
                   stringsAsFactors = FALSE)
  df <- df[!is.na(df$OS_time) & !is.na(df$OS_event) & df$OS_time > 0, ]
  expr <- expr[, df$sample, drop = FALSE]

  # 5) Kiểm gene signature có đủ không
  miss <- setdiff(sel_genes, rownames(expr))
  if (length(miss) > 0)
    warning(gse_id, ": THIẾU gene signature trên platform này: ",
            paste(miss, collapse = ", "),
            " -> risk score dùng phần gene còn lại (nêu rõ trong limitations).")
  present <- intersect(sel_genes, rownames(expr))

  list(gse = gse_id, expr = expr[present, , drop = FALSE],
       surv = df, genes_present = present, genes_missing = miss)
}

# ---- A2. Nạp 2 cohort. ⚠ time_field/event_field bên dưới là PLACEHOLDER ------
#   Em BẮT BUỘC mở pData ra kiểm tra tên cột thật rồi thay lại:
#     g <- getGEO("GSE31210")[[1]]; View(pData(g))  (tìm cột OS time / status)
#   Ví dụ GSE31210: thường có "days before death/censor" & "death"; GSE72094:
#   "survival_time_in_days" & "vital_status". KIỂM TRA rồi mới chạy tiếp.
coh1 <- load_geo_cohort("GSE31210",
                        time_field  = "days before death/censor:ch1",   # <-- kiểm tra lại
                        event_field = "death:ch1",                      # <-- kiểm tra lại
                        time_unit   = "day")

coh2 <- load_geo_cohort("GSE72094",
                        time_field  = "survival_time_in_days:ch1",      # <-- kiểm tra lại
                        event_field = "vital_status:ch1",               # <-- kiểm tra lại
                        time_unit   = "day")

# ---- A3. Áp risk score + KM + timeROC cho 1 cohort --------------------------
validate_cohort <- function(coh, tag) {
  rs <- compute_risk(coh$expr, beta_vec)
  df <- coh$surv
  df$risk_score <- rs

  # A5. Chia high/low ở MEDIAN NỘI COHORT (khai báo minh bạch trong methods).
  #     Lý do: risk score = β·z, mà z chuẩn hoá riêng từng cohort -> median
  #     training không cùng thang. Median nội cohort là lựa chọn pre-specified,
  #     tái lập được, tránh optimal-cutpoint bias.
  cut_v <- median(df$risk_score)
  df$risk_group <- factor(ifelse(df$risk_score > cut_v, "High", "Low"),
                          levels = c("Low", "High"))

  # KM
  fit_km <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = df)
  p_km <- ggsurvplot(
    fit_km, data = df, pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    palette = c("#2c7bb6", "#c0392b"),
    legend.title = "Risk", legend.labs = c("Low", "High"),
    xlab = "Time (days)", ylab = "Overall survival",
    title = sprintf("Figure 3 — KM validation (%s, n=%d)", tag, nrow(df)),
    risk.table.height = 0.28)
  ggsave(sprintf("outputs/Figures/Fig3_KM_%s.png", tag),
         p_km$plot, width = 7, height = 5.5, dpi = 300)
  ggsave(sprintf("outputs/Figures/Fig3_KM_%s_table.png", tag),
         p_km$table, width = 7, height = 2, dpi = 300)

  # time-ROC (bỏ time điểm nào vượt max follow-up)
  valid_times <- ROC_TIMES[ROC_TIMES < max(df$OS_time)]
  roc <- timeROC(T = df$OS_time, delta = df$OS_event,
                 marker = df$risk_score, cause = 1,
                 times = valid_times, iid = TRUE)
  list(df = df, cut = cut_v, roc = roc, times = valid_times,
       km_pval = survdiff(Surv(OS_time, OS_event) ~ risk_group, data = df))
}

res1 <- validate_cohort(coh1, "GSE31210")
res2 <- validate_cohort(coh2, "GSE72094")

# ---- A4. Figure 3C — timeROC gộp 2 cohort trên 1 khung ----------------------
png("outputs/Figures/Fig3C_timeROC_validation.png",
    width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
for (r in list(list(res1, "GSE31210"), list(res2, "GSE72094"))) {
  rr <- r[[1]]; tg <- r[[2]]
  cols <- c("#2c7bb6", "#f39c12", "#c0392b")[seq_along(rr$times)]
  plot(rr$roc, time = rr$times[1], col = cols[1], title = FALSE)
  if (length(rr$times) > 1)
    for (k in 2:length(rr$times)) plot(rr$roc, time = rr$times[k], col = cols[k], add = TRUE)
  legend("bottomright", bty = "n",
         legend = sprintf("%d-yr AUC = %.3f", round(rr$times/365), rr$roc$AUC),
         col = cols, lwd = 2)
  title(sprintf("Figure 3C — time-ROC (%s)", tg))
}
dev.off()

# Bảng tổng hợp AUC validation
val_auc <- bind_rows(
  data.frame(cohort = "GSE31210", years = round(res1$times/365), AUC = res1$roc$AUC),
  data.frame(cohort = "GSE72094", years = round(res2$times/365), AUC = res2$roc$AUC))
write.csv(val_auc, "outputs/Tables/Table_S4_validation_AUC.csv", row.names = FALSE)
cat("Validation AUC lưu ở Table_S4.\n")

# =============================================================================
# PHẦN B — DRUG SENSITIVITY (oncoPredict, Figure 3D)
# =============================================================================
# Ý tưởng: dùng oncoPredict huấn luyện trên GDSC2 (cell line + IC50 đã biết),
# dự đoán IC50 cho từng bệnh nhân TCGA-LUAD, rồi so IC50 giữa High vs Low risk.
# Thuốc mà nhóm High risk NHẠY hơn (IC50 thấp hơn) -> gợi ý điều trị.
#
# CHUẨN BỊ (làm 1 lần): tải training data GDSC2 của oncoPredict, đặt vào
#   outputs/oncoPredict_training/ :  GDSC2_Expr.rds  &  GDSC2_Res.rds
#   (hướng dẫn: vignette oncoPredict / DepMap portal). KHÔNG commit file lớn này.
GDSC_DIR <- "outputs/oncoPredict_training"

run_drug_sensitivity <- function() {
  expr_path <- file.path(GDSC_DIR, "GDSC2_Expr.rds")
  res_path  <- file.path(GDSC_DIR, "GDSC2_Res.rds")
  if (!file.exists(expr_path) || !file.exists(res_path)) {
    message("Chưa có training GDSC2 -> BỎ QUA Part B. Tải file rồi chạy lại.")
    return(invisible(NULL))
  }
  GDSC2_Expr <- readRDS(expr_path)
  GDSC2_Res  <- readRDS(res_path)

  # Ma trận biểu hiện TCGA-LUAD (gene × sample) từ STEP 1
  tcga_expr <- step1$expr_mat[, step2$cox_df |> rownames() %||% colnames(step1$expr_mat)]
  tcga_expr <- as.matrix(step1$expr_mat)   # gene × sample (log2-CPM)

  calcPhenotype(
    trainingExprData = GDSC2_Expr, trainingPtype = GDSC2_Res,
    testExprData = tcga_expr, batchCorrect = "eb",
    powerTransformPhenotype = TRUE, removeLowVaryingGenes = 0.2,
    minNumSamples = 10, printOutput = TRUE, removeLowVaringGenesFrom = "homogenizeData")
  # calcPhenotype ghi ra ./calcPhenotype_Output/DrugPredictions.csv
  pred <- read.csv("calcPhenotype_Output/DrugPredictions.csv", row.names = 1,
                   check.names = FALSE)
  pred
}

drug_pred <- run_drug_sensitivity()

if (!is.null(drug_pred)) {
  # Gắn risk_group từ TRAINING (STEP 2) cho các bệnh nhân TCGA
  rg <- step2$cox_df[, "risk_group", drop = FALSE]
  common <- intersect(rownames(drug_pred), rownames(rg))
  dp <- drug_pred[common, , drop = FALSE]
  grp <- rg[common, "risk_group"]

  # Wilcoxon mỗi thuốc: IC50 High vs Low + BH-FDR
  drug_stat <- lapply(colnames(dp), function(d) {
    x <- dp[[d]]
    w <- suppressWarnings(wilcox.test(x ~ grp))
    data.frame(drug = d,
               median_High = median(x[grp == "High"], na.rm = TRUE),
               median_Low  = median(x[grp == "Low"],  na.rm = TRUE),
               p = w$p.value)
  }) |> bind_rows() |>
    mutate(padj = p.adjust(p, method = "BH"),
           dir  = ifelse(median_High < median_Low,
                         "High risk nhạy hơn", "Low risk nhạy hơn")) |>
    arrange(padj)
  write.csv(drug_stat, "outputs/Tables/Table_S5_drug_sensitivity.csv", row.names = FALSE)

  # Figure 3D — boxplot IC50 cho DRUG_TOP_N thuốc khác biệt nhất
  top_drugs <- head(drug_stat$drug, DRUG_TOP_N)
  plot_df <- dp[, top_drugs, drop = FALSE] |>
    rownames_to_column("sample") |>
    mutate(risk = grp) |>
    pivot_longer(all_of(top_drugs), names_to = "drug", values_to = "IC50")
  # dán FDR vào tên facet để figure tự sinh số
  lab_map <- setNames(sprintf("%s\n(FDR=%.1e)", drug_stat$drug, drug_stat$padj),
                      drug_stat$drug)
  plot_df$drug <- factor(lab_map[plot_df$drug], levels = lab_map[top_drugs])

  p_drug <- ggplot(plot_df, aes(risk, IC50, fill = risk)) +
    geom_boxplot(outlier.size = 0.4, width = 0.6) +
    facet_wrap(~ drug, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c("Low" = "#2c7bb6", "High" = "#c0392b")) +
    theme_bw(base_size = 11) +
    theme(legend.position = "top", strip.text = element_text(size = 8)) +
    labs(title = "Figure 3D — Predicted drug sensitivity (GDSC2) by risk group",
         subtitle = "Estimated IC50 (oncoPredict); lower = more sensitive; Wilcoxon BH-FDR",
         x = NULL, y = "Estimated log IC50")
  ggsave("outputs/Figures/Fig3D_drug_sensitivity_boxplots.png",
         p_drug, width = 12, height = 8, dpi = 300)
  cat("Drug sensitivity: Table_S5 + Fig3D đã lưu.\n")
}

# =============================================================================
# 7. Lưu objects + sessionInfo
# =============================================================================
saveRDS(
  list(beta_vec = beta_vec,
       coh1 = coh1, coh2 = coh2,
       res_GSE31210 = res1, res_GSE72094 = res2,
       val_auc = val_auc,
       drug_stat = if (exists("drug_stat")) drug_stat else NULL,
       params = list(ROC_TIMES = ROC_TIMES, FDR_ALPHA = FDR_ALPHA)),
  "outputs/Rdata/STEP3_objects.rds")

writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_STEP3.txt")

cat("\n========================================\n")
cat(" STEP 3 completed.\n")
cat(" Validation cohorts: GSE31210, GSE72094\n")
cat(" 5-yr AUC:", paste(sprintf("%s=%.3f", val_auc$cohort,
      val_auc$AUC)[val_auc$years == 5], collapse = " | "), "\n")
cat(" Figures -> outputs/Figures/  | Tables -> outputs/Tables/\n")
cat("========================================\n")

# `%||%` helper (nếu chưa có từ package khác)
`%||%` <- function(a, b) if (is.null(a)) b else a
