# =============================================================================
# STEP 5 FINAL — Genomic landscape + Immunotherapy prediction + Single-cell
# Cuproptosis-related gene (CRG) prognostic signature in TCGA-LUAD
#
# Bản FINAL (2026-10-28), dựa trên STEP5_pipeline.R của Thanh (2026-10-24).
# ĐÃ CHẠY ĐƯỢC rồi, nên bản này KHÔNG viết lại logic — chỉ sửa 10 điểm
# làm sai số/sai nhãn trên figure và bổ sung bảng số để audit được.
#
# NHẬT KÝ SỬA (xem chi tiết trong phần FIX#n ở từng chỗ):
#  FIX#1  Đường dẫn: 1 biến PROJ_ROOT thay vì 2 absolute path của máy Thanh.
#  FIX#2  Kiểm tra risk_group đúng level "Low"/"High" (trước đây sai tên cột
#         là ra NA âm thầm, toàn bộ Phần A/B sẽ sai mà không báo lỗi).
#  FIX#3  Chọn aliquot cho mỗi bệnh nhân theo SỐ VARIANT NHIỀU NHẤT thay vì
#         thứ tự alphabet (alphabet là ngẫu nhiên về mặt sinh học).
#  FIX#4  Fig7C: nhãn in ODDS RATIO THẬT (Inf thì ghi "Inf"), KHÔNG in giá trị
#         đã bị clamp = 50. Bản cũ in "OR 50.00" cho 4 gene mà OR thật là Inf.
#  FIX#5  Fig7C: FDR in bằng %.2g (bản cũ %.2f -> STK11 hiện "FDR = 0.00").
#  FIX#6  Fig7C: subtitle ghi rõ SỐ GENE ĐƯỢC TEST + số dương tính kỳ vọng do
#         ngẫu nhiên + số gene qua BH-FDR. (353 gene p<0.05 là gần mức ngẫu
#         nhiên, chỉ 1 gene qua FDR -> không được viết như 353 phát hiện.)
#  FIX#7  Thêm phân tích GIẢ THUYẾT-DẪN trên panel 20 driver gene LUAD, BH chỉ
#         trên 20 test -> FDR có nghĩa, dùng cho Results; scan toàn genome
#         giữ lại nhưng dán nhãn exploratory.
#  FIX#8  TIDE: xuất MEDIAN + IQR + CHIỀU (High cao hơn hay thấp hơn) cho từng
#         score; đưa Fisher test của Responder vào CÙNG họ BH (4 test).
#  FIX#9  Mọi con số trước đây chỉ in ra console (Wilcoxon p, Spearman rho,
#         Cox HR, số mẫu khớp) nay ghi vào Table_S9d + STEP5_run_log.txt.
#  FIX#10 Fig8: tách tên file Fig8A (CancerSEA) vs Fig8B (TISCH2) — bản cũ hai
#         hình khác nhau cùng mang tiêu đề "Figure 8A".
# =============================================================================

# ---- 0. Môi trường ----------------------------------------------------------
stopifnot(getRversion() >= "4.2.0")

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(tidyr)
})

# Chạy script từ BÊN TRONG thư mục Step_5/ (giống Step_2/3/4)
dir.create("outputs/Figures/vector", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Tables",          recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Rdata",           recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/TIDE",            recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/CancerSEA",       recursive = TRUE, showWarnings = FALSE)

set.seed(123)
FDR_ALPHA <- 0.05
PALETTE   <- c(Low = "#2c7bb6", High = "#c0392b")   # khớp STEP2/3/4

# ---- CONFIG -----------------------------------------------------------------
# FIX#1: khai báo 1 chỗ duy nhất. Đổi PROJ_ROOT là xong, không phải sửa 2 dòng
# absolute path (bản cũ hard-code /Users/aithanh/... -> người khác chạy là lỗi).
PROJ_ROOT <- "/Users/aithanh/Cuproptosis_LUAD"

STEP1_RDS <- file.path(PROJ_ROOT, "Step1/outputs_step1_new/Rdata/STEP1_objects.rds")
STEP2_RDS <- file.path(PROJ_ROOT, "Step_2/outputs/Rdata/STEP2_objects.rds")

MAF_RDS       <- "outputs/Rdata/maf_raw.rds"
GDC_DIR       <- "GDCdata"
TIDE_IN       <- "outputs/TIDE/TIDE_input_expression.txt"
TIDE_OUT      <- "outputs/TIDE/TIDE_output.txt"
CANCERSEA_IN  <- "outputs/CancerSEA/sel_genes_for_CancerSEA.txt"
CANCERSEA_CSV <- "outputs/CancerSEA/CancerSEA_gene_state_correlation.csv"
TISCH_CSV     <- "outputs/CancerSEA/TISCH2_gene_celltype_expression.csv"
SCRNA_RDS     <- NA_character_
CELLTYPE_COL  <- "cell_type"

TMB_EXOME_MB <- 38    # Mb vùng exome dùng để chuẩn hoá TMB.
# Ghi trong Methods: "TMB = số biến thể non-synonymous / 38 Mb" và trích nguồn
# (Chalmers et al. Genome Med 2017). Đừng để trống con số này khi reviewer hỏi.
NONSYN_CLASSES <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
                    "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
                    "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")

# FIX#7: panel driver LUAD cho phân tích giả thuyết-dẫn (BH trên 20 test).
LUAD_DRIVERS <- c("TP53", "KRAS", "STK11", "KEAP1", "EGFR", "BRAF", "MET",
                  "NF1", "RIT1", "RBM10", "SMARCA4", "SETD2", "ARID1A",
                  "U2AF1", "RB1", "CDKN2A", "ATM", "PIK3CA", "ERBB2", "MGA")

# ---- Log: mọi con số đều phải ra file, không chỉ ra console ------------------
# FIX#9
.LOG <- character(0)
say <- function(...) {
  txt <- paste0(...)
  .LOG <<- c(.LOG, txt)
  cat(txt, "\n", sep = "")
}
.STATS <- list()
stat_add <- function(metric, value, note = "") {
  .STATS[[length(.STATS) + 1]] <<- data.frame(metric = metric,
                                              value = as.character(value),
                                              note = note,
                                              stringsAsFactors = FALSE)
  invisible(NULL)
}

# ---- Helper: lưu figure ra CẢ png (300 dpi) LẪN pdf vector ------------------
pdf_dev <- if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf else grDevices::pdf

save_gg <- function(p, name, w, h) {
  ggsave(file.path("outputs/Figures", paste0(name, ".png")), p,
         width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path("outputs/Figures/vector", paste0(name, ".pdf")), p,
         width = w, height = h, device = pdf_dev)
  invisible(NULL)
}

save_base <- function(name, w, h, draw) {
  png(file.path("outputs/Figures", paste0(name, ".png")),
      width = w, height = h, units = "in", res = 300, bg = "white")
  ok <- tryCatch({ draw(); TRUE }, error = function(e) { message("png: ", conditionMessage(e)); FALSE })
  dev.off()
  if (!ok) return(invisible(NULL))
  pdf_dev(file.path("outputs/Figures/vector", paste0(name, ".pdf")), width = w, height = h)
  tryCatch(draw(), error = function(e) message("pdf: ", conditionMessage(e)))
  dev.off()
  invisible(NULL)
}

theme_step5 <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5),
          panel.grid.minor = element_blank())
}

fmt_p <- function(p) ifelse(is.na(p), "NA", formatC(p, format = "g", digits = 3))

# ---- 1. Nạp objects đóng băng từ STEP 1 + STEP 2 ----------------------------
for (f in c(STEP1_RDS, STEP2_RDS))
  if (!file.exists(f)) stop("Không tìm thấy: ", f, "\n  -> sửa PROJ_ROOT ở phần CONFIG.")

step1 <- readRDS(STEP1_RDS)
step2 <- readRDS(STEP2_RDS)

expr_mat  <- as.matrix(step1$expr_mat)   # gene x sample, log2-CPM
cox_df    <- step2$cox_df                # rownames = sample (full barcode)
sel_genes <- step2$sel_genes

stopifnot(all(c("risk_score", "risk_group", "OS_time", "OS_event") %in% colnames(cox_df)))

# FIX#2: bản cũ gọi factor(..., levels = c("Low","High")) mà KHÔNG kiểm tra.
# Nếu STEP2 ghi "low"/"high" thì toàn bộ risk_group thành NA, Wilcoxon/Fisher
# vẫn chạy trên dữ liệu rỗng và không ai biết. Giờ fail-fast.
rg_lv <- sort(unique(as.character(cox_df$risk_group)))
if (!identical(rg_lv, c("High", "Low")))
  stop("risk_group trong STEP2 có giá trị: ", paste(rg_lv, collapse = "/"),
       " — phải là 'Low'/'High'. Sửa STEP2 hoặc map lại ở đây.")

say("== STEP 5 FINAL ==  ", format(Sys.time(), "%Y-%m-%d %H:%M"))
say("Samples (STEP2): ", nrow(cox_df),
    " | Low = ", sum(cox_df$risk_group == "Low"),
    " | High = ", sum(cox_df$risk_group == "High"))
say("Signature genes (", length(sel_genes), "): ", paste(sel_genes, collapse = ", "))
stat_add("n_samples_step2", nrow(cox_df))
stat_add("n_signature_genes", length(sel_genes), paste(sel_genes, collapse = "|"))

# sample (full barcode) -> risk_group  (dùng cho TIDE, giống STEP 4)
risk_lookup <- cox_df |> rownames_to_column("sample") |> select(sample, risk_group)

# patient-level (barcode 12 ký tự) -> dùng để khớp với MAF
dup_n <- sum(duplicated(substr(rownames(cox_df), 1, 12)))
if (dup_n > 0)
  say("CHÚ Ý: ", dup_n, " mẫu trùng bệnh nhân bị loại khi gộp về barcode12 ",
      "(STEP1 đã quy ước 1 BN/1 mẫu u nên con số này phải = 0).")
stat_add("n_duplicate_patient_samples_dropped", dup_n)

risk_pat <- cox_df |>
  rownames_to_column("sample") |>
  mutate(barcode12 = substr(sample, 1, 12)) |>
  distinct(barcode12, .keep_all = TRUE) |>
  select(barcode12, risk_group, risk_score, OS_time, OS_event) |>
  mutate(risk_group = factor(as.character(risk_group), levels = c("Low", "High")))

# =============================================================================
# PHẦN A — Somatic mutation landscape + TMB (Figure 7A–7C, 7E)
# =============================================================================
get_maf <- function() {
  if (file.exists(MAF_RDS)) {
    say(">> Dùng MAF đã cache: ", MAF_RDS); return(readRDS(MAF_RDS))
  }
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    message("Chưa cài 'TCGAbiolinks' và chưa có ", MAF_RDS, " -> BỎ QUA Phần A."); return(NULL)
  }
  maf <- tryCatch({
    q <- TCGAbiolinks::GDCquery(
      project = "TCGA-LUAD",
      data.category = "Simple Nucleotide Variation",
      data.type = "Masked Somatic Mutation",
      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
    TCGAbiolinks::GDCdownload(q, directory = GDC_DIR)
    TCGAbiolinks::GDCprepare(q, directory = GDC_DIR)
  }, error = function(e) {
    message("Tải MAF lỗi (", conditionMessage(e), "). Chạy trên máy em, ",
            "hoặc đặt sẵn maf_raw.rds vào outputs/Rdata/ -> BỎ QUA Phần A."); NULL
  })
  if (!is.null(maf)) saveRDS(maf, MAF_RDS)
  maf
}

run_part_A <- function() {
  if (!requireNamespace("maftools", quietly = TRUE)) {
    message("Chưa cài 'maftools' -> BỎ QUA Phần A."); return(NULL) }
  maf_raw <- get_maf()
  if (is.null(maf_raw)) return(NULL)

  maf_df <- as.data.frame(maf_raw)
  need <- c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")
  if (!all(need %in% colnames(maf_df))) {
    message("MAF thiếu cột: ", paste(setdiff(need, colnames(maf_df)), collapse = ", "),
            " -> BỎ QUA Phần A."); return(NULL)
  }

  # --- A1. Khớp MAF <-> risk_group bằng barcode 12 ký tự ---------------------
  maf_df$barcode12 <- substr(maf_df$Tumor_Sample_Barcode, 1, 12)

  # FIX#3: bệnh nhân có >1 aliquot -> bản cũ lấy aliquot đầu theo ALPHABET.
  # Giờ lấy aliquot có NHIỀU VARIANT NHẤT (đại diện tốt hơn cho khối u, và
  # vẫn deterministic nhờ tie-break bằng tên).
  aliquot_rank <- maf_df |>
    count(barcode12, Tumor_Sample_Barcode, name = "n_var") |>
    arrange(barcode12, desc(n_var), Tumor_Sample_Barcode)
  keep_tsb <- aliquot_rank |> distinct(barcode12, .keep_all = TRUE)
  n_multi <- sum(table(aliquot_rank$barcode12) > 1)
  say(">> Bệnh nhân có nhiều aliquot trong MAF: ", n_multi,
      " (đã giữ aliquot nhiều variant nhất)")
  stat_add("n_patients_multi_aliquot_maf", n_multi)

  maf_df <- maf_df |>
    filter(Tumor_Sample_Barcode %in% keep_tsb$Tumor_Sample_Barcode,
           barcode12 %in% risk_pat$barcode12)

  n_match <- length(unique(maf_df$barcode12))
  say(">> Khớp MAF: ", n_match, " / ", nrow(risk_pat),
      " bệnh nhân (STEP 2) có dữ liệu đột biến (",
      sprintf("%.1f%%", 100 * n_match / nrow(risk_pat)), ")")
  stat_add("n_patients_with_maf", n_match, paste0("of ", nrow(risk_pat)))
  if (n_match < 0.5 * nrow(risk_pat))
    warning("MAF: <50% bệnh nhân khớp — kiểm tra barcode / project TCGA-LUAD.")
  stopifnot(n_match > 0)

  maf_df$Tumor_Sample_Barcode <- maf_df$barcode12       # patient-level
  clin_maf <- risk_pat |>
    filter(barcode12 %in% maf_df$Tumor_Sample_Barcode) |>
    transmute(Tumor_Sample_Barcode = barcode12,
              risk_group = as.character(risk_group),
              risk_score, OS_time, OS_event)

  maf_obj  <- maftools::read.maf(maf = maf_df, clinicalData = clin_maf, verbose = FALSE)
  high_tsb <- clin_maf$Tumor_Sample_Barcode[clin_maf$risk_group == "High"]
  low_tsb  <- clin_maf$Tumor_Sample_Barcode[clin_maf$risk_group == "Low"]
  maf_high <- maftools::subsetMaf(maf_obj, tsb = high_tsb, verbose = FALSE)
  maf_low  <- maftools::subsetMaf(maf_obj, tsb = low_tsb,  verbose = FALSE)
  say(">> Sau khi tách: High = ", length(high_tsb), " | Low = ", length(low_tsb), " bệnh nhân")
  stat_add("n_high_maf", length(high_tsb)); stat_add("n_low_maf", length(low_tsb))

  # --- A2. TMB (non-synonymous / 38 Mb) -> Figure 7A ------------------------
  tmb_df <- maf_df |>
    filter(Variant_Classification %in% NONSYN_CLASSES) |>
    count(Tumor_Sample_Barcode, name = "n_nonsyn") |>
    right_join(clin_maf, by = "Tumor_Sample_Barcode") |>
    mutate(n_nonsyn = ifelse(is.na(n_nonsyn), 0L, n_nonsyn),
           TMB = n_nonsyn / TMB_EXOME_MB,
           risk_group = factor(risk_group, levels = c("Low", "High")))
  write.csv(tmb_df, "outputs/Tables/Table_S9b_TMB_by_patient.csv", row.names = FALSE)

  wp  <- wilcox.test(TMB ~ risk_group, data = tmb_df)$p.value
  rho <- suppressWarnings(cor.test(tmb_df$risk_score, tmb_df$TMB,
                                   method = "spearman", exact = FALSE))

  # FIX#9: median/IQR theo nhóm — con số phải đưa vào Results, không chỉ boxplot
  tmb_sum <- tmb_df |>
    group_by(risk_group) |>
    summarise(n = n(),
              median_TMB = median(TMB),
              Q1 = quantile(TMB, .25), Q3 = quantile(TMB, .75), .groups = "drop")
  print(as.data.frame(tmb_sum))
  for (i in seq_len(nrow(tmb_sum)))
    say(sprintf("   TMB %-4s n = %d, median = %.2f (IQR %.2f-%.2f) mut/Mb",
                as.character(tmb_sum$risk_group[i]), tmb_sum$n[i],
                tmb_sum$median_TMB[i], tmb_sum$Q1[i], tmb_sum$Q3[i]))
  say(sprintf(">> TMB High vs Low: Wilcoxon p = %s | Spearman(risk_score, TMB) rho = %.3f (p = %s)",
              fmt_p(wp), rho$estimate, fmt_p(rho$p.value)))
  stat_add("TMB_wilcoxon_p", fmt_p(wp))
  stat_add("TMB_median_Low",  sprintf("%.2f", tmb_sum$median_TMB[tmb_sum$risk_group == "Low"]))
  stat_add("TMB_median_High", sprintf("%.2f", tmb_sum$median_TMB[tmb_sum$risk_group == "High"]))
  stat_add("TMB_spearman_rho", sprintf("%.3f", rho$estimate), paste0("p = ", fmt_p(rho$p.value)))

  n_tab <- table(tmb_df$risk_group)
  p7a <- ggplot(tmb_df, aes(risk_group, log2(TMB + 1), fill = risk_group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.35) +
    scale_fill_manual(values = PALETTE) +
    theme_step5() + theme(legend.position = "none") +
    labs(title = "Figure 7A — Tumor mutation burden by risk group",
         subtitle = sprintf("Wilcoxon p = %s | n(Low) = %d, n(High) = %d | Spearman rho = %.2f (p = %s)",
                            fmt_p(wp), n_tab[["Low"]], n_tab[["High"]], rho$estimate, fmt_p(rho$p.value)),
         x = NULL, y = "log2(TMB + 1)  [non-synonymous / Mb]")
  save_gg(p7a, "Fig7A_TMB_by_risk", 5, 5)
  say(">> Figure 7A saved.")

  # --- A3. Oncoplot đồng thời High vs Low -> Figure 7B ----------------------
  top_genes <- head(maftools::getGeneSummary(maf_obj)$Hugo_Symbol, 20)
  save_base("Fig7B_coOncoplot_high_vs_low", 12, 8, function() {
    maftools::coOncoplot(m1 = maf_high, m2 = maf_low, genes = top_genes,
                         m1Name = "High risk", m2Name = "Low risk",
                         removeNonMutated = FALSE)
  })
  say(">> Figure 7B saved.")

  # --- A4. mafCompare -> Figure 7C + Table S9 ------------------------------
  mc_res <- NULL; panel_res <- NULL
  tryCatch({
    mc <- maftools::mafCompare(m1 = maf_high, m2 = maf_low,
                               m1Name = "High", m2Name = "Low", minMut = 5)
    mc_res <- as.data.frame(mc$results)
    mc_res$padj_BH <- p.adjust(mc_res$pval, method = "BH")
    mc_res <- mc_res[order(mc_res$pval), ]
    write.csv(mc_res, "outputs/Tables/Table_S9_mafCompare_high_vs_low.csv", row.names = FALSE)

    n_tested   <- nrow(mc_res)
    n_p05      <- sum(mc_res$pval < 0.05, na.rm = TRUE)
    n_expected <- round(0.05 * n_tested)
    n_fdr      <- sum(mc_res$padj_BH < FDR_ALPHA, na.rm = TRUE)
    # FIX#6: đây là con số PHẢI viết trong Results. 353/5000 gene p<0.05 khi
    # kỳ vọng ngẫu nhiên ~250 KHÔNG phải "353 gene khác biệt".
    say(sprintf(">> mafCompare: %d gene được test (minMut >= 5) | p < 0.05: %d (kỳ vọng ngẫu nhiên ~%d) | BH-FDR < %.2f: %d",
                n_tested, n_p05, n_expected, FDR_ALPHA, n_fdr))
    stat_add("mafCompare_n_genes_tested", n_tested)
    stat_add("mafCompare_n_p_lt_0.05", n_p05, paste0("expected by chance ~", n_expected))
    stat_add("mafCompare_n_FDR_sig", n_fdr)
    if (n_fdr > 0) {
      say("   Gene qua BH-FDR: ",
          paste(mc_res$Hugo_Symbol[which(mc_res$padj_BH < FDR_ALPHA)], collapse = ", "))
      stat_add("mafCompare_FDR_sig_genes",
               paste(mc_res$Hugo_Symbol[which(mc_res$padj_BH < FDR_ALPHA)], collapse = "|"))
    }

    # FIX#7: panel giả thuyết-dẫn — BH chỉ trên 20 driver LUAD.
    panel_res <- mc_res[mc_res$Hugo_Symbol %in% LUAD_DRIVERS, ]
    if (nrow(panel_res) > 0) {
      panel_res$padj_BH_panel <- p.adjust(panel_res$pval, method = "BH")
      panel_res <- panel_res[order(panel_res$pval), ]
      write.csv(panel_res, "outputs/Tables/Table_S9e_LUAD_driver_panel.csv", row.names = FALSE)
      say(">> Panel 20 driver LUAD (BH trên ", nrow(panel_res), " test):")
      show_cols <- intersect(c("Hugo_Symbol", "High", "Low", "or", "pval", "padj_BH_panel"),
                             colnames(panel_res))
      print(head(panel_res[, show_cols], 10), row.names = FALSE)
      sig_panel <- panel_res$Hugo_Symbol[panel_res$padj_BH_panel < FDR_ALPHA]
      say("   Qua FDR trong panel: ",
          if (length(sig_panel)) paste(sig_panel, collapse = ", ") else "(không có)")
      stat_add("driver_panel_FDR_sig_genes",
               if (length(sig_panel)) paste(sig_panel, collapse = "|") else "none",
               paste0("BH over ", nrow(panel_res), " driver genes"))
    }

    # Forest plot: vị trí điểm bị clamp để vẽ được, NHƯNG nhãn in số thật.
    plot_df <- head(mc_res, 15) |>
      mutate(or_plot = pmin(pmax(ifelse(is.finite(or), or, 50), 0.05), 50),
             lo_plot = pmin(pmax(ifelse(is.finite(ci.low), ci.low, 0.05), 0.05), 50),
             hi_plot = pmin(pmax(ifelse(is.finite(ci.up),  ci.up,  50), 0.05), 50),
             # FIX#4 + FIX#5
             or_lab  = ifelse(is.finite(or), sprintf("%.2f", or), "Inf"),
             label   = sprintf("%s  (OR %s, p = %s, FDR = %s)",
                               Hugo_Symbol, or_lab, fmt_p(pval), fmt_p(padj_BH)),
             label   = factor(label, levels = rev(label)))
    p7c <- ggplot(plot_df, aes(or_plot, label)) +
      geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
      geom_errorbarh(aes(xmin = lo_plot, xmax = hi_plot), height = 0.25) +
      geom_point(aes(colour = or_plot > 1), shape = 18, size = 3.4) +
      scale_colour_manual(values = c(`TRUE` = "#c0392b", `FALSE` = "#2c7bb6"),
                          labels = c(`TRUE` = "Enriched in High", `FALSE` = "Enriched in Low"),
                          name = NULL) +
      scale_x_log10() +
      theme_step5() +
      labs(title = "Figure 7C — Differentially mutated genes, High vs Low risk",
           subtitle = sprintf("Fisher's exact test, %d genes tested; top 15 by p. p<0.05: %d (≈%d expected by chance); BH-FDR<%.2f: %d.\nPoints at the axis limits are capped for display; labels show the true OR.",
                              n_tested, n_p05, n_expected, FDR_ALPHA, n_fdr),
           x = "Odds ratio (log scale; High / Low)", y = NULL)
    save_gg(p7c, "Fig7C_mafCompare_forest", 11, 6.5)
    say(">> Figure 7C + Table S9 saved.")
  }, error = function(e) message("mafCompare lỗi (", conditionMessage(e), ") -> bỏ Fig7C/Table S9."))

  # --- A5. KM 4 nhóm TMB x risk -> Figure 7E -------------------------------
  cox_tab <- NULL
  tryCatch({
    km_df <- tmb_df |>
      mutate(tmb_group = factor(ifelse(TMB >= median(TMB), "TMB-high", "TMB-low"),
                                levels = c("TMB-low", "TMB-high")),
             grp4 = factor(paste(tmb_group, risk_group, sep = " & "),
                           levels = c("TMB-low & Low", "TMB-low & High",
                                      "TMB-high & Low", "TMB-high & High")))
    say(">> Số bệnh nhân mỗi nhóm TMB x risk:")
    print(table(km_df$grp4))
    for (g in levels(km_df$grp4))
      say("   ", g, " = ", sum(km_df$grp4 == g))

    fit4 <- survfit(Surv(OS_time, OS_event) ~ grp4, data = km_df)
    pv4  <- surv_pvalue(fit4, data = km_df)$pval

    # risk_group còn giá trị bổ sung sau khi hiệu chỉnh TMB?
    cx <- summary(coxph(Surv(OS_time, OS_event) ~ risk_group + tmb_group, data = km_df))$coefficients
    cox_tab <- data.frame(term = rownames(cx),
                          HR = cx[, "exp(coef)"],
                          p  = cx[, "Pr(>|z|)"],
                          row.names = NULL)
    write.csv(cox_tab, "outputs/Tables/Table_S9f_Cox_risk_plus_TMB.csv", row.names = FALSE)
    say(">> Cox (risk_group + tmb_group) — risk_group có độc lập với TMB?")
    for (i in seq_len(nrow(cox_tab))) {
      say(sprintf("   %-22s HR = %.2f, p = %s", cox_tab$term[i], cox_tab$HR[i], fmt_p(cox_tab$p[i])))
      stat_add(paste0("Cox_", cox_tab$term[i], "_HR"), sprintf("%.2f", cox_tab$HR[i]),
               paste0("p = ", fmt_p(cox_tab$p[i])))
    }
    stat_add("KM4_logrank_p", fmt_p(pv4))

    km4 <- ggsurvplot(fit4, data = km_df, pval = TRUE, risk.table = TRUE,
                      palette = c("#9ecae1", "#fb6a4a", "#3182bd", "#a50f15"),
                      legend.title = "", legend.labs = levels(km_df$grp4),
                      xlab = "Time (days)", ggtheme = theme_step5(),
                      title = sprintf("Figure 7E — OS by TMB x risk group (log-rank p = %s)", fmt_p(pv4)))
    png("outputs/Figures/Fig7E_KM_TMB_x_risk.png", width = 8, height = 7,
        units = "in", res = 300, bg = "white"); print(km4); dev.off()
    pdf_dev("outputs/Figures/vector/Fig7E_KM_TMB_x_risk.pdf", width = 8, height = 7)
    print(km4); dev.off()
    say(">> Figure 7E saved.")
  }, error = function(e) { if (dev.cur() > 1) dev.off()
    message("KM TMB x risk lỗi (", conditionMessage(e), ") -> bỏ Fig7E.") })

  list(tmb_df = tmb_df, tmb_summary = tmb_sum, maf_compare = mc_res,
       driver_panel = panel_res, cox_risk_tmb = cox_tab,
       n_high = length(high_tsb), n_low = length(low_tsb))
}
res_A <- run_part_A()

# =============================================================================
# PHẦN B — Immunotherapy response (Figure 7D)
#   B1. TMB-as-biomarker = Fig7A (đã làm ở Phần A)
#   B2. TIDE (web tool): XUẤT input -> em upload -> NHẬP output -> Fig7D
# =============================================================================
run_tide <- function() {
  common <- intersect(colnames(expr_mat), risk_lookup$sample)
  say(">> TIDE: ", length(common), "/", nrow(risk_lookup), " mẫu có trong expr_mat")
  stat_add("n_samples_TIDE_input", length(common), paste0("of ", nrow(risk_lookup)))
  stopifnot(length(common) > 0)

  if (!file.exists(TIDE_IN)) {
    m <- expr_mat[, common, drop = FALSE]
    m <- m[apply(m, 1, sd) > 0, , drop = FALSE]
    m <- sweep(m, 1, rowMeans(m))      # average-centering theo cohort (yêu cầu của TIDE)
    write.table(m, TIDE_IN, sep = "\t", quote = FALSE, col.names = NA)
    say(">> Đã XUẤT input TIDE: ", TIDE_IN, " (", nrow(m), " gene x ", ncol(m), " mẫu)")
  }

  if (!file.exists(TIDE_OUT)) {
    message("Chưa có ", TIDE_OUT, ".\n",
            "  -> Upload ", TIDE_IN, " lên http://tide.dfci.harvard.edu (cancer type = NSCLC),\n",
            "     hoặc chạy:  tidepy ", TIDE_IN, " -o ", TIDE_OUT, " -c NSCLC\n",
            "     rồi chạy lại script -> hiện BỎ QUA Fig7D.")
    return(invisible(NULL))
  }

  tide_sep <- if (grepl("\\.csv$", TIDE_OUT, ignore.case = TRUE)) "," else "\t"
  tide <- read.table(TIDE_OUT, sep = tide_sep, header = TRUE, row.names = 1,
                     check.names = FALSE, stringsAsFactors = FALSE,
                     quote = "\"", comment.char = "")
  tide$sample <- rownames(tide)
  if (!any(tide$sample %in% risk_lookup$sample)) tide$sample <- gsub("\\.", "-", tide$sample)

  merged <- inner_join(risk_lookup, tide, by = "sample") |>
    mutate(risk_group = factor(as.character(risk_group), levels = c("Low", "High")))
  say(">> Merge TIDE: ", nrow(merged), "/", nrow(risk_lookup), " mẫu khớp")
  stat_add("n_samples_TIDE_merged", nrow(merged), paste0("of ", nrow(risk_lookup)))
  if (nrow(merged) < 0.5 * nrow(risk_lookup))
    warning("TIDE: <50% mẫu khớp — kiểm tra tên mẫu trong TIDE_output.")
  stopifnot(nrow(merged) > 0)

  score_cols <- intersect(c("TIDE", "Dysfunction", "Exclusion"), colnames(merged))
  if (length(score_cols) == 0) {
    message("TIDE_output thiếu cột TIDE/Dysfunction/Exclusion -> bỏ Fig7D."); return(NULL) }

  long <- merged |> pivot_longer(all_of(score_cols), names_to = "score", values_to = "value")

  # FIX#8: median + IQR + CHIỀU. Bản cũ chỉ có p -> đọc figure không biết nhóm
  # nào cao hơn, mà chiều mới là điều người đọc paper cần.
  stat <- long |>
    group_by(score) |>
    summarise(n_low      = sum(risk_group == "Low"),
              n_high     = sum(risk_group == "High"),
              median_low  = median(value[risk_group == "Low"],  na.rm = TRUE),
              median_high = median(value[risk_group == "High"], na.rm = TRUE),
              IQR_low     = IQR(value[risk_group == "Low"],  na.rm = TRUE),
              IQR_high    = IQR(value[risk_group == "High"], na.rm = TRUE),
              p = wilcox.test(value ~ risk_group)$p.value, .groups = "drop") |>
    mutate(direction = ifelse(median_high > median_low, "higher in High risk",
                              "higher in Low risk"),
           test = "Wilcoxon rank-sum")

  # Responder: đưa VÀO CÙNG họ BH với 3 score (4 test), không để riêng lẻ.
  fisher_p <- NA_real_; p7d2 <- NULL; resp_tab <- NULL
  if ("Responder" %in% colnames(merged)) {
    merged$responder <- toupper(as.character(merged$Responder)) %in% c("TRUE", "T", "1", "YES")
    # factor cố định 2 level: nếu cả cohort cùng 1 nhãn thì cột "TRUE"/"FALSE"
    # vẫn tồn tại (bản cũ sẽ lỗi "subscript out of bounds").
    merged$responder <- factor(merged$responder, levels = c(FALSE, TRUE),
                               labels = c("FALSE", "TRUE"))
    resp_tab <- table(risk = merged$risk_group, responder = merged$responder)
    fisher_p <- tryCatch(fisher.test(resp_tab)$p.value, error = function(e) NA_real_)
    say(">> Bảng 2x2 Responder (TIDE) theo risk group:")
    print(resp_tab)
    for (r in rownames(resp_tab))
      say(sprintf("   %-5s responder = %d / %d (%.1f%%)", r,
                  resp_tab[r, "TRUE"], sum(resp_tab[r, ]),
                  100 * resp_tab[r, "TRUE"] / sum(resp_tab[r, ])))
    stat <- bind_rows(stat, data.frame(
      score = "Responder (proportion)", n_low = sum(resp_tab["Low", ]),
      n_high = sum(resp_tab["High", ]),
      median_low = NA, median_high = NA, IQR_low = NA, IQR_high = NA,
      p = fisher_p,
      direction = ifelse(is.na(fisher_p), "NA",
                         ifelse(resp_tab["High", "TRUE"] / sum(resp_tab["High", ]) >
                                resp_tab["Low", "TRUE"] / sum(resp_tab["Low", ]),
                                "more responders in High risk", "more responders in Low risk")),
      test = "Fisher's exact"))
  }
  stat$padj <- p.adjust(stat$p, method = "BH")
  write.csv(stat, "outputs/Tables/Table_S9c_TIDE_score_stat.csv", row.names = FALSE)
  for (i in seq_len(nrow(stat)))
    say(sprintf("   TIDE %-22s p = %-10s BH-FDR = %-10s  %s",
                stat$score[i], fmt_p(stat$p[i]), fmt_p(stat$padj[i]), stat$direction[i]))
  for (i in seq_len(nrow(stat)))
    stat_add(paste0("TIDE_", gsub("[^A-Za-z]", "_", stat$score[i])),
             paste0("BH-FDR = ", fmt_p(stat$padj[i])), stat$direction[i])

  ann <- stat |> filter(score %in% score_cols) |>
    mutate(label = sprintf("BH-FDR = %s\n%s", fmt_p(padj), direction))
  p7d1 <- ggplot(long, aes(risk_group, value, fill = risk_group)) +
    geom_boxplot(outlier.size = 0.5) +
    geom_text(data = ann, aes(x = 1.5, y = Inf, label = label),
              vjust = 1.2, size = 3, inherit.aes = FALSE) +
    facet_wrap(~ score, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = PALETTE) + theme_step5() + theme(legend.position = "none") +
    labs(title = "Figure 7D — Predicted immunotherapy response (TIDE)",
         subtitle = "Wilcoxon rank-sum, BH-FDR across 4 tests. TIDE scores are transcriptome-derived,\nlike the risk score itself — association is hypothesis-generating, not independent validation.",
         x = NULL, y = "Score")

  if (!is.null(resp_tab)) {
    bar_df <- merged |> count(risk_group, responder) |> group_by(risk_group) |>
      mutate(pct = 100 * n / sum(n)) |> ungroup()
    p7d2 <- ggplot(bar_df, aes(risk_group, pct, fill = responder)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = c(`TRUE` = "#27ae60", `FALSE` = "grey70"),
                        labels = c(`TRUE` = "Responder", `FALSE` = "Non-responder"), name = NULL) +
      theme_step5() +
      labs(subtitle = sprintf("Responder proportion; Fisher p = %s (BH-FDR = %s)",
                              fmt_p(fisher_p),
                              fmt_p(stat$padj[stat$score == "Responder (proportion)"])),
           x = NULL, y = "% of patients")
  }

  if (!is.null(p7d2) && requireNamespace("patchwork", quietly = TRUE)) {
    save_gg(patchwork::wrap_plots(p7d1, p7d2, widths = c(length(score_cols), 1.2)),
            "Fig7D_TIDE_immunotherapy", 4 * length(score_cols) + 3, 5.5)
  } else {
    save_gg(p7d1, "Fig7D_TIDE_scores", 4 * length(score_cols), 5.5)
    if (!is.null(p7d2)) save_gg(p7d2, "Fig7D_TIDE_responder", 4.5, 5)
  }
  say(">> Figure 7D saved.")
  list(merged = merged, stat = stat, fisher_p = fisher_p, resp_tab = resp_tab)
}
res_B <- tryCatch(run_tide(), error = function(e) { message("TIDE lỗi: ", conditionMessage(e)); NULL })

# =============================================================================
# PHẦN C — Single-cell confirmation (Figure 8)
#   C1. CancerSEA (web, nhẹ): XUẤT gene list -> tra tay -> NHẬP bảng -> Fig8A
#   C1b. TISCH2 (phương án dự phòng nếu CancerSEA không vào được)  -> Fig8B
#   C2. (tuỳ chọn) scRNA thật bằng Seurat AddModuleScore           -> Fig8C/8D
# =============================================================================
writeLines(sel_genes, CANCERSEA_IN)

# FIX#10 + tiện cho Thanh: nếu chưa có bảng, script tạo luôn FILE MẪU đúng
# định dạng để em chỉ việc điền số, khỏi đoán tên cột.
make_template <- function(path, cols, example) {
  if (file.exists(path)) return(invisible(NULL))
  tmpl <- sub("\\.csv$", "_TEMPLATE.csv", path)
  if (!file.exists(tmpl)) {
    writeLines(c(paste(cols, collapse = ","), example), tmpl)
    message("Đã tạo file mẫu: ", tmpl, " — điền số rồi đổi tên thành ", basename(path))
  }
  invisible(NULL)
}

run_cancersea <- function() {
  if (!file.exists(CANCERSEA_CSV)) {
    make_template(CANCERSEA_CSV, c("gene", "state", "correlation", "fdr"),
                  c("DLD,Cell cycle,0.42,0.001", "DLD,Hypoxia,-0.10,0.30"))
    message("Chưa có ", CANCERSEA_CSV, ".\n",
            "  -> Tra từng gene trong ", CANCERSEA_IN,
            " tại http://biocc.hrbmu.edu.cn/CancerSEA/home.jsp\n",
            "     (404 thì thử http://202.97.205.69/CancerSEA/ ; vào không được -> dùng TISCH2)\n",
            "     cancer = LUAD/NSCLC, long format (gene,state,correlation[,fdr]) -> hiện BỎ QUA Fig8A.")
    return(invisible(NULL))
  }
  cs <- read.csv(CANCERSEA_CSV, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(c("gene", "state", "correlation") %in% colnames(cs)))
  miss <- setdiff(sel_genes, unique(cs$gene))
  if (length(miss) > 0) say("CancerSEA: chưa có dữ liệu cho gene: ", paste(miss, collapse = ", "))
  say(">> CancerSEA: ", length(unique(cs$gene)), " gene x ", length(unique(cs$state)), " functional state")

  state_order <- cs |> group_by(state) |>
    summarise(m = mean(correlation, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(m)) |> pull(state)
  cs$gene  <- factor(cs$gene, levels = rev(intersect(sel_genes, unique(cs$gene))))
  cs$state <- factor(cs$state, levels = state_order)
  cs$star  <- if ("fdr" %in% colnames(cs)) ifelse(!is.na(cs$fdr) & cs$fdr < FDR_ALPHA, "*", "") else ""

  p8a <- ggplot(cs, aes(state, gene, fill = correlation)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = paste0(sprintf("%.2f", correlation), star)), size = 3) +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#c0392b",
                         limits = c(-1, 1), name = "Correlation") +
    theme_step5() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Figure 8A — Signature genes vs. functional states (CancerSEA)",
         subtitle = "Single-cell correlation, LUAD; * FDR < 0.05 (nếu có cột fdr)", x = NULL, y = NULL)
  save_gg(p8a, "Fig8A_CancerSEA_heatmap", 10, 2 + 0.6 * nlevels(cs$gene))
  write.csv(cs, "outputs/Tables/Table_S10_CancerSEA_correlation.csv", row.names = FALSE)
  say(">> Figure 8A saved. (t-SNE của gene tiêu biểu chụp từ web CancerSEA, ghép tay khi ráp figure.)")
  cs
}
res_C1 <- tryCatch(run_cancersea(), error = function(e) { message("CancerSEA lỗi: ", conditionMessage(e)); NULL })

run_tisch2 <- function() {
  if (!file.exists(TISCH_CSV)) {
    make_template(TISCH_CSV, c("gene", "cell_type", "expression"),
                  c("DLD,Malignant,1.85", "DLD,CD8T,0.60"))
    message("Chưa có ", TISCH_CSV, " -> BỎ QUA heatmap TISCH2 (phương án dự phòng).")
    return(invisible(NULL))
  }
  ts <- read.csv(TISCH_CSV, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(c("gene", "cell_type", "expression") %in% colnames(ts)))
  miss <- setdiff(sel_genes, unique(ts$gene))
  if (length(miss) > 0) say("TISCH2: chưa có dữ liệu cho gene: ", paste(miss, collapse = ", "))
  say(">> TISCH2: ", length(unique(ts$gene)), " gene x ", length(unique(ts$cell_type)), " cell type")

  ts$gene <- factor(ts$gene, levels = rev(intersect(sel_genes, unique(ts$gene))))
  p8b <- ggplot(ts, aes(cell_type, gene, fill = expression)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", expression)), size = 3) +
    scale_fill_gradient(low = "white", high = "#c0392b", name = "Expression") +
    theme_step5() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Figure 8B — Signature gene expression by cell type (TISCH2, scRNA-seq)",
         subtitle = "Supportive evidence of cellular expression; hypothesis-generating", x = NULL, y = NULL)
  # FIX#10: tên file + tiêu đề là 8B, không trùng 8A nữa.
  save_gg(p8b, "Fig8B_TISCH2_celltype_heatmap", 10, 2 + 0.6 * nlevels(ts$gene))
  write.csv(ts, "outputs/Tables/Table_S10b_TISCH2_expression.csv", row.names = FALSE)
  say(">> Figure 8B (TISCH2) saved.")
  ts
}
res_C1b <- tryCatch(run_tisch2(), error = function(e) { message("TISCH2 lỗi: ", conditionMessage(e)); NULL })

run_scrna <- function() {
  if (is.na(SCRNA_RDS) || !file.exists(SCRNA_RDS)) {
    message("SCRNA_RDS chưa đặt/không tồn tại -> BỎ QUA scRNA (Fig8C/8D, tuỳ chọn)."); return(invisible(NULL)) }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    message("Chưa cài 'Seurat' -> BỎ QUA scRNA."); return(invisible(NULL)) }
  obj <- readRDS(SCRNA_RDS)
  if (!CELLTYPE_COL %in% colnames(obj@meta.data)) {
    message("Không thấy cột '", CELLTYPE_COL, "' trong meta.data (có: ",
            paste(head(colnames(obj@meta.data), 10), collapse = ", "), ") -> sửa CELLTYPE_COL.")
    return(NULL) }
  genes_ok <- intersect(sel_genes, rownames(obj))
  say(">> scRNA: ", length(genes_ok), "/", length(sel_genes), " gene signature có trong dataset")
  if (length(genes_ok) < 2) { message("Quá ít gene -> BỎ QUA scRNA."); return(NULL) }

  # AddModuleScore = điểm biểu hiện TẬP GENE, KHÔNG phải risk score sống còn
  obj <- Seurat::AddModuleScore(obj, features = list(genes_ok), name = "CuproSig")
  f8c <- Seurat::FeaturePlot(obj, features = "CuproSig1", reduction = "umap") +
    ggtitle("Figure 8C — Signature module score (UMAP)")
  f8d <- Seurat::VlnPlot(obj, features = "CuproSig1", group.by = CELLTYPE_COL, pt.size = 0) +
    ggtitle("Figure 8D — Signature module score by cell type") + theme(legend.position = "none")
  save_gg(f8c, "Fig8C_scRNA_module_umap", 7, 6)
  save_gg(f8d, "Fig8D_scRNA_module_violin", 9, 5)
  say(">> Figure 8C/8D saved.")
  obj@meta.data[, c(CELLTYPE_COL, "CuproSig1")]
}
res_C2 <- tryCatch(run_scrna(), error = function(e) { message("scRNA lỗi: ", conditionMessage(e)); NULL })

# =============================================================================
# 7. Lưu objects + bảng số + log + sessionInfo
# =============================================================================
saveRDS(list(part_A = res_A, tide = res_B, cancersea = res_C1,
             tisch2 = res_C1b, scrna_scores = res_C2),
        "outputs/Rdata/STEP5_objects.rds")

# FIX#9: bảng mọi con số dùng để viết Results — không phải đọc lại console.
if (length(.STATS) > 0)
  write.csv(do.call(rbind, .STATS), "outputs/Tables/Table_S9d_STEP5_key_statistics.csv",
            row.names = FALSE)

say("")
say("==========================================")
say(" STEP 5 completed (phần thiếu dữ liệu/package đã SKIP có message).")
say(paste(" Phần A:", ifelse(is.null(res_A), "SKIP", "OK"),
          "| TIDE:", ifelse(is.null(res_B), "chờ output web", "OK"),
          "| CancerSEA:", ifelse(is.null(res_C1), "chờ bảng tay", "OK"),
          "| TISCH2:", ifelse(is.null(res_C1b), "chờ bảng tay (dự phòng)", "OK")))
say(" Figures -> outputs/Figures/ (png) + outputs/Figures/vector/ (pdf)")
say(" Số liệu  -> outputs/Tables/Table_S9d_STEP5_key_statistics.csv")
say("==========================================")

writeLines(.LOG, "outputs/STEP5_run_log.txt")
writeLines(capture.output(sessionInfo()), "outputs/sessionInfo_STEP5.txt")
