# =============================================================================
# STEP 6 — ĐÓNG GÓI NỘP BÀI (assemble submission package)
# Cuproptosis-related gene prognostic signature in TCGA-LUAD
#
# STEP 6 KHÔNG có phân tích mới. Script này chỉ:
#   Part A — kiểm kê toàn bộ figure/table của STEP 1-5  -> MANIFEST
#   Part B — ghép panel thành Figure 1-8 (bản nháp)     -> composite + TODO
#   Part C — gộp Supplementary Tables thành 1 file xlsx -> đánh số lại S1..Sn
#   Part D — đóng gói reproducibility pack              -> .zip + README
#   Part E — checklist nộp bài                          -> báo cáo CSV
#
# CÁCH CHẠY: đặt file này trong thư mục Step_6/ rồi
#   setwd(".../Step_6"); source("STEP6_pipeline.R", echo = TRUE)
#
# BA VIỆC PHẢI XONG TRƯỚC KHI CHẠY (xem STEP6_GUIDANCE_and_ROADMAP.md mục 1):
#   1. Đã chạy lại STEP5_pipeline_FINAL.R (bản cũ in sai odds ratio ở Fig7C)
#   2. Đã chốt Figure 8 hoặc quyết định bỏ
#   3. Đã chuyển code retire (Step_3.R, STEP5_pipeline.R) vào _archive/
# =============================================================================

stopifnot(getRversion() >= "4.2.0")
options(stringsAsFactors = FALSE)

# ---- CONFIG: chỉ sửa khối này -----------------------------------------------
PROJ_ROOT <- "/Users/aithanh/Cuproptosis_LUAD"

# Tên thư mục từng step, đúng như trên ổ đĩa. Sửa nếu khác.
STEP_DIRS <- c(
  STEP1 = file.path(PROJ_ROOT, "Step1/outputs_step1_new"),
  STEP2 = file.path(PROJ_ROOT, "Step_2/outputs"),
  STEP3 = file.path(PROJ_ROOT, "Step_3/outputs"),
  STEP4 = file.path(PROJ_ROOT, "Step_4/outputs"),
  STEP5 = file.path(PROJ_ROOT, "Step_5/outputs")
)

# Script gốc từng step -> tên mới đánh số theo THỨ TỰ CHẠY trong repro pack.
# Tên trái phải trỏ tới file THẬT. Bỏ dòng nào chưa có file.
SCRIPT_MAP <- c(
  "Step1/fetch_luad_metadata.R"      = "01_STEP1_fetch_and_normalize.R",
  "Step_2/STEP2_pipeline_REVIEWED.R" = "02_STEP2_signature_construction.R",
  "Step_3/STEP3_pipeline_FINAL.R"    = "03_STEP3_external_validation.R",
  "Step_4/STEP4_pipeline_REVIEWED.R" = "04_STEP4_downstream_characterization.R",
  "Step_5/STEP5_pipeline_FINAL.R"    = "05_STEP5_genomics_immunotherapy.R",
  "Step_6/STEP6_pipeline.R"          = "06_STEP6_assemble_submission.R"
)

# Dữ liệu dùng trong bài — đi vào data_accessions.csv của repro pack.
DATA_ACCESSIONS <- data.frame(
  dataset   = c("TCGA-LUAD", "TCGA-LUAD", "GSE31210", "GSE72094", "GDSC2", "MSigDB"),
  role      = c("training (expression)", "training (mutation)",
                "external validation", "external validation",
                "drug sensitivity reference", "gene sets"),
  platform  = c("Illumina RNA-seq (STAR-Counts)", "WXS Masked Somatic Mutation",
                "Affymetrix microarray", "Affymetrix microarray", "-", "-"),
  n_samples = c("504 tumours (after QC)", "492 patients", "226", "398", "-", "-"),
  source    = c("GDC Data Portal", "GDC Data Portal", "NCBI GEO", "NCBI GEO",
                "oncoPredict / GDSC2", "msigdbr"),
  date_accessed = c("<ĐIỀN>", "<ĐIỀN>", "<ĐIỀN>", "<ĐIỀN>", "<ĐIỀN>", "<ĐIỀN>"),
  stringsAsFactors = FALSE
)
# TRIPOD bắt buộc ghi ngày truy cập — TCGA và GEO cập nhật liên tục.

SEED_USED      <- 123
TARGET_JOURNAL <- "<ĐIỀN tên tạp chí đích>"
MAX_WIDTH_MM   <- 180    # khổ 2 cột rộng nhất mà tạp chí chấp nhận
MIN_DPI        <- 300

dir.create("outputs/Figures_final", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/Supplementary", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/reproducibility_pack", recursive = TRUE, showWarnings = FALSE)

.LOG <- character(0)
say <- function(...) { txt <- paste0(...); .LOG <<- c(.LOG, txt); cat(txt, "\n", sep = "") }

# =============================================================================
# PART A — KIỂM KÊ
# =============================================================================

# Đọc kích thước + dpi của PNG bằng base R (không cần package).
# Width/height nằm trong chunk IHDR; dpi nằm trong chunk pHYs (pixels per metre).
png_info <- function(path) {
  con <- file(path, "rb"); on.exit(close(con), add = TRUE)
  sig <- as.integer(readBin(con, "raw", 8))
  if (length(sig) < 8 || !identical(sig, c(137L, 80L, 78L, 71L, 13L, 10L, 26L, 10L)))
    return(list(width_px = NA_integer_, height_px = NA_integer_, dpi = NA_real_))
  be32 <- function(r) sum(as.integer(r) * c(256^3, 256^2, 256, 1))
  w <- h <- NA_integer_; dpi <- NA_real_
  repeat {
    len_raw <- readBin(con, "raw", 4)
    if (length(len_raw) < 4) break
    len  <- be32(len_raw)
    type <- rawToChar(readBin(con, "raw", 4))
    dat  <- if (len > 0) readBin(con, "raw", len) else raw(0)
    invisible(readBin(con, "raw", 4))                      # CRC
    if (type == "IHDR" && length(dat) >= 8) {
      w <- be32(dat[1:4]); h <- be32(dat[5:8])
    } else if (type == "pHYs" && length(dat) >= 9) {
      if (as.integer(dat[9]) == 1L) dpi <- be32(dat[1:4]) * 0.0254   # px/m -> px/inch
    } else if (type == "IEND") break
  }
  list(width_px = w, height_px = h, dpi = dpi)
}

csv_shape <- function(path) {
  tryCatch({
    nf <- utils::count.fields(path, sep = ",", quote = "\"")
    list(n_row = max(0L, length(nf) - 1L),            # trừ dòng header
         n_col = as.integer(stats::median(nf, na.rm = TRUE)))
  }, error = function(e) list(n_row = NA_integer_, n_col = NA_integer_))
}

inventory <- function() {
  rows <- list()
  for (st in names(STEP_DIRS)) {
    base <- STEP_DIRS[[st]]
    if (!dir.exists(base)) { say("THIẾU thư mục ", st, ": ", base); next }

    figs <- list.files(file.path(base, "Figures"), pattern = "\\.png$",
                       full.names = TRUE, ignore.case = TRUE)
    for (f in figs) {
      info <- png_info(f)
      dpi  <- if (is.na(info$dpi)) MIN_DPI else info$dpi   # ggsave ghi pHYs; nếu thiếu thì giả định
      vec  <- file.path(base, "Figures", "vector",
                        sub("\\.png$", ".pdf", basename(f), ignore.case = TRUE))
      rows[[length(rows) + 1]] <- data.frame(
        step = st, kind = "figure", name = basename(f), path = f,
        width_px = info$width_px, height_px = info$height_px,
        dpi_declared = info$dpi,
        print_width_mm  = round(info$width_px  / dpi * 25.4, 1),
        print_height_mm = round(info$height_px / dpi * 25.4, 1),
        has_vector = file.exists(vec),
        n_row = NA_integer_, n_col = NA_integer_,
        size_kb = round(file.size(f) / 1024, 1))
    }

    tbls <- list.files(file.path(base, "Tables"), pattern = "\\.csv$",
                       full.names = TRUE, ignore.case = TRUE)
    for (t in tbls) {
      sh <- csv_shape(t)
      rows[[length(rows) + 1]] <- data.frame(
        step = st, kind = "table", name = basename(t), path = t,
        width_px = NA_integer_, height_px = NA_integer_, dpi_declared = NA_real_,
        print_width_mm = NA_real_, print_height_mm = NA_real_, has_vector = NA,
        n_row = sh$n_row, n_col = sh$n_col,
        size_kb = round(file.size(t) / 1024, 1))
    }
  }
  if (!length(rows)) stop("Không tìm thấy output nào — kiểm tra PROJ_ROOT và STEP_DIRS.")
  man <- do.call(rbind, rows)
  man <- man[order(man$step, man$kind, man$name), ]
  write.csv(man, "outputs/MANIFEST_figures_tables.csv", row.names = FALSE)

  say("== PART A — KIỂM KÊ ==")
  say("Figure: ", sum(man$kind == "figure"), " | Table: ", sum(man$kind == "table"))

  # Ba cảnh báo phải xử lý trước khi chốt figure
  no_vec <- man[man$kind == "figure" & !is.na(man$has_vector) & !man$has_vector, ]
  if (nrow(no_vec)) {
    say("\n[1] Figure chỉ có PNG, THIẾU PDF vector (", nrow(no_vec), "):")
    for (n in no_vec$name) say("    - ", n)
    say("    -> xuất lại bằng save_gg() của chính step đó.")
  }
  too_wide <- man[man$kind == "figure" & !is.na(man$print_width_mm) &
                    man$print_width_mm > MAX_WIDTH_MM, ]
  if (nrow(too_wide)) {
    say("\n[2] Figure rộng hơn ", MAX_WIDTH_MM, " mm khi in ở ", MIN_DPI,
        " dpi -> tạp chí sẽ co lại, chữ tụt dưới 7 pt (", nrow(too_wide), "):")
    for (i in seq_len(nrow(too_wide)))
      say(sprintf("    - %-42s %.0f mm  (scale %.2fx -> font 12pt thành %.1fpt)",
                  too_wide$name[i], too_wide$print_width_mm[i],
                  MAX_WIDTH_MM / too_wide$print_width_mm[i],
                  12 * MAX_WIDTH_MM / too_wide$print_width_mm[i]))
    say("    -> XUẤT LẠI ở đúng khổ in, KHÔNG resize. Xem GUIDANCE mục 3.1.")
  }
  empty_tbl <- man[man$kind == "table" & !is.na(man$n_row) & man$n_row <= 0, ]
  if (nrow(empty_tbl)) {
    say("\n[3] Table rỗng hoặc chỉ có header (", nrow(empty_tbl), "):")
    for (n in empty_tbl$name) say("    - ", n)
  }
  if (!nrow(no_vec) && !nrow(too_wide) && !nrow(empty_tbl))
    say("\nMANIFEST sạch — chuyển sang Part B.")

  say("\n-> outputs/MANIFEST_figures_tables.csv")
  man
}
manifest <- inventory()

# =============================================================================
# PART B — GHÉP FIGURE 1-8
#
# QUAN TRỌNG: composite mà Part B tạo ra là BẢN NHÁP để nhóm xem và để gửi kèm
# bản thảo vòng đầu. Nó ghép từ PNG đã xuất, nên vẫn là raster đã resize.
# Figure NỘP CHO TẠP CHÍ phải xuất lại từ script gốc ở đúng khổ in cuối cùng
# (ggsave(width = 6.69, ...) + base_size phù hợp). Part B in ra TODO cho việc đó.
#
# Panel có tỉ lệ khác nhau sẽ để lại khoảng trắng trong bản nháp (vì ghép theo
# cùng chiều rộng). Chuyện bình thường, không cần sửa ở bước này — khi xuất lại
# ở khổ in cuối thì đặt height cho từng panel sao cho mỗi hàng cao bằng nhau.
# =============================================================================

# Thanh điền: mỗi dòng một panel. 'file' lấy đúng tên trong MANIFEST.
# ncol = số panel trên một hàng khi ghép. width_mm = khổ in ĐÍCH của figure.
FIGURE_SPEC <- data.frame(
  figure  = c(7, 7, 7, 7, 7),
  panel   = c("A", "B", "C", "D", "E"),
  file    = c("Fig7A_TMB_by_risk.png",
              "Fig7B_coOncoplot_high_vs_low.png",
              "Fig7C_mafCompare_forest.png",
              "Fig7D_TIDE_immunotherapy.png",
              "Fig7E_KM_TMB_x_risk.png"),
  step    = c("STEP5", "STEP5", "STEP5", "STEP5", "STEP5"),
  ncol    = c(2, 2, 2, 2, 2),
  width_mm = c(180, 180, 180, 180, 180),
  stringsAsFactors = FALSE
)
# TODO Thanh: thêm các dòng cho Figure 1-6 (và 8 nếu giữ) theo bảng ở
# GUIDANCE mục 3.3, rồi chạy lại Part B.

assemble_figures <- function(spec, man) {
  say("\n== PART B — GHÉP FIGURE ==")
  if (!nrow(spec)) { say("FIGURE_SPEC còn rỗng — Thanh điền trước."); return(invisible(NULL)) }

  # Panel nào khai trong spec mà không có trong MANIFEST?
  miss <- spec[!spec$file %in% man$name, ]
  if (nrow(miss)) {
    say("Panel khai trong FIGURE_SPEC nhưng không tìm thấy file:")
    for (i in seq_len(nrow(miss))) say("    - Figure ", miss$figure[i], miss$panel[i],
                                      ": ", miss$file[i])
  }

  # TODO xuất lại ở khổ in: panel nào đang rộng hơn khổ đích
  todo <- merge(spec, man[man$kind == "figure", c("name", "print_width_mm", "has_vector")],
                by.x = "file", by.y = "name", all.x = TRUE)
  todo$needs_reexport <- is.na(todo$print_width_mm) |
    todo$print_width_mm > todo$width_mm * 1.05 | !isTRUE(todo$has_vector)
  write.csv(todo[order(todo$figure, todo$panel), ],
            "outputs/Figures_final/FIGURE_REEXPORT_TODO.csv", row.names = FALSE)
  say("Panel cần xuất lại ở khổ in: ", sum(todo$needs_reexport, na.rm = TRUE),
      "/", nrow(todo), " -> outputs/Figures_final/FIGURE_REEXPORT_TODO.csv")

  if (!requireNamespace("magick", quietly = TRUE)) {
    say("Chưa cài 'magick' -> BỎ QUA việc ghép composite.")
    say("  install.packages('magick')  (macOS: brew install imagemagick trước)")
    return(invisible(NULL))
  }

  for (fig in sort(unique(spec$figure))) {
    sp <- spec[spec$figure == fig, ]
    sp <- sp[order(sp$panel), ]
    paths <- vapply(seq_len(nrow(sp)), function(i) {
      hit <- man$path[man$kind == "figure" & man$name == sp$file[i]]
      if (length(hit)) hit[1] else NA_character_
    }, character(1))
    if (anyNA(paths)) { say("Figure ", fig, ": thiếu panel -> bỏ qua."); next }

    target_px <- round(sp$width_mm[1] / 25.4 * MIN_DPI)
    ncol      <- max(1L, as.integer(sp$ncol[1]))
    panel_px  <- round(target_px / ncol)

    imgs <- lapply(seq_along(paths), function(i) {
      im <- magick::image_read(paths[i])
      im <- magick::image_resize(im, paste0(panel_px, "x"))
      magick::image_annotate(im, sp$panel[i], size = round(panel_px * 0.055),
                             weight = 700, gravity = "northwest",
                             location = "+8+6", color = "black", font = "Helvetica")
    })
    rows <- split(imgs, ceiling(seq_along(imgs) / ncol))
    row_imgs <- lapply(rows, function(r) magick::image_append(do.call(c, r)))
    comp <- magick::image_append(do.call(c, row_imgs), stack = TRUE)
    comp <- magick::image_background(comp, "white", flatten = TRUE)

    out_png <- sprintf("outputs/Figures_final/Figure_%d_DRAFT.png", fig)
    magick::image_write(comp, out_png, format = "png", density = MIN_DPI)
    say(sprintf("Figure %d: ghép %d panel, %d cột, khổ đích %.0f mm -> %s (BẢN NHÁP)",
                fig, nrow(sp), ncol, sp$width_mm[1], basename(out_png)))
  }
  say("Nhắc lại: bản nháp chỉ để nhóm xem. Figure nộp phải xuất lại từ script gốc.")
  invisible(NULL)
}
assemble_figures(FIGURE_SPEC, manifest)

# =============================================================================
# PART C — SUPPLEMENTARY TABLES: đánh số lại S1..Sn, gộp vào 1 file xlsx
# =============================================================================

# Thanh + Tâm điền: THỨ TỰ dòng ở đây = thứ tự S1, S2, S3... trong bản thảo.
# Chỉ đưa bảng ĐƯỢC TRÍCH DẪN trong bài. Bảng không trích dẫn thì bỏ.
SUPP_SPEC <- data.frame(
  file = c("Table_S9b_TMB_by_patient.csv",
           "Table_S9_mafCompare_high_vs_low.csv",
           "Table_S9e_LUAD_driver_panel.csv",
           "Table_S9c_TIDE_score_stat.csv",
           "Table_S9f_Cox_risk_plus_TMB.csv",
           "Table_S9d_STEP5_key_statistics.csv"),
  caption = c(
    "Tumour mutation burden per patient with CRG risk group assignment.",
    "Differentially mutated genes between high- and low-risk groups (Fisher's exact test, genome-wide exploratory scan).",
    "Mutation frequency of 20 established LUAD driver genes by risk group (BH-FDR over the 20-gene panel).",
    "TIDE scores and predicted responder proportion by risk group.",
    "Multivariable Cox model of overall survival including CRG risk group and TMB group.",
    "Key statistics reported in the main text for the genomic and immunotherapy analyses."),
  cited_in = c("Results 3.5", "Results 3.5", "Results 3.5",
               "Results 3.6", "Results 3.5", "Results 3.5-3.6"),
  stringsAsFactors = FALSE
)
# TODO: thêm bảng của STEP1-STEP4 (Table1 clinical, Table_S4 log-rank validation,
# Table_S6 multivariate Cox...) vào ĐÚNG thứ tự xuất hiện trong bản thảo.

build_supplementary <- function(spec, man) {
  say("\n== PART C — SUPPLEMENTARY ==")
  if (!nrow(spec)) { say("SUPP_SPEC rỗng — điền trước."); return(invisible(NULL)) }

  spec$new_id  <- paste0("S", seq_len(nrow(spec)))
  spec$path    <- vapply(spec$file, function(f) {
    hit <- man$path[man$kind == "table" & man$name == f]
    if (length(hit)) hit[1] else NA_character_
  }, character(1))
  spec$old_id  <- sub("^Table_(S[0-9a-z]+)_.*$", "\\1", spec$file)
  spec$found   <- !is.na(spec$path)

  write.csv(spec[, c("old_id", "new_id", "file", "caption", "cited_in", "found")],
            "outputs/Supplementary/Supplementary_Table_renumbering_map.csv",
            row.names = FALSE)
  say("Bản đồ đánh số lại -> outputs/Supplementary/Supplementary_Table_renumbering_map.csv")
  say("  Tâm dùng file này để sửa trích dẫn trong bản thảo, ĐỪNG sửa bằng mắt.")
  for (i in seq_len(nrow(spec)))
    say(sprintf("  %-6s -> %-4s %s%s", spec$old_id[i], spec$new_id[i], spec$file[i],
                ifelse(spec$found[i], "", "   [KHÔNG TÌM THẤY FILE]")))

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    say("Chưa cài 'openxlsx' -> không gộp được xlsx.  install.packages('openxlsx')")
    return(invisible(spec))
  }
  ok <- spec[spec$found, ]
  if (!nrow(ok)) { say("Không có bảng nào đọc được."); return(invisible(spec)) }

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Contents")
  contents <- data.frame(Table = ok$new_id, Title = ok$caption, `Cited in` = ok$cited_in,
                         check.names = FALSE)
  openxlsx::writeData(wb, "Contents", "Supplementary Tables — Contents", startRow = 1)
  openxlsx::writeData(wb, "Contents", contents, startRow = 3)
  openxlsx::addStyle(wb, "Contents", openxlsx::createStyle(textDecoration = "bold"),
                     rows = c(1, 3), cols = 1:3, gridExpand = TRUE)
  openxlsx::setColWidths(wb, "Contents", cols = 1:3, widths = c(10, 90, 20))

  for (i in seq_len(nrow(ok))) {
    sh <- ok$new_id[i]
    dat <- tryCatch(read.csv(ok$path[i], check.names = FALSE), error = function(e) NULL)
    if (is.null(dat)) { say("  Đọc lỗi: ", ok$file[i]); next }
    openxlsx::addWorksheet(wb, sh)
    openxlsx::writeData(wb, sh, paste0("Table ", sh, ". ", ok$caption[i]), startRow = 1)
    openxlsx::addStyle(wb, sh, openxlsx::createStyle(textDecoration = "bold"),
                       rows = 1, cols = 1)
    openxlsx::writeData(wb, sh, dat, startRow = 3)
    openxlsx::freezePane(wb, sh, firstActiveRow = 4)
  }
  out <- "outputs/Supplementary/Supplementary_Tables.xlsx"
  openxlsx::saveWorkbook(wb, out, overwrite = TRUE)
  say("Đã gộp ", nrow(ok), " bảng + sheet Contents -> ", out)
  invisible(spec)
}
supp <- build_supplementary(SUPP_SPEC, manifest)

# =============================================================================
# PART D — REPRODUCIBILITY PACK
# =============================================================================
build_repro_pack <- function() {
  say("\n== PART D — REPRODUCIBILITY PACK ==")
  root <- "outputs/reproducibility_pack"
  unlink(root, recursive = TRUE)
  for (d in c("scripts", "session_info", "derived")) dir.create(file.path(root, d), recursive = TRUE)

  # 1. Script, đánh số theo thứ tự chạy
  copied <- character(0); missing <- character(0)
  for (src in names(SCRIPT_MAP)) {
    from <- file.path(PROJ_ROOT, src)
    if (file.exists(from)) {
      file.copy(from, file.path(root, "scripts", SCRIPT_MAP[[src]]), overwrite = TRUE)
      copied <- c(copied, SCRIPT_MAP[[src]])
    } else missing <- c(missing, src)
  }
  say("Script đã đóng gói: ", length(copied))
  if (length(missing)) {
    say("THIẾU script (sửa SCRIPT_MAP hoặc kiểm lại tên file):")
    for (m in missing) say("    - ", m)
  }

  # 2. sessionInfo của từng step
  n_si <- 0
  for (st in names(STEP_DIRS)) {
    f <- file.path(STEP_DIRS[[st]], paste0("sessionInfo_", st, ".txt"))
    if (!file.exists(f)) f <- file.path(STEP_DIRS[[st]], "..",
                                        paste0("sessionInfo_", st, ".txt"))
    if (file.exists(f)) {
      file.copy(f, file.path(root, "session_info", paste0("sessionInfo_", st, ".txt")),
                overwrite = TRUE); n_si <- n_si + 1
    } else say("Thiếu sessionInfo cho ", st, " -> chạy lại step đó để sinh ra.")
  }
  say("sessionInfo đã đóng gói: ", n_si, "/", length(STEP_DIRS))
  writeLines(capture.output(sessionInfo()), file.path(root, "session_info", "sessionInfo_STEP6.txt"))

  # 3. Accession
  write.csv(DATA_ACCESSIONS, file.path(root, "data_accessions.csv"), row.names = FALSE)
  if (any(DATA_ACCESSIONS$date_accessed == "<ĐIỀN>"))
    say("CHƯA ĐIỀN ngày truy cập dữ liệu — TRIPOD bắt buộc có. Sửa DATA_ACCESSIONS.")

  # 4. README
  readme <- c(
    "# Reproducibility pack — Cuproptosis-related gene prognostic signature in LUAD",
    "",
    "## 1. Chạy gì, theo thứ tự nào",
    "Chạy các script trong `scripts/` theo đúng số thứ tự. Mỗi script đọc file .rds",
    "do script trước sinh ra, nên không bỏ bước được.",
    "",
    paste0("  ", sort(copied), collapse = "\n"),
    "",
    "## 2. Môi trường",
    paste0("R >= 4.2.0. Phiên bản package chính xác của từng bước nằm trong ",
           "`session_info/sessionInfo_STEP*.txt`."),
    "Package Bioconductor cần: TCGAbiolinks, edgeR, limma, maftools, GSVA,",
    "clusterProfiler, msigdbr. Package CRAN cần: survival, survminer, glmnet,",
    "timeROC, rms, ggplot2, dplyr, tidyr, tibble, openxlsx, magick.",
    "",
    "## 3. Seed",
    paste0("Mọi bước có yếu tố ngẫu nhiên đều đặt `set.seed(", SEED_USED, ")` trước khi chạy."),
    "Quan trọng nhất là `cv.glmnet` 10-fold ở STEP 2: đổi seed là ra bộ gene khác.",
    "Chia nhóm nguy cơ ở MEDIAN của risk score, không dùng optimal cutpoint.",
    "",
    "## 4. Dữ liệu",
    "Xem `data_accessions.csv`. Pack này KHÔNG chứa dữ liệu bệnh nhân thô:",
    "tải lại từ GDC Data Portal và NCBI GEO bằng accession trong file đó.",
    "`derived/` chỉ chứa bảng dẫn xuất không định danh được.",
    "",
    "## 5. Liên hệ",
    "<ĐIỀN tên tác giả liên hệ và email>"
  )
  writeLines(readme, file.path(root, "README.md"))

  # 5. Bảng dẫn xuất — không định danh, an toàn để công bố
  for (f in c("STEP2_coefficients.csv", "risk_score_per_sample.csv")) {
    hit <- manifest$path[manifest$kind == "table" & manifest$name == f]
    if (length(hit)) file.copy(hit[1], file.path(root, "derived", f), overwrite = TRUE)
    else say("Chưa có ", f, " -> xuất từ STEP2 (hệ số beta và risk score từng mẫu).")
  }

  # 6. Zip
  zipfile <- "outputs/reproducibility_pack.zip"
  old <- setwd("outputs")
  ok <- tryCatch({ utils::zip("reproducibility_pack.zip", "reproducibility_pack",
                              flags = "-r9Xq"); TRUE },
                 error = function(e) { message("zip lỗi: ", conditionMessage(e)); FALSE })
  setwd(old)
  if (ok && file.exists(zipfile))
    say("-> ", zipfile, " (", round(file.size(zipfile) / 1024, 1), " KB)")
  else
    say("Không zip được (thiếu lệnh zip). Nén tay thư mục outputs/reproducibility_pack/")
  say("Bước tiếp: upload zip lên Zenodo lấy DOI, rồi trích DOI đó vào Data Availability Statement.")
  invisible(NULL)
}
build_repro_pack()

# =============================================================================
# PART E — CHECKLIST NỘP BÀI
# =============================================================================
submission_checklist <- function(man, spec) {
  say("\n== PART E — CHECKLIST ==")
  fig_names <- man$name[man$kind == "figure"]
  chk <- function(item, pass, note = "") data.frame(item = item,
    status = ifelse(isTRUE(pass), "OK", "CHƯA"), note = note)

  rows <- rbind(
    chk("Figure 7C đã xuất lại từ bản STEP5 FINAL (odds ratio thật)",
        FALSE, "kiểm bằng mắt: nhãn không được có 'OR 50.00'"),
    chk("Figure 8 đã chốt hoặc đã quyết định bỏ",
        any(grepl("^Fig8", fig_names)), "deadline 07/11/2026, xem GUIDANCE 1.2"),
    chk("Mọi figure có PDF vector",
        all(man$has_vector[man$kind == "figure"] %in% TRUE)),
    chk("Không figure nào rộng quá 180 mm ở 300 dpi",
        all(is.na(man$print_width_mm) | man$print_width_mm <= MAX_WIDTH_MM)),
    chk("FIGURE_SPEC đã khai đủ Figure 1-8",
        length(unique(spec$figure)) >= 7,
        paste("hiện khai: Figure", paste(sort(unique(spec$figure)), collapse = ", "))),
    chk("Supplementary_Tables.xlsx đã tạo",
        file.exists("outputs/Supplementary/Supplementary_Tables.xlsx")),
    chk("Bảng đánh số liên tục S1..Sn",
        file.exists("outputs/Supplementary/Supplementary_Table_renumbering_map.csv")),
    chk("reproducibility_pack.zip đã tạo",
        file.exists("outputs/reproducibility_pack.zip")),
    chk("Ngày truy cập dữ liệu đã điền",
        !any(DATA_ACCESSIONS$date_accessed == "<ĐIỀN>"), "TRIPOD bắt buộc"),
    chk("Tạp chí đích đã chốt", TARGET_JOURNAL != "<ĐIỀN tên tạp chí đích>"),
    chk("TRIPOD checklist đã điền số trang", FALSE, "tải từ equator-network.org"),
    chk("Limitations có nêu TIDE cùng nguồn biểu hiện với risk score", FALSE,
        "điểm yếu thật — tự nêu, đừng để reviewer nêu"),
    chk("Limitations có nêu chưa so sánh với signature đã công bố", FALSE),
    chk("Đã báo cáo số bệnh nhân bị loại do thiếu stage/smoking ở STEP4", FALSE),
    chk("Mọi HR và C-index có 95% CI", FALSE),
    chk("Code retire đã chuyển khỏi pack", !file.exists(file.path(PROJ_ROOT, "Step_3/Step_3.R")))
  )
  write.csv(rows, "outputs/SUBMISSION_CHECKLIST.csv", row.names = FALSE)
  for (i in seq_len(nrow(rows)))
    say(sprintf("  [%-4s] %s%s", rows$status[i], rows$item[i],
                ifelse(nzchar(rows$note[i]), paste0("  — ", rows$note[i]), "")))
  say("\nCòn ", sum(rows$status == "CHƯA"), "/", nrow(rows), " mục chưa xong.")
  say("-> outputs/SUBMISSION_CHECKLIST.csv")
  rows
}
checklist <- submission_checklist(manifest, FIGURE_SPEC)

# =============================================================================
writeLines(.LOG, "outputs/STEP6_run_log.txt")
say("\n-> outputs/STEP6_run_log.txt")
