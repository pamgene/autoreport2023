# Tests for the "Visual assessment of overall signal" QC heatmap (R/01_BasicProcessing.R
# "QC HEATMAP" section): render_qc_heatmap()/render_qc_heatmap_to_file() (drawing +
# graceful degradation), select_qc_heatmap_source()/render_qc_heatmap_for_assay()
# (TR-preferred/BR-fallback data selection, "latest normalization" convention, and
# per-assay caption/filename generation).

qc_heatmap_base_df <- function() {
  tibble(
    ID = rep(c("p1", "p2", "p3"), times = 4),
    Supergroup = rep(c("Sg1", "Sg1", "Sg2", "Sg2"), each = 3),
    `Test Condition` = rep(c("C", "T", "C", "T"), each = 3),
    Barcode = rep(c("B1", "B2", "B3", "B4"), each = 3),
    Row = rep(1, 12),
    logTransformed = rnorm(12, mean = 10, sd = 2)
  )
}

test_that("rainbow_scale maps the data range from blue to red", {
  col_fun <- rainbow_scale(c(0, 5, 10))
  expect_equal(col_fun(0), "#0000FFFF")
  expect_equal(col_fun(10), "#FF0000FF")
})

test_that("contiguous_runs finds runs of identical consecutive values", {
  runs <- contiguous_runs(c("A", "A", "B", "B", "B", "A"))
  expect_equal(runs$values, c("A", "B", "A"))
  expect_equal(runs$starts, c(1, 3, 6))
  expect_equal(runs$ends, c(2, 5, 6))
})

test_that("format_sample_component drops trailing .0 from whole-number Row values but keeps decimals", {
  expect_equal(format_sample_component(c(3, 4)), c("3", "4"))
  expect_equal(format_sample_component(3.0), "3")
  expect_equal(format_sample_component(3.5), "3.5")
  expect_equal(format_sample_component("A"), "A")
})

test_that("render_qc_heatmap draws without error for a well-formed 2-level QC data frame", {
  df <- qc_heatmap_base_df()
  cls <- classify_qc_columns(df, "logTransformed")
  out <- tempfile(fileext = ".png")
  grDevices::png(out, width = 800, height = 600)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(render_qc_heatmap(df, "logTransformed", cls, "test"))
})

test_that("render_qc_heatmap errors on a non-numeric value column", {
  df <- qc_heatmap_base_df() %>% mutate(logTransformed = as.character(logTransformed))
  cls <- classify_qc_columns(qc_heatmap_base_df(), "logTransformed")
  expect_error(render_qc_heatmap(df, "logTransformed", cls, "test"), "not numeric")
})

test_that("render_qc_heatmap errors when there are no condition columns", {
  df <- qc_heatmap_base_df() %>% select(-Supergroup, -`Test Condition`)
  cls <- list(peptide_col = "ID", sample_cols = c("Barcode", "Row"), condition_cols = character(0))
  expect_error(render_qc_heatmap(df, "logTransformed", cls, "test"), "no condition columns")
})

test_that("render_qc_heatmap errors on NA in a condition or sample column", {
  df_na_cond <- qc_heatmap_base_df()
  df_na_cond$Supergroup[1:3] <- NA
  cls <- classify_qc_columns(qc_heatmap_base_df(), "logTransformed")
  expect_error(render_qc_heatmap(df_na_cond, "logTransformed", cls, "test"), "missing \\(NA\\) values in condition")

  df_na_sample <- qc_heatmap_base_df()
  df_na_sample$Barcode[1:3] <- NA
  expect_error(render_qc_heatmap(df_na_sample, "logTransformed", cls, "test"), "missing \\(NA\\) values in sample")
})

test_that("render_qc_heatmap errors when a sample maps to more than one condition combination", {
  df <- qc_heatmap_base_df()
  df$`Test Condition`[df$Barcode == "B1" & df$ID == "p1"] <- "T"
  cls <- classify_qc_columns(qc_heatmap_base_df(), "logTransformed")
  expect_error(render_qc_heatmap(df, "logTransformed", cls, "test"), "more than one condition combination")
})

test_that("render_qc_heatmap errors when fewer than 2 samples or 2 peptides remain", {
  df <- qc_heatmap_base_df()
  cls <- classify_qc_columns(df, "logTransformed")
  expect_error(render_qc_heatmap(df %>% filter(Barcode == "B1"), "logTransformed", cls, "test"), "fewer than 2 distinct samples")
  expect_error(render_qc_heatmap(df %>% filter(ID == "p1"), "logTransformed", cls, "test"), "fewer than 2 distinct peptides")
})

test_that("render_qc_heatmap handles a single condition column (no outer grouping) without error", {
  df <- qc_heatmap_base_df() %>% select(-Supergroup)
  cls <- list(peptide_col = "ID", sample_cols = c("Barcode", "Row"), condition_cols = "Test Condition")
  out <- tempfile(fileext = ".png")
  grDevices::png(out, width = 800, height = 600)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(render_qc_heatmap(df, "logTransformed", cls, "test"))
})

test_that("render_qc_heatmap handles a single-sample-column (BR-style) data frame without error", {
  df <- qc_heatmap_base_df() %>% mutate(Biol_Rep = Barcode) %>% select(-Barcode, -Row)
  cls <- list(peptide_col = "ID", sample_cols = "Biol_Rep", condition_cols = c("Supergroup", "Test Condition"))
  out <- tempfile(fileext = ".png")
  grDevices::png(out, width = 800, height = 600)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(render_qc_heatmap(df, "logTransformed", cls, "test"))
})

test_that("render_qc_heatmap_to_file returns FALSE and leaves no file behind on failure", {
  df <- qc_heatmap_base_df() %>% select(-Supergroup, -`Test Condition`)
  cls <- list(peptide_col = "ID", sample_cols = c("Barcode", "Row"), condition_cols = character(0))
  out <- tempfile(fileext = ".png")
  ok <- render_qc_heatmap_to_file(df, "logTransformed", cls, "test", out)
  expect_false(ok)
  expect_false(file.exists(out))
})

test_that("render_qc_heatmap_to_file returns TRUE and writes a file on success", {
  df <- qc_heatmap_base_df()
  cls <- classify_qc_columns(df, "logTransformed")
  out <- tempfile(fileext = ".png")
  ok <- render_qc_heatmap_to_file(df, "logTransformed", cls, "test", out)
  expect_true(ok)
  expect_true(file.exists(out))
})

# --- select_qc_heatmap_source() / render_qc_heatmap_for_assay() ---

test_that("select_qc_heatmap_source prefers TR over BR for the same assay type", {
  qc_files <- tibble(
    Assay_Type = c("PTK", "PTK"),
    qc_file = c(test_input_path("QC_PTK_01_TR.csv"), test_input_path("QC_PTK_01_BR.csv")),
    Replicate_Type = c("TR", "BR")
  )
  src <- select_qc_heatmap_source(qc_files, "tercen", "PTK")
  expect_equal(src$replicate_type, "TR")
})

test_that("select_qc_heatmap_source picks the latest present normalization", {
  qc_files <- tibble(
    Assay_Type = "STK",
    qc_file = test_input_path("QC_STK_02_BR.csv"), # VSN + ComBat fixture
    Replicate_Type = "BR"
  )
  src <- select_qc_heatmap_source(qc_files, "tercen", "STK")
  expect_equal(src$label, "VSN + ComBat")
})

test_that("select_qc_heatmap_source returns NULL for bionav, missing assay type, or no QC files", {
  qc_files <- tibble(
    Assay_Type = "PTK", qc_file = test_input_path("QC_PTK_01_TR.csv"), Replicate_Type = "TR"
  )
  expect_null(select_qc_heatmap_source(qc_files, "bionav", "PTK"))
  expect_null(select_qc_heatmap_source(qc_files, "tercen", "STK"))
  expect_null(select_qc_heatmap_source(data.frame(), "tercen", "PTK"))
})

test_that("select_qc_heatmap_source skips an unreadable/malformed file but still uses a good one for the same assay", {
  tmp <- tempfile("qcmix_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  writeLines('"Foo","Bar","ID","logTransformed"\n"x","y","p1","1.0"', file.path(tmp, "QC_PTK_02_TR_malformed.csv"))
  file.copy(test_input_path("QC_PTK_01_TR.csv"), file.path(tmp, "QC_PTK_01_TR.csv"))

  qc_files <- tibble(
    Assay_Type = c("PTK", "PTK"),
    qc_file = c(file.path(tmp, "QC_PTK_02_TR_malformed.csv"), file.path(tmp, "QC_PTK_01_TR.csv")),
    Replicate_Type = c("TR", "TR")
  )
  src <- select_qc_heatmap_source(qc_files, "tercen", "PTK")
  expect_false(is.null(src))
})

test_that("select_qc_heatmap_source falls back to BR when every TR file for the assay is unreadable", {
  tmp <- tempfile("qcfallback_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  writeLines('"Foo","Bar","ID","logTransformed"\n"x","y","p1","1.0"', file.path(tmp, "QC_PTK_01_TR_malformed.csv"))
  file.copy(test_input_path("QC_PTK_01_BR.csv"), file.path(tmp, "QC_PTK_01_BR.csv"))

  qc_files <- tibble(
    Assay_Type = c("PTK", "PTK"),
    qc_file = c(file.path(tmp, "QC_PTK_01_TR_malformed.csv"), file.path(tmp, "QC_PTK_01_BR.csv")),
    Replicate_Type = c("TR", "BR")
  )
  src <- select_qc_heatmap_source(qc_files, "tercen", "PTK")
  expect_false(is.null(src))
  expect_equal(src$replicate_type, "BR")
})

test_that("render_qc_heatmap_for_assay writes QC_Heatmap_<assay>_<TR|BR>.png with a matching caption", {
  qc_files <- tibble(
    Assay_Type = "PTK", qc_file = test_input_path("QC_PTK_01_TR.csv"), Replicate_Type = "TR"
  )
  out_dir <- tempfile("figs_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- render_qc_heatmap_for_assay(qc_files, "tercen", "PTK", out_dir = out_dir)
  expect_false(is.null(result))
  expect_equal(basename(result$path), "QC_Heatmap_PTK_TR.png")
  expect_true(file.exists(result$path))
  expect_match(result$caption, "PamChip array")
  expect_match(result$caption, "Log2-transformed")
})

test_that("render_qc_heatmap_for_assay uses 'biological replicate' in the caption for BR-sourced data", {
  qc_files <- tibble(
    Assay_Type = "PTK", qc_file = test_input_path("QC_PTK_01_BR.csv"), Replicate_Type = "BR"
  )
  out_dir <- tempfile("figs_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)

  result <- render_qc_heatmap_for_assay(qc_files, "tercen", "PTK", out_dir = out_dir)
  expect_false(is.null(result))
  expect_equal(basename(result$path), "QC_Heatmap_PTK_BR.png")
  expect_match(result$caption, "biological replicate")
})

test_that("render_qc_heatmap_for_assay returns NULL (not an error) when the assay type has no usable QC data", {
  qc_files <- tibble(
    Assay_Type = "PTK", qc_file = test_input_path("QC_PTK_01_TR.csv"), Replicate_Type = "TR"
  )
  out_dir <- tempfile("figs_")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  expect_null(render_qc_heatmap_for_assay(qc_files, "tercen", "STK", out_dir = out_dir))
})
