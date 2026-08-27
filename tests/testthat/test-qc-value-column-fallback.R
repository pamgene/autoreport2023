# Tests for the "value"-column fallback: some QC files (e.g. a biological-replicate file
# that only exists because it was aggregated up from technical-replicate rows) lose the
# semantic logTransformed/identity/CmbCor column name and end up with a generic "value"
# column instead. In that case the normalization is read from a tag in the filename
# (Log/LogCmb/VSN/VSNCmb, case-insensitive) instead of the column name.

test_that("extract_normalization_hint finds a tag as an exact token, case-insensitively", {
  expect_equal(extract_normalization_hint(c("QC", "STK", "03", "BR", "VSNCmb")), "VSN + ComBat")
  expect_equal(extract_normalization_hint(c("QC", "STK", "03", "BR", "vsncmb")), "VSN + ComBat")
  expect_equal(extract_normalization_hint(c("QC", "STK", "03", "BR", "VsnCmb")), "VSN + ComBat")
  expect_equal(extract_normalization_hint(c("QC", "PTK", "01", "BR", "Log")), "Log")
  expect_equal(extract_normalization_hint(c("QC", "PTK", "01", "BR", "LOG")), "Log")
  expect_equal(extract_normalization_hint(c("QC", "PTK", "01", "BR", "LogCmb")), "Log + ComBat")
  expect_equal(extract_normalization_hint(c("QC", "PTK", "01", "BR", "VSN")), "VSN")
})

test_that("extract_normalization_hint returns NA when there's no tag", {
  expect_true(is.na(extract_normalization_hint(c("QC", "STK", "01", "BR"))))
  expect_true(is.na(extract_normalization_hint(c("QC", "STK", "LogCombined")))) # not an exact token match
})

test_that("read_qc_dir parses the normalization tag from the filename, case-insensitively", {
  tmp <- tempfile("qc_dir_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  file.create(file.path(tmp, "QC_STK_03_BR_VSNCmb.csv"))
  file.create(file.path(tmp, "QC_STK_04_BR_vsncmb.csv"))
  file.create(file.path(tmp, "QC_PTK_01_BR.csv")) # no tag -> NA

  qc_files <- read_qc_dir(paste0(tmp, "/"))

  expect_equal(
    qc_files$Normalization_Hint[grepl("03_BR_VSNCmb", qc_files$qc_file)], "VSN + ComBat"
  )
  expect_equal(
    qc_files$Normalization_Hint[grepl("04_BR_vsncmb", qc_files$qc_file)], "VSN + ComBat"
  )
  expect_true(is.na(qc_files$Normalization_Hint[grepl("QC_PTK_01_BR.csv$", qc_files$qc_file)]))
})

test_that("identify_value_columns falls back to a generic 'value' column when given a filename hint", {
  df <- tibble(ID = "p1", Barcode = "b1", Row = 1, value = 4.2)
  norms <- identify_value_columns(df, "VSN + ComBat")
  expect_named(norms, "VSN + ComBat")
  expect_equal(norms[["VSN + ComBat"]], "value")

  # Case-insensitive column match too (e.g. "Value" instead of "value").
  df2 <- tibble(ID = "p1", Barcode = "b1", Row = 1, Value = 4.2)
  norms2 <- identify_value_columns(df2, "Log")
  expect_named(norms2, "Log")
  expect_equal(norms2[["Log"]], "Value")
})

test_that("identify_value_columns ignores the 'value' fallback without a filename hint", {
  df <- tibble(ID = "p1", Barcode = "b1", Row = 1, value = 4.2)
  expect_length(identify_value_columns(df, NA_character_), 0)
})

test_that("identify_value_columns prefers logTransformed/identity/CmbCor over the value fallback", {
  # A file with a real semantic column name should never need (or use) the hint.
  df <- tibble(ID = "p1", Barcode = "b1", Row = 1, logTransformed = 4.2)
  norms <- identify_value_columns(df, "VSN + ComBat") # hint present but irrelevant here
  expect_named(norms, "Log")
  expect_equal(norms[["Log"]], "logTransformed")
})

test_that("end-to-end: a synthetic BR 'value' file, correctly tagged, parses and computes SD", {
  dir <- make_input_dir(c(
    "QC_STK_03_BR_VSNCmb.csv" = file.path("fixtures", "QC_STK_synthetic_value_column.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_equal(qc_files$Normalization_Hint, "VSN + ComBat")
  expect_equal(qc_files$Replicate_Type, "BR")

  parsed <- parse_qc(qc_files, "tercen")

  expect_equal(nrow(parsed$qc_table), 1) # BR -> single row
  expect_equal(nrow(parsed$variability_table), 1)
  expect_equal(parsed$variability_table$Normalization, "VSN + ComBat")

  # Reference calculation, independent of the function under test: peptide SD per
  # condition (only 1 condition here), median across the 2 peptides (with exactly 2
  # peptides the median equals their mean, so this doesn't distinguish the two - see
  # test-qc-variability.R for cases where they actually differ).
  # p1: sd(4.0, 4.2, 4.4); p2: sd(5.0, 5.4, 5.8)
  expected_median_sd <- median(c(sd(c(4.0, 4.2, 4.4)), sd(c(5.0, 5.4, 5.8))))
  expect_equal(parsed$variability_table$Variability_Value, expected_median_sd, tolerance = 1e-6)
})

test_that("end-to-end: an untagged file with only a 'value' column is skipped gracefully (no crash)", {
  dir <- make_input_dir(c(
    "QC_STK_03_BR.csv" = file.path("fixtures", "QC_STK_synthetic_value_column.csv") # no tag
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_true(is.na(qc_files$Normalization_Hint))

  parsed <- parse_qc(qc_files, "tercen")
  expect_equal(nrow(parsed$qc_table), 0) # no recognized normalization -> nothing to report
  expect_equal(nrow(parsed$variability_table), 0)
})

test_that("detect_normalizations_from_qc also resolves the value-column fallback via filename tag", {
  dir <- make_input_dir(c(
    "QC_STK_03_BR_VSNCmb.csv" = file.path("fixtures", "QC_STK_synthetic_value_column.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_setequal(detect_normalizations_from_qc(qc_files, "tercen"), c("vsn", "combat"))
})

test_that("a simple TR-only project needs no tag - real files, untagged, still fully parse", {
  # Real fixtures with proper logTransformed/CmbCor and identity/CmbCor columns (not a
  # generic "value" column), renamed to a plain untagged TR filename - exactly what a
  # study with only technical replicates (no BR at all) would look like.
  dir <- make_input_dir(c(
    "QC_STK_01_TR.csv" = test_input_path("QC_STK_01_BR.csv"),
    "QC_STK_02_TR.csv" = test_input_path("QC_STK_02_BR.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_true(all(is.na(qc_files$Normalization_Hint))) # no tag needed or present
  expect_true(all(qc_files$Replicate_Type == "TR"))

  parsed <- parse_qc(qc_files, "tercen")
  # TR -> only the latest normalization in Table 2 (tag absence is irrelevant - detected
  # entirely from the files' own columns).
  expect_equal(nrow(parsed$qc_table), 1)
  expect_equal(parsed$qc_table$Normalization, "VSN + ComBat")
  expect_true(!is.na(parsed$qc_table$Variability_Flag))
  # But all 4 normalizations still show up in the full detail table (Table 4).
  expect_setequal(parsed$variability_table$Normalization, c("Log", "Log + ComBat", "VSN", "VSN + ComBat"))
})

test_that("a simple BR-only project needs no tag - real files, untagged, still fully parse", {
  dir <- make_input_dir(c(
    "QC_STK_01_BR.csv" = test_input_path("QC_STK_01_BR.csv"),
    "QC_STK_02_BR.csv" = test_input_path("QC_STK_02_BR.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_true(all(is.na(qc_files$Normalization_Hint)))
  expect_true(all(qc_files$Replicate_Type == "BR"))

  parsed <- parse_qc(qc_files, "tercen")
  expect_equal(nrow(parsed$qc_table), 1) # BR -> single collapsed row
  expect_equal(nrow(parsed$variability_table), 4) # all 4 normalizations still detailed here
  expect_setequal(parsed$variability_table$Normalization, c("Log", "Log + ComBat", "VSN", "VSN + ComBat"))
})
