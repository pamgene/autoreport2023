test_that("identify_value_columns finds Log + ComBat from a logTransformed file", {
  df <- tibble(ID = "p1", logTransformed = 1, CmbCor = 2)
  norms <- identify_value_columns(df)
  expect_named(norms, c("Log", "Log + ComBat"))
  expect_equal(norms[["Log"]], "logTransformed")
  expect_equal(norms[["Log + ComBat"]], "CmbCor")
})

test_that("identify_value_columns finds VSN + ComBat from an identity file", {
  df <- tibble(ID = "p1", identity = 1, CmbCor = 2)
  norms <- identify_value_columns(df)
  expect_named(norms, c("VSN", "VSN + ComBat"))
})

test_that("identify_value_columns handles a file with no ComBat correction", {
  df <- tibble(ID = "p1", logTransformed = 1)
  norms <- identify_value_columns(df)
  expect_named(norms, "Log")
})

test_that("identify_value_columns returns empty for a file with neither base column", {
  df <- tibble(ID = "p1", CmbCor = 2)
  norms <- identify_value_columns(df)
  expect_length(norms, 0)
})

test_that("classify_qc_columns uses Barcode+Row when both are present", {
  df <- tibble(ID = "p1", Barcode = "b1", Row = 1, Supergroup = "s1", `Test Condition` = "Test", logTransformed = 1)
  cls <- classify_qc_columns(df, "logTransformed")
  expect_equal(cls$peptide_col, "ID")
  expect_setequal(cls$sample_cols, c("Barcode", "Row"))
  expect_setequal(cls$condition_cols, c("Supergroup", "Test Condition"))
})

test_that("classify_qc_columns falls back to a single 'sample'-matching column", {
  df <- tibble(ID = "p1", Sample_no = "s1", `Test Condition` = "Test", logTransformed = 1)
  cls <- classify_qc_columns(df, "logTransformed")
  expect_equal(cls$sample_cols, "Sample_no")
  expect_equal(cls$condition_cols, "Test Condition")
})

test_that("classify_qc_columns recognizes case-insensitive 'sample' variants", {
  for (col in c("sample_id", "SampleName", "Sample Name")) {
    df <- tibble(ID = "p1", `Test Condition` = "Test", logTransformed = 1)
    df[[col]] <- "s1"
    cls <- classify_qc_columns(df, "logTransformed")
    expect_equal(cls$sample_cols, col)
  }
})

test_that("classify_qc_columns errors when sample columns can't be identified", {
  df <- tibble(ID = "p1", Foo = "x", Bar = "y", logTransformed = 1)
  expect_error(classify_qc_columns(df, "logTransformed"), "could not identify sample column")
})

test_that("classify_qc_columns errors when 'sample' matches more than one column ambiguously", {
  df <- tibble(ID = "p1", Sample_A = "x", Sample_B = "y", logTransformed = 1)
  expect_error(classify_qc_columns(df, "logTransformed"), "could not identify sample column")
})

test_that("read_qc_dir parses assay type and replicate type from the new filename convention", {
  tmp <- tempfile("qc_dir_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  file.create(file.path(tmp, "QC_PTK_01_TR.csv"))
  file.create(file.path(tmp, "QC_PTK_TR.csv"))
  file.create(file.path(tmp, "QC_STK_01_BR.csv"))
  file.create(file.path(tmp, "QC_STK_LogCmb.csv")) # legacy filename, no TR/BR suffix

  qc_files <- read_qc_dir(paste0(tmp, "/"))

  expect_equal(nrow(qc_files), 4)
  ptk_rows <- qc_files %>% filter(Assay_Type == "PTK")
  expect_true(all(ptk_rows$Replicate_Type == "TR"))

  stk_new <- qc_files %>% filter(Assay_Type == "STK", grepl("01_BR", qc_file))
  expect_equal(stk_new$Replicate_Type, "BR")

  stk_legacy <- qc_files %>% filter(Assay_Type == "STK", grepl("LogCmb", qc_file))
  expect_equal(stk_legacy$Replicate_Type, "BR") # defaults to BR when suffix is missing
})

test_that("read_qc_dir finds TR/BR even when there's extra text after it in the filename", {
  tmp <- tempfile("qc_dir_")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  # TR/BR doesn't have to be the last underscore-separated element.
  file.create(file.path(tmp, "QC_PTK_01_BR_extra_description.csv"))
  file.create(file.path(tmp, "QC_PTK_02_BR.csv")) # no extra text - still works
  file.create(file.path(tmp, "QC_STK_01_TR_repeat1.csv"))

  qc_files <- read_qc_dir(paste0(tmp, "/"))

  ptk_rows <- qc_files %>% filter(Assay_Type == "PTK")
  expect_true(all(ptk_rows$Replicate_Type == "BR"))

  stk_rows <- qc_files %>% filter(Assay_Type == "STK")
  expect_equal(stk_rows$Replicate_Type, "TR")
})

test_that("read_qc_dir on the real QC_PTK.csv / QC_STK.csv fixtures defaults to BR", {
  qc_files <- read_qc_dir(paste0(qc_input_path(), "/"))
  legacy <- qc_files %>% filter(qc_file %in% c(qc_input_path("QC_PTK.csv"), qc_input_path("QC_STK.csv")))
  expect_true(all(legacy$Replicate_Type == "BR"))
})
