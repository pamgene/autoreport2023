# Integration tests for UKA kinase-analysis inputs, per README's "use data/test input
# files to test all inputs" convention. Fixtures: data/test inputs/UKA_PTK_01_
# csUKA-180307Breast.csv, UKA_STK_01_csUKA-180307Breast.csv - real "all-vs-all" UKA
# output (identified by a `contrast` column). The older UKA_MTvC/UKA_TGC app formats
# (filenames containing "_ukam-"/"_ukat-", handled by process_ukat_ukam()) are no longer
# used in practice, so they're intentionally not covered here.
#
# Scope note: only read_kinase_dir() (R/00_GeneralFunctions.R) is exercised. The
# downstream plotting/table functions in R/03_KinaseAnalysis.R pull in a heavy,
# partly non-CRAN dependency chain (CORALcli, Cairo, magick, rsvg, ...) unrelated to
# input handling, so they're out of scope here.

make_kinase_input_dir <- function(files) {
  dir <- make_input_dir(files, subfolder = "03_Kinase Analysis")
  dir.create(file.path(dir, "02_DATA"))
  dir
}

test_that("read_kinase_dir splits a real all-vs-all UKA file by comparison", {
  dir <- make_kinase_input_dir(c(
    "UKA_PTK_01_csUKA-180307Breast.csv" = test_input_path("UKA_PTK_01_csUKA-180307Breast.csv")
  ))
  old_wd <- setwd(dir)
  on.exit({
    setwd(old_wd)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)

  kinase_files <- read_kinase_dir(folder = "03_Kinase Analysis/", csUKA = TRUE)

  # 5 distinct Sgroup_contrast values in the fixture -> 5 split files.
  expect_equal(nrow(kinase_files), 5)
  expect_true(all(kinase_files$Assay_Type == "PTK"))
  expect_equal(kinase_files$Order, 1:5)
  expect_true(all(file.exists(kinase_files$UKA_file)))
  expect_true(all(grepl("^Sgroup1 - ", kinase_files$Comparison)))

  # Each split file has one comparison's worth of ranked kinase rows.
  split_df <- read_delim(kinase_files$UKA_file[1], show_col_types = FALSE)
  expect_true(nrow(split_df) > 0)
  expect_true(all(c("Kinase Name", "Final score", "Specificity Score") %in% colnames(split_df)))

  # The cleaned key-column table (same name as the upload) is archived to 02_DATA/;
  # the full "_all_columns" copy goes to the 02_DATA/Extended data/ subfolder.
  archived <- list.files(file.path(dir, "02_DATA"))
  expect_true("UKA_PTK_01_csUKA-180307Breast.csv" %in% archived)
  extended <- list.files(file.path(dir, "02_DATA", "Extended data"))
  expect_true("UKA_PTK_01_csUKA-180307Breast_all_columns.csv" %in% extended)
})

test_that("read_kinase_dir numbers PTK and STK files with one shared counter, not reset per assay", {
  dir <- make_kinase_input_dir(c(
    "UKA_PTK_01_csUKA-180307Breast.csv" = test_input_path("UKA_PTK_01_csUKA-180307Breast.csv"),
    "UKA_STK_01_csUKA-180307Breast.csv" = test_input_path("UKA_STK_01_csUKA-180307Breast.csv")
  ))
  old_wd <- setwd(dir)
  on.exit({
    setwd(old_wd)
    unlink(dir, recursive = TRUE)
  }, add = TRUE)

  kinase_files <- read_kinase_dir(folder = "03_Kinase Analysis/", csUKA = TRUE)

  expect_equal(nrow(kinase_files), 10)
  expect_equal(kinase_files$Order, 1:10)
  # This documents actual current behavior (the split-file counter in
  # process_uka_allvsall() is shared across all input files in one call, not reset per
  # assay type) rather than asserting it's the ideal behavior - PTK ends up 01-05 and
  # STK continues at 06-10 instead of also starting at 01.
  ptk_orders <- kinase_files$Order[kinase_files$Assay_Type == "PTK"]
  stk_orders <- kinase_files$Order[kinase_files$Assay_Type == "STK"]
  expect_equal(ptk_orders, 1:5)
  expect_equal(stk_orders, 6:10)
})

test_that("read_kinase_dir returns an empty data frame (with a warning) when no files are present", {
  dir <- make_kinase_input_dir(character(0))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  expect_warning(result <- read_kinase_dir(folder = file.path(dir, "03_Kinase Analysis/"), csUKA = TRUE), "No kinase files")
  expect_equal(nrow(result), 0)
})
