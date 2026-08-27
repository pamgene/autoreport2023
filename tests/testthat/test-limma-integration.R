# Integration tests for Limma phosphosite-analysis inputs, per README's "use data/test
# input files to test all inputs" convention. Fixtures: data/test inputs/Limma_PTK_01_
# Supergroup.csv, Limma_STK_01_Supergroup.csv (tercen format).

test_that("read_phosphosite_dir parses Order/Comparison/Stats/Assay_Type/Group from real Limma files", {
  dir <- make_input_dir(c(
    "Limma_PTK_01_Supergroup.csv" = test_input_path("Limma_PTK_01_Supergroup.csv"),
    "Limma_STK_01_Supergroup.csv" = test_input_path("Limma_STK_01_Supergroup.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  result <- read_phosphosite_dir(paste0(dir, "/"), datatype = "tercen")

  expect_equal(nrow(result), 2)
  expect_setequal(result$Assay_Type, c("PTK", "STK"))
  expect_true(all(result$Stats == "Limma"))
  expect_true(all(result$Comparison == "Limma"))
  expect_true(all(result$Group == "Supergroup"))
  expect_true(all(result$Order == 1))
  expect_true(all(file.exists(result$File)))
})

test_that("Limma fixture content is readable and has the expected tercen columns after cleaning", {
  raw <- read_delim(test_input_path("Limma_PTK_01_Supergroup.csv"), show_col_types = FALSE)
  raw <- clean_tercen_columns(raw)

  expect_true(all(c("Supergroup_annot", "contrast", "ID", "logFC", "pvalue", "FDR") %in% colnames(raw)))
  expect_true(nrow(raw) > 0)
})

test_that("read_phosphosite_dir returns an empty data frame (with a warning) when no files are present", {
  dir <- make_input_dir(character(0))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  expect_warning(result <- read_phosphosite_dir(paste0(dir, "/"), datatype = "tercen"), "No stats files")
  expect_equal(nrow(result), 0)
})
