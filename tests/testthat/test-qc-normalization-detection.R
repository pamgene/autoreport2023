# Tests for auto-detecting which normalization was applied, from the QC files
# themselves, so the "Normalizations" UI checkbox can't silently drift out of sync with
# what was actually uploaded (see R/01_BasicProcessing.R detect_normalizations_from_qc()
# / get_effective_normalizations()).

test_that("detect_normalizations_from_qc reads the real Log+ComBat and VSN+ComBat fixtures", {
  qc_files <- tibble(
    Assay_Type = factor(c("STK", "STK")),
    qc_file = c(test_input_path("QC_STK_01_BR.csv"), test_input_path("QC_STK_02_BR.csv")),
    Replicate_Type = c("BR", "BR")
  )

  # Both a Log-based and a VSN-based file are present -> "latest" per normalization_order
  # is VSN + ComBat -> both vsn and combat detected.
  expect_setequal(detect_normalizations_from_qc(qc_files, "tercen"), c("vsn", "combat"))
})

test_that("detect_normalizations_from_qc picks up ComBat without VSN when only a Log+ComBat file is present", {
  qc_files <- tibble(
    Assay_Type = factor("STK"),
    qc_file = test_input_path("QC_STK_01_BR.csv"),
    Replicate_Type = "BR"
  )
  expect_equal(detect_normalizations_from_qc(qc_files, "tercen"), "combat")
})

test_that("detect_normalizations_from_qc picks up VSN+ComBat when only the VSN file is present", {
  qc_files <- tibble(
    Assay_Type = factor("STK"),
    qc_file = test_input_path("QC_STK_02_BR.csv"),
    Replicate_Type = "BR"
  )
  expect_setequal(detect_normalizations_from_qc(qc_files, "tercen"), c("vsn", "combat"))
})

test_that("detect_normalizations_from_qc returns empty (not NULL) for a plain log-only legacy file", {
  qc_files <- tibble(
    Assay_Type = factor("PTK"),
    qc_file = test_input_path("QC_PTK.csv"),
    Replicate_Type = "BR"
  )
  result <- detect_normalizations_from_qc(qc_files, "tercen")
  expect_length(result, 0)
  expect_false(is.null(result)) # must stay distinguishable from "detection impossible"
})

test_that("detect_normalizations_from_qc returns NULL (detection impossible) for bionav, or no QC files", {
  qc_files <- tibble(
    Assay_Type = factor("STK"),
    qc_file = test_input_path("QC_STK_01_BR.csv"),
    Replicate_Type = "BR"
  )
  expect_null(detect_normalizations_from_qc(qc_files, "bionav"))
  expect_null(detect_normalizations_from_qc(data.frame(), "tercen"))
})

test_that("get_effective_normalizations prefers detection over the manual param when both are available", {
  qc_files <- tibble(
    Assay_Type = factor("STK"),
    qc_file = test_input_path("QC_STK_02_BR.csv"),
    Replicate_Type = "BR"
  )
  # Manual param says "nothing" but the QC file is VSN+ComBat - detection should win.
  expect_setequal(get_effective_normalizations(qc_files, "tercen", NULL), c("vsn", "combat"))
})

test_that("get_effective_normalizations falls back to the manual param for bionav", {
  qc_files <- tibble(
    Assay_Type = factor("STK"),
    qc_file = test_input_path("QC_STK_02_BR.csv"),
    Replicate_Type = "BR"
  )
  expect_equal(get_effective_normalizations(qc_files, "bionav", c("combat")), c("combat"))
})

test_that("get_effective_normalizations falls back to the manual param when there are no QC files yet", {
  expect_equal(get_effective_normalizations(data.frame(), "tercen", c("vsn")), c("vsn"))
  expect_null(get_effective_normalizations(data.frame(), "tercen", NULL))
})
