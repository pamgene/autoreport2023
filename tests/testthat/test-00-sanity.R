test_that("QC fixture files are reachable via qc_input_path()", {
  expect_true(file.exists(qc_input_path("QC_PTK.csv")))
  expect_true(file.exists(qc_input_path("QC_STK_01_BR.csv")))
})
