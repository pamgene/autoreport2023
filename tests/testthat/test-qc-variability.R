test_that("compute_condition_sd matches a hand-computed reference for a small dataset", {
  # Note: only 2 peptides per condition here, so median == mean for this particular
  # dataset - this test alone wouldn't catch a regression back to mean(). See the
  # dedicated "uses the median, not the mean" test below for that.
  df <- tribble(
    ~Condition, ~ID, ~Value,
    "Test",    "p1", 4.0,
    "Test",    "p1", 4.2,
    "Test",    "p2", 5.0,
    "Test",    "p2", 5.6,
    "Control", "p1", 3.0,
    "Control", "p1", 3.6,
    "Control", "p2", 2.0,
    "Control", "p2", 2.1,
  )

  result <- compute_condition_sd(df, "Condition", "ID", "Value")

  # Reference computed directly, independent of the function under test.
  expected <- df %>%
    group_by(Condition, ID) %>%
    summarise(peptide_sd = sd(Value), .groups = "drop") %>%
    group_by(Condition) %>%
    summarise(median_sd = median(peptide_sd), .groups = "drop") %>%
    arrange(Condition) %>%
    pull(median_sd)

  expect_equal(sort(result), sort(expected))
  expect_length(result, 2)
})

test_that("compute_condition_sd uses the median across peptides, not the mean", {
  # 3 peptides with deliberately skewed per-peptide SDs, so mean and median genuinely
  # differ - this is what actually guards against a regression back to mean().
  df <- tribble(
    ~Condition, ~ID, ~Value,
    "Test", "p1", 1,
    "Test", "p1", 2,   # p1 sd = 0.7071068
    "Test", "p2", 1,
    "Test", "p2", 3,   # p2 sd = 1.4142136
    "Test", "p3", 1,
    "Test", "p3", 10,  # p3 sd = 6.3639610
  )

  result <- compute_condition_sd(df, "Condition", "ID", "Value")

  peptide_sds <- c(sd(c(1, 2)), sd(c(1, 3)), sd(c(1, 10)))
  expect_equal(result, median(peptide_sds))
  expect_false(isTRUE(all.equal(result, mean(peptide_sds)))) # mean would give a different (wrong) answer
})

test_that("compute_condition_sd supports multiple condition columns (Supergroup + Test Condition)", {
  df <- tribble(
    ~Supergroup, ~`Test Condition`, ~ID, ~Value,
    "Sgroup1", "Test",    "p1", 4.0,
    "Sgroup1", "Test",    "p1", 4.4,
    "Sgroup1", "Control", "p1", 3.0,
    "Sgroup1", "Control", "p1", 3.2,
  )
  result <- compute_condition_sd(df, c("Supergroup", "Test Condition"), "ID", "Value")
  expect_length(result, 2)
  expect_equal(sort(round(result, 4)), sort(round(c(sd(c(4.0, 4.4)), sd(c(3.0, 3.2))), 4)))
})

test_that("stringify_sd formats a min-max range rounded to 2 decimals", {
  expect_equal(stringify_sd(c(0.1234, 0.5678, 0.3)), "0.12 - 0.57")
  expect_equal(stringify_sd(c(0.5, 0.5)), "0.5 - 0.5")
})

test_that("compute_condition_sd on the real QC_STK_01_BR.csv fixture (Log) matches a captured baseline", {
  raw <- read_delim(qc_input_path("QC_STK_01_BR.csv"), show_col_types = FALSE)
  raw <- clean_tercen_columns(raw)
  norms <- identify_value_columns(raw)
  cls <- classify_qc_columns(raw, unname(unlist(norms)))

  sd_log <- compute_condition_sd(raw, cls$condition_cols, cls$peptide_col, norms[["Log"]])
  # Baseline captured from a verified run of the implementation - a regression guard,
  # not an independently derived reference value.
  expect_equal(stringify_sd(sd_log), "0.08 - 1.67")

  sd_logcmb <- compute_condition_sd(raw, cls$condition_cols, cls$peptide_col, norms[["Log + ComBat"]])
  expect_equal(stringify_sd(sd_logcmb), "0.3 - 1.45")
})

test_that("compute_condition_sd on the real QC_STK_02_BR.csv fixture (VSN) matches a captured baseline", {
  raw <- read_delim(qc_input_path("QC_STK_02_BR.csv"), show_col_types = FALSE)
  raw <- clean_tercen_columns(raw)
  norms <- identify_value_columns(raw)
  cls <- classify_qc_columns(raw, unname(unlist(norms)))

  sd_vsn <- compute_condition_sd(raw, cls$condition_cols, cls$peptide_col, norms[["VSN"]])
  expect_equal(stringify_sd(sd_vsn), "0.07 - 0.74")

  sd_vsncmb <- compute_condition_sd(raw, cls$condition_cols, cls$peptide_col, norms[["VSN + ComBat"]])
  expect_equal(stringify_sd(sd_vsncmb), "0.07 - 0.68")
})
