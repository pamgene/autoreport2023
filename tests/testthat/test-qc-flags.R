test_that("combine_flags matches the original sum-based formula for every 3-criteria combination", {
  # Original formula (R/01_BasicProcessing.R pre-redesign): overall_val <- sum(a,b,c);
  # overall_val == 9 -> good(3); overall_val > 4 -> fair(2); else poor(1).
  for (a in 1:3) for (b in 1:3) for (c in 1:3) {
    overall_val <- a + b + c
    expected <- if (overall_val == 9) 3 else if (overall_val > 4) 2 else 1
    expect_equal(combine_flags(c(a, b, c)), expected, info = sprintf("a=%d b=%d c=%d", a, b, c))
  }
})

test_that("combine_flags 2-criteria rule: good iff both good, poor iff both poor, else fair", {
  expect_equal(combine_flags(c(3, 3)), 3)
  expect_equal(combine_flags(c(1, 1)), 1)
  expect_equal(combine_flags(c(3, 1)), 2)
  expect_equal(combine_flags(c(3, 2)), 2)
  expect_equal(combine_flags(c(2, 1)), 2)
  expect_equal(combine_flags(c(2, 2)), 2)
})

test_that("combine_flags errors on an unsupported number of criteria", {
  expect_error(combine_flags(c(1)), "expected 2 or 3")
  expect_error(combine_flags(c(1, 1, 1, 1)), "expected 2 or 3")
})

test_that("flag_low_is_bad / flag_high_is_bad boundaries", {
  expect_equal(flag_low_is_bad(999, c(1000, 2000)), 1)
  expect_equal(flag_low_is_bad(1000, c(1000, 2000)), 2)
  expect_equal(flag_low_is_bad(1999, c(1000, 2000)), 2)
  expect_equal(flag_low_is_bad(2000, c(1000, 2000)), 3)

  expect_equal(flag_high_is_bad(0.36, c(0.35, 0.20)), 1)
  expect_equal(flag_high_is_bad(0.35, c(0.35, 0.20)), 2)
  expect_equal(flag_high_is_bad(0.21, c(0.35, 0.20)), 2)
  expect_equal(flag_high_is_bad(0.20, c(0.35, 0.20)), 3)
})

test_that("determine_flag reproduces the documented PTK/STK thresholds (2-criteria, no variability)", {
  good <- determine_flag(list(Signal = 5000, num_peptides = 150), "PTK")
  expect_equal(good$Signal_Flag, 3)
  expect_equal(good$Pep_Flag, 3)
  expect_true(is.na(good$Variability_Flag))
  expect_equal(good$Overall_Flag, 3)

  poor <- determine_flag(list(Signal = 500, num_peptides = 50), "PTK")
  expect_equal(poor$Overall_Flag, 1)

  # STK has different (lower) peptide thresholds (56/90) than PTK (78/123): 80
  # peptides is still "fair" for both here, but 95 clears STK's good threshold (90)
  # while still falling short of PTK's (123).
  ptk_95 <- determine_flag(list(Signal = 5000, num_peptides = 95), "PTK")
  stk_95 <- determine_flag(list(Signal = 5000, num_peptides = 95), "STK")
  expect_equal(ptk_95$Pep_Flag, 2)
  expect_equal(stk_95$Pep_Flag, 3)
})

test_that("determine_flag includes variability as a 3rd criterion when supplied", {
  flags <- determine_flag(list(Signal = 5000, num_peptides = 150, Variability = 0.5), "PTK")
  expect_equal(flags$Variability_Flag, 1)
  # signal=good, peptides=good, variability=poor -> not all-good, not (0-good & >=2-poor)
  # -> fair, per the fully explicit 3-criteria rule.
  expect_equal(flags$Overall_Flag, 2)

  all_good <- determine_flag(list(Signal = 5000, num_peptides = 150, Variability = 0.1), "PTK")
  expect_equal(all_good$Overall_Flag, 3)
})

test_that("classify_variability_tier(..., 'TR') matches determine_flag's variability thresholds", {
  # determine_flag always uses the TR band (0.30/0.20) - a BR row's variability is never a
  # flag, so it never reaches determine_flag at all.
  expect_equal(classify_variability_tier(0.10, "TR"), "Good")
  expect_equal(classify_variability_tier(0.20, "TR"), "Good")
  expect_equal(classify_variability_tier(0.25, "TR"), "Fair")
  expect_equal(classify_variability_tier(0.36, "TR"), "Poor")
})

test_that("classify_variability_tier uses a different (wider) band for 'BR'", {
  expect_equal(classify_variability_tier(0.25, "BR"), "Good")  # would be Fair under TR
  expect_equal(classify_variability_tier(0.30, "BR"), "Good")
  expect_equal(classify_variability_tier(0.35, "BR"), "Fair")
  expect_equal(classify_variability_tier(0.45, "BR"), "Poor")
})
