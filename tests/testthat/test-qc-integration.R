# Integration tests: exercise the full pipeline (read_qc_dir -> parse_qc -> render_* ->
# build_qc_result_text) against real QC input files, per README's "use data/test input
# files to test all inputs" convention. Real "main input" fixtures used:
#   - data/test inputs/QC_PTK.csv, QC_STK.csv       (legacy tercen format, no Supergroup)
#   - data/test inputs/QC_STK_01_BR.csv (Log)        (new format, biological replicates)
#   - data/test inputs/QC_STK_02_BR.csv (VSN)        (new format, biological replicates)
# No real technical-replicate ("_TR") fixture exists yet, so the TR and mixed scenarios
# copy the real BR files into a temp dir under TR-suffixed names - the file *content* is
# still a real input, only the replicate-type label is synthesized.
#
# The Data Variability Indicator (variability_table / render_variability_boxplot()) is
# BR-only, always: an assay type with no BR file gets no entry in variability_table at
# all, regardless of how many normalizations its TR file(s) have.

make_qc_dir_from <- function(mapping) {
  # mapping: named vector, new_filename -> source fixture path
  tmp <- tempfile("qc_dir_")
  dir.create(tmp)
  for (new_name in names(mapping)) {
    file.copy(mapping[[new_name]], file.path(tmp, new_name))
  }
  tmp
}

test_that("BR scenario: real QC_STK_01_BR.csv + QC_STK_02_BR.csv end-to-end", {
  dir <- make_qc_dir_from(c(
    "QC_STK_01_BR.csv" = qc_input_path("QC_STK_01_BR.csv"),
    "QC_STK_02_BR.csv" = qc_input_path("QC_STK_02_BR.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_equal(nrow(qc_files), 2)
  expect_true(all(qc_files$Replicate_Type == "BR"))

  parsed <- parse_qc(qc_files, "tercen")

  # BR-only assay (no TR file at all) -> still gets a QC-results row, from its BR data,
  # using only the 2 criteria a BR flag can support (signal, peptides) - no per-replicate
  # noise here to support a 3rd variability criterion.
  expect_equal(nrow(parsed$qc_table), 1)
  expect_equal(parsed$qc_table$Assay_Label, "STK (VSN + ComBat)") # latest of the 4 normalizations present
  expect_equal(parsed$qc_table$Signal_Flag, 3)
  expect_equal(parsed$qc_table$Pep_Flag, 2)
  expect_true(is.na(parsed$qc_table$Variability_Flag)) # not part of a BR-fallback flag
  expect_equal(parsed$qc_table$Variability_String, "—")
  expect_equal(parsed$qc_table$Overall_Flag, 2) # combine_flags_br(c(good, fair)) -> fair

  # BR -> variability table has all 4 normalization approaches (both files cover Log,
  # Log+ComBat, VSN, VSN+ComBat between them) - the Data Variability Indicator is
  # unaffected by there also being a QC flag now.
  expect_equal(nrow(parsed$variability_table), 4)
  expect_setequal(parsed$variability_table$Normalization, c("Log", "Log + ComBat", "VSN", "VSN + ComBat"))
  expect_true(all(parsed$variability_table$Assay_Type == "STK"))
  expect_true(should_show_variability_section(qc_files))

  # Render functions run cleanly and produce the expected column shape.
  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_s3_class(t1, "flextable")
  expect_false(any(grepl("^Tech_Variability", colnames(t1$body$dataset)))) # no TR columns

  # Table 2 itself: since every row's Variability_Flag is NA here (all-BR study), the
  # Technical Variability column is omitted entirely, not just blank.
  t2 <- render_qc_results_table(parsed$qc_table)
  expect_false("Variability_String" %in% t2$col_keys)

  # Data Variability Indicator: one boxplot, STK panel only (the only assay type with BR
  # data here), reading the precomputed historical distribution from
  # data/qc_variability_distribution.rds.
  p <- with_repo_root(function() render_variability_boxplot(parsed$variability_table))
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_setequal(unique(built$plot$data$Assay_Type), "STK")

  qc_result <- build_qc_result_text(parsed$qc_table)
  expect_match(qc_result, 'the STK QC flag was "fair"')
})

test_that("TR scenario: real STK fixture data, relabeled as technical replicates", {
  dir <- make_qc_dir_from(c(
    "QC_STK_01_TR.csv" = qc_input_path("QC_STK_01_BR.csv"),
    "QC_STK_02_TR.csv" = qc_input_path("QC_STK_02_BR.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_true(all(qc_files$Replicate_Type == "TR"))

  parsed <- parse_qc(qc_files, "tercen")

  # TR -> only the "latest" normalization approach appears in Table 2 (4 are present
  # across the two files: Log, Log + ComBat, VSN, VSN + ComBat - latest is VSN + ComBat).
  expect_equal(nrow(parsed$qc_table), 1)
  expect_equal(parsed$qc_table$Normalization, "VSN + ComBat")
  expect_equal(parsed$qc_table$Assay_Label, "STK (VSN + ComBat)")
  expect_true(!is.na(parsed$qc_table$Variability_Flag))
  # Table 2 shows the single median across conditions (matching the real fixture baseline
  # from test-qc-variability.R's sd_vsncmb vector), not a min-max range - and the flag
  # itself is now driven by that same median, not the worst-case max.
  expect_equal(parsed$qc_table$Variability_String, "0.25")

  # No BR file anywhere in this study -> the Data Variability Indicator is BR-only, so it
  # gets nothing here at all, regardless of Table 2 having multiple normalizations.
  expect_equal(nrow(parsed$variability_table), 0)
  expect_false(should_show_variability_section(qc_files))

  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_true("Tech_Variability_STK" %in% colnames(t1$body$dataset))

  t2 <- render_qc_results_table(parsed$qc_table)
  expect_true("Variability_String" %in% t2$col_keys)
  # htmltools_value() (not information_data_chunk(), which isn't available across all
  # flextable versions this project has been run under) - not the abbreviated "Tech
  # Variability".
  t2_html <- with_repo_root(function() as.character(flextable::htmltools_value(t2)))
  expect_true(grepl("Technical Variability", t2_html))

  qc_result <- build_qc_result_text(parsed$qc_table)
  expect_match(qc_result, "the STK QC flag was")
  expect_no_match(qc_result, "Biological variability") # no BR assay type at all
})

test_that("Mixed scenario: PTK=TR, STK=BR, both from real fixture data", {
  dir <- make_qc_dir_from(c(
    "QC_PTK_01_TR.csv" = qc_input_path("QC_STK_01_BR.csv"), # content reused, assay relabeled
    "QC_STK_01_BR.csv" = qc_input_path("QC_STK_01_BR.csv"),
    "QC_STK_02_BR.csv" = qc_input_path("QC_STK_02_BR.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_equal(qc_files$Replicate_Type[qc_files$Assay_Type == "PTK"], "TR")
  expect_true(all(qc_files$Replicate_Type[qc_files$Assay_Type == "STK"] == "BR"))

  parsed <- parse_qc(qc_files, "tercen")

  ptk_rows <- parsed$qc_table %>% filter(Assay_Type == "PTK")
  stk_rows <- parsed$qc_table %>% filter(Assay_Type == "STK")
  expect_equal(nrow(ptk_rows), 1) # TR -> only the latest normalization (Log + ComBat)
  expect_equal(ptk_rows$Normalization, "Log + ComBat")
  expect_true(!is.na(ptk_rows$Variability_Flag)) # TR -> real 3rd criterion

  # STK (BR-only, no TR file) still gets its own row now, from the BR-fallback flag - just
  # 2 criteria, no variability flag.
  expect_equal(nrow(stk_rows), 1)
  expect_equal(stk_rows$Normalization, "VSN + ComBat")
  expect_true(is.na(stk_rows$Variability_Flag))

  # Variability table (BR-only) includes STK only - PTK has no BR file at all here, so it
  # gets no entry, even though its TR file has >1 normalization present.
  expect_setequal(parsed$variability_table$Assay_Type, "STK")
  expect_equal(nrow(parsed$variability_table), 4)

  # Table 1 gets exactly one Tech Variability column (PTK), not two.
  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_true("Tech_Variability_PTK" %in% colnames(t1$body$dataset))
  expect_false("Tech_Variability_STK" %in% colnames(t1$body$dataset))

  # Boxplot only ever shows the STK panel here - PTK has no BR data to plot.
  p <- with_repo_root(function() render_variability_boxplot(parsed$variability_table))
  built <- ggplot2::ggplot_build(p)
  expect_setequal(unique(built$plot$data$Assay_Type), "STK")

  qc_result <- build_qc_result_text(parsed$qc_table)
  expect_match(qc_result, "the PTK QC flag was") # TR -> flagged
  expect_match(qc_result, "the STK QC flag was") # BR-fallback -> also flagged now
  expect_no_match(qc_result, "Biological variability") # sentence removed entirely
})

test_that("legacy tercen fixtures (QC_PTK.csv, QC_STK.csv, no Supergroup/CmbCor) parse end-to-end", {
  # Tagged _TR (rather than the legacy files' natural BR default) so a QC-results row
  # actually gets built here - this test's purpose is exercising legacy-format column
  # parsing (no Supergroup/CmbCor), not the TR/BR flagging rule itself.
  dir <- make_qc_dir_from(c(
    "QC_PTK_01_TR.csv" = qc_input_path("QC_PTK.csv"),
    "QC_STK_01_TR.csv" = qc_input_path("QC_STK.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  parsed <- parse_qc(qc_files, "tercen")

  expect_equal(nrow(parsed$qc_table), 2) # one row per assay
  expect_true(all(parsed$qc_table$Signal > 0))
  expect_true(all(parsed$qc_table$num_peptides > 0))
  # No BR file anywhere -> variability_table (BR-only) is empty.
  expect_equal(nrow(parsed$variability_table), 0)

  qc_result <- build_qc_result_text(parsed$qc_table)
  expect_match(qc_result, "the PTK QC flag was")
  expect_match(qc_result, "the STK QC flag was")

  # Table 1: PTK and STK are both TR here, and share identical Technical Variability
  # thresholds (variability_criteria is not assay-specific) - one merged column, not two.
  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_true("Tech_Variability" %in% colnames(t1$body$dataset))
  expect_false(any(c("Tech_Variability_PTK", "Tech_Variability_STK") %in% colnames(t1$body$dataset)))
})

test_that("a pure-TR study (no BR at all) does NOT trigger the variability section, regardless of normalization count", {
  dir <- make_qc_dir_from(c("QC_PTK_01_TR.csv" = qc_input_path("QC_PTK.csv")))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  parsed <- parse_qc(qc_files, "tercen")

  expect_equal(nrow(parsed$qc_table), 1)
  expect_equal(parsed$qc_table$Normalization, "Log")
  expect_equal(nrow(parsed$variability_table), 0) # BR-only -> nothing here, no BR file exists
  expect_false(should_show_variability_section(qc_files))
})

test_that("one assay type with BOTH a TR file and a BR file: Table 2 uses TR, the boxplot uses BR", {
  # Per ref/HowToUse.md's documented workflow: export raw per-array TR data, then also
  # export a biological-replicate-level file averaged up from it. QC_STK_01_BR.csv (has
  # Log, Log + ComBat) stands in for the TR export; QC_STK_02_BR.csv (has VSN,
  # VSN + ComBat) stands in for the separately-computed BR export - same assay type
  # (STK), two different replicate-type files.
  dir <- make_qc_dir_from(c(
    "QC_STK_01_TR.csv" = qc_input_path("QC_STK_01_BR.csv"),
    "QC_STK_02_BR.csv" = qc_input_path("QC_STK_02_BR.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  expect_setequal(qc_files$Replicate_Type, c("TR", "BR"))
  expect_true(all(qc_files$Assay_Type == "STK"))

  parsed <- parse_qc(qc_files, "tercen")

  # Table 2: TR-driven - one row, latest of the TR file's normalizations (Log, Log +
  # ComBat), with a real variability flag and its median (not range) displayed.
  expect_equal(nrow(parsed$qc_table), 1)
  expect_equal(parsed$qc_table$Normalization, "Log + ComBat")
  expect_equal(parsed$qc_table$Assay_Label, "STK (Log + ComBat)")
  expect_true(!is.na(parsed$qc_table$Variability_Flag))
  expect_equal(parsed$qc_table$Variability_String, "0.55")

  # The Data Variability Indicator is BR-only, always - it uses the BR file's
  # normalizations (VSN, VSN + ComBat), never the TR file's (Log, Log + ComBat), with no
  # per-row "which replicate type" field needed since every row is BR by construction.
  expect_setequal(parsed$variability_table$Normalization, c("VSN", "VSN + ComBat"))
  expect_true(should_show_variability_section(qc_files))

  # Table 1 still gets a Tech Variability column (TR data exists for this assay).
  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_true("Tech_Variability_STK" %in% colnames(t1$body$dataset))

  p <- with_repo_root(function() render_variability_boxplot(parsed$variability_table))
  built <- ggplot2::ggplot_build(p)
  expect_setequal(unique(built$plot$data$Assay_Type), "STK")

  qc_result <- build_qc_result_text(parsed$qc_table)
  expect_match(qc_result, "the STK QC flag was")
  expect_no_match(qc_result, "Biological variability") # sentence removed entirely
})
