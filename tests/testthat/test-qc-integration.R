# Integration tests: exercise the full pipeline (read_qc_dir -> parse_qc -> render_* ->
# build_qc_result_text) against real QC input files, per README's "use data/test input
# files to test all inputs" convention. Real "main input" fixtures used:
#   - data/test inputs/QC_PTK.csv, QC_STK.csv       (legacy tercen format, no Supergroup)
#   - data/test inputs/QC_STK_01_BR.csv (Log)        (new format, biological replicates)
#   - data/test inputs/QC_STK_02_BR.csv (VSN)        (new format, biological replicates)
# No real technical-replicate ("_TR") fixture exists yet, so the TR and mixed scenarios
# copy the real BR files into a temp dir under TR-suffixed names - the file *content* is
# still a real input, only the replicate-type label is synthesized.

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

  # BR -> exactly one QC-results row for STK, no per-normalization breakdown in it.
  expect_equal(nrow(parsed$qc_table), 1)
  expect_equal(parsed$qc_table$Assay_Type, "STK")
  expect_equal(parsed$qc_table$Variability_String, "—")
  expect_true(is.na(parsed$qc_table$Variability_Flag))
  expect_true(parsed$qc_table$Overall_Flag %in% 1:3)

  # BR -> variability table has all 4 normalization approaches (both files cover Log,
  # Log+ComBat, VSN, VSN+ComBat between them).
  expect_equal(nrow(parsed$variability_table), 4)
  expect_setequal(parsed$variability_table$Normalization, c("Log", "Log + ComBat", "VSN", "VSN + ComBat"))
  expect_true(all(parsed$variability_table$Assay_Type == "STK"))
  expect_true(should_show_variability_section(qc_files, parsed$variability_table))

  # Render functions run cleanly and produce the expected column shape.
  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_s3_class(t1, "flextable")
  expect_false(any(grepl("^Tech_Variability", colnames(t1$body$dataset)))) # no TR columns

  t2 <- render_qc_results_table(parsed$qc_table)
  expect_s3_class(t2, "flextable")
  expect_false("Variability_String" %in% t2$col_keys) # no Tech Variability column in BR-only table

  t3 <- render_variability_ref_table(parsed$variability_table)
  expect_s3_class(t3, "flextable")
  # BR-only study -> a single, unsuffixed "Variability" column, using the BR (biological)
  # band, not the TR (technical) one.
  expect_true("Variability" %in% colnames(t3$body$dataset))
  expect_false(any(grepl("_BR$|_TR$", colnames(t3$body$dataset))))
  expect_equal(t3$body$dataset$Variability, c("< 0.3", "0.3 - 0.4", "> 0.4"))

  t4 <- render_variability_table(parsed$variability_table)
  expect_s3_class(t4, "flextable")
  expect_true("Variability_STK" %in% colnames(t4$body$dataset))
  # Generic caption - no "Biological"/"Technical" suffix, since Table 4 can now cover
  # either replicate type (or both).
  expect_equal(t4$caption$value, "Data Variability Indicator")

  qc_result <- build_qc_result_text(parsed$qc_table, parsed$variability_table, qc_files)
  expect_match(qc_result, "the STK QC flag was")
  expect_match(qc_result, "Biological variability was")
  expect_no_match(qc_result, "\\(STK\\)") # single assay type in the whole study -> no parenthetical
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

  # TR with >1 normalization -> the full per-normalization detail still goes into
  # variability_table, for Table 4 (even though there's no BR assay at all here).
  expect_equal(nrow(parsed$variability_table), 4)
  expect_setequal(parsed$variability_table$Normalization, c("Log", "Log + ComBat", "VSN", "VSN + ComBat"))
  expect_true(all(parsed$variability_table$Assay_Type == "STK"))
  # Pure TR, but this assay has >1 normalization present -> section still shows, even
  # though there's no BR assay anywhere in the study.
  expect_true(should_show_variability_section(qc_files, parsed$variability_table))

  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_true("Tech_Variability_STK" %in% colnames(t1$body$dataset))

  t2 <- render_qc_results_table(parsed$qc_table)
  expect_true("Variability_String" %in% t2$col_keys)

  t3 <- render_variability_ref_table(parsed$variability_table)
  # Pure TR study -> a single, unsuffixed "Variability" column, using the TR (technical)
  # band, tighter than BR's.
  expect_true("Variability" %in% colnames(t3$body$dataset))
  expect_equal(t3$body$dataset$Variability, c("< 0.2", "0.2 - 0.3", "> 0.3"))

  qc_result <- build_qc_result_text(parsed$qc_table, parsed$variability_table, qc_files)
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
  expect_equal(nrow(stk_rows), 1) # BR -> single row
  expect_true(!is.na(ptk_rows$Variability_Flag))
  expect_true(is.na(stk_rows$Variability_Flag))

  # Variability table now includes both: STK (BR, always) and PTK (TR, because it has
  # >1 normalization present - Log and Log + ComBat).
  expect_setequal(parsed$variability_table$Assay_Type, c("PTK", "STK"))
  expect_setequal(
    parsed$variability_table$Normalization[parsed$variability_table$Assay_Type == "PTK"],
    c("Log", "Log + ComBat")
  )

  # Table 1 gets exactly one Tech Variability column (PTK), not two.
  t1 <- with_repo_root(function() render_ref_qc_table(qc_files))
  expect_true("Tech_Variability_PTK" %in% colnames(t1$body$dataset))
  expect_false("Tech_Variability_STK" %in% colnames(t1$body$dataset))

  # Table 3 gets both bands side by side (STK contributes BR rows, PTK contributes TR
  # rows), each using its own threshold set.
  t3 <- render_variability_ref_table(parsed$variability_table)
  expect_true(all(c("Variability_BR", "Variability_TR") %in% colnames(t3$body$dataset)))
  expect_equal(t3$body$dataset$Variability_BR, c("< 0.3", "0.3 - 0.4", "> 0.4"))
  expect_equal(t3$body$dataset$Variability_TR, c("< 0.2", "0.2 - 0.3", "> 0.3"))

  qc_result <- build_qc_result_text(parsed$qc_table, parsed$variability_table, qc_files)
  expect_match(qc_result, "the PTK QC flag was")
  expect_match(qc_result, "the STK QC flag was")
  expect_match(qc_result, "Biological variability was")
  expect_match(qc_result, "\\(STK\\)") # 2 assay types total -> parenthetical needed
})

test_that("legacy tercen fixtures (QC_PTK.csv, QC_STK.csv, no Supergroup/CmbCor) parse end-to-end", {
  dir <- make_qc_dir_from(c(
    "QC_PTK_01_BR.csv" = qc_input_path("QC_PTK.csv"),
    "QC_STK_01_BR.csv" = qc_input_path("QC_STK.csv")
  ))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  parsed <- parse_qc(qc_files, "tercen")

  expect_equal(nrow(parsed$qc_table), 2) # one row per assay, both BR
  expect_true(all(parsed$qc_table$Signal > 0))
  expect_true(all(parsed$qc_table$num_peptides > 0))
  # Only "Log" normalization is available (no CmbCor in these legacy files).
  expect_setequal(parsed$variability_table$Normalization, "Log")

  qc_result <- build_qc_result_text(parsed$qc_table, parsed$variability_table, qc_files)
  expect_match(qc_result, "the PTK QC flag was")
  expect_match(qc_result, "the STK QC flag was")
})

test_that("a pure-TR study with only one normalization does NOT trigger the variability section", {
  # QC_PTK.csv has only "Log" (no CmbCor) - a single normalization, no BR anywhere in
  # the study - Table 2's one row already says everything there is to say.
  dir <- make_qc_dir_from(c("QC_PTK_01_TR.csv" = qc_input_path("QC_PTK.csv")))
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  qc_files <- read_qc_dir(paste0(dir, "/"))
  parsed <- parse_qc(qc_files, "tercen")

  expect_equal(nrow(parsed$qc_table), 1)
  expect_equal(parsed$qc_table$Normalization, "Log")
  expect_equal(nrow(parsed$variability_table), 1) # still populated...
  expect_false(should_show_variability_section(qc_files, parsed$variability_table)) # ...but not shown
})
