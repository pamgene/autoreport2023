# Standalone script - NOT sourced by app.R or the Rmds. Run once (or whenever
# data/QC_for_SD_distribution/ changes) from the repo root:
#   Rscript R/calculate_variability_distribution.R
#
# Computes the historical Data Variability Indicator distribution from every QC file in
# data/QC_for_SD_distribution/ and persists it to data/qc_variability_distribution.rds,
# so the Supplement can load it at render time without re-reading the source CSVs.
#
# Every file in that folder is treated as BR (see ref/median_sd_variability_notes.md and
# plan discussion: the Data Variability Indicator is BR-only, and this corpus has no
# TR/BR filename tags to say otherwise). Reuses clean_tercen_columns(), identify_value_
# columns(), classify_qc_columns(), compute_condition_sd() from R/01_BasicProcessing.R
# unchanged - only sourced here, never edited.

source("R/00_GeneralFunctions.R")
source("R/01_BasicProcessing.R")

# Loads one historical QC file and returns its per-condition median-peptide-SD values
# (one row per condition, per normalization present in the file). Two data-cleaning steps
# are needed that read_qc_dir()/load_qc_file() don't do, because they're specific to this
# ad hoc historical corpus (see ref/median_sd_variability_notes.md investigation):
#   1. Drop readr's auto-named blank trailing column (many files have a trailing comma in
#      their header), before clean_tercen_columns() can turn "...7" into a bare "7" that
#      would otherwise be picked up as a spurious condition column.
#   2. Drop rows with an empty/NA peptide ID - almost every file has a garbage first data
#      row (empty Barcode/ID/value, only Supergroup populated), an export artifact.
load_historical_file <- function(path) {
  df <- read_delim(path, show_col_types = FALSE)
  df <- df[, !grepl("^\\.\\.\\.[0-9]+$", colnames(df))]
  df <- clean_tercen_columns(df)
  df <- df[, colnames(df) != ""]

  norms <- identify_value_columns(df)
  if (length(norms) == 0) {
    return(list(rows = NULL, empty_id_dropped = 0L, singleton_dropped = 0L))
  }

  n_before <- nrow(df)
  df <- df %>% filter(!is.na(ID), ID != "")
  empty_id_dropped <- n_before - nrow(df)

  value_cols <- unname(unlist(norms))
  classification <- classify_qc_columns(df, value_cols)
  assay_type <- if (grepl("_PTK_", basename(path))) "PTK" else "STK"

  rows <- list()
  singleton_dropped <- 0L
  for (label in names(norms)) {
    sds <- compute_condition_sd(df, classification$condition_cols, classification$peptide_col, norms[[label]])
    n_na <- sum(is.na(sds))
    singleton_dropped <- singleton_dropped + n_na
    sds <- sds[!is.na(sds)]
    if (length(sds) > 0) {
      rows[[length(rows) + 1]] <- tibble(
        Assay_Type = assay_type, Source_File = basename(path), Median_SD = sds
      )
    }
  }

  list(rows = bind_rows(rows), empty_id_dropped = empty_id_dropped, singleton_dropped = singleton_dropped)
}

hist_files <- list.files("data/QC_for_SD_distribution", pattern = "[.]csv$", full.names = TRUE)
results <- lapply(hist_files, load_historical_file)

variability_distribution <- bind_rows(lapply(results, `[[`, "rows"))
total_empty_id_dropped <- sum(sapply(results, `[[`, "empty_id_dropped"))
total_singleton_dropped <- sum(sapply(results, `[[`, "singleton_dropped"))

cat("Historical Data Variability Indicator distribution:\n")
cat("  Files processed:", length(hist_files), "\n")
cat("  Garbage rows dropped (empty/NA peptide ID):", total_empty_id_dropped, "\n")
cat("  Single-replicate conditions dropped (SD undefined):", total_singleton_dropped, "\n")
cat("  Rows per assay type:\n")
print(table(variability_distribution$Assay_Type))

saveRDS(variability_distribution, "data/qc_variability_distribution.rds")
cat("\nSaved: data/qc_variability_distribution.rds\n")
