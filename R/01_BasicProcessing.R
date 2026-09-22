require(tidyverse)
require(flextable)
source("R/report_theme.R")

#########################################################
################### QC CRITERIA CONFIG ##################
#########################################################

# Signal / QC_peptides: c(poor_below, fair_below) - value >= fair_below is "good".
qc_criteria <- list(
  "PTK" = list(
    "Signal" = c(1000, 2000),
    "QC_peptides" = c(78, 123)
  ),
  "STK" = list(
    "Signal" = c(1000, 2000),
    "QC_peptides" = c(56, 90)
  )
)

# Variability: c(poor_above, fair_above) - value <= fair_above is "good". Not assay-specific
# (PTK and STK share the same bands, matching the original pre-redesign code, where bio_cv
# was also identical for both assays). Differs by replicate type: BR is informational only
# (Low/Medium/High tiers, Table 3), TR is a real flagging criterion (Good/Fair/Poor). Values
# reuse the original pre-redesign CV% thresholds directly as SD thresholds (BR <- old
# bio_cv, TR <- old tech_cv, per the reference table's documented bands - tech_cv itself was
# never actually computed pre-redesign, only ever hardcoded to "n.a."). CV and SD are
# different units, but per explicit product decision these are indicative bands only, not a
# rigorous conversion. TODO (future work, no data yet): replace these with bands derived
# from real reference data - compute quartiles across a large sample of past studies and
# classify each new report against them; the same applies to the Signal and QC_peptides
# bands above, which are likewise illustrative placeholders rather than data-derived.
variability_criteria <- list(
  "BR" = c(0.40, 0.30),
  "TR" = c(0.30, 0.20)
)

variability_tier_labels <- list("Good" = "Low", "Fair" = "Medium", "Poor" = "High")

normalization_order <- c("Log", "Log + ComBat", "VSN", "VSN + ComBat")

# Maps a normalization tag as it appears in a QC filename (case-insensitive) to its
# canonical label. Used as a fallback when a file's value column has been reduced to a
# generic "value" (e.g. a BR file that only exists because it was aggregated to the
# biological-replicate level from technical-replicate rows, which drops the semantic
# logTransformed/identity/CmbCor column name).
normalization_filename_tags <- list(
  "log" = "Log", "logcmb" = "Log + ComBat", "vsn" = "VSN", "vsncmb" = "VSN + ComBat"
)

#########################################################
##################### LOAD QC FILE ######################
#########################################################

# Reads and cleans a single QC file. For tercen data, column names are normalized via
# clean_tercen_columns() so downstream logic can rely on fixed names (ID, Barcode, Row,
# Supergroup, Test Condition, logTransformed/identity, CmbCor).
load_qc_file <- function(qc_file, datatype) {
  if (datatype == "bionav") {
    qc_df <- read_delim(qc_file, skip = 1, show_col_types = FALSE)
    colnames(qc_df)[length(colnames(qc_df))] <- "Value"
  } else if (datatype == "tercen") {
    qc_df <- read_delim(qc_file, show_col_types = FALSE)
    qc_df <- clean_tercen_columns(qc_df)
  }
  return(qc_df)
}

# Identifies which normalization approach(es) a cleaned tercen QC file provides, and the
# column holding each one's values. A file has at most one of logTransformed/identity
# (never both), and optionally a CmbCor column belonging to whichever base transform is
# present - so CmbCor is resolved to "Log + ComBat" or "VSN + ComBat" per-file, with no
# ambiguity and no cross-file lookup needed.
#
# Fallback: if none of those are present but the file has a generic "value" column
# instead (matched case-insensitively) - e.g. a BR file aggregated to the
# biological-replicate level, which drops the semantic column name - the normalization
# is instead read from a tag in the filename (see normalization_filename_tags),
# passed in via filename_hint.
identify_value_columns <- function(df, filename_hint = NA_character_) {
  cols <- colnames(df)
  norms <- list()

  if ("logTransformed" %in% cols) {
    norms[["Log"]] <- "logTransformed"
    if ("CmbCor" %in% cols) {
      norms[["Log + ComBat"]] <- "CmbCor"
    }
  } else if ("identity" %in% cols) {
    norms[["VSN"]] <- "identity"
    if ("CmbCor" %in% cols) {
      norms[["VSN + ComBat"]] <- "CmbCor"
    }
  } else {
    value_col <- cols[tolower(cols) == "value"]
    if (length(value_col) == 1 && !is.na(filename_hint)) {
      norms[[filename_hint]] <- value_col[1]
    }
  }

  return(norms)
}

# Extracts a normalization tag (Log/LogCmb/VSN/VSNCmb, case-insensitive) from a QC
# filename, as an exact underscore-separated token - same mechanism as the TR/BR
# suffix. Returns NA if no such tag is present (the normal case, when the file's own
# columns already self-describe the normalization).
extract_normalization_hint <- function(file_elements) {
  matches <- tolower(file_elements) %in% names(normalization_filename_tags)
  if (!any(matches)) {
    return(NA_character_)
  }
  normalization_filename_tags[[tolower(file_elements[matches][1])]]
}

# Classifies the remaining columns of a cleaned tercen QC file (after removing ID and the
# value column(s)) into sample columns (identify an individual replicate) and condition
# columns (identify an experimental condition). This matters because the SD calculation
# must group rows by condition, using all sample-replicate rows within each
# condition x peptide group.
classify_qc_columns <- function(df, value_cols) {
  cols <- colnames(df)
  remaining <- setdiff(cols, c("ID", value_cols))

  if (all(c("Barcode", "Row") %in% remaining)) {
    sample_cols <- c("Barcode", "Row")
  } else {
    sample_candidates <- remaining[grepl("sample", remaining, ignore.case = TRUE) | remaining == "Biol_Rep"]
    if (length(sample_candidates) == 1) {
      sample_cols <- sample_candidates
    } else {
      stop(
        "classify_qc_columns: could not identify sample column(s) - expected both ",
        "'Barcode' and 'Row', or exactly one column matching 'sample' (case-insensitive) ",
        "or named exactly 'Biol_Rep'. Remaining columns: ", paste(remaining, collapse = ", ")
      )
    }
  }

  condition_cols <- setdiff(remaining, sample_cols)

  list(peptide_col = "ID", sample_cols = sample_cols, condition_cols = condition_cols)
}

#########################################################
############# CALCULATE VARIABILITY (SD) ################
#########################################################

# For a given normalization value column, computes the peptide signal SD within each
# condition, then takes the median across all peptides - one number per condition.
compute_condition_sd <- function(df, condition_cols, peptide_col, value_col) {
  df %>%
    group_by(across(all_of(c(condition_cols, peptide_col)))) %>%
    summarise(peptide_sd = sd(.data[[value_col]], na.rm = TRUE), .groups = "drop") %>%
    group_by(across(all_of(condition_cols))) %>%
    summarise(median_sd = median(peptide_sd, na.rm = TRUE), .groups = "drop") %>%
    pull(median_sd)
}

stringify_sd <- function(sd_vec) {
  min_sd <- round(min(sd_vec, na.rm = TRUE), 2)
  max_sd <- round(max(sd_vec, na.rm = TRUE), 2)
  paste0(min_sd, " - ", max_sd)
}

# Table 2's Technical Variability cell: a single median value across conditions, not the
# min-max range stringify_sd() produces.
stringify_median_sd <- function(sd_vec) {
  as.character(round(median(sd_vec, na.rm = TRUE), 2))
}

classify_variability_tier <- function(value, replicate_type) {
  thresholds <- variability_criteria[[replicate_type]]
  if (value > thresholds[1]) {
    return("Poor")
  } else if (value > thresholds[2]) {
    return("Fair")
  } else {
    return("Good")
  }
}

format_pep_string <- function(assay_type, num_peptides) {
  if (assay_type == "PTK") {
    paste(num_peptides, "/ 195")
  } else if (assay_type == "STK") {
    paste(num_peptides, "/ 142")
  } else {
    as.character(num_peptides)
  }
}

#########################################################
######################## QC FLAG ########################
#########################################################

flag_low_is_bad <- function(value, thresholds) {
  if (value < thresholds[1]) {
    1
  } else if (value < thresholds[2]) {
    2
  } else {
    3
  }
}

flag_high_is_bad <- function(value, thresholds) {
  if (value > thresholds[1]) {
    1
  } else if (value > thresholds[2]) {
    2
  } else {
    3
  }
}

# Combines 2 or 3 individual criterion flags (each 1=poor/2=fair/3=good) into one overall
# flag, per the fully explicit rule (equivalent to the original sum-based formula for the
# 3-criteria case, restated without requiring arithmetic to read it):
#   2 criteria: good iff both good; poor iff both poor; else fair.
#   3 criteria: good iff all three good; poor iff no criterion is good AND at least two
#               criteria are poor; else fair.
combine_flags <- function(flags) {
  n <- length(flags)
  n_good <- sum(flags == 3)
  n_poor <- sum(flags == 1)

  if (n == 2) {
    if (n_good == 2) {
      return(3)
    } else if (n_poor == 2) {
      return(1)
    } else {
      return(2)
    }
  } else if (n == 3) {
    if (n_good == 3) {
      return(3)
    } else if (n_good == 0 && n_poor >= 2) {
      return(1)
    } else {
      return(2)
    }
  } else {
    stop("combine_flags: expected 2 or 3 criteria, got ", n)
  }
}

# Table 2's BR-fallback row (no TR file for this assay type - see parse_qc()): a
# 2-criteria rule distinct from combine_flags()'s (bionav's) 2-criteria rule - poor
# whenever NEITHER criterion is good (not just when both are poor), since there's no
# 3rd (variability) criterion here to soften a single weak signal.
combine_flags_br <- function(flags) {
  n_good <- sum(flags == 3)
  if (n_good == 2) {
    3
  } else if (n_good == 0) {
    1
  } else {
    2
  }
}

# values must contain Signal, num_peptides, and optionally Variability (present iff this
# is a technical-replicate row and variability is a flagging criterion for it).
determine_flag <- function(values, assay_type) {
  criteria <- qc_criteria[[assay_type]]

  signal_flag <- flag_low_is_bad(values$Signal, criteria$Signal)
  pep_flag <- flag_low_is_bad(values$num_peptides, criteria$QC_peptides)

  if (is.null(values$Variability) || is.na(values$Variability)) {
    overall_flag <- combine_flags(c(signal_flag, pep_flag))
    return(tibble(
      "Signal_Flag" = signal_flag, "Pep_Flag" = pep_flag,
      "Variability_Flag" = NA_integer_, "Overall_Flag" = overall_flag
    ))
  }

  # Variability is only ever passed for a TR row (a BR row's variability is informational,
  # not a flag - see the early-return above), so the TR band always applies here.
  variability_flag <- flag_high_is_bad(values$Variability, variability_criteria$TR)
  overall_flag <- combine_flags(c(signal_flag, pep_flag, variability_flag))
  tibble(
    "Signal_Flag" = signal_flag, "Pep_Flag" = pep_flag,
    "Variability_Flag" = variability_flag, "Overall_Flag" = overall_flag
  )
}

#########################################################
####################### QC TABLE ########################
#########################################################

# Parses all QC files into: (1) qc_table - one row per assay type that has any QC data
# at all (TR-preferred, BR-fallback when no TR file exists - see below), using only the
# "latest" normalization approach present (see normalization_order), each with its own
# QC flag; (2) variability_table - one row per assay-type x normalization-approach
# combination, for every assay type that has a BR file, used for the Data Variability
# Indicator boxplot (BR-only - an assay type with no BR file gets no row here at all).
# datatype == "bionav" files only ever have a single value column, so they're always
# treated as a simple 2-criteria BR row with no variability breakdown.
parse_qc <- function(qc_files, datatype) {
  qc_rows <- list()
  variability_rows <- list()

  for (assay_type in levels(qc_files$Assay_Type)) {
    assay_files <- qc_files %>% filter(Assay_Type == assay_type)

    if (datatype == "bionav") {
      qc_df <- load_qc_file(assay_files$qc_file[1], datatype)
      signal <- as.vector(round(quantile(2 ^ qc_df$Value, .99, na.rm = TRUE)))
      num_peptides <- qc_df %>% distinct(ID) %>% nrow()
      pep_string <- format_pep_string(assay_type, num_peptides)

      flags <- determine_flag(list(Signal = signal, num_peptides = num_peptides), assay_type)
      qc_rows[[length(qc_rows) + 1]] <- bind_cols(
        tibble(
          "Assay_Type" = assay_type, "Assay_Label" = assay_type, "Normalization" = NA_character_,
          "Signal" = signal, "num_peptides" = num_peptides, "pep_string" = pep_string,
          "Variability_String" = "—"
        ),
        flags
      )
      next
    }

    # datatype == "tercen": gather every normalization approach across TR files and BR
    # files for this assay type *separately* - a single assay type can legitimately have
    # both (a raw per-array TR export, plus a biological-replicate-level export averaged
    # up from it - see ref/HowToUse.md).
    #
    # Table 2 (flagging) is TR-only: a QC flag requires real per-replicate noise, and a BR
    # file is (per the documented workflow) already averaged up from TR rows, which washes
    # that noise out - it can't meaningfully support a flag. An assay type with no TR file
    # at all gets NO Table 2 row (it can still get a Data Variability Indicator entry, see
    # below - that's informational, not a flag, so BR data is fine there).
    #
    # The Data Variability Indicator (Table 3/4) prefers BR data when present (biological
    # variability is the informational number worth reporting), falling back to TR only
    # when no BR file exists for this assay at all (e.g. a TR-only study).
    has_hint_col <- "Normalization_Hint" %in% colnames(assay_files)

    gather_normalizations <- function(files_df) {
      normalizations <- list()
      for (i in seq_len(nrow(files_df))) {
        f <- files_df$qc_file[i]
        hint <- if (has_hint_col) files_df$Normalization_Hint[i] else NA_character_
        raw <- load_qc_file(f, datatype)
        file_norms <- identify_value_columns(raw, hint)
        if (length(file_norms) == 0) {
          next
        }
        value_cols <- unname(unlist(file_norms))
        classification <- classify_qc_columns(raw, value_cols)
        for (label in names(file_norms)) {
          if (!label %in% names(normalizations)) {
            normalizations[[label]] <- list(
              df = raw, column = file_norms[[label]], classification = classification
            )
          }
        }
      }
      normalizations
    }

    per_norm_sd_for <- function(normalizations, present_norms) {
      per_norm_sd <- list()
      for (label in present_norms) {
        entry <- normalizations[[label]]
        per_norm_sd[[label]] <- compute_condition_sd(
          entry$df, entry$classification$condition_cols, entry$classification$peptide_col, entry$column
        )
      }
      per_norm_sd
    }

    tr_normalizations <- gather_normalizations(assay_files %>% filter(Replicate_Type == "TR"))
    br_normalizations <- gather_normalizations(assay_files %>% filter(Replicate_Type == "BR"))

    if (length(tr_normalizations) == 0 && length(br_normalizations) == 0) {
      next
    }

    # --- Table 2 / Table 1 (flagging): TR-preferred. An assay type with a TR file gets a
    # 3-criteria flag from it (signal, peptides, technical variability). An assay type
    # with no TR file at all still gets a Table 2 row, but from its BR data and only a
    # 2-criteria flag (signal, peptides) - biological variability isn't folded into this
    # flag; it's reported separately, informationally, via the Data Variability Indicator
    # below.
    if (length(tr_normalizations) > 0) {
      present_norms <- intersect(normalization_order, names(tr_normalizations))

      # Signal strength and QC-passed phosphosites, from the "primary" TR normalization
      # (Log preferred, since 2^Value assumes a log2 scale; falls back to VSN if no
      # Log-based file is present).
      primary_label <- if ("Log" %in% present_norms) "Log" else present_norms[1]
      primary_df <- tr_normalizations[[primary_label]]$df
      primary_col <- tr_normalizations[[primary_label]]$column

      signal <- as.vector(round(quantile(2 ^ primary_df[[primary_col]], .99, na.rm = TRUE)))
      num_peptides <- primary_df %>% distinct(ID) %>% nrow()
      pep_string <- format_pep_string(assay_type, num_peptides)

      per_norm_sd <- per_norm_sd_for(tr_normalizations, present_norms)

      # Only the "latest" normalization approach appears in Table 2 (same convention used
      # elsewhere - the Main Report's per-assay flag) - not one row per normalization.
      label <- present_norms[length(present_norms)]
      # Both the flag and the displayed string are driven by the median across conditions -
      # a single representative number, not the worst-case max or a min-max range.
      variability_value <- median(per_norm_sd[[label]], na.rm = TRUE)
      flags <- determine_flag(
        list(Signal = signal, num_peptides = num_peptides, Variability = variability_value),
        assay_type
      )
      qc_rows[[length(qc_rows) + 1]] <- bind_cols(
        tibble(
          "Assay_Type" = assay_type, "Assay_Label" = paste0(assay_type, " (", label, ")"),
          "Normalization" = label, "Signal" = signal, "num_peptides" = num_peptides,
          "pep_string" = pep_string, "Variability_String" = stringify_median_sd(per_norm_sd[[label]])
        ),
        flags
      )
    } else if (length(br_normalizations) > 0) {
      present_norms <- intersect(normalization_order, names(br_normalizations))

      primary_label <- if ("Log" %in% present_norms) "Log" else present_norms[1]
      primary_df <- br_normalizations[[primary_label]]$df
      primary_col <- br_normalizations[[primary_label]]$column

      signal <- as.vector(round(quantile(2 ^ primary_df[[primary_col]], .99, na.rm = TRUE)))
      num_peptides <- primary_df %>% distinct(ID) %>% nrow()
      pep_string <- format_pep_string(assay_type, num_peptides)
      label <- present_norms[length(present_norms)]

      criteria <- qc_criteria[[assay_type]]
      signal_flag <- flag_low_is_bad(signal, criteria$Signal)
      pep_flag <- flag_low_is_bad(num_peptides, criteria$QC_peptides)

      qc_rows[[length(qc_rows) + 1]] <- tibble(
        "Assay_Type" = assay_type, "Assay_Label" = paste0(assay_type, " (", label, ")"),
        "Normalization" = label, "Signal" = signal, "num_peptides" = num_peptides,
        "pep_string" = pep_string, "Variability_String" = "—",
        "Signal_Flag" = signal_flag, "Pep_Flag" = pep_flag,
        "Variability_Flag" = NA_integer_, "Overall_Flag" = combine_flags_br(c(signal_flag, pep_flag))
      )
    }

    # --- Data Variability Indicator (boxplot): BR-only, always - no TR fallback. An assay
    # type with no BR file gets no entry here at all (should_show_variability_section()
    # decides the section's overall visibility from qc_files directly).
    if (length(br_normalizations) > 0) {
      var_present_norms <- intersect(normalization_order, names(br_normalizations))
      var_per_norm_sd <- per_norm_sd_for(br_normalizations, var_present_norms)

      for (label in var_present_norms) {
        variability_rows[[length(variability_rows) + 1]] <- tibble(
          "Assay_Type" = assay_type, "Normalization" = label,
          "Variability_Value" = median(var_per_norm_sd[[label]], na.rm = TRUE)
        )
      }
    }
  }

  qc_table <- bind_rows(qc_rows)
  if (nrow(qc_table) > 0) {
    qc_table <- qc_table %>% mutate(Flag_img = case_when(
      Overall_Flag == 1 ~ "img/red_25.png",
      Overall_Flag == 2 ~ "img/orange_25.png",
      Overall_Flag == 3 ~ "img/green_25.png"
    ))
  }

  variability_table <- bind_rows(variability_rows)

  list(qc_table = qc_table, variability_table = variability_table)
}

render_qc_results_table <- function(qc_table) {
  bg_red <- "#FF7F81"
  bg_orange <- "#FFCE78"
  bg_green <- "#67AD67"

  has_variability <- any(!is.na(qc_table$Variability_Flag))

  col_keys <- c("Assay_Label", "Flag_img", "Signal", "pep_string")
  header_labels <- list(
    "Assay_Label" = "Assay Type", "Flag_img" = "Flag",
    "Signal" = "Signal Strength (AU)", "pep_string" = "QC passed phosphosites"
  )
  if (has_variability) {
    col_keys <- c(col_keys, "Variability_String")
    header_labels[["Variability_String"]] <- "Technical Variability"
  }

  ft <- qc_table %>%
    flextable(col_keys = col_keys) %>%
    colformat_image(j = "Flag_img", width = .1, height = .1) %>%
    bg(~Signal_Flag == 3, ~ Signal, bg = bg_green) %>%
    bg(~Signal_Flag == 2, ~ Signal, bg = bg_orange) %>%
    bg(~Signal_Flag == 1, ~ Signal, bg = bg_red) %>%
    bg(~Pep_Flag == 3, ~ pep_string, bg = bg_green) %>%
    bg(~Pep_Flag == 2, ~ pep_string, bg = bg_orange) %>%
    bg(~Pep_Flag == 1, ~ pep_string, bg = bg_red)

  if (has_variability) {
    ft <- ft %>%
      bg(~Variability_Flag %in% 3, ~ Variability_String, bg = bg_green) %>%
      bg(~Variability_Flag %in% 2, ~ Variability_String, bg = bg_orange) %>%
      bg(~Variability_Flag %in% 1, ~ Variability_String, bg = bg_red)
  }

  ft %>%
    set_header_labels(values = header_labels) %>%
    autofit() %>%
    theme_box() %>%
    fontsize(size = 10, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption("QC Results")
}

#########################################################
################## REFERENCE QC TABLE  ##################
#########################################################

render_ref_qc_table <- function(qc_files) {
  base <- read_delim("ref/ref_qc_table.txt", show_col_types = FALSE)

  header_labels <- list(
    "Level" = "Level", "Flag" = "Flag", "Signal" = "Signal Strength (AU)",
    "PTK_QC" = "QC-passed phosphosites (PTK) of 195", "STK_QC" = "QC-passed phosphosites (STK) of 142"
  )

  tr_assays <- character(0)
  if (nrow(qc_files) > 0 && "Replicate_Type" %in% colnames(qc_files)) {
    tr_assays <- qc_files %>%
      filter(Replicate_Type == "TR") %>%
      pull(Assay_Type) %>%
      unique() %>%
      as.character()
    tr_assays <- intersect(c("PTK", "STK"), tr_assays)
  }

  # The threshold values are identical for PTK and STK (variability_criteria is not
  # assay-specific) - when both are TR, show one combined column instead of two
  # identical ones.
  if (length(tr_assays) == 2) {
    thresholds <- variability_criteria$TR
    base[["Tech_Variability"]] <- c(
      paste0("< ", thresholds[2]), paste0(thresholds[2], " - ", thresholds[1]), paste0("> ", thresholds[1])
    )
    header_labels[["Tech_Variability"]] <- "Technical Variability (PTK & STK)"
  } else {
    for (assay in tr_assays) {
      thresholds <- variability_criteria$TR
      colname <- paste0("Tech_Variability_", assay)
      base[[colname]] <- c(
        paste0("< ", thresholds[2]), paste0(thresholds[2], " - ", thresholds[1]), paste0("> ", thresholds[1])
      )
      header_labels[[colname]] <- paste0("Technical Variability (", assay, ")")
    }
  }

  base %>%
    flextable() %>%
    colformat_image(j = "Flag", width = .1, height = .1) %>%
    set_header_labels(values = header_labels) %>%
    footnote(
      j = "Flag", part = "header",
      value = as_paragraph("The QC flag is based on the criteria shown. All criteria must be met to get a green QC flag.")
    ) %>%
    autofit() %>%
    theme_box() %>%
    fontsize(size = 10, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption("Flag system used to assess data quality")
}

#########################################################
############### DATA VARIABILITY INDICATOR ##############
#########################################################

# The Data Variability Indicator is BR-only, always (see parse_qc()) - shown iff at least
# one assay type has a BR QC file, regardless of normalization count.
should_show_variability_section <- function(qc_files) {
  if (is.null(qc_files) || nrow(qc_files) == 0) {
    return(FALSE)
  }
  any(qc_files$Replicate_Type == "BR")
}

# Normalization display labels/colors for the boxplot's current-experiment reference
# lines. Fixed color per normalization (never shifts with how many are present): brown /
# orange for Log, blue / lime green for VSN - chosen so the base and ComBat variant of
# each normalization contrast strongly with each other.
variability_norm_labels <- c(
  "Log" = "Log2", "Log + ComBat" = "Log2-Cmb", "VSN" = "VSN", "VSN + ComBat" = "VSN-Cmb"
)
variability_line_colors <- c(
  "Log2" = "#7A2E1D", "Log2-Cmb" = "#FF9F1C", "VSN" = "#4B60D2", "VSN-Cmb" = "#A3D24B"
)

# Renders the Data Variability Indicator: one boxplot panel per assay type present in
# variability_table (i.e. present as BR in the current study - see parse_qc()), showing
# the historical per-condition median-SD distribution for that assay type (from the
# one-time precomputed data/qc_variability_distribution.rds - see
# R/calculate_variability_distribution.R), with one vertical reference line per
# normalization the current experiment has, at its own median-across-conditions value.
render_variability_boxplot <- function(variability_table) {
  reference_data <- readRDS("data/qc_variability_distribution.rds")

  current_lines <- variability_table %>%
    mutate(Label = variability_norm_labels[Normalization])

  plot_data <- reference_data %>% filter(Assay_Type %in% unique(current_lines$Assay_Type))

  ggplot(plot_data, aes(x = Median_SD, y = "")) +
    geom_boxplot(width = 0.4, outlier.alpha = 0.4) +
    geom_vline(data = current_lines, aes(xintercept = Variability_Value, color = Label), linewidth = 1) +
    facet_wrap(~Assay_Type, ncol = 1) +
    scale_color_manual(values = variability_line_colors, name = "Current experiment\n(median)") +
    labs(
      title = "Data Variability Indicator",
      subtitle = paste0(
        "Distribution of variability across a variety of PamDx experiments.\n",
        "Lines mark the current experiment's median variability of all conditions."
      ),
      x = "Variability (median peptide SD per condition)", y = NULL
    ) +
    report_theme +
    theme(
      strip.placement = "outside", strip.background = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
}

#########################################################
############# MAIN REPORT QC SUMMARY TEXT ################
#########################################################

# Builds the Main Report's one-line QC summary: one flag per assay type (for TR assay
# types with multiple normalization rows, uses the "latest" one present, per
# normalization_order).
build_qc_result_text <- function(qc_table) {
  flag_string_for <- function(overall_flag) {
    case_when(
      overall_flag == 1 ~ "\"poor\", ![](img/red_10.png)",
      overall_flag == 2 ~ "\"fair\", ![](img/orange_10.png)",
      overall_flag == 3 ~ "\"good\", ![](img/green_10.png)"
    )
  }

  # qc_table is TR-only (see parse_qc()) - an all-BR study (no TR file anywhere) has no
  # rows here at all, i.e. no QC flag sentence at all.
  qc_result <- ""
  if (nrow(qc_table) > 0) {
    latest_per_assay <- qc_table %>%
      mutate(.norm_rank = match(Normalization, normalization_order)) %>%
      group_by(Assay_Type) %>%
      arrange(.norm_rank, .by_group = TRUE) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      mutate(flag_string = flag_string_for(Overall_Flag))

    assay_types <- unique(qc_table$Assay_Type)
    if (length(assay_types) == 2) {
      ptk_result <- latest_per_assay %>% filter(Assay_Type == "PTK") %>% pull(flag_string)
      stk_result <- latest_per_assay %>% filter(Assay_Type == "STK") %>% pull(flag_string)
      qc_result <- paste0(
        "In this study, the PTK QC flag was ", ptk_result, " and the STK QC flag was ", stk_result, "."
      )
    } else {
      assay_type <- latest_per_assay %>% pull(Assay_Type)
      flag_result <- latest_per_assay %>% pull(flag_string)
      qc_result <- paste0("In this study, the ", assay_type, " QC flag was ", flag_result, ".")
    }
  }

  qc_result
}

#########################################################
############ AUTO-DETECT NORMALIZATIONS (QC) ############
#########################################################

# Looks at every uploaded QC file (tercen only - bionav QC files carry no normalization
# metadata) and returns which normalization approach was used, in the same c("vsn",
# "combat") shape as the manual "Normalizations" UI param, so callers can use it as a
# drop-in replacement. When multiple normalization approaches are present across the QC
# files (e.g. both a Log and a VSN file were uploaded), the "latest" one is used - per
# normalization_order (1 log, 2 log+ComBat, 3 VSN, 4 VSN+ComBat) - the same convention
# already used elsewhere (build_qc_result_text(), Table 4) for picking a single
# representative normalization out of several available ones.
# Returns NULL when detection isn't possible (bionav, no QC files, or no recognized
# normalization columns in any of them) so the caller can fall back to the manual param.
detect_normalizations_from_qc <- function(qc_files, datatype) {
  if (datatype != "tercen" || is.null(qc_files) || nrow(qc_files) == 0) {
    return(NULL)
  }

  has_hint_col <- "Normalization_Hint" %in% colnames(qc_files)
  present_norms <- character(0)
  for (i in seq_len(nrow(qc_files))) {
    f <- qc_files$qc_file[i]
    hint <- if (has_hint_col) qc_files$Normalization_Hint[i] else NA_character_
    raw <- tryCatch(load_qc_file(f, datatype), error = function(e) NULL)
    if (is.null(raw)) {
      next
    }
    present_norms <- union(present_norms, names(identify_value_columns(raw, hint)))
  }

  present_norms <- intersect(normalization_order, present_norms)
  if (length(present_norms) == 0) {
    return(NULL)
  }

  latest <- present_norms[length(present_norms)]

  result <- character(0)
  if (latest %in% c("VSN", "VSN + ComBat")) {
    result <- c(result, "vsn")
  }
  if (latest %in% c("Log + ComBat", "VSN + ComBat")) {
    result <- c(result, "combat")
  }
  result
}

# The effective normalization state to describe in the report: auto-detected from the QC
# files when possible (tercen), otherwise the manually-set "Normalizations" UI param
# (bionav, or no usable QC files yet).
get_effective_normalizations <- function(qc_files, datatype, manual_normalizations) {
  detected <- detect_normalizations_from_qc(qc_files, datatype)
  if (!is.null(detected)) {
    return(detected)
  }
  manual_normalizations
}
