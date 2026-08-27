require(tidyverse)
require(flextable)

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
    sample_candidates <- remaining[grepl("sample", remaining, ignore.case = TRUE)]
    if (length(sample_candidates) == 1) {
      sample_cols <- sample_candidates
    } else {
      stop(
        "classify_qc_columns: could not identify sample column(s) - expected both ",
        "'Barcode' and 'Row', or exactly one column matching 'sample' (case-insensitive). ",
        "Remaining columns: ", paste(remaining, collapse = ", ")
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

# Parses all QC files into: (1) qc_table - one row per assay type (BR: single row;
# TR: single row too, using only the "latest" normalization approach present - see
# normalization_order), each with its own QC flag; (2) variability_table - one row per
# assay-type x normalization-approach combination, for every assay type (BR or TR) that
# has normalization data, used for the Data Variability Indicator (Table 3/4). datatype
# == "bionav" files only ever have a single value column, so they're always treated as a
# simple 2-criteria BR row with no variability breakdown.
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

    # datatype == "tercen": gather every normalization approach across all files for
    # this assay type.
    replicate_type <- unique(assay_files$Replicate_Type)[1]
    normalizations <- list()

    has_hint_col <- "Normalization_Hint" %in% colnames(assay_files)
    for (i in seq_len(nrow(assay_files))) {
      f <- assay_files$qc_file[i]
      hint <- if (has_hint_col) assay_files$Normalization_Hint[i] else NA_character_
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

    if (length(normalizations) == 0) {
      next
    }

    present_norms <- intersect(normalization_order, names(normalizations))

    # Signal strength and QC-passed phosphosites are computed once per assay type, from
    # the "primary" normalization (Log preferred, since 2^Value assumes a log2 scale;
    # falls back to VSN if no Log-based file is present).
    primary_label <- if ("Log" %in% present_norms) "Log" else present_norms[1]
    primary_df <- normalizations[[primary_label]]$df
    primary_col <- normalizations[[primary_label]]$column

    signal <- as.vector(round(quantile(2 ^ primary_df[[primary_col]], .99, na.rm = TRUE)))
    num_peptides <- primary_df %>% distinct(ID) %>% nrow()
    pep_string <- format_pep_string(assay_type, num_peptides)

    # Variability, per normalization approach present for this assay type.
    per_norm_sd <- list()
    for (label in present_norms) {
      entry <- normalizations[[label]]
      sd_vals <- compute_condition_sd(
        entry$df, entry$classification$condition_cols, entry$classification$peptide_col, entry$column
      )
      per_norm_sd[[label]] <- sd_vals
    }

    # Full per-normalization detail (Table 4), for every assay type regardless of
    # replicate type - visibility of the whole Data Variability Indicator section is
    # decided separately (see the Rmd: shown if any assay is BR, or any assay - BR or TR
    # - has more than one normalization present, since Table 2 already fully covers a
    # single-normalization assay on its own).
    for (label in present_norms) {
      max_sd <- max(per_norm_sd[[label]], na.rm = TRUE)
      variability_rows[[length(variability_rows) + 1]] <- tibble(
        "Assay_Type" = assay_type, "Replicate_Type" = replicate_type, "Normalization" = label,
        "Variability_String" = stringify_sd(per_norm_sd[[label]]),
        "Variability_Value" = max_sd,
        "Variability_Tier" = variability_tier_labels[[classify_variability_tier(max_sd, replicate_type)]]
      )
    }

    if (replicate_type == "BR") {
      flags <- determine_flag(list(Signal = signal, num_peptides = num_peptides), assay_type)
      qc_rows[[length(qc_rows) + 1]] <- bind_cols(
        tibble(
          "Assay_Type" = assay_type, "Assay_Label" = assay_type, "Normalization" = NA_character_,
          "Signal" = signal, "num_peptides" = num_peptides, "pep_string" = pep_string,
          "Variability_String" = "—"
        ),
        flags
      )
    } else {
      # TR: only the "latest" normalization approach appears in Table 2 (same convention
      # used elsewhere - the Main Report's per-assay flag, BR's variability sentence) -
      # not one row per normalization. The full per-normalization breakdown lives in
      # Table 4 instead (populated above) when there's more than one to show.
      label <- present_norms[length(present_norms)]
      variability_value <- max(per_norm_sd[[label]], na.rm = TRUE)
      flags <- determine_flag(
        list(Signal = signal, num_peptides = num_peptides, Variability = variability_value),
        assay_type
      )
      qc_rows[[length(qc_rows) + 1]] <- bind_cols(
        tibble(
          "Assay_Type" = assay_type, "Assay_Label" = paste0(assay_type, " (", label, ")"),
          "Normalization" = label, "Signal" = signal, "num_peptides" = num_peptides,
          "pep_string" = pep_string, "Variability_String" = as.character(round(variability_value, 2))
        ),
        flags
      )
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
    header_labels[["Variability_String"]] <- "Tech Variability"
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

  for (assay in tr_assays) {
    thresholds <- variability_criteria$TR
    colname <- paste0("Tech_Variability_", assay)
    base[[colname]] <- c(
      paste0("< ", thresholds[2]), paste0(thresholds[2], " - ", thresholds[1]), paste0("> ", thresholds[1])
    )
    header_labels[[colname]] <- paste0("Technical Variability (", assay, ")")
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

# Table 3/4 (Data Variability Indicator) are shown if any assay is BR (always
# informational there), or any assay - BR or TR - has more than one normalization
# approach present. Table 2 only ever shows the "latest" normalization for a TR assay
# type, so anything beyond that would otherwise be invisible without Table 4.
should_show_variability_section <- function(qc_files, variability_table) {
  if (is.null(qc_files) || nrow(qc_files) == 0) {
    return(FALSE)
  }
  any(qc_files$Replicate_Type == "BR") || any(table(variability_table$Assay_Type) > 1)
}

# Renders the tier legend for whichever replicate type(s) actually contribute rows to the
# Data Variability Indicator (BR always informational; TR only when it has more than one
# normalization present - see should_show_variability_section()). BR and TR reuse
# different bands (see variability_criteria), so a study mixing both gets two columns.
render_variability_ref_table <- function(variability_table) {
  band_col <- function(thresholds) {
    c(paste0("< ", thresholds[2]), paste0(thresholds[2], " - ", thresholds[1]), paste0("> ", thresholds[1]))
  }

  replicate_types <- intersect(c("BR", "TR"), unique(variability_table$Replicate_Type))
  if (length(replicate_types) == 0) {
    replicate_types <- "BR" # defensive fallback, shouldn't happen if the section is shown at all
  }

  ref <- tibble("Level" = c("Low", "Medium", "High"))
  header_labels <- list("Level" = "Level")
  for (rt in replicate_types) {
    colname <- if (length(replicate_types) == 2) paste0("Variability_", rt) else "Variability"
    ref[[colname]] <- band_col(variability_criteria[[rt]])
    header_labels[[colname]] <- if (length(replicate_types) == 2) {
      paste0(if (rt == "BR") "Biological" else "Technical", " Variability (SD)")
    } else {
      "Variability (SD)"
    }
  }

  ref %>%
    flextable() %>%
    set_header_labels(values = header_labels) %>%
    footnote(
      j = "Level", part = "header",
      value = as_paragraph("Informational only — not part of the QC flag.")
    ) %>%
    autofit() %>%
    theme_box() %>%
    fontsize(size = 10, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption("Data Variability Indicator tiers")
}

render_variability_table <- function(variability_table) {
  wide <- variability_table %>%
    mutate(Normalization = factor(Normalization, levels = normalization_order)) %>%
    arrange(Normalization) %>%
    select(Assay_Type, Normalization, Variability_String) %>%
    pivot_wider(names_from = Assay_Type, values_from = Variability_String, names_prefix = "Variability_")

  header_labels <- list("Normalization" = "Normalization")
  assay_cols <- setdiff(colnames(wide), "Normalization")
  for (col in assay_cols) {
    assay <- sub("^Variability_", "", col)
    header_labels[[col]] <- paste0("Variability (", assay, ")")
  }

  wide %>%
    flextable() %>%
    set_header_labels(values = header_labels) %>%
    autofit() %>%
    theme_box() %>%
    fontsize(size = 10, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption("Data Variability Indicator")
}

#########################################################
############# MAIN REPORT QC SUMMARY TEXT ################
#########################################################

# Builds the Main Report's one-line QC summary: one flag per assay type (for TR assay
# types with multiple normalization rows, uses the "latest" one present, per
# normalization_order), plus, for BR assay type(s) only, a "Biological variability was
# '<tier>'" clause using the latest normalization's tier for each.
build_qc_result_text <- function(qc_table, variability_table, qc_files) {
  flag_string_for <- function(overall_flag) {
    case_when(
      overall_flag == 1 ~ "\"poor\", ![](img/red_10.png)",
      overall_flag == 2 ~ "\"fair\", ![](img/orange_10.png)",
      overall_flag == 3 ~ "\"good\", ![](img/green_10.png)"
    )
  }

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

  br_assays <- qc_files %>%
    filter(Replicate_Type == "BR") %>%
    pull(Assay_Type) %>%
    unique() %>%
    as.character()

  if (length(br_assays) > 0 && nrow(variability_table) > 0) {
    br_variability <- variability_table %>%
      filter(Assay_Type %in% br_assays) %>%
      mutate(.norm_rank = match(Normalization, normalization_order)) %>%
      group_by(Assay_Type) %>%
      arrange(.norm_rank, .by_group = TRUE) %>%
      slice_tail(n = 1) %>%
      ungroup()

    if (nrow(br_variability) > 0) {
      n_assay_total <- length(assay_types)
      if (n_assay_total == 1) {
        variability_text <- paste0(" Biological variability was \"", tolower(br_variability$Variability_Tier[1]), "\".")
      } else {
        parts <- paste0(
          "\"", tolower(br_variability$Variability_Tier), "\" (", br_variability$Assay_Type, ")"
        )
        variability_text <- paste0(" Biological variability was ", paste(parts, collapse = " and "), ".")
      }
      qc_result <- paste0(qc_result, variability_text)
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
