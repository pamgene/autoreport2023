library(tidyverse)
library(flextable)


extract_cv_values <- function(df, group_name) {
  idx <- which(colnames(df) %in% group_name)

  cv_vals <- df %>%
    group_by(across(idx), ID) %>%
    summarise(cv = sd(Value) / mean(Value)) %>%
    ungroup() %>%
    group_by(across(1)) %>%
    summarise(mean_cv = mean(cv)) %>%
    pull()

  cv_string <- stringify_cv(cv_vals)
  max_cv <- max(cv_vals, na.rm = TRUE)

  output <- list(cv_string = cv_string, max = max_cv)

  return(output)
}

stringify_cv <- function(cv_vec) {
  max_cv <- max(cv_vec, na.rm = TRUE)
  min_cv <- min(cv_vec, na.rm = TRUE)

  max_cv <- round(max_cv * 100)
  min_cv <- round(min_cv * 100)

  cv_string <- paste0(min_cv, "-", max_cv, "%")
  return(cv_string)
}

determine_flag <- function(values, assay_type) {
  ptk_criteria <- list(
    "Signal" = c(1000, 2000, 4000),
    "QC_peptides" = c(78, 123, 192),
    "bio_cv" = c(0.40, 0.30, 0)
  )

  stk_criteria <- list(
    "Signal" = c(1000, 2000, 4000),
    "QC_peptides" = c(56, 90, 144),
    "bio_cv" = c(0.40, 0.30, 0)
  )

  if (assay_type == "PTK") {
    criteria <- ptk_criteria
  } else if (assay_type == "STK") {
    criteria <- stk_criteria
  }

  if (values$Signal < criteria$Signal[1]) {
    signal_flag <- 1
  } else if (values$Signal < criteria$Signal[2]) {
    signal_flag <- 2
  } else {
    signal_flag <- 3
  }

  if (values$num_peptides < criteria$QC_peptides[1]) {
    pep_flag <- 1
  } else if (values$num_peptides < criteria$QC_peptides[2]) {
    pep_flag <- 2
  } else {
    pep_flag <- 3
  }

  if (values$cv_max > criteria$bio_cv[1]) {
    cv_flag <- 1
  } else if (values$cv_max > criteria$bio_cv[2]) {
    cv_flag <- 2
  } else {
    cv_flag <- 3
  }

  overall_val <- sum(signal_flag, pep_flag, cv_flag)

  if (overall_val == 9) {
    overall_flag <- 3
  } else if (overall_val > 4) {
    overall_flag <- 2
  } else {
    overall_flag <- 1
  }

  flag_df <- tibble("Signal_Flag" = signal_flag, "Pep_Flag" = pep_flag, "CV_Flag" = cv_flag, "Overall_Flag" = overall_flag)
  return(flag_df)
}

extract_qc <- function(qc_file, assay_type, group_name, datatype = "bionav") {
  if (datatype == "bionav") {
    qc_df <- read_delim(qc_file, skip = 1, show_col_types = FALSE)
  } else if (datatype == "tercen") {
    qc_df <- read_delim(qc_file, show_col_types = FALSE)
    qc_df <- clean_tercen_columns(qc_df)
  }
  
  colnames(qc_df)[length(colnames(qc_df))] <- "Value"
  signal <- round(quantile(2^qc_df$Value, .99, na.rm = TRUE))
  num_peptides <- qc_df %>%
    distinct(ID) %>%
    summarise(n = n()) %>%
    pull()
  cv_vals <- extract_cv_values(qc_df, group_name)

  if (assay_type == "PTK") {
    pep_string <- paste(num_peptides, "/ 195")
  } else if (assay_type == "STK") {
    pep_string <- paste(num_peptides, "/ 142")
  }

  out_df <- data.frame("Assay_Type" = assay_type, "Signal" = as.vector(signal), "num_peptides" = num_peptides, "pep_string" = pep_string, "tech_cv" = "n.a.", "cv_max" = cv_vals$max, "cv_string" = cv_vals$cv_string)
  flags <- determine_flag(out_df, assay_type)
  out_df <- bind_cols(out_df, flags)
  return(out_df)
}


parse_qc <- function(qc_files, group_name, datatype) {
  dfs <- list()
  for (assay_type in levels(qc_files$Assay_Type)) {
    df <- qc_files %>%
      filter(Assay_Type == assay_type) %>%
      do(extract_qc(.$qc_file, .$Assay_Type, group_name, datatype))
    dfs <- append(dfs, list(df))
  }

  output <- bind_rows(dfs)
  output <- output %>% relocate(Overall_Flag, .after = Assay_Type)
  output <- output %>% mutate(Flag_img = case_when(
    Overall_Flag == 1 ~ "img/red_25.png",
    Overall_Flag == 2 ~ "img/orange_25.png",
    Overall_Flag == 3 ~ "img/green_25.png"
  ))
  return(output)
}

make_qc_table <- function(qc_files, group_name, datatype) {
  qc_table <- parse_qc(qc_files, group_name, datatype)

  bg_red <- "#FF7F81"
  bg_orange <- "#FFCE78"
  bg_green <- "#67AD67"
  qc_table %>%
    flextable(col_keys = c("Assay_Type", "Flag_img", "Signal", "pep_string", "tech_cv", "cv_string")) %>%
    colformat_image(j = "Flag_img", width = .1, height = .1) %>%
    bg(~ Signal_Flag == 3, ~Signal, bg = bg_green) %>%
    bg(~ Signal_Flag == 2, ~Signal, bg = bg_orange) %>%
    bg(~ Signal_Flag == 1, ~Signal, bg = bg_red) %>%
    bg(~ Pep_Flag == 3, ~pep_string, bg = bg_green) %>%
    bg(~ Pep_Flag == 2, ~pep_string, bg = bg_orange) %>%
    bg(~ Pep_Flag == 1, ~pep_string, bg = bg_red) %>%
    bg(~ CV_Flag == 3, ~cv_string, bg = bg_green) %>%
    bg(~ CV_Flag == 2, ~cv_string, bg = bg_orange) %>%
    bg(~ CV_Flag == 1, ~cv_string, bg = bg_red) %>%
    set_header_labels(Assay_Type = "Assay Type", Flag_img = "Flag", Signal = "Signal Strength (AU)", pep_string = "QC passed phosphosites", tech_cv = "Technical CV", cv_string = "Biological CV") %>%
    autofit() %>%
    theme_box() %>%
    fontsize(size = 10, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption("QC Results")
}