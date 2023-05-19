library(flextable)
library(officer)
library(tidyverse)

source("R/report_theme.R")

extract_kinase_table <- function(kinase_file, assay_type, kinase_number = 20) {
  df <- read_delim(kinase_file, show_col_types = FALSE)
  df <- df %>%
    select(`Kinase Name`, `Median Final score`, `Median Kinase Statistic`) %>%
    slice_head(n = 20) %>%
    mutate(
      `Median Final score` = round(`Median Final score`, digits = 2),
      `Median Kinase Statistic` = round(`Median Kinase Statistic`, digits = 2)
    )
  df <- rownames_to_column(df, var = "Rank")
  colnames(df) <- c("Rank", paste(assay_type, "Name", sep = "_"), paste(assay_type, "Score", sep = "_"), paste(assay_type, "Statistic", sep = "_"))

  return(df)
}

render_kinase_table <- function(res_table, assay_types, caption) {
  small_border <- fp_border(color = "gray", width = 1)
  header_row <- c("", assay_types)

  if (length(assay_types) == 2) {
    header_widths <- c(1, 3, 3)
    ft <- flextable(res_table) %>%
      set_header_labels(Rank = "Rank", PTK_Name = "Kinase Name", PTK_Score = "Kinase Score", PTK_Statistic = "Kinase Statistic", STK_Name = "Kinase Name", STK_Score = "Kinase Score", STK_Statistic = "Kinase Statistic") %>%
      add_header_row(values = header_row, colwidths = header_widths)
  } else if (length(assay_types) == 1) {
    header_widths <- c(1, 3)
    if (assay_types[1] == "PTK") {
      ft <- flextable(res_table) %>%
        set_header_labels(Rank = "Rank", PTK_Name = "Kinase Name", PTK_Score = "Kinase Score", PTK_Statistic = "Kinase Statistic") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    } else if (assay_types[1] == "STK") {
      ft <- flextable(res_table) %>%
        set_header_labels(Rank = "Rank", STK_Name = "Kinase Name", STK_Score = "Kinase Score", STK_Statistic = "Kinase Statistic") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    }
  }
  ft <- ft %>%
    autofit() %>%
    theme_kin_table()

  if (length(assay_types) == 2) {
    ft <- ft %>%
      vline(j = "PTK_Statistic", border = small_border)
  }
  ft <- ft %>%
    vline(j = "Rank", border = small_border) %>%
    fontsize(size = 7, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption(caption)
  return(ft)
}

make_kinase_tables <- function(kinase_files, kinase_number = 20) {
  # determine first whether study has both PTK and STK
  assay_types <- kinase_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()
  # determine the number of unique comparisons
  comparisons <- kinase_files %>%
    distinct(Comparison) %>%
    pull()
  # then iterate over comparisons, make table for each
  fts <- list()
  for (comparison in comparisons) {
    c_rows <- kinase_files %>% filter(Comparison == comparison)
    if (length(assay_types) == 2) {
      c_ptk <- c_rows %>%
        filter(Assay_Type == "PTK") %>%
        do(extract_kinase_table(.$UKA_file, .$Assay_Type, kinase_number))
      c_stk <- c_rows %>%
        filter(Assay_Type == "STK") %>%
        do(extract_kinase_table(.$UKA_file, .$Assay_Type, kinase_number))
      df <- left_join(c_ptk, c_stk, by = "Rank")
    } else if (length(assay_types) == 1) {
      df <- c_rows %>% do(extract_kinase_table(.$UKA_file, .$Assay_Type, kinase_number))
    }
    ft <- render_kinase_table(df, assay_types, paste("Kinase Score Table", comparison))
    fts <- append(fts, list(ft))
  }
  return(fts)
}