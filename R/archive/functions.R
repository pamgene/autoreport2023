source("R/load_data.R")
source("R/read_qc_dir.R")
source("R/make_ref_qc_table.R")
source("R/make_qc_table.R")
source("R/read_phosphosite_dir.R")
source("R/make_volcano_plots.R")
source("R/make_peptide_heatmaps.R")
source("R/read_kinase_dir.R")
source("R/make_kinase_tables.R")
source("R/make_kinase_plots.R")
source("R/output_kinase_analysis.R")

make_report_table <- function(res_table, caption) {
  # identify first which assay types are included
  assay_types <- str_extract(colnames(res_table), ".TK") %>%
    na.omit() %>%
    unique()
  header_row <- c("Assay Type", assay_types)
  if (length(assay_types) == 2) {
    header_widths <- c(1, 2, 2)
    ft <- flextable(res_table) %>%
      set_header_labels(comparison = "Comparisons", up_PTK = "Up", down_PTK = "Down", up_STK = "Up", down_STK = "Down") %>%
      color(j = c("up_PTK", "up_STK"), color = "red", part = "body") %>%
      color(j = c("down_PTK", "down_STK"), color = "blue", part = "body") %>%
      add_header_row(values = header_row, colwidths = header_widths)
  } else if (length(assay_types) == 1) {
    header_widths <- c(1, 2)

    if (assay_types[1] == "PTK") {
      ft <- flextable(res_table) %>%
        set_header_labels(comparison = "Comparisons", up_PTK = "Up", down_PTK = "Down") %>%
        color(j = c("up_PTK"), color = "red", part = "body") %>%
        color(j = c("down_PTK"), color = "blue", part = "body") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    } else if (assay_types[1] == "STK") {
      ft <- flextable(res_table) %>%
        set_header_labels(comparison = "Comparisons", up_STK = "Up", down_STK = "Down") %>%
        color(j = c("up_STK"), color = "red", part = "body") %>%
        color(j = c("down_STK"), color = "blue", part = "body") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    }
  }
  ft %>%
    autofit() %>%
    theme_box() %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    set_caption(caption)
}

stats_footnotes <- function(ft, res_table) {
  stats_types <- res_table %>%
    ungroup() %>%
    distinct(stats) %>%
    pull()
  if (length(stats_types) == 2) {
    ft <- ft %>%
      footnote(~ stats == "mtvc", ~comparison, value = as_paragraph("Significance was obtained using a one-way ANOVA followed by a post-hoc Dunnett\u0027s test, p<0.05"), ref_symbols = c("a"), inline = TRUE) %>%
      footnote(~ stats == "tt", ~comparison, value = as_paragraph("Significance was obtained using a two-sided unpaired Student\u0027s T-test, p<0.05."), ref_symbols = c("b"), inline = TRUE)
  } else if (stats_types == "mtvc") {
    ft <- ft %>%
      footnote(~ stats == "mtvc", ~comparison, value = as_paragraph("Significance was obtained using a one-way ANOVA followed by a post-hoc Dunnett\u0027s test, p<0.05"), ref_symbols = c("a"), inline = TRUE)
  } else if (stats_types == "tt") {
    ft <- ft %>%
      footnote(~ stats == "tt", ~comparison, value = as_paragraph("Significance was obtained using a two-sided unpaired Student\u0027s T-test, p<0.05."), ref_symbols = c("a"), inline = TRUE)
  }
  return(ft)
}

make_report_table_stats <- function(res_table, caption) {
  # identify first which assay types are included
  assay_types <- str_extract(colnames(res_table), ".TK") %>%
    na.omit() %>%
    unique()
  header_row <- c("Assay Type", assay_types)
  if (length(assay_types) == 2) {
    header_widths <- c(1, 2, 2)
    ft <- flextable(res_table, col_keys = c("comparison", "up_PTK", "down_PTK", "up_STK", "down_STK")) %>%
      set_header_labels(comparison = "Comparisons", up_PTK = "Up", down_PTK = "Down", up_STK = "Up", down_STK = "Down") %>%
      color(j = c("up_PTK", "up_STK"), color = "red", part = "body") %>%
      color(j = c("down_PTK", "down_STK"), color = "blue", part = "body") %>%
      add_header_row(values = header_row, colwidths = header_widths)
  } else if (length(assay_types) == 1) {
    header_widths <- c(1, 2)

    if (assay_types[1] == "PTK") {
      ft <- flextable(res_table, col_keys = c("comparison", "up_PTK", "down_PTK")) %>%
        set_header_labels(comparison = "Comparisons", up_PTK = "Up", down_PTK = "Down") %>%
        color(j = c("up_PTK"), color = "red", part = "body") %>%
        color(j = c("down_PTK"), color = "blue", part = "body") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    } else if (assay_types[1] == "STK") {
      ft <- flextable(res_table, col_keys = c("comparison", "up_STK", "down_STK")) %>%
        set_header_labels(comparison = "Comparisons", up_STK = "Up", down_STK = "Down") %>%
        color(j = c("up_STK"), color = "red", part = "body") %>%
        color(j = c("down_STK"), color = "blue", part = "body") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    }
  }
  ft %>%
    stats_footnotes(., res_table = res_table) %>%
    fontsize(size = 9, part = "footer") %>%
    font(part = "all", fontname = "Arial") %>%
    bold(bold = TRUE, part = "header") %>%
    bold(bold = TRUE, j = "comparison") %>%
    theme_box() %>%
    align(align = "center", part = "all") %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    set_caption(caption)
}

clean_tercen_columns <- function(df) {
  cols <- colnames(df)
  split <- str_split(cols, pattern = "\\.")
  cols <- sapply(split, tail, 1)
  colnames(df) <- cols
  return(df)
}
