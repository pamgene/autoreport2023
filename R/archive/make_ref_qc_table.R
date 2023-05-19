library(tidyverse)
library(flextable)

render_ref_qc_table <- function() {
  ref_qc_table <- read_delim("ref/ref_qc_table.txt", show_col_types = FALSE)

  ft <- ref_qc_table %>%
    flextable() %>%
    colformat_image(j = "Flag", width = .1, height = .1) %>%
    set_header_labels(Level = "Level", Flag = "Flag", Signal = "Signal Strength (AU)", PTK_QC = "QC-passed phosphosites (PTK) of 195", STK_QC = "QC-passed phosphosites (STK) of 142", Tech_cv = "Technical CV", Bio_cv = "Biological CV") %>%
    footnote(j = "Flag", part = "header", value = as_paragraph("The QC flag is based on the 3 individual criteria. All criteria must be met to get a green QC flag.")) %>%
    autofit() %>%
    theme_box() %>%
    fontsize(size = 10, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption("Flag system used to assess data quality")

  return(ft)
}