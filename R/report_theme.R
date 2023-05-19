library(tidyverse)
# General plotting theme
report_theme <-
  theme_classic() +
  theme(
    line = element_line(linewidth = 1),
    text = element_text(size = 10),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.key.size = unit(8, "pt"),
    legend.box.spacing = unit(0.1, "mm"),
    legend.text = element_text(size = 8),
    axis.line.x = element_line(linewidth = 0.25),
    axis.title.x = element_text(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.line.y = element_line(linewidth = 0.25),
    axis.ticks = element_line(linewidth = 0.25),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = NA, color = NA),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") # trbl (top right bottom left)
  )

library(flextable)
theme_kin_table <- function(x) {
  if (!inherits(x, "flextable")) {
    stop("theme_kin_table supports only flextable objects.")
  }

  fp_bdr <- fp_border(
    width = 1,
    color = "#000000"
  )

  x <- border_remove(x)
  x <- bg(x, bg = "transparent", part = "all")
  x <- color(x, color = "#000000", part = "all")
  x <- bold(x = x, bold = FALSE, part = "all")
  x <- italic(x = x, italic = FALSE, part = "all")
  x <- padding(x = x, padding = 3, part = "all")

  x <- align_text_col(x, align = "left", header = TRUE)
  x <- align_nottext_col(x, align = "right", header = TRUE)
  x <- hline_bottom(x, part = "header", border = fp_bdr)
  x <- hline_top(x, part = "body", border = fp_bdr)
  fix_border_issues(x)
}