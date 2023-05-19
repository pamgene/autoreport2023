source("R/make_kinase_plots.R")
source("R/make_kinase_tables.R")
source("R/make_coral_trees.R")

library(flextable)

output_kinase_analysis <- function(kinase_files, kin_params = params$`kinase_analysis`, color_type = params$`score_plot_color`, xax_scale = params$`xax_scale`) {
  temp <- list()
  if ("table" %in% kin_params) {
    tables <- make_kinase_tables(kinase_files)
    temp["table"] <- list(tables)
  }
  if ("plot" %in% kin_params) {
    temp["plot"] <- list(make_kinase_plots(kinase_files, kinase_number = 20, color_type = color_type, xax_scale = xax_scale))
  }

  if ("tree" %in% kin_params) {
    temp["tree"] <- list(make_coral_trees(kinase_files))
  }

  comparisons <- kinase_files %>%
    distinct(Comparison) %>%
    pull()
  for (i in seq_along(comparisons)) {
    for (vis in kin_params) {
      item <- temp[vis][[1]][[i]]
      if (vis == "table") {
        flextable_to_rmd(item)
      } else if (vis == "plot") {
        cat(paste0("**Kinase Score Plot ", comparisons[i], "** \n"))
        print(item)
      } else {
        cat(paste0("**Kinase Coral Tree ", comparisons[i], "** \n\n"))
        cat(paste0("![](", item, ")"))
      }
    }
    cat("\n \n")
    cat("#########")
    cat("\n \n")
  }
}