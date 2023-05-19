library(tidyverse)
library(cowplot)
library(EnhancedVolcano)
source("R/report_theme.R")

render_volcano_plot_tt <- function(df, lfc_range, p_range) {
  p <- EnhancedVolcano(df,
    lab = NA,
    x = "LogFC",
    y = "P",
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    pCutoff = 0.05,
    cutoffLineType = "blank",
    FCcutoff = 0,
    xlim = lfc_range,
    ylim = p_range,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    col = c("black", "black", "black", "red3"),
    legendPosition = "none"
  )

  p <- p +
    facet_wrap(~Assay_type, scales = "free", nrow = 1) +
    report_theme +
    theme(legend.position = "none")
  return(p)
}

render_volcano_plot_mtvc <- function(df, lfc_range, p_range) {
  p <- EnhancedVolcano(df,
    lab = NA,
    x = "LogFC",
    y = "P",
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    pCutoff = 0.05,
    cutoffLineType = "blank",
    FCcutoff = 0,
    xlim = lfc_range,
    ylim = p_range,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    col = c("black", "black", "black", "red3"),
    legendPosition = "none"
  )

  p <- p +
    facet_wrap(~Comparison, scales = "free", nrow = 1) +
    report_theme +
    theme(legend.position = "none")
  return(p)
}



make_volcano_plots_mtvc <- function(stats_files, datatype = "bionav") {
  # determine range of data
  stats_range <- get_stats_data_range(stats_files, datatype = datatype)

  # determine first whether study has both PTK and STK
  assay_types <- stats_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()

  plots <- list()
  if (length(assay_types) == 2) {
    mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
    groups <- stats_files %>%
      filter(Stats == "MTvC") %>%
      distinct(Group) %>%
      pull()

    if (!is.na(groups)) {
      for (group in groups) {
        c_rows <- mtvc_rows %>% filter(Group == group)
        if (datatype == "bionav") {
          mtvc_ptk <- c_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
          mtvc_stk <- c_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
        } else if (datatype == "tercen") {
          mtvc_ptk <- c_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
          mtvc_stk <- c_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
        }
       
        p_ptk <- render_volcano_plot_mtvc(mtvc_ptk, stats_range$lfc, stats_range$p)
        p_stk <- render_volcano_plot_mtvc(mtvc_stk, stats_range$lfc, stats_range$p)
        pg <- plot_grid(p_ptk, p_stk, ncol = 1, labels = c(paste("PTK"), paste("STK")))
        title <- ggdraw() + draw_label(paste(group), fontface = "bold", x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))
        pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
        save_plot(paste0("99_Saved Plots/", group, "_MTvC_Volcano.pdf"), pg, dpi = 300)
        plots[[length(plots) + 1]] <- pg
      }
    } else {
      if (datatype == "bionav") {
        mtvc_ptk <- mtvc_rows %>%
          filter(Assay_Type == "PTK") %>%
          do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
        mtvc_stk <- mtvc_rows %>%
          filter(Assay_Type == "STK") %>%
          do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      } else if (datatype == "tercen") {
        mtvc_ptk <- mtvc_rows %>%
          filter(Assay_Type == "PTK") %>%
          do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
        mtvc_stk <- mtvc_rows %>%
          filter(Assay_Type == "STK") %>%
          do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
      }
      
      p_ptk <- render_volcano_plot_mtvc(mtvc_ptk, stats_range$lfc, stats_range$p)
      p_stk <- render_volcano_plot_mtvc(mtvc_stk, stats_range$lfc, stats_range$p)
      pg <- plot_grid(p_ptk, p_stk, ncol = 1, labels = c(paste("PTK"), paste("STK")))
      save_plot(paste0("99_Saved Plots/", "MTvC_Volcano.pdf"), pg, dpi = 300)
      plots[[length(plots) + 1]] <- pg
    }
  } else if (length(assay_types) == 1) {
    mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
    groups <- stats_files %>%
      filter(Stats == "MTvC") %>%
      distinct(Group) %>%
      pull()
    if (!is.na(groups)) {
      for (group in groups) {
        c_rows <- mtvc_rows %>% filter(Group == group)
        if (datatype == "bionav") {
          df <- c_rows %>% do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
        } else if (datatype == "tercen") {
          df <- c_rows %>% do(extract_phosphosite_data_mtvc_tercen(.$File, .$P_file, .$Assay_Type))
        }
        p <- render_volcano_plot_mtvc(df, stats_range$lfc, stats_range$p)
        title <- ggdraw() + draw_label(paste(group), fontface = "bold", x = 0, hjust = 0) + theme(plot.margin = margin(0, 0, 0, 7))
        pg <- plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
        save_plot(paste0("99_Saved Plots/", group, "_MTvC_Volcano.pdf"), pg, dpi = 300)
        plots[[length(plots) + 1]] <- pg
      }
    } else {
      if (datatype == "bionav") {
        df <- mtvc_rows %>% do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      } else if (datatype == "tercen") {
        df <- mtvc_rows %>% do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
        
      }
      p <- render_volcano_plot_mtvc(df, stats_range$lfc, stats_range$p)
      save_plot(paste0("99_Saved Plots/", "MTvC_Volcano.pdf"), p, dpi = 300)
      plots[[length(plots) + 1]] <- p
    }
  }
  return(plots)
}

make_volcano_plots_tt <- function(stats_files, datatype = "bionav") {
  # determine range of data
  stats_range <- get_stats_data_range(stats_files, datatype = datatype)

  # determine first whether study has both PTK and STK
  assay_types <- stats_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()

  plots <- list()
  if (length(assay_types) == 2) {
    ttest_rows <- stats_files %>% filter(Stats != "MTvC")
    comparisons <- stats_files %>%
      filter(Stats != "MTvC") %>%
      distinct(Comparison) %>%
      pull()
    for (comparison in comparisons) {
      c_rows <- ttest_rows %>% filter(Comparison == comparison)
      
      if (datatype == "tercen") {
        c_ptk <- c_rows %>%
          filter(Assay_Type == "PTK") %>%
          do(extract_phosphosite_data_tt_tercen(.$File, .$Assay_Type))
        c_stk <- c_rows %>%
          filter(Assay_Type == "STK") %>%
          do(extract_phosphosite_data_tt_tercen(.$File, .$Assay_Type))
      } else if (datatype == "bionav") {
        c_ptk <- c_rows %>%
          filter(Assay_Type == "PTK") %>%
          do(extract_phosphosite_data_tt_bionav(.$LFC_file, .$P_file, .$Assay_Type))
        c_stk <- c_rows %>%
          filter(Assay_Type == "STK") %>%
          do(extract_phosphosite_data_tt_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      }
      df <- bind_rows(c_ptk, c_stk)
      if (length(levels(df$Comparison)) > 1) {
        super_df_ptk <- df %>% filter(Assay_type == "PTK")
        super_df_stk <- df %>% filter(Assay_type == "STK")
        p_ptk <- render_volcano_plot_mtvc(super_df_ptk, stats_range$lfc, stats_range$p)
        p_stk <- render_volcano_plot_mtvc(super_df_stk, stats_range$lfc, stats_range$p)
        pg <- plot_grid(p_ptk, p_stk, ncol = 1, labels = c(paste("PTK"), paste("STK")))
      } else {
        p <- render_volcano_plot_tt(df, stats_range$lfc, stats_range$p)
        pg <- plot_grid(p, labels = comparison)
      }
      save_plot(paste0("99_Saved Plots/", comparison, "_Volcano.pdf"), pg, dpi = 300)
      plots[[length(plots) + 1]] <- pg
    }
  } else if (length(assay_types) == 1) {
    ttest_rows <- stats_files %>% filter(Stats != "MTvC")
    comparisons <- stats_files %>%
      filter(Stats != "MTvC") %>%
      distinct(Comparison) %>%
      pull()
    for (comparison in comparisons) {
      c_rows <- ttest_rows %>% filter(Comparison == comparison)
      if (datatype == "tercen") {
        df <- c_rows %>% do(extract_phosphosite_data_tt_tercen(.$File, .$Assay_Type))
      } else if (datatype == "bionav") {
        df <- c_rows %>% do(extract_phosphosite_data_tt_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      }
      if (length(levels(df$Comparison)) > 1) {
        pg <- render_volcano_plot_mtvc(df, stats_range$lfc, stats_range$p)
      } else {
        p <- render_volcano_plot_tt(df, stats_range$lfc, stats_range$p)
        pg <- plot_grid(p, labels = comparison)
      }
      save_plot(paste0("99_Saved Plots/", comparison, "_Volcano.pdf"), pg, dpi = 300)
      plots[[length(plots) + 1]] <- pg
    }
  }
  return(plots)
}

