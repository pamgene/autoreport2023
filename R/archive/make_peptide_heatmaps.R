library(tidyverse)
library(cowplot)
source("R/report_theme.R")

render_heatmap_tt <- function(df, lfc_range, comparison = NULL) {
  assay_types <- df %>% pull(Assay_type)
  if (!"Comparison" %in% colnames(df)) {
    df$Comparison <- comparison
  }
  df <- df %>%
    filter(P < 0.05) %>%
    select(ID, LogFC, Comparison, Assay_type)
  
  clust_mat <- df %>%
    pivot_wider(names_from = Comparison, values_from = LogFC, values_fill = 0) %>%
    column_to_rownames("ID") %>%
    as.matrix(.)
  
  row_fontsize <- case_when(
    # terrible hack
    nrow(clust_mat) / length(assay_types) > 100 ~ 3,
    nrow(clust_mat) / length(assay_types) > 80 ~ 4,
    nrow(clust_mat) / length(assay_types) > 60 ~ 5,
    nrow(clust_mat) / length(assay_types) > 40 ~ 6,
    nrow(clust_mat) / length(assay_types) > 20 ~ 7,
    nrow(clust_mat) / length(assay_types) > 10 ~ 8,
    nrow(clust_mat) / length(assay_types) > 0 ~ 9
  )
  
  # ord <- hclust(dist(clust_mat, method = "euclidean"))$order
    df$ID <- factor(df$ID)

  p <- ggplot(df, aes(x = Comparison, y = fct_reorder(ID, desc(LogFC)))) +
    geom_tile(aes(fill = LogFC), colour = "lightgray") +
    scale_fill_gradient2(limits = c(lfc_range[1], lfc_range[2]), low = "blue", mid = "white", high = "red") +
    # coord_fixed(ratio = 0.5) +
    facet_wrap(~Assay_type, scales = "free", nrow=1) +
    report_theme +
    theme(
      axis.text.y = element_text(size = row_fontsize),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6)
    )
  
  return(p)
}




render_heatmap_mtvc <- function(df, lfc_range, comparison = NULL) {
  if (!"Comparison" %in% colnames(df)) {
    df$Comparison <- comparison
  }
  df <- df %>%
    filter(P < 0.05) %>%
    select(ID, LogFC, Comparison)
  
  clust_mat <- df %>%
    pivot_wider(names_from = Comparison, values_from = LogFC, values_fill = 0) %>%
    column_to_rownames("ID") %>%
    as.matrix(.)
  
  row_fontsize <- case_when(
    nrow(clust_mat) > 100 ~ 3,
    nrow(clust_mat) > 80 ~ 4,
    nrow(clust_mat) > 60 ~ 5,
    nrow(clust_mat) > 40 ~ 6,
    nrow(clust_mat) > 20 ~ 7,
    nrow(clust_mat) > 10 ~ 8,
    nrow(clust_mat) > 0 ~ 9
  )
  
  num_comparisons <- df %>% distinct(Comparison) %>% pull() %>% length(.)
  
  if (num_comparisons > 1) {
    ord <- hclust(dist(clust_mat, method = "euclidean"))$order
    df$ID <- factor(df$ID, levels = rownames(clust_mat)[ord])
  }

  p <- ggplot(df, aes(x = Comparison, y = ID)) +
    geom_tile(aes(fill = LogFC), colour = "lightgray") +
    scale_fill_gradient2(limits = c(lfc_range[1], lfc_range[2]), low = "blue", mid = "white", high = "red") +
    # coord_fixed(ratio = 0.5) +
    report_theme +
    theme(
      axis.text.y = element_text(size = row_fontsize),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6)
    )
  
  return(p)
}


make_heatmaps_mtvc <- function(stats_files, datatype = "bionav") {
  # determine range of data
  stats_range <- get_stats_data_range(stats_files, datatype)
  
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

        heatmap_ptk <- render_heatmap_mtvc(mtvc_ptk, stats_range$lfc)
        heatmap_stk <- render_heatmap_mtvc(mtvc_stk, stats_range$lfc)
        pg <- plot_grid(heatmap_ptk, heatmap_stk, nrow = 1, labels = c(paste("PTK"), paste("STK"), label_x = 0.1, label_y = -4))
        pg <- plot_grid(title, pg, ncol = 1, rel_heights = c(0.1, 1))
        save_plot(paste0("99_Saved Plots/", group, "_MTvC_Heatmap.pdf"), pg, dpi = 300)
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
      heatmap_ptk <- render_heatmap_mtvc(mtvc_ptk, stats_range$lfc)
      heatmap_stk <- render_heatmap_mtvc(mtvc_stk, stats_range$lfc)
      pg <- plot_grid(heatmap_ptk, heatmap_stk, nrow = 1, labels = c(paste("PTK"), paste("STK")))
      save_plot(paste0("99_Saved Plots/", "MTvC_Heatmap.pdf"), pg, dpi = 300)
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
          df <- c_rows %>% do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
        }
        heatmap <- render_heatmap_mtvc(df, stats_range$lfc)
        pg <- plot_grid(title, heatmap, ncol = 1, rel_heights = c(0.1, 1))
        save_plot(paste0("99_Saved Plots/", group, "_MTvC_Heatmap.pdf"), pg, dpi = 300)
        plots[[length(plots) + 1]] <- pg
      }
    } else {
      if (datatype == "bionav") {
        df <- mtvc_rows %>% do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      } else if (datatype == "tercen") {
        df <- mtvc_rows %>% do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
      }      
      heatmap <- render_heatmap_mtvc(df, stats_range$lfc)
      save_plot(paste0("99_Saved Plots/", "MTvC_Heatmap.pdf"), heatmap, dpi = 300)
      plots[[length(plots) + 1]] <- heatmap
    }
  }
  return(plots)
}

make_heatmaps_tt <- function(stats_files, datatype = "bionav") {
  # determine range of data
  stats_range <- get_stats_data_range(stats_files, datatype)
  
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
        h_ptk <- render_heatmap_mtvc(c_ptk, stats_range$lfc)
        h_stk <- render_heatmap_mtvc(c_stk, stats_range$lfc)
        pg <- plot_grid(h_ptk, h_stk, nrow = 1, labels = c(paste("PTK"), paste("STK")))
      } else {
        p <- render_heatmap_tt(df, stats_range$lfc, comparison)
        pg <- plot_grid(p, labels = comparison)
      }
      save_plot(paste0("99_Saved Plots/", comparison, "_Heatmap.pdf"), pg, dpi = 300)
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
      if(length(levels(df$Comparison)) > 1) {
        pg <- render_heatmap_mtvc(df, stats_range$lfc)
      } else {
        p <- render_heatmap_tt(df, stats_range$lfc)
        pg <- plot_grid(p, labels = comparison)
      }
      save_plot(paste0("99_Saved Plots/", comparison, "_Heatmap.pdf"), pg, dpi = 300)
      plots[[length(plots) + 1]] <- pg
    }
  }
  return(plots)
}