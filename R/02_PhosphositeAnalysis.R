require(tidyverse)
require(cowplot)
require(EnhancedVolcano)
require(egg)
library(gridGraphics)
library(gridExtra)
source("R/report_theme.R")





#########################################################
################# LOAD PHOSPHOSITE DATA #################
#########################################################

extract_phosphosite_data_tt_tercen <- function(filepath, assay_type) {
  df <- read_delim(filepath, show_col_types = FALSE)
  df <- clean_tercen_columns(df)
  if ("UniprotAccession" %in% colnames(df)){
    df <- df %>% select(-UniprotAccession) %>% distinct()
  }
  if (length(colnames(df)) == 4) {
    # supergroup
    df <- df %>%
      rename(Comparison = 1, LogFC = delta, P = p) %>%
      mutate(Comparison = as.factor(Comparison))
  } else if (length(colnames(df)) == 3) {
    # no supergroup
    df <- df %>%
      rename(LogFC = delta, P = p)
  }
  df$Assay_type <- assay_type
  return(df)
}

extract_phosphosite_data_tt_bionav <- function(lfc_file, p_file, assay_type) {
  # first identify whether it's a supergroup file or not
  titles <- read_lines(p_file, n_max = 2)
  split_titles <- unlist(str_split(titles[2], "\t"))
  if (length(split_titles) == 1) {
    lfc_cols <- c("ID", "Uniprot", "Sequence", "Empty", "LogFC")
    p_cols <- c("ID", "Uniprot", "Sequence", "Empty", "P")

    logfc <- read_delim(lfc_file, skip = 3, col_names = lfc_cols, show_col_types = FALSE)
    p <- read_delim(p_file, skip = 3, col_names = p_cols, show_col_types = FALSE)
    df <- logfc %>% select(ID, LogFC)
    df$P <- p %>%
      select(P) %>%
      pull()
    df <- df %>% mutate(
      LogFC = as.numeric(LogFC),
      P = as.numeric(P)
    )
  } else {
    p_names <- c("ID", "Uniprot", "Sequence", "Empty", split_titles[5:length(split_titles)])
    p_names <- p_names[p_names != ""]

    logfc <- read_delim(lfc_file, skip = 3, col_names = p_names, show_col_types = FALSE)
    pvalues <- read_delim(p_file, skip = 3, col_names = p_names, show_col_types = FALSE)

    df <- list()
    for (title in split_titles[5:length(split_titles)]) {
      if ((title != "")) {
        comp_df <- logfc %>% select(ID, title)

        comp_df$pvalue <- pvalues %>%
          select(title) %>%
          pull()

        colnames(comp_df) <- c("ID", "LogFC", "P")
        comp_df$LogFC <- as.numeric(comp_df$LogFC)
        comp_df$P <- as.numeric(comp_df$P)

        comparison <- paste(title)
        comp_df$Comparison <- comparison
        df[[length(df) + 1]] <- comp_df
      }
    }
    df <- bind_rows(df)
    df$Comparison <- as.factor(df$Comparison)
  }
  df$Assay_type <- assay_type
  return(df)

}

extract_phosphosite_data_mtvc_bionav <- function(lfc_file, p_file, assay_type) {
  titles <- read_lines(p_file, n_max = 2)
  split_titles <- unlist(str_split(titles[2], "\t"))
  p_names <- c("c_cluster", "ID", "Uniprot", "Empty", split_titles[5:length(split_titles)])
  p_names <- p_names[p_names != ""]

  logfc <- read_delim(lfc_file, skip = 3, col_names = p_names, show_col_types = FALSE)
  pvalues <- read_delim(p_file, skip = 3, col_names = p_names, show_col_types = FALSE)

  dfs <- list()
  ctrl_title <- ""

  # first find which column is control, e.g., all NAs
  for (title in split_titles[5:length(split_titles)]) {
    if (title != "") {
      df <- logfc %>% select(ID, title)
      colnames(df) <- c("ID", "LogFC")
      df$LogFC <- as.numeric(df$LogFC)

      if (suppressWarnings(is.na(any(df$LogFC)))) {
        ctrl_title <- title
      }
    }
  }

  for (title in split_titles[5:length(split_titles)]) {
    if ((title != ctrl_title) & (title != "")) {
      df <- logfc %>% select(ID, title)
      df$pvalue <- pvalues %>%
        select(title) %>%
        pull()
      colnames(df) <- c("ID", "LogFC", "P")
      df$LogFC <- as.numeric(df$LogFC)
      df$P <- as.numeric(df$P)

      comparison <- paste(title, "vs", ctrl_title)
      df$Comparison <- comparison
      dfs[[length(dfs) + 1]] <- df
    }
  }

  dfs <- bind_rows(dfs)
  dfs$Comparison <- as.factor(dfs$Comparison)
  dfs$Assay_type <- assay_type

  return(dfs)
}

extract_phosphosite_data_mtvc_tercen <- function(filepath, assay_type) {
  df <- read_delim(filepath, show_col_types = FALSE)
  df <- clean_tercen_columns(df)
  df <- df %>%
    rename(Comparison = 1) %>%
    mutate(variable = case_when(grepl("\\.LogFC", variable) ~ "LogFC", grepl("\\.pvalue", variable) ~ "P"))
  ctrl <- df %>%
    filter(is.na(value)) %>%
    distinct(Comparison) %>%
    pull()

  df <- df %>% filter(Comparison != ctrl) %>%
    mutate(Comparison = paste(Comparison, "vs", ctrl),
           Comparison = as.factor(Comparison),
           Assay_type = assay_type) %>%
    pivot_wider(names_from = variable, values_from = value)

  return(df)
}

#########################################################
############ STRINGIFY PHOSPHOSITE DIRECTION ############
#########################################################

extract_direction_tt <- function(filepath, comparison, assay_type, datatype, p_file = NULL) {
  if (datatype == "tercen") {
    df <- extract_phosphosite_data_tt_tercen(filepath, assay_type)
    if (length(colnames(df)) == 4){
      supergroup_file <- F
    } else if (length(colnames(df)) == 5){
      supergroup_file <- T
    }
  } else if (datatype == "bionav") {
    df <- extract_phosphosite_data_tt_bionav(filepath, p_file, assay_type)
    if (length(colnames(df)) == 4){
      supergroup_file <- F
    } else if (length(colnames(df)) == 5){
      supergroup_file <- T
    }
  }
  
  

  if (!supergroup_file) {
    if (assay_type == "PTK") {
      df <- df %>%
        drop_na() %>%
        summarise(
          comparison = comparison,
          up_PTK = sum((LogFC > 0) & (P < 0.05)),
          down_PTK = sum((LogFC < 0) & (P < 0.05)),
        )
    } else if (assay_type == "STK") {
      df <- df %>%
        drop_na() %>%
        summarise(
          comparison = comparison,
          up_STK = sum((LogFC > 0) & (P < 0.05)),
          down_STK = sum((LogFC < 0) & (P < 0.05)),
        )
    }
  } else if (supergroup_file) {
    output <- list()
    for (group in levels(df$Comparison)) {
      comparison_sg <- paste(comparison, group)
      direction <- df %>%
        filter(Comparison == group) %>%
        drop_na() %>%
        summarise(
          comparison = comparison_sg,
          up = sum((LogFC > 0) & (P < 0.05)),
          down = sum((LogFC < 0) & (P < 0.05))
        )
      direction$assay_type <- assay_type
      direction$stats <- "tt"
      direction <- direction %>% pivot_wider(names_from = assay_type, values_from = c(up, down))
      output[[length(output) + 1]] <- direction
    }
    df <- bind_rows(output)
  }
  df$stats <- "tt"

  return(df)
}

extract_direction_mtvc <- function(filepath, assay_type, group, datatype, p_file = NULL) {
  if (datatype == "tercen") {
    df <- extract_phosphosite_data_mtvc_tercen(filepath, assay_type)
  } else if (datatype == "bionav") {
    df <- extract_phosphosite_data_mtvc_bionav(filepath, p_file, assay_type)
  }

  output <- df %>%
              drop_na() %>%
              group_by(Comparison, Assay_type) %>%
              summarise(up = sum((LogFC > 0) & (P < 0.05)),
                        down = sum(LogFC < 0 & (P < 0.05)))

  if (!is.na(group)) {
    output <- output %>% mutate(Comparison = paste(group, Comparison))
  }

  output <- output %>% rename(comparison = Comparison,
                              assay_type = Assay_type) %>%
                      mutate(stats = "mtvc")

  output <- output %>% pivot_wider(names_from = assay_type, values_from = c(up, down))

  return(output)
}



parse_mtvc <- function(stat_files, datatype, assay_types){
  dfs <- list()
  mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
  groups <- stats_files %>%
    filter(Stats == "MTvC") %>%
    distinct(Group) %>%
    pull()
  if (!is.na(groups)) {
    for (group in groups) {
      group_rows <- mtvc_rows %>% filter(Group == group)
      if (datatype == "bionav") {
        if (length(assay_types) == 2){
          mtvc_ptk <- group_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_direction_mtvc(.$LFC_file, .$Assay_Type, .$Group, datatype, .$P_file))
          mtvc_stk <- group_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_direction_mtvc(.$LFC_file, .$Assay_Type, .$Group, datatype, .$P_file))
          df <- left_join(mtvc_ptk, mtvc_stk, by = c("comparison", "stats"))
        } else if (length(assay_types) == 1){
          df <- group_rows %>% do(extract_direction_mtvc(.$LFC_file, .$Assay_Type, .$Group, datatype, .$P_file))
        } 
      } else if (datatype == "tercen") {
        if (length(assay_types) == 2){
          mtvc_ptk <- group_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_direction_mtvc(.$File, .$Assay_Type, .$Group, datatype))
          mtvc_stk <- group_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_direction_mtvc(.$File, .$Assay_Type, .$Group, datatype))
          df <- left_join(mtvc_ptk, mtvc_stk, by = c("comparison", "stats"))
        } else if (length(assay_types) == 1){
          df <- group_rows %>% do(extract_direction_mtvc(.$File, .$Assay_Type, .$Group, datatype))
        }
      }
      dfs <- append(dfs, list(df))
    }
  } 
  # else {
  #   # this functionality should be deprecated
  #   if (datatype == "bionav") {
  #     if (length(assay_types) == 2){
  #       mtvc_ptk <- mtvc_rows %>%
  #         filter(Assay_Type == "PTK") %>%
  #         do(extract_direction_mtvc(.$LFC_file, .$Assay_Type, .$Group, datatype, .$P_file))
  #       mtvc_stk <- mtvc_rows %>%
  #         filter(Assay_Type == "STK") %>%
  #         do(extract_direction_mtvc(.$LFC_file, .$Assay_Type, .$Group, datatype, .$P_file))
  #       df <- left_join(mtvc_ptk, mtvc_stk, by = c("comparison", "stats"))
  #     } else if (length(assay_types) == 1) {
  #       df <- mtvc_rows %>% do(extract_direction_mtvc(.$LFC_file, .$Assay_Type, .$Group, datatype, .$P_file))
  #     } 
  #   } else if (datatype == "tercen") {
  #     if (length(assay_types) == 2){
  #       mtvc_ptk <- mtvc_rows %>%
  #         filter(Assay_Type == "PTK") %>%
  #         do(extract_direction_mtvc(.$File, .$Assay_Type, .$Group, datatype))
  #       mtvc_stk <- mtvc_rows %>%
  #         filter(Assay_Type == "STK") %>%
  #         do(extract_direction_mtvc(.$File, .$Assay_Type, .$Group, datatype))
  #       df <- left_join(mtvc_ptk, mtvc_stk, by = c("comparison", "stats"))
  #     } else if (length(assay_types) == 1){
  #       df <- mtvc_rows %>% do(extract_direction_mtvc(.$File, .$Assay_Type, .$Group, datatype))
  #     }
  #   }
  #   dfs <- append(dfs, list(df))
  # }
  return(dfs)
}



parse_tt <- function(stats_files, datatype, assay_types){
  dfs <- list()
  ttest_rows <- stats_files %>% filter(Stats != "MTvC")
  comparisons <- stats_files %>%
    filter(Stats != "MTvC") %>%
    distinct(Comparison) %>%
    pull()
  for (comparison in comparisons) {
    c_rows <- ttest_rows %>% filter(Comparison == comparison)
    if (datatype == "bionav") {
      if (length(assay_types) == 2){
        c_ptk <- c_rows %>%
          filter(Assay_Type == "PTK") %>%
          do(extract_direction_tt(.$LFC_file, .$Comparison, .$Assay_Type, datatype, .$P_file))
        c_stk <- c_rows %>%
          filter(Assay_Type == "STK") %>%
          do(extract_direction_tt(.$LFC_file, .$Comparison, .$Assay_Type, datatype, .$P_file))
        df <- left_join(c_ptk, c_stk, by = c("comparison", "stats"))
      } else if (length(assay_types) == 1) {
        df <- c_rows %>% do(extract_direction_tt(.$LFC_file, .$Comparison, .$Assay_Type, datatype, .$P_file))
      }
    } else if (datatype == "tercen") {
      if (length(assay_types) == 2){
        c_ptk <- c_rows %>%
          filter(Assay_Type == "PTK") %>%
          do(extract_direction_tt(.$File, .$Comparison, .$Assay_Type, datatype))
        c_stk <- c_rows %>%
          filter(Assay_Type == "STK") %>%
          do(extract_direction_tt(.$File, .$Comparison, .$Assay_Type, datatype))
        df <- left_join(c_ptk, c_stk, by = c("comparison", "stats"))
      } else if (length(assay_types) == 1){
        df <- c_rows %>% do(extract_direction_tt(.$File, .$Comparison, .$Assay_Type, datatype))
      }
    }
    dfs <- append(dfs, list(df))
  }
  return(dfs)
}


enrich_results <- function(stats_files){
  e <- read_csv("data/enrichment_86402_86412_87102_arrays.csv") %>% select(-family, -SeqMatch)
  
  if ("TT" %in% stats_files$Stats){
    ttest_rows <- stats_files %>% filter(Stats != "MTvC")
    for (i in 1:nrow(ttest_rows)){
      assaytype <- ttest_rows[i,]$Assay_Type %>% as.character()
      df <- ttest_rows[i,] %>% 
        do(extract_phosphosite_data_tt_tercen(filepath = .$File, assay_type = .$Assay_Type))
      df_e <- left_join(df, e, by = "ID")
      first_cols <- colnames(df_e)[! colnames(df_e) %in% c("LogFC", "P", "Assay_type")]
      df_e <- df_e[c(first_cols, "LogFC", "P", "Assay_type")]
      write_csv(df_e, paste0("99_Saved Plots/TT_", assaytype, "_", ttest_rows[i,]$Comparison, "_enriched", ".csv"))
    }
  }
  
  if ("MTvC" %in% stats_files$Stats){
    mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
    for (i in 1:nrow(mtvc_rows)){
      assaytype <- mtvc_rows[i,]$Assay_Type %>% as.character()
      df <- mtvc_rows[i,] %>% 
        do(extract_phosphosite_data_mtvc_tercen(filepath = .$File, assay_type = .$Assay_Type))
      df_e <- left_join(df, e, by = "ID")
      first_cols <- colnames(df_e)[! colnames(df_e) %in% c("LogFC", "P", "Assay_type")]
      df_e <- df_e[c(first_cols, "LogFC", "P", "Assay_type")]
      write_csv(df_e, paste0("99_Saved Plots/MTvC_", assaytype, "_", mtvc_rows[i,]$Group, "_enriched", ".csv"))
    }
  }
}


parse_stats_files <- function(stats_files, datatype = "bionav") {
  if (datatype == "tercen"){
    enrich_results(stats_files)
  }
  
  # determine whether study has both PTK and STK
  assay_types <- stats_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()
  
  if ("MTvC" %in% stats_files$Stats) {
    dfs_m <- parse_mtvc(stats_files, datatype, assay_types)
  }
  if ("TT" %in% stats_files$Stats) {
    dfs_t <- parse_tt(stats_files, datatype, assay_types)
  }
  all_dfs <- append(dfs_m, dfs_t)
  return(bind_rows(all_dfs))
}


#########################################################
################ OUTPUT PHOSPHOSITE TABLE ###############
#########################################################

stats_footnotes <- function(ft, res_table) {
  stats_types <- res_table %>%
    ungroup() %>%
    distinct(stats) %>%
    pull()
  if (length(stats_types) == 2) {
    ft <- ft %>%
      footnote(~stats == "mtvc", ~ comparison, value = as_paragraph("Significance was obtained using a one-way ANOVA followed by a post-hoc Dunnett\u0027s test, p<0.05"), ref_symbols = c("a"), inline = TRUE) %>%
      footnote(~stats == "tt", ~ comparison, value = as_paragraph("Significance was obtained using a two-sided unpaired Student\u0027s T-test, p<0.05."), ref_symbols = c("b"), inline = TRUE)
  } else if (stats_types == "mtvc") {
    ft <- ft %>%
      footnote(~stats == "mtvc", ~ comparison, value = as_paragraph("Significance was obtained using a one-way ANOVA followed by a post-hoc Dunnett\u0027s test, p<0.05"), ref_symbols = c("a"), inline = TRUE)
  } else if (stats_types == "tt") {
    ft <- ft %>%
      footnote(~stats == "tt", ~ comparison, value = as_paragraph("Significance was obtained using a two-sided unpaired Student\u0027s T-test, p<0.05."), ref_symbols = c("a"), inline = TRUE)
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

#########################################################
############## CALCULATE PHOSPHOSITE RANGE ##############
#########################################################

get_stats_data_range <- function(stats_files, datatype = "bionav") {
  lfc <- c()
  p <- c()
  if ("MTvC" %in% stats_files$Stats) {
    mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
    if (datatype == "bionav") {
      mtvc <- bind_rows(apply(mtvc_rows, 1, function(x) extract_phosphosite_data_mtvc_bionav(lfc_file = x[6], p_file = x[7], assay_type = x[4])))
    } else if (datatype == "tercen") {
      mtvc <- bind_rows(apply(mtvc_rows, 1, function(x) extract_phosphosite_data_mtvc_tercen(filepath = x[6], assay_type = x[4])))
    }
    lfc <- c(lfc, mtvc %>% pull(LogFC))
    p <- c(p, mtvc %>% pull(P))
  }
  if ("TT" %in% stats_files$Stats) {
    ttest_rows <- stats_files %>% filter(Stats != "MTvC")
    if (datatype == "bionav") {
      tt <- bind_rows(apply(ttest_rows, 1, function(x) extract_phosphosite_data_tt_bionav(lfc_file = x[6], p_file = x[7], assay_type = x[4])))
    } else if (datatype == "tercen") {
      tt <- bind_rows(apply(ttest_rows, 1, function(x) extract_phosphosite_data_tt_tercen(filepath = x[6], assay_type = x[4])))
    }
    lfc <- c(lfc, tt %>% pull(LogFC))
    p <- c(p, tt %>% pull(P))
  }

  max_lfc <- max(lfc, na.rm = T)
  max_lfc <- ceiling(max_lfc / 0.5) * 0.5
  min_lfc <- min(lfc, na.rm = T)
  min_lfc <- floor(min_lfc / 0.5) * 0.5

  if (abs(max_lfc) > abs(min_lfc)) {
    lfc_range <- c(-max_lfc, max_lfc)
  } else {
    lfc_range <- c(min_lfc, - min_lfc)
  }

  # now p
  max_p <- (ceiling(-log10(min(p, na.rm = TRUE))))
  min_p <- 0
  p_range <- c(min_p, max_p)

  return(list(lfc = lfc_range, p = p_range))
}

#########################################################
##################### VOLCANO PLOTS #####################
#########################################################

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
        save_plot(paste0("99_Saved Plots/", "MTvC_", group, "_Volcano.pdf"), pg, dpi = 300)
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
          df <- c_rows %>% do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
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
      save_plot(paste0("99_Saved Plots/TT_", comparison, "_Volcano.pdf"), pg, dpi = 300)
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
      save_plot(paste0("99_Saved Plots/TT_", comparison, "_Volcano.pdf"), pg, dpi = 300)
      plots[[length(plots) + 1]] <- pg
    }
  }
  return(plots)
}

#########################################################
####################### HEATMAPS ########################
#########################################################

get_plotparams <- function(dflist, comparison = NULL){
  myvect <- c()
  
  for (i in seq_along(dflist)){
    df <- dflist[[i]] 
    if (!"Comparison" %in% colnames(df)) {
      df <- df %>% mutate(Comparison = comparison)
    }
    df <- df %>%
      filter(P < 0.05) %>%
      select(ID, LogFC, Comparison) %>%
      pivot_wider(names_from = Comparison, values_from = LogFC, values_fill = 0)
    
    myvect <- append(myvect, nrow(df))
  }
  
  dfmax <- max(myvect)
  # which plot will be bigger?
  bigger <- ifelse(match(dfmax, myvect) == 1, "PTK", "STK")
  
  nrows <- 1
  
  # Nrows could be 2 in some cases --> to test
  # nrows <- case_when(
  #   myvect[1] > 30 & myvect[2] > 30 ~ 1
  #   myvect[1] < 10 & myvect[2] > 30 ~ 2)
  
  dfsum <- sum(myvect)
  
  h <- case_when(
    nrows == 1 & dfmax > 140 ~ 19,
    nrows == 1 & dfmax > 120 ~ 18,
    nrows == 1 & dfmax > 100 ~ 17,
    nrows == 1 & dfmax > 90 ~ 16.5,
    nrows == 1 & dfmax > 80 ~ 16,
    nrows == 1 & dfmax > 70 ~ 15.5,
    nrows == 1 & dfmax > 60 ~ 15,
    nrows == 1 & dfmax > 50 ~ 14.5,
    nrows == 1 & dfmax > 40 ~ 14,
    nrows == 1 & dfmax > 30 ~ 13.5,
    nrows == 1 & dfmax > 20 ~ 13,
    nrows == 1 & dfmax > 10 ~ 12.5,
    nrows == 1 & dfmax > 0 ~ 12
    # nrows == 2 & dfsum > 120 ~ 16,
    # nrows == 2 & dfsum > 110 ~ 15.5,
    # nrows == 2 & dfsum > 100 ~ 15,
    # nrows == 2 & dfsum > 90 ~ 14.5,
    # nrows == 2 & dfsum > 80 ~ 14
  )
  
  
  return(list(h = h, n = nrows, b = bigger))
  
}



# render_heatmap_tt <- function(df, lfc_range, comparison = NULL) {
#   assay_types <- df %>% pull(Assay_type)
#   if (!"Comparison" %in% colnames(df)) {
#     df$Comparison <- comparison
#   }
#   df <- df %>%
#     filter(P < 0.05) %>%
#     select(ID, LogFC, Comparison, Assay_type)
# 
#   clust_mat <- df %>%
#     pivot_wider(names_from = Comparison, values_from = LogFC, values_fill = 0) %>%
#     column_to_rownames("ID") %>%
#     as.matrix(.)
# 
#   row_fontsize <- case_when(
#   # terrible hack
#     nrow(clust_mat) / length(assay_types) > 100 ~ 3,
#     nrow(clust_mat) / length(assay_types) > 90 ~ 3.5,
#     nrow(clust_mat) / length(assay_types) > 80 ~ 4,
#     nrow(clust_mat) / length(assay_types) > 70 ~ 4.5,
#     nrow(clust_mat) / length(assay_types) > 60 ~ 5,
#     nrow(clust_mat) / length(assay_types) > 50 ~ 5.5,
#     nrow(clust_mat) / length(assay_types) > 40 ~ 6,
#     nrow(clust_mat) / length(assay_types) > 30 ~ 6.5,
#     nrow(clust_mat) / length(assay_types) > 20 ~ 7,
#     nrow(clust_mat) / length(assay_types) > 10 ~ 8,
#     nrow(clust_mat) / length(assay_types) > 0 ~ 9
#   )
# 
#   # ord <- hclust(dist(clust_mat, method = "euclidean"))$order
#   df$ID <- factor(df$ID)
# 
#   p <- ggplot(df, aes(x = Comparison, y = fct_reorder(ID, desc(LogFC)))) +
#     geom_tile(aes(fill = LogFC), colour = "lightgray") +
#     scale_fill_gradient2(limits = c(lfc_range[1], lfc_range[2]), low = "blue", mid = "white", high = "red") +
#     facet_wrap(~Assay_type, scales = "free", nrow = 1) +
#     report_theme +
#     theme(
#       #axis.text.y = element_text(size = row_fontsize),
#       axis.title.y = element_blank(),
#       axis.text.x = element_text(angle = 90, vjust = 0.5),
#       axis.title.x = element_blank(),
#       legend.key.size = unit(0.4, "cm"),
#       legend.title = element_text(size = 8),
#       legend.text = element_text(size = 6)
#     )
#   
#   # p1 <- egg::set_panel_size(p, height=unit(nrow(clust_mat)/2.5, "cm"),
#   #                           width=unit(ncol(clust_mat)*1.5, "cm") )
#   
#   return(p)
# }


render_heatmap_mtvc <- function(df, lfc_range, comparison = NULL, assay_type) {
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
  
  num_comparisons <- df %>% distinct(Comparison) %>% pull() %>% length(.)
  
  if (num_comparisons > 1) {
    ord <- hclust(dist(clust_mat, method = "euclidean"))$order
    df$ID <- factor(df$ID, levels = rownames(clust_mat)[ord])
  }
  p <- ggplot(df, aes(x = Comparison, y = ID)) +
    geom_tile(aes(fill = LogFC), colour = "lightgray") +
    scale_fill_gradient2(limits = c(lfc_range[1], lfc_range[2]), low = "blue", mid = "white", high = "red") +
    #coord_fixed(clip = 'off') +
    ggtitle(assay_type) +
    report_theme +
    theme(
      #axis.text.y = element_text(size = row_fontsize),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      axis.title.x = element_blank(),
      legend.key.size = unit(0.4, "cm"),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 6),
      #plot.margin = margin(10,2,10,2, "cm")
    )
  
  
  tile_h_scale <- ifelse(nrow(clust_mat) > 70, 3.5, 3)
  tile_w_scale <- ifelse(nrow(clust_mat) > 70, 1.2, 1.5)
  #print(tile_h_scale)
  p1 <- egg::set_panel_size(p, height=unit(nrow(clust_mat)/tile_h_scale, "cm"),
                            width=unit(ncol(clust_mat)*tile_w_scale, "cm") )
  return(p1)
}




make_heatmaps_mtvc <- function(stats_files, datatype = "bionav") {
  # determine range of data
  stats_range <- get_stats_data_range(stats_files, datatype)

  # determine first whether study has both PTK and STK
  assay_types <- stats_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()

  mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
  groups <- stats_files %>%
    filter(Stats == "MTvC") %>%
    distinct(Group) %>%
    pull()
  
  plots <- list()
  if (length(assay_types) == 2) {
    if (!is.na(groups)) {
      # If there are groups
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

        heatmap_ptk <- render_heatmap_mtvc(mtvc_ptk, stats_range$lfc, assay_type = "PTK")
        heatmap_stk <- render_heatmap_mtvc(mtvc_stk, stats_range$lfc, assay_type = "STK")
        # get plotting size parameters
        ppars <- get_plotparams(list(mtvc_ptk, mtvc_stk))
        pg <- cowplot::plot_grid(heatmap_ptk, heatmap_stk, nrow = ppars$n)
        ggsave(paste0("99_Saved Plots/", "MTvC_", group, "_Heatmap.pdf"), pg, width = 10, height = ppars$h)
        plots[[length(plots) + 1]] <- pg
      }
    } else {
      # else: there are no groups
      # DO NOT USE THIS OPTION!
      # if (datatype == "bionav") {
      #   mtvc_ptk <- mtvc_rows %>%
      #     filter(Assay_Type == "PTK") %>%
      #     do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      #   mtvc_stk <- mtvc_rows %>%
      #     filter(Assay_Type == "STK") %>%
      #     do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      # } else if (datatype == "tercen") {
      #   mtvc_ptk <- mtvc_rows %>%
      #     filter(Assay_Type == "PTK") %>%
      #     do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
      #   mtvc_stk <- mtvc_rows %>%
      #     filter(Assay_Type == "STK") %>%
      #     do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
      # }
      # heatmap_ptk <- render_heatmap_mtvc(mtvc_ptk, stats_range$lfc, assay_type = "PTK")
      # heatmap_stk <- render_heatmap_mtvc(mtvc_stk, stats_range$lfc, assay_type = "STK")
      # ppars <- get_plotparams(list(mtvc_ptk, mtvc_stk))
      # pg <- cowplot::plot_grid(heatmap_ptk, heatmap_stk, nrow=ppars$n)
      # ggsave(paste0("99_Saved Plots/", group, "_MTvC_Heatmap.pdf"), pg, width = 10, height = ppars$h)
      # plots[[length(plots) + 1]] <- pg
    }
  } else if (length(assay_types) == 1) {
    mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
    groups <- stats_files %>%
      filter(Stats == "MTvC") %>%
      distinct(Group) %>%
      pull()
    assay_type = mtvc_rows %>% pull(Assay_Type)
    if (!is.na(groups)) {
      for (group in groups) {
        c_rows <- mtvc_rows %>% filter(Group == group)
        if (datatype == "bionav") {
          df <- c_rows %>% do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
        } else if (datatype == "tercen") {
          df <- c_rows %>% do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
        }
        heatmap <- render_heatmap_mtvc(df, stats_range$lfc, assay_type = assay_type)
        ppars <- get_plotparams(list(df))
        ggsave(paste0("99_Saved Plots/", "MTvC_", group, "_Heatmap.pdf"), heatmap, width = 10, height = ppars$h)
        plots[[length(plots) + 1]] <- heatmap
      }
    } else {
      if (datatype == "bionav") {
        df <- mtvc_rows %>% do(extract_phosphosite_data_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      } else if (datatype == "tercen") {
        df <- mtvc_rows %>% do(extract_phosphosite_data_mtvc_tercen(.$File, .$Assay_Type))
      }
      
      heatmap <- render_heatmap_mtvc(df, stats_range$lfc, assay_type = assay_type)
      ppars <- get_plotparams(list(df))
      ggsave(paste0("99_Saved Plots/", "MTvC_Heatmap.pdf"), heatmap, width = 10, height = ppars$h)
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
  
  ttest_rows <- stats_files %>% filter(Stats != "MTvC")
  comparisons <- stats_files %>%
    filter(Stats != "MTvC") %>%
    distinct(Comparison) %>%
    pull()
  
  plots <- list()
  if (length(assay_types) == 2) {
    for (comparison in comparisons) {
      c_rows <- ttest_rows %>% filter(Comparison == comparison)
      c_ptk <- c_rows %>% filter(Assay_Type == "PTK")
      c_stk <- c_rows %>% filter(Assay_Type == "STK")
      if (datatype == "tercen") {
        c_ptk <- c_ptk %>%
          do(extract_phosphosite_data_tt_tercen(.$File, .$Assay_Type))
        c_stk <- c_stk %>%
          do(extract_phosphosite_data_tt_tercen(.$File, .$Assay_Type))
      } else if (datatype == "bionav") {
        c_ptk <- c_ptk %>%
          do(extract_phosphosite_data_tt_bionav(.$LFC_file, .$P_file, .$Assay_Type))
        c_stk <- c_stk %>%
          do(extract_phosphosite_data_tt_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      }
      ppars <- get_plotparams(list(c_ptk, c_stk), comparison = comparison)
      df <- bind_rows(c_ptk, c_stk)
      if (length(levels(df$Comparison)) > 1) {
        # TT with supergroup of min 2 levels
        h_ptk <- render_heatmap_mtvc(c_ptk, stats_range$lfc, assay_type = "PTK")
        h_stk <- render_heatmap_mtvc(c_stk, stats_range$lfc, assay_type = "STK")
        pg <- cowplot::plot_grid(h_ptk, h_stk, nrow = 1)
      } else {
        h_ptk <- render_heatmap_mtvc(c_ptk, stats_range$lfc, comparison = comparison, assay_type = "PTK")
        h_stk <- render_heatmap_mtvc(c_stk, stats_range$lfc, comparison = comparison, assay_type = "STK")
        pg <- cowplot::plot_grid(h_ptk, h_stk, nrow = 1)
      }
      ggsave(paste0("99_Saved Plots/", "TT_", comparison, "_Heatmap.pdf"), pg, width = 10, height = ppars$h)
      plots[[length(plots) + 1]] <- pg
    }
  } else if (length(assay_types) == 1) {
    for (comparison in comparisons) {
      c_rows <- ttest_rows %>% filter(Comparison == comparison)
      if (datatype == "tercen") {
        df <- c_rows %>% do(extract_phosphosite_data_tt_tercen(.$File, .$Assay_Type))
      } else if (datatype == "bionav") {
        df <- c_rows %>% do(extract_phosphosite_data_tt_bionav(.$LFC_file, .$P_file, .$Assay_Type))
      }
      ppars <- get_plotparams(list(df), comparison = comparison)
      if (length(levels(df$Comparison)) > 1) {
        p <- render_heatmap_mtvc(df, stats_range$lfc, assay_type = df[1, "Assay_type"])
      } else {
        p <- render_heatmap_mtvc(df, stats_range$lfc, comparison = comparison, assay_type = df[1, "Assay_type"])
      }
      pg <- plot_grid(p)
      ggsave(paste0("99_Saved Plots/", "TT_", comparison, "_Heatmap.pdf"), pg, width = 10, height = ppars$h)
      plots[[length(plots) + 1]] <- pg
    }
  }
  return(plots)
}

