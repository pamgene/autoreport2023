library(tidyverse)

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

  max_lfc <- max(lfc, na.rm = TRUE)
  max_lfc <- ceiling(max_lfc / 0.5) * 0.5
  min_lfc <- min(lfc, na.rm = TRUE)
  min_lfc <- floor(min_lfc / 0.5) * 0.5

  if (abs(max_lfc) > abs(min_lfc)) {
    lfc_range <- c(-max_lfc, max_lfc)
  } else {
    lfc_range <- c(min_lfc, -min_lfc)
  }

  # now p
  max_p <- (ceiling(-log10(min(p, na.rm = TRUE))))
  min_p <- 0
  p_range <- c(min_p, max_p)

  return(list(lfc = lfc_range, p = p_range))
}

extract_phosphosite_data_tt_tercen <- function(filepath, assay_type) {
  df <- read_delim(filepath, show_col_types = FALSE)
  df <- clean_tercen_columns(df)
  if (length(colnames(df)) == 5) {
    # supergroup
    df <- df %>% 
      rename(Comparison = 1) %>%
      mutate(variable = case_when(grepl("\\.delta", variable) ~ "LogFC", grepl("\\.p$", variable) ~ "P"),
             Comparison = as.factor(Comparison)) %>%
      distinct(Comparison, ID, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value)
  } else {
    # no supergroup
    df <- df %>% 
      mutate(variable = case_when(grepl("\\.delta", variable) ~ "LogFC", grepl("\\.p$", variable) ~ "P")) %>%
      distinct(ID, variable, value) %>%
      pivot_wider(names_from = variable, values_from = value)
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