read_phosphosite_dir <- function(folder = "02_Phosphosite Analysis/", datatype = "bionav") {
  # will return df with files for easy processing
  files <- list.files(folder, pattern = ".txt$|.csv$", full.names = TRUE, include.dirs = FALSE)
  if (length(files) == 0) {
    warning("No stats files")
    return(data.frame())
  }
  
  if (datatype == "bionav") {
    files <- files[grepl("LogFC", files)]
  }

  dfs <- list()
  for (i in seq_along(files)) {
    file <- files[i]

    file_base <- basename(file)
    file_elements <- str_split(tools::file_path_sans_ext(file_base), pattern = "_")
    stats_type <- file_elements[[1]][1]
    assay_type <- file_elements[[1]][2]
    
    if (stats_type == "MTvC") {
      if (length(file_elements[[1]]) > 3) {
        order <- as.integer(file_elements[[1]][3])
        group <- file_elements[[1]][4]
      } else {
        order <- 1
        group <- NA
      }
      comparison <- "MTvC"
    } else if (stats_type == "TT") {
      order <- as.integer(file_elements[[1]][3]) + 1
      comparison <- file_elements[[1]][4]
      comparison <- str_replace(comparison, "\\Bvs\\B", " vs ")
      group <- NA
    }
    
    if (datatype == "bionav") {
      if (str_detect(file_elements[[1]][length(file_elements[[1]])], "LogFC")) {
        p_val_file <- str_replace(file, "LogFC", "p")
        df <- tibble("Order" = order, "Comparison" = comparison, "Stats" = stats_type, "Assay_Type" = assay_type, "Group" = group, "LFC_file" = file, "P_file" = p_val_file)
      } else if (str_detect(file_elements[[1]][length(file_elements[[1]])], "_p.txt")) {
        break
      }
    } else if (datatype == "tercen") {
      df <- tibble("Order" = order, "Comparison" = comparison, "Stats" = stats_type, "Assay_Type" = assay_type, "Group" = group, "File" = file)
    }
    dfs[[i]] <- df
  }
  output <- bind_rows(dfs)
  return(output %>% arrange(Assay_Type, Order))
}

extract_direction_tt_tercen <- function(filepath, comparison, assay_type) {
  df <- extract_phosphosite_data_tt_tercen(filepath, assay_type)
  
  if (length(colnames(df)) == 4) {
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
  } else if (length(colnames(df)) == 5) {
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

extract_direction_tt_bionav <- function(lfc_file, p_file, comparison, assay_type) {
  titles <- read_lines(p_file, n_max = 2)
  split_titles <- unlist(str_split(titles[2], "\t"))

if(length(split_titles) == 1) {
  lfc_cols <- c("ID", "Uniprot", "Sequence", "Empty", "LogFC")
    p_cols <- c("ID", "Uniprot", "Sequence", "Empty", "P")

    logfc <- read_delim(lfc_file, skip = 3, col_names = lfc_cols, show_col_types = FALSE)
    p <- read_delim(p_file, skip = 3, col_names = p_cols, show_col_types = FALSE)
    df <- logfc %>% select(ID, LogFC)
    df$p <- p %>%
      select(P) %>%
      pull()
    df <- df %>% mutate(
      LogFC = as.numeric(LogFC),
      p = as.numeric(p)
    )

    if (assay_type == "PTK") {
      df <- df %>%
        drop_na() %>%
        summarise(
          comparison = comparison,
          up_PTK = sum((LogFC > 0) & (p < 0.05)),
          down_PTK = sum((LogFC < 0) & (p < 0.05)),
        )
    } else if (assay_type == "STK") {
      df <- df %>%
        drop_na() %>%
        summarise(
          comparison = comparison,
          up_STK = sum((LogFC > 0) & (p < 0.05)),
          down_STK = sum((LogFC < 0) & (p < 0.05)),
        )
    }
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
        
        comparison_sg <- paste(comparison, title)
         direction <- comp_df %>%
           drop_na() %>%
            summarise(
              comparison = comparison_sg,
              up = sum((LogFC > 0) & (P < 0.05)),
              down = sum((LogFC < 0) & (P < 0.05))
            )
        direction$assay_type <- assay_type
        direction$stats <- "mtvc"
        direction <- direction %>% pivot_wider(names_from = assay_type, values_from = c(up, down))
        df[[length(df) + 1]] <- direction
      }
    }
    df <- bind_rows(df)
  }

  
  df$stats <- "tt"

  return(df)
}

extract_direction_mtvc_tercen <- function(filepath, assay_type, group) {
  df <- extract_phosphosite_data_mtvc_tercen(filepath, assay_type)
  
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

extract_direction_mtvc_bionav <- function(lfc_file, p_file, assay_type, group) {
  titles <- read_lines(p_file, n_max = 2)
  split_titles <- unlist(str_split(titles[2], "\t"))
  p_names <- c("c_cluster", "ID", "Uniprot", "Empty", split_titles[5:length(split_titles)])
  p_names <- p_names[p_names != ""]

  logfc <- read_delim(lfc_file, skip = 3, col_names = p_names, show_col_types = FALSE)
  pvalues <- read_delim(p_file, skip = 3, col_names = p_names, show_col_types = FALSE)

  directions <- list()
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
      colnames(df) <- c("ID", "Log2FoldChange", "pvalue")
      df$Log2FoldChange <- as.numeric(df$Log2FoldChange)
      df$pvalue <- as.numeric(df$pvalue)

      if (is.na(group)) {
        comparison <- paste(title, "vs", ctrl_title)
      } else {
        comparison <- paste(group, title, "vs", ctrl_title)
      }

      direction <- df %>%
        drop_na() %>%
        summarise(
          comparison = comparison,
          up = sum((Log2FoldChange > 0) & (pvalue < 0.05)),
          down = sum((Log2FoldChange < 0) & (pvalue < 0.05))
        )
      direction$assay_type <- assay_type
      direction$stats <- "mtvc"
      direction <- direction %>% pivot_wider(names_from = assay_type, values_from = c(up, down))
      directions[[length(directions) + 1]] <- direction
    }
  }
  directions <- bind_rows(directions)

  return(directions)
}

parse_stats_files <- function(stats_files, datatype = "bionav") {
  # determine first whether study has both PTK and STK
  assay_types <- stats_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()

  dfs <- list()
  if (length(assay_types) == 2) {
    # MTvC first
    if ("MTvC" %in% stats_files$Stats) {
      mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
      groups <- stats_files %>%
        filter(Stats == "MTvC") %>%
        distinct(Group) %>%
        pull()
      if (!is.na(groups)) {
        for (group in groups) {
          group_rows <- mtvc_rows %>% filter(Group == group)
          if (datatype == "bionav") {
            mtvc_ptk <- group_rows %>%
              filter(Assay_Type == "PTK") %>%
              do(extract_direction_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type, .$Group))
            mtvc_stk <- group_rows %>%
              filter(Assay_Type == "STK") %>%
              do(extract_direction_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type, .$Group))
          } else if (datatype == "tercen") {
            mtvc_ptk <- group_rows %>%
              filter(Assay_Type == "PTK") %>%
              do(extract_direction_mtvc_tercen(.$File, .$Assay_Type, .$Group))
            mtvc_stk <- group_rows %>%
              filter(Assay_Type == "STK") %>%
              do(extract_direction_mtvc_tercen(.$File, .$Assay_Type, .$Group))
          }

          df <- left_join(mtvc_ptk, mtvc_stk, by = c("comparison", "stats"))
          dfs <- append(dfs, list(df))
        }
      } else {
        if (datatype == "bionav") {
          mtvc_ptk <- mtvc_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_direction_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type, .$Group))
          mtvc_stk <- mtvc_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_direction_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type, .$Group))
        } else if (datatype == "tercen") {
          mtvc_ptk <- mtvc_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_direction_mtvc_tercen(.$File, .$Assay_Type, .$Group))
          mtvc_stk <- mtvc_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_direction_mtvc_tercen(.$File, .$Assay_Type, .$Group))
        }
        df <- left_join(mtvc_ptk, mtvc_stk, by = c("comparison", "stats"))
        dfs <- append(dfs, list(df))
      }
    }
    if ("TT" %in% stats_files$Stats) {
      ttest_rows <- stats_files %>% filter(Stats != "MTvC")
      comparisons <- stats_files %>%
        filter(Stats != "MTvC") %>%
        distinct(Comparison) %>%
        pull()
      for (comparison in comparisons) {
        c_rows <- ttest_rows %>% filter(Comparison == comparison)
        if (datatype == "bionav") {
          c_ptk <- c_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_direction_tt_bionav(.$LFC_file, .$P_file, .$Comparison, .$Assay_Type))
          c_stk <- c_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_direction_tt_bionav(.$LFC_file, .$P_file, .$Comparison, .$Assay_Type))
        } else if (datatype == "tercen") {
          c_ptk <- c_rows %>%
            filter(Assay_Type == "PTK") %>%
            do(extract_direction_tt_tercen(.$File, .$Comparison, .$Assay_Type))
          c_stk <- c_rows %>%
            filter(Assay_Type == "STK") %>%
            do(extract_direction_tt_tercen(.$File, .$Comparison, .$Assay_Type))
        }

        df <- left_join(c_ptk, c_stk, by = c("comparison", "stats"))
        dfs <- append(dfs, list(df))
      }
    }
  } else if (length(assay_types) == 1) {
    if ("MTvC" %in% stats_files$Stats) {
      mtvc_rows <- stats_files %>% filter(Stats == "MTvC")
      groups <- stats_files %>%
        filter(Stats == "MTvC") %>%
        distinct(Group) %>%
        pull()
      if (!is.na(groups)) {
        for (group in groups) {
          group_rows <- mtvc_rows %>% filter(Group == group)
          if (datatype == "bionav") {
            df <- group_rows %>% do(extract_direction_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type, .$Group))
          } else if (datatype == "tercen") {
            df <- group_rows %>% do(extract_direction_mtvc_tercen(.$File, .$Assay_Type, .$Group))
          }
          dfs <- append(dfs, list(df))
        }
      } else {
        if (datatype == "bionav") {
          df <- mtvc_rows %>% do(extract_direction_mtvc_bionav(.$LFC_file, .$P_file, .$Assay_Type, .$Group))
        } else if (datatype == "tercen") {
          df <- mtvc_rows %>% do(extract_direction_mtvc_tercen(.$File, .$Assay_Type, .$Group))
        }        
        dfs <- append(dfs, list(df))
      }
    }
    if ("TT" %in% stats_files$Stats) {
      ttest_rows <- stats_files %>% filter(Stats != "MTvC")
      comparisons <- stats_files %>%
        filter(Stats != "MTvC") %>%
        distinct(Comparison) %>%
        pull()
      for (comparison in comparisons) {
        c_rows <- ttest_rows %>% filter(Comparison == comparison)
        if (datatype == "bionav") {
          df <- c_rows %>% do(extract_direction_tt_bionav(.$LFC_file, .$P_file, .$Comparison, .$Assay_Type))
        } else if (datatype == "tercen") {
          df <- c_rows %>% do(extract_direction_tt_tercen(.$File, .$Comparison, .$Assay_Type))
        }
        
        dfs <- append(dfs, list(df))
      }
    }
  }


  return(bind_rows(dfs))
}