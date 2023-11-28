require(tidyverse)

clean_tercen_columns <- function(df) {
  cols <- colnames(df)
  split <- str_split(cols, pattern = "\\.")
  cols <- sapply(split, tail, 1)
  colnames(df) <- cols
  return(df)
}

read_qc_dir <- function(folder = "01_Basic Processing/") {
  # will return df with files for easy processing
  files <- list.files(folder, pattern = ".txt$|.csv$", full.names = TRUE, include.dirs = FALSE)

  if (length(files) == 0) {
    warning("No QC files")
    return(data.frame())
  }

  dfs <- list()
  for (i in seq_along(files)) {
    file <- files[i]

    file_base <- basename(file)
    file_base <- tools::file_path_sans_ext(file_base)
    file_elements <- str_split(file_base, pattern = "_")

    assay_type <- file_elements[[1]][2]

    df <- tibble("Assay_Type" = assay_type, "qc_file" = file)
    dfs[[i]] <- df
  }
  output <- bind_rows(dfs)
  output$Assay_Type <- as.factor(output$Assay_Type)
  return(output)
}

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

read_kinase_dir <- function(folder = "03_Kinase Analysis/") {
  # will return df with files for easy processing
  files <- list.files(folder, pattern = ".txt$|.csv$", full.names = TRUE, include.dirs = FALSE)
  
  # decide if it is new uka
  # if it is new uka, make the same dataformat as old uka, save and read the files again
  if (any(grepl("_ukam-|_ukat-", files))){
    files_to_process <- list.files(folder, pattern = "_uka[m|t]-", full.names = TRUE, include.dirs = FALSE)
    counter <- 0
    for (f in files_to_process){
      uka <- read_delim(f, show_col_types = FALSE) %>% clean_tercen_columns()
      colnames(uka)[1] <- "Comparison"
      uka <- uka %>% arrange(-`Median Final score`)
      ctrl <- sub(".*vs(.+)\\.csv", "\\1", f) %>% str_trim()
      test <- sub(".*uka[m|t]-(.+)vs.*", "\\1", f) %>% str_trim()
      comparison <- sub(".*uka[m|t]-(.+ vs .+)\\.csv", "\\1", f)
      assay_type <- sub(".*UKA_(.TK)_.*", "\\1", f)
      ukatype <- sub(".*[0-9]{2,}_(uka.)-.*", "\\1", f)
      
      split_uka <- split(uka, uka[,1])
      # pivot wider & clean columns
      #split_uka_w <- lapply(split_uka, function(x) pivot_wider(x, names_from = variable, values_from = value))
      split_uka_w <- lapply(split_uka, function(x) clean_tercen_columns(x) %>% select(-Comparison))
      # save
      sapply(names(split_uka_w), 
             function (x){
               counter <<- counter + 1
               filename <- ifelse(ukatype == "ukam", 
                                  paste0(folder, "/UKA_", assay_type, "_", sprintf("%02d",counter), "_", test, " - ", x, " vs ", ctrl, ".csv"),
                                  paste0(folder, "/UKA_", assay_type, "_", sprintf("%02d",counter), "_", x, " - ", comparison, ".csv"))
               # save to the 03_Kinase Analysis folder to read again
               write_csv(split_uka_w[[x]], filename)
               filename_to_report <- ifelse(ukatype == "ukam", 
                                            paste0("99_Saved Plots/UKA_", assay_type, "_", sprintf("%02d",counter), "_", test, " - ", x, " vs ", ctrl, ".csv"),
                                            paste0("99_Saved Plots/UKA_", assay_type, "_", sprintf("%02d",counter), "_", x, " - ", comparison, ".csv"))
               # save to output folder
               write_csv(split_uka_w[[x]], filename_to_report)
             })
      counter <- counter
    }
    files <- list.files(folder, pattern = ".txt$|.csv$", full.names = TRUE, include.dirs = FALSE)
    files <- files[!files %in% files_to_process]
  }

  
  
  
  if (length(files) == 0) {
    warning("No kinase files")
    return(data.frame())
  }

  dfs <- list()
  for (i in seq_along(files)) {
    file <- files[i]

    file_base <- basename(file)
    file_base <- tools::file_path_sans_ext(file_base)
    file_elements <- str_split(file_base, pattern = "_")

    assay_type <- file_elements[[1]][2]
    order <- as.integer(file_elements[[1]][3])
    comparison <- file_elements[[1]][4]
    comparison <- str_replace(comparison, "\\Bvs\\B", " vs ")

    df <- tibble("Order" = order, "Comparison" = comparison, "Assay_Type" = assay_type, "UKA_file" = file)
    dfs[[i]] <- df
  }
  output <- bind_rows(dfs)
  return(output %>% arrange(Order))
}