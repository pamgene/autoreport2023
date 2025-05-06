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
      # only multiple MTvCs are allowed
      order <- as.integer(file_elements[[1]][3]) + 1
      group <- file_elements[[1]][4]
      comparison <- "MTvC"
    } else if (stats_type == "TT") {
      order <- as.integer(file_elements[[1]][3]) + 2
      comparison <- file_elements[[1]][4]
      comparison <- str_replace(comparison, "\\Bvs\\B", " vs ")
      group <- NA
    } else if (stats_type == "Limma"){
      order <- as.integer(file_elements[[1]][3])
      comparison <- "Limma"
      group <- file_elements[[1]][4]
    }
    
    if (datatype == "bionav") {
      if (str_detect(file_elements[[1]][length(file_elements[[1]])], "LogFC")) {
        p_val_file <- str_replace(file, "LogFC", "p")
        df <- tibble("Order" = order, "Comparison" = comparison, "Stats" = stats_type, 
                     "Assay_Type" = assay_type, "Group" = group, "LFC_file" = file, "P_file" = p_val_file)
      } else if (str_detect(file_elements[[1]][length(file_elements[[1]])], "_p.txt")) {
        break
      }
    } else if (datatype == "tercen") {
      df <- tibble("Order" = order, "Comparison" = comparison, "Stats" = stats_type, 
                   "Assay_Type" = assay_type, "Group" = group, "File" = file)
    }
    dfs[[i]] <- df
  }
  output <- bind_rows(dfs)
  return(output %>% arrange(Assay_Type, Order))
}

process_ukat_ukam <- function(folder){
  # To process UKA files from UKA_MTvC app or UKA_TGC_app
  
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
  files <- list.files(folder, pattern = ".csv$", full.names = TRUE, include.dirs = FALSE)
  files <- files[!files %in% files_to_process]
  return(files)
}


process_uka_allvsall <- function(files, folder, counter = 0) {
  # To process UKA files from UKA_app that gives a Sgroup_contrast column
  files_to_process <- files
  
  for (f in files_to_process) {
    uka <- read_delim(f, show_col_types = FALSE) %>% clean_tercen_columns()
    # write this cleaned version to the output of 99_Saved_plots
    write_csv(uka, paste0("99_Saved Plots/", basename(f) ))
    if ("Sgroup_contrast" %in% colnames(uka)) {
      uka <- uka %>% dplyr::rename('Comparison' = 'Sgroup_contrast')
    }
    
    uka <- uka %>% arrange(Comparison, -`Median Final score`)
    uka$Comparison <- sub("_", " - ", uka$Comparison)
    assay_type <- sub(".*UKA_(.TK)_.*", "\\1", f)
    
    split_uka <- split(uka, uka[, 'Comparison'])
    split_uka_clean <- lapply(split_uka, clean_tercen_columns)
    
    
    
    # Save each split file
    for (x in names(split_uka_clean)) {
      counter <- counter + 1
      filename <- paste0("03_Kinase Analysis/UKA_", assay_type, "_", sprintf("%02d", counter), "_", x, ".csv")
      write_csv(split_uka_clean[[x]], filename)
    }
  }
  
  # List all CSV files in the folder, excluding the ones that were processed
  files <- list.files(folder, pattern = ".csv$", full.names = TRUE, include.dirs = FALSE)
  original_files <- files[files %in% files_to_process] 
  files <- files[!files %in% original_files]
  # delete original files
  file.remove(original_files)
  
  return(files)
}




read_kinase_dir <- function(folder = "03_Kinase Analysis/") {
  # will return df with files for easy processing
  files <- list.files(folder, pattern = ".txt$|.csv$", full.names = TRUE, include.dirs = FALSE)
  if (length(files) == 0) {
    warning("No kinase files")
    return(data.frame())
  }
  
  # if it is UKA from the UKA_MTvC or TGC app, make the same dataformat as old uka, save and read the files again
  if (any(grepl("_ukam-|_ukat-", files))){
    files <- process_ukat_ukam(folder = folder)
  }
  
  # if it is uka all vs all comparison, make the same dataformat as old uka, save and read the files again
  uka_example <- read_delim(files[1], show_col_types = FALSE) %>% clean_tercen_columns()
  if ("contrast" %in% colnames(uka_example)){
    files <- process_uka_allvsall(files = files, folder = folder)
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