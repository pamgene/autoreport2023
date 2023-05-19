library(tidyverse)

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
