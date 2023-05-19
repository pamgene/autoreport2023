read_kinase_dir <- function(folder = "03_Kinase Analysis/") {
  # will return df with files for easy processing
  files <- list.files(folder, pattern = ".txt$|.csv$", full.names = TRUE, include.dirs = FALSE)

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

extract_top_kinases <- function(uka_file, assay_type, comparison, number_kinases) {
  uka <- read_delim(uka_file, show_col_types = FALSE)
  uka <- uka %>%
    filter(`Median Final score` > 1.2) %>%
    slice_head(n = number_kinases)
  uka <- uka %>% select(`Kinase Name`, `Median Kinase Statistic`)
  down <- uka %>%
    filter(`Median Kinase Statistic` < 0) %>%
    pull(`Kinase Name`)
  up <- uka %>%
    filter(`Median Kinase Statistic` > 0) %>%
    pull(`Kinase Name`)

  down <- paste(down, collapse = ", ")
  up <- paste(up, collapse = ", ")

  if (assay_type == "PTK") {
    df <- tibble("comparison" = comparison, "up_PTK" = up, "down_PTK" = down)
  } else if (assay_type == "STK") {
    df <- tibble("comparison" = comparison, "up_STK" = up, "down_STK" = down)
  }
  return(df)
}

parse_kinase_files <- function(kinase_files, kinase_number = 10) {
  # determine first whether study has both PTK and STK
  assay_types <- kinase_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()

  # then get all unique comparisons
  comparisons <- kinase_files %>%
    distinct(Comparison) %>%
    pull()

  dfs <- list()
  for (comparison in comparisons) {
    if (length(assay_types) == 2) {
      c_rows <- kinase_files %>% filter(Comparison == comparison)
      c_ptk <- c_rows %>%
        filter(Assay_Type == "PTK") %>%
        do(extract_top_kinases(.$UKA_file, .$Assay_Type, .$Comparison, kinase_number))
      c_stk <- c_rows %>%
        filter(Assay_Type == "STK") %>%
        do(extract_top_kinases(.$UKA_file, .$Assay_Type, .$Comparison, kinase_number))
      df <- left_join(c_ptk, c_stk, by = "comparison")
      dfs <- append(dfs, list(df))
    } else if (length(assay_types) == 1) {
      c_rows <- kinase_files %>% filter(Comparison == comparison)
      df <- c_rows %>% do(extract_top_kinases(.$UKA_file, .$Assay_Type, .$Comparison, kinase_number))
      dfs <- append(dfs, list(df))
    }
  }
  return(bind_rows(dfs))
}