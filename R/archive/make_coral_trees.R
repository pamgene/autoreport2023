require(CORALcli)
require(rsvg)
require(magick)

extract_coral_data <- function(kinase_file, assay_type) {
  df <- read_delim(kinase_file, show_col_types = FALSE)
  df$Assay_type <- assay_type
  return(df)
}

render_single_coral <- function(df, comparison, tree_dir = "99_Saved Plots", ks_min, ks_max) {
  coral_file <- CORALcli::plot_tree(df, comparison, tree_dir = tree_dir, min_col = ks_min, max_col = ks_max)
  return(coral_file)
}

get_uka_data_range <- function(kinase_files) {
  df <- bind_rows(lapply(kinase_files$UKA_file, read_uka_file))
  
  df <- df %>% 
    group_by(Assay_type) %>% 
    summarise(max_mks = max(`Median Kinase Statistic`), min_mks = min(`Median Kinase Statistic`)) %>%
    mutate(max_mks = ceiling(max_mks / 0.5) * 0.5, min_mks = floor(min_mks / 0.5) * 0.5) %>% 
    mutate(max_mks = case_when(abs(min_mks) > max_mks ~ -min_mks,
                               TRUE ~ max_mks),
           min_mks = case_when(abs(min_mks) < max_mks ~ -max_mks,
                               TRUE ~ min_mks))
  return(df)
}

calc_crop <- function(width) {
  new_width <- round(0.85 * width)
  new_height <- round(0.65 * width)
  w_offset <- round(0.05 * width)
  crop_string <- paste0(new_width, "x", new_height, "+", w_offset)

  return(crop_string)
}

crop_coral_tree <- function(file_location, dimension) {
  img_loc <- tools::file_path_sans_ext(basename(file_location))
  img <- magick::image_read_svg(file_location, width = dimension, height = dimension)
  crop_string <- calc_crop(dimension)
  img <- magick::image_crop(img, crop_string)
  img <- magick::image_convert(img, "png")
  png_name <- paste0("temp/", img_loc, ".png")
  magick::image_write(img, path = png_name, format = "png")
  return(png_name)
}

make_coral_trees <- function(kinase_files, kinase_number = 20) {
  # determine range of data
  range_df <- get_uka_data_range(kinase_files)
  ks_max <- range_df %>% pull(max_mks) %>% max(.)
  ks_min <- range_df %>% pull(min_mks) %>% min(.)
  # determine first whether study has both PTK and STK
  assay_types <- kinase_files %>%
    select(Assay_Type) %>%
    unique() %>%
    pull()
  # determine the number of unique comparisons
  comparisons <- kinase_files %>%
    distinct(Comparison) %>%
    pull()
  # then iterate over comparisons, make table for each
  plots <- list()
  for (comparison in comparisons) {
    c_rows <- kinase_files %>% filter(Comparison == comparison)
    if (length(assay_types) == 2) {
      c_ptk <- c_rows %>%
        filter(Assay_Type == "PTK") %>%
        do(extract_coral_data(.$UKA_file, .$Assay_Type))
      c_stk <- c_rows %>%
        filter(Assay_Type == "STK") %>%
        do(extract_coral_data(.$UKA_file, .$Assay_Type))
      df <- bind_rows(c_ptk, c_stk)
    } else if (length(assay_types) == 1) {
      df <- c_rows %>% do(extract_coral_data(.$UKA_file, .$Assay_Type))
    }
    coral_file <- render_single_coral(df, comparison, ks_min = ks_min, ks_max = ks_max)
    plot <- crop_coral_tree(coral_file, 1880)
    plots[[length(plots) + 1]] <- plot
  }
  return(plots)
}