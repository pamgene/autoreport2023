require(tidyverse)
source("R/report_theme.R")
require(RColorBrewer)
require(cowplot)
require(Cairo)
require(extrafont)
require(flextable)
require(officer)
require(CORALcli)
require(rsvg)
require(magick)
set_null_device("cairo")

#########################################################
################### LOAD KINASE DATA ####################
#########################################################

read_uka_file <- function(uka_file, assay_type) {
  df <- read_delim(uka_file, show_col_types = FALSE)

  file_base <- basename(uka_file)
  file_base <- tools::file_path_sans_ext(file_base)
  file_elements <- str_split(file_base, pattern = "_")

  assay_type <- file_elements[[1]][2]
  df$Assay_type <- assay_type
  return(df)
}

extract_top_kinases <- function(uka_file, assay_type, comparison) {
  uka <- read_uka_file(uka_file)
  uka <- uka %>%
    filter(`Median Final score` > 1.3) %>%
    slice_head(n = 10)
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

extract_kinase_plot_data <- function(kinase_file, assay_type, comparison) {
  df <- read_uka_file(kinase_file)
  df <- df %>%
    select(`Kinase Name`, `Mean Specificity Score`, `Median Kinase Statistic`, `Kinase Family`) %>%
    slice_head(n = 20)
  df$`Mean Specificity Score`[df$`Mean Specificity Score` > 2] <- 2
  df$`Mean Specificity Score` <- 0.5 * df$`Mean Specificity Score`
  df <- rownames_to_column(df, var = "Rank")
  df$Assay_type <- assay_type
  return(df)
}

extract_kinase_table <- function(kinase_file, assay_type) {
  df <- read_uka_file(kinase_file)
  df <- df %>%
    select(`Kinase Name`, `Median Final score`, `Median Kinase Statistic`) %>%
    #filter(`Median Final score` > 1.3) %>%
    slice_head(n = 10) %>%
    mutate(
      `Median Final score` = round(`Median Final score`, digits = 2),
      `Median Kinase Statistic` = round(`Median Kinase Statistic`, digits = 2)
    )
  df <- rownames_to_column(df, var = "Rank")
  colnames(df) <- c("Rank", paste(assay_type, "Name", sep = "_"), paste(assay_type, "Score", sep = "_"), paste(assay_type, "Statistic", sep = "_"))

  return(df)
}

extract_coral_data <- function(kinase_file, assay_type) {
  df <- read_uka_file(kinase_file)
  return(df)
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

#########################################################
############## STRINGIFY KINASE DIRECTION ###############
#########################################################


parse_kinase_files <- function(kinase_files) {
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
        do(extract_top_kinases(.$UKA_file, .$Assay_Type, .$Comparison))
      c_stk <- c_rows %>%
        filter(Assay_Type == "STK") %>%
        do(extract_top_kinases(.$UKA_file, .$Assay_Type, .$Comparison))
      df <- left_join(c_ptk, c_stk, by = "comparison")
      dfs <- append(dfs, list(df))
    } else if (length(assay_types) == 1) {
      c_rows <- kinase_files %>% filter(Comparison == comparison)
      df <- c_rows %>% do(extract_top_kinases(.$UKA_file, .$Assay_Type, .$Comparison))
      dfs <- append(dfs, list(df))
    }
  }
  return(bind_rows(dfs))
}

make_report_table <- function(res_table, caption) {
  # identify first which assay types are included
  assay_types <- str_extract(colnames(res_table), ".TK") %>%
    na.omit() %>%
    unique()
  header_row <- c("Assay Type", assay_types)
  if (length(assay_types) == 2) {
    header_widths <- c(1, 2, 2)
    ft <- flextable(res_table) %>%
      set_header_labels(comparison = "Comparisons", up_PTK = "Up", down_PTK = "Down", up_STK = "Up", down_STK = "Down") %>%
      color(j = c("up_PTK", "up_STK"), color = "red", part = "body") %>%
      color(j = c("down_PTK", "down_STK"), color = "blue", part = "body") %>%
      add_header_row(values = header_row, colwidths = header_widths)
  } else if (length(assay_types) == 1) {
    header_widths <- c(1, 2)

    if (assay_types[1] == "PTK") {
      ft <- flextable(res_table) %>%
        set_header_labels(comparison = "Comparisons", up_PTK = "Up", down_PTK = "Down") %>%
        color(j = c("up_PTK"), color = "red", part = "body") %>%
        color(j = c("down_PTK"), color = "blue", part = "body") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    } else if (assay_types[1] == "STK") {
      ft <- flextable(res_table) %>%
        set_header_labels(comparison = "Comparisons", up_STK = "Up", down_STK = "Down") %>%
        color(j = c("up_STK"), color = "red", part = "body") %>%
        color(j = c("down_STK"), color = "blue", part = "body") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    }
  }
  ft %>%
    autofit() %>%
    theme_box() %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(width = 1, layout = "autofit") %>%
    set_caption(caption)
}

#########################################################
################## KINASE SCORE PLOTS ###################
#########################################################

render_kinase_plot <- function(df, range_df, color_type = "specificity", xax_scale = "yes") {
  df$Assay_type <- as.factor(df$Assay_type)

  assay_types <- df %>%
    select(Assay_type) %>%
    unique() %>%
    pull()

  if (color_type == "family") {
    df$KinaseFamilyClr <- df$`Kinase Family`
    for (family in unique(df$`Kinase Family`)) {
      if (sum(str_detect(df$`Kinase Family`, family)) == 1) {
        df <- df %>% mutate(KinaseFamilyClr = replace(KinaseFamilyClr, `Kinase Family` == family, NA))
      }
    }
    df$KinaseFamilyClr <- as.factor(df$KinaseFamilyClr)



    plotlist <- list()
    for (i in seq_along(assay_types)) {
      assay_type <- assay_types[i]

      ks_min <- range_df %>% filter(Assay_type == assay_type) %>% pull(min_mks)
      ks_max <- range_df %>% filter(Assay_type == assay_type) %>% pull(max_mks)

      p <- df %>%
        filter(Assay_type == assay_type) %>%
        ggplot(., aes(x = `Median Kinase Statistic`, y = reorder(`Kinase Name`, desc(as.integer(`Rank`))))) +
        geom_bar(stat = "identity", aes(fill = KinaseFamilyClr)) +
        scale_color_brewer(palette = "Paired", na.value = brewer.pal(8, "Dark2")[c(8)], breaks = levels(df$KinaseFamilyClr)) +
        scale_fill_brewer(palette = "Paired", na.value = brewer.pal(8, "Dark2")[c(8)], breaks = levels(df$KinaseFamilyClr)) +
        guides(fill = guide_legend(title = "Kinase Family")) +
        labs(y = "Kinase Top List")

      if (xax_scale == "yes") {
        p <- p + xlim(ks_min, ks_max)
      }

      p <- p + report_theme
      plotlist[[i]] <- p
    }
    pg <- plot_grid(plotlist = plotlist, labels = assay_types, nrow = 1)
  }

  if (color_type == "specificity") {
    plotlist <- list()

    for (i in seq_along(assay_types)) {
      assay_type <- assay_types[i]

      ks_min <- range_df %>% filter(Assay_type == assay_type) %>% pull(min_mks)
      ks_max <- range_df %>% filter(Assay_type == assay_type) %>% pull(max_mks)

      p <- df %>%
        filter(Assay_type == assay_type) %>%
        ggplot(., aes(x = `Median Kinase Statistic`, y = reorder(`Kinase Name`, desc(as.integer(`Rank`))))) +
        geom_bar(stat = "identity", aes(fill = `Mean Specificity Score`)) +
        scale_fill_gradientn(
          name = "Mean\nSpecificity\nScore", space = "Lab",
          colours = c("black", "grey", "red"),
          values = c(0, 0.65, 1),
          breaks = c(0, 0.65, 1),
          limits = c(0, 1),
          labels = c(0, 0.65, 1) * 2
        ) +
        labs(y = "Kinase Top List")


      if (xax_scale == "yes") {
        p <- p + xlim(ks_min, ks_max)
      }

      p <- p + report_theme

      plotlist[[i]] <- p
    }
    pg <- plot_grid(plotlist = plotlist, labels = assay_types, nrow = 1)
  }
  return(pg)
}

make_kinase_plots <- function(kinase_files, color_type = "specificity", xax_scale = "yes") {
  # determine range of data
  range_df <- get_uka_data_range(kinase_files)
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
        do(extract_kinase_plot_data(.$UKA_file, .$Assay_Type, .$Comparison))
      c_stk <- c_rows %>%
        filter(Assay_Type == "STK") %>%
        do(extract_kinase_plot_data(.$UKA_file, .$Assay_Type, .$Comparison))
      df <- bind_rows(c_ptk, c_stk)
    } else if (length(assay_types) == 1) {
      df <- c_rows %>% do(extract_kinase_plot_data(.$UKA_file, .$Assay_Type, .$Comparison))
    }
    plot <- render_kinase_plot(df, range_df, color_type, xax_scale)
    save_plot(paste0("99_Saved Plots/", "UKA_Score Plot_", color_type, "_", comparison, ".pdf"), plot, dpi = 300, base_height = 4)
    plots[[length(plots) + 1]] <- plot
  }
  return(plots)
}

#########################################################
################## KINASE SCORE TABLES ##################
#########################################################

render_kinase_table <- function(res_table, assay_types, caption) {
  small_border <- fp_border(color = "gray", width = 1)
  header_row <- c("", assay_types)

  if (length(assay_types) == 2) {
    header_widths <- c(1, 3, 3)
    ft <- flextable(res_table) %>%
      set_header_labels(Rank = "Rank", PTK_Name = "Kinase Name", PTK_Score = "Kinase Score", PTK_Statistic = "Kinase Statistic", STK_Name = "Kinase Name", STK_Score = "Kinase Score", STK_Statistic = "Kinase Statistic") %>%
      add_header_row(values = header_row, colwidths = header_widths)
  } else if (length(assay_types) == 1) {
    header_widths <- c(1, 3)
    if (assay_types[1] == "PTK") {
      ft <- flextable(res_table) %>%
        set_header_labels(Rank = "Rank", PTK_Name = "Kinase Name", PTK_Score = "Kinase Score", PTK_Statistic = "Kinase Statistic") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    } else if (assay_types[1] == "STK") {
      ft <- flextable(res_table) %>%
        set_header_labels(Rank = "Rank", STK_Name = "Kinase Name", STK_Score = "Kinase Score", STK_Statistic = "Kinase Statistic") %>%
        add_header_row(values = header_row, colwidths = header_widths)
    }
  }
  ft <- ft %>%
    autofit() %>%
    theme_kin_table()

  if (length(assay_types) == 2) {
    ft <- ft %>%
      vline(j = "PTK_Statistic", border = small_border)
  }
  ft <- ft %>%
    vline(j = "Rank", border = small_border) %>%
    fontsize(size = 7, part = "all") %>%
    font(part = "all", fontname = "Arial") %>%
    align(align = "center", part = "all") %>%
    set_table_properties(layout = "autofit") %>%
    set_caption(caption)
  return(ft)
}

make_kinase_tables <- function(kinase_files) {
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
  fts <- list()
  for (comparison in comparisons) {
    c_rows <- kinase_files %>% filter(Comparison == comparison)
    if (length(assay_types) == 2) {
      c_ptk <- c_rows %>%
        filter(Assay_Type == "PTK") %>%
        do(extract_kinase_table(.$UKA_file, .$Assay_Type))
      c_stk <- c_rows %>%
        filter(Assay_Type == "STK") %>%
        do(extract_kinase_table(.$UKA_file, .$Assay_Type))
      df <- left_join(c_ptk, c_stk, by = "Rank")
    } else if (length(assay_types) == 1) {
      df <- c_rows %>% do(extract_kinase_table(.$UKA_file, .$Assay_Type))
    }
    ft <- render_kinase_table(df, assay_types, paste("Kinase Score Table", comparison))
    fts <- append(fts, list(ft))
  }
  return(fts)
}

#########################################################
################## KINASE CORAL TREES ###################
#########################################################

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

make_coral_trees <- function(kinase_files, coral_ks_thrs, coral_min, coral_max) {
  if (coral_ks_thrs == "coral_auto"){
    # determine range of data
    range_df <- get_uka_data_range(kinase_files)
    ks_max <- range_df %>% pull(max_mks) %>% max(.)
    ks_min <- range_df %>% pull(min_mks) %>% min(.)
  } else if (coral_ks_thrs == "coral_man"){
    ks_max <- coral_max
    ks_min <- coral_min
  }
  
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
    coral_filename <- paste0("UKA_", comparison)
    coral_file <- render_single_coral(df, comparison = coral_filename, ks_min = ks_min, ks_max = ks_max)
    plot <- crop_coral_tree(coral_file, 1880)
    plots[[length(plots) + 1]] <- plot
  }
  return(plots)
}

#########################################################
################ KINASE ANALYSIS OUTPUT #################
#########################################################

output_kinase_analysis <- function(kinase_files, kin_params, xax_scale, coral_ks_thrs, coral_min, coral_max) {
  temp <- list()
  if ("table" %in% kin_params) {
    tables <- make_kinase_tables(kinase_files)
    temp["table"] <- list(tables)
  }
  if ("splotf" %in% kin_params) {
    temp["splotf"] <- list(make_kinase_plots(kinase_files, color_type = "family", xax_scale = xax_scale))
  }
  if ("splots" %in% kin_params) {
    temp["splots"] <- list(make_kinase_plots(kinase_files, color_type = "specificity", xax_scale = xax_scale))
  }

  if ("tree" %in% kin_params) {
    temp["tree"] <- list(make_coral_trees(kinase_files, coral_ks_thrs, coral_min, coral_max))
  }

  comparisons <- kinase_files %>%
    distinct(Comparison) %>%
    pull()
  for (i in seq_along(comparisons)) {
    for (vis in kin_params) {
      item <- temp[vis][[1]][[i]]
      if (vis == "table") {
        flextable_to_rmd(item)
      } else if (vis == "splotf") {
        cat(paste0("**Kinase Score Plot (family) - ", comparisons[i], "** \n"))
        print(item)
      } else if (vis == "splots") {
        cat(paste0("**Kinase Score Plot (specificity) - ", comparisons[i], "** \n"))
        print(item)
      } else {
        cat(paste0("**Kinase Coral Tree - ", comparisons[i], "** \n\n"))
        cat(paste0("![](", item, ")"))
      }
    }
    cat("\n \n")
    cat("#########")
    cat("\n \n")
  }
}