source("R/report_theme.R")
library(RColorBrewer)
library(cowplot)
library(Cairo)
library(extrafont)
set_null_device("cairo")

extract_kinase_plot_data <- function(kinase_file, assay_type, comparison, kinase_number = 20) {
  df <- read_delim(kinase_file, show_col_types = FALSE)
  df <- df %>%
    select(`Kinase Name`, `Mean Specificity Score`, `Median Kinase Statistic`, `Kinase Family`) %>%
    slice_head(n = 20)
  df$`Mean Specificity Score`[df$`Mean Specificity Score` > 2] <- 2
  df$`Mean Specificity Score` <- 0.5 * df$`Mean Specificity Score`
  df <- rownames_to_column(df, var = "Rank")
  df$Assay_type <- assay_type
  return(df)
}

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
    
        p<- p + report_theme
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

read_uka_file <- function(uka_file, assay_type) {
  df <- read_delim(uka_file, show_col_types = FALSE)
  
  file_base <- basename(uka_file)
  file_base <- tools::file_path_sans_ext(file_base)
  file_elements <- str_split(file_base, pattern = "_")
  
  assay_type <- file_elements[[1]][2]
  df$Assay_type <- assay_type
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



make_kinase_plots <- function(kinase_files, kinase_number = 20, color_type = "specificity", xax_scale = "yes") {
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
        do(extract_kinase_plot_data(.$UKA_file, .$Assay_Type, .$Comparison, kinase_number))
      c_stk <- c_rows %>%
        filter(Assay_Type == "STK") %>%
        do(extract_kinase_plot_data(.$UKA_file, .$Assay_Type, .$Comparison, kinase_number))
      df <- bind_rows(c_ptk, c_stk)
    } else if (length(assay_types) == 1) {
      df <- c_rows %>% do(extract_kinase_plot_data(.$UKA_file, .$Assay_Type, .$Comparison, kinase_number))
    }
    plot <- render_kinase_plot(df, range_df, color_type, xax_scale)
    save_plot(paste0("99_Saved Plots/", comparison, "_Score Plot.pdf"), plot, dpi = 300, base_height = 4)
    plots[[length(plots) + 1]] <- plot
  }
  return(plots)
}