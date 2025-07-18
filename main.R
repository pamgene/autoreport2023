renv::restore()

if (!dir.exists("99_Saved Plots")){
  dir.create("99_Saved Plots")
}

if (!dir.exists("temp")){
  dir.create("temp")
}

date_str <- format(Sys.Date(), "%y%m%d")
main_report_file <- paste0("01_MainReport_", date_str, ".docx")
supplement_file <- paste0("02_ReportSupplement_", date_str, ".docx")

rmarkdown::render("01_MainReport.Rmd", 
                  params = yaml::read_yaml("./params.yml"),
                  output_file = main_report_file)
rmarkdown::render("02_ReportSupplement.Rmd", 
                  params = yaml::read_yaml("./params.yml"),
                  output_file = supplement_file)

unlink("temp", recursive = TRUE)

