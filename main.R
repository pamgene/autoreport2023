renv::restore()

if (!dir.exists("01_REPORTS")){
  dir.create("01_REPORTS")
}

if (!dir.exists("02_DATA")){
  dir.create("02_DATA")
}

if (!dir.exists("03_FIGURES")){
  dir.create("03_FIGURES")
}

if (!dir.exists("temp")){
  dir.create("temp")
}

date_str <- format(Sys.Date(), "%y%m%d")
main_report_file <- paste0("01_REPORTS/01_MainReport_PamDx_", date_str, ".docx")
supplement_file <- paste0("01_REPORTS/02_ReportSupplement_PamDx_", date_str, ".docx")

rmarkdown::render("01_MainReport.Rmd", 
                  params = yaml::read_yaml("./params.yml"),
                  output_file = main_report_file)
rmarkdown::render("02_ReportSupplement.Rmd", 
                  params = yaml::read_yaml("./params.yml"),
                  output_file = supplement_file)

unlink("temp", recursive = TRUE)

