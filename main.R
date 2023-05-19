renv::restore()

if (!dir.exists("99_Saved Plots")){
  dir.create("99_Saved Plots")
}

if (!dir.exists("temp")){
  dir.create("temp")
}

rmarkdown::render("01_MainReport.Rmd", params = yaml::read_yaml("./params.yml"))
rmarkdown::render("02_ReportSupplement.Rmd", params = yaml::read_yaml("./params.yml"))

unlink("temp", recursive = TRUE)

