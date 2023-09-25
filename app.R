APP_VERSION <- "v1.5.10"

library(shiny)
library(sortable)
library(yaml)
library(lubridate)
library(shinyjs)

source("R/00_GeneralFunctions.R")

ui <- fluidPage(
    useShinyjs(),
    fluidPage(
        sidebarLayout(
            sidebarPanel(
                h1("Report Parameters"),
                checkboxInput("checkparams", "Load Previous Report Parameters"),
                conditionalPanel(
                    condition = "input.checkparams == true",
                    fileInput("params", "Upload params.yml File to Restore Report Parameters"),
                ),
                textInput("author", "Author Name"),
                dateInput("date", "Report Date"),
                textAreaInput("aim", "Project Aim", rows = 3),
                textAreaInput("comparisons", "Experiment Comparisons", rows = 3),
                checkboxGroupInput("normalizations", "Normalizations", choices = c("VSN" = "vsn", "ComBat Correction" = "combat")),
                textInput("qc_cv_factor", "Factor Used for CV Calculation"),
                radioButtons("signal_heatmap", "Include Overall Signal Heatmap Text", choices = c("Yes" = "yes", "No" = "no")),
                checkboxGroupInput("heatmap", "Significant Peptide Heatmap", choices = c("Yes" = "heatmap")),
                checkboxGroupInput("kinase_analysis", "Kinase Analysis", 
                                   choices = c("Score Plot - family" = "splotf", "Score Plot - specificity" = "splots", "Score Table" = "table", "Coral Tree" = "tree")),
                radioButtons(
                  "coral_ks_thrs", "Coral Kinase Statistics thresholds", choices = c("Automatic" = "coral_auto", "Manual" = "coral_man")),
                conditionalPanel("input.coral_ks_thrs == 'coral_man'",
                                 numericInput("coral_min", "Coral KS min", -5, -30, 30),
                                 numericInput("coral_max", "Coral KS max", 5, -30, 30)
                ),

                radioButtons("xax_scale", "Same X axis for all score plots", 
                             choices = c("No" = "no", "Yes" = "yes")),
                fluidRow(
                  actionButton("save", "Save Parameters", class = "btn-lg btn-primary"),
                  disabled(actionButton("knit", "Knit Report", class = "btn-lg btn-success")),
                  br(),
                  uiOutput("downloadParams"),
                  uiOutput("download")
                )
            ),
            mainPanel(
                br(),
                h1("Report Files"),
                radioButtons("datatype", "Report Data Type", choices = c("BioNavigator" = "bionav", "Tercen" = "tercen"), inline = TRUE),
                fileInput("reportFiles", "Upload all the files for the report separately or as a .zip.", multiple = TRUE),
                h2("QC Files"),
                tableOutput("qc_table"),
                h2("Phosphosite Analysis Files"),
                tableOutput("phosphosite_table"),
                h2("Kinase Analysis Files"),
                tableOutput("kinase_table")
                # textOutput("kin_analysis"),
                # uiOutput("kin_analysis_chosen")
            )
        ),
        hr(),
        sprintf("PamGene Automated Report Version %s", APP_VERSION),
        br(),
        br()
    )
)

server <- function(input, output, session) {
  make_data_folders <- function() {
    folders <- c("01_Basic Processing", "02_Phosphosite Analysis", "03_Kinase Analysis", "99_Saved Plots", "unzipped", "output")
    for (folder in folders) {
      if (!dir.exists(folder)) {
        dir.create(folder)
      }
    }
  }

  remove_data_folders <- function() {
    folders <- c("01_Basic Processing", "02_Phosphosite Analysis", "03_Kinase Analysis", "99_Saved Plots", "unzipped", "output")
    for (folder in folders) {
      if (dir.exists(folder)) {
        unlink(folder, recursive = TRUE)
      }
    }
  }

  load_saved_params <- function(params_list) {
    # params_list <- read_yaml("params.yml")

    updateTextInput(session, "author", value = params_list$author)
    updateDateInput(session, "date", value = mdy(params_list$date))
    updateTextAreaInput(session, "aim", value = params_list$aim)
    updateTextAreaInput(session, "comparisons", value = params_list$comparisons)
    updateCheckboxGroupInput(session, "normalizations", selected = params_list$normalizations)
    updateTextInput(session, "qc_cv_factor", value = params_list$`qc_cv_factor`)
    updateCheckboxGroupInput(session, "heatmap", selected = ifelse(is.null(params_list$`phosphosite_heatmap`), character(0), "heatmap"))
    #updateCheckboxGroupInput(session, "kinase_analysis", selected = params_list$`kinase_analysis`)
    updateRadioButtons(session, 'coral_ks_thrs', selected = param_list$`coral_ks_thrs`)
    updateRadioButtons(session, "datatype", selected = params_list$datatype)
  }

  save_params <- function() {
    params <- list(
            "author" = input$author,
            "date" = format(input$date, "%B %d, %Y"),
            "aim" = input$aim,
            "comparisons" = input$comparisons,
            "normalizations" = input$normalizations,
            "qc_cv_factor" = input$`qc_cv_factor`,
            "signal_heatmap" = input$`signal_heatmap`,
            "phosphosite_heatmap" = input$heatmap,
            "kinase_analysis" = input$`kinase_analysis`,
            "coral_ks_thrs" = input$`coral_ks_thrs`,
            "coral_min" = input$`coral_min`,
            "coral_max" = input$`coral_max`,
            "xax_scale" = input$`xax_scale`,
            "datatype" = input$datatype
        )
    write_yaml(params, "params.yml")
  }

  strip_duplicate_extensions <- function(filename) {
    pattern <- "(.csv|.txt){2,}$"
    stripped_name <- str_replace(filename, pattern, "\\1")
    return(stripped_name)
  }

  move_uploaded_files <- function(input_files) {
    for (i in 1:nrow(input_files)) {
      aFile <- input_files[i,]

      stripped_name <- strip_duplicate_extensions(aFile$name)

      if (grepl("^QC", stripped_name)) {
        file.copy(from = aFile$datapath, to = file.path("01_Basic Processing", stripped_name))
      } else if (grepl("^MTvC|TT", stripped_name)) {
        file.copy(from = aFile$datapath, to = file.path("02_Phosphosite Analysis", stripped_name))
      } else if (grepl("^UKA", stripped_name)) {
        file.copy(from = aFile$datapath, to = file.path("03_Kinase Analysis", stripped_name))
      } else {
        print(paste("File not recognized:", aFile$name))
      }
    }
  }

  make_report_zip <- function() {
    file_name <- paste0(format(input$date, "%y%m%d"), "_Report.zip")
    zip_contents <- c("01_MainReport.docx", "02_ReportSupplement.docx", 
                      "99_Saved Plots/", "params.yml")
    zip(paste0("output/", file_name), zip_contents)
    return(file_name)
  }

  observe({
    if (is.null(input$reportFiles)) return()
    remove_data_folders()
    make_data_folders()

    filetypes <- input$reportFiles %>% distinct(type) %>% pull()

    if (filetypes == "text/csv") {
      showNotification("Tercen output detected!", type = "message")
      updateRadioButtons(session, "datatype", selected = "tercen")
    } else if (filetypes == "text/plain") {
      showNotification("BioNavigator output detected!", type = "message")
      updateRadioButtons(session, "datatype", selected = "bionav")
    }

    if (grepl(".zip$", input$reportFiles$name)) {
      files <- unzip(input$reportFiles$datapath, exdir = "unzipped")
      filenames <- basename(files)
      file_df <- data.frame(name = filenames, datapath = files)
      move_uploaded_files(file_df)
    } else {
      move_uploaded_files(input$reportFiles)
    }

    output$download <- renderUI({ NULL })
    output$`qc_table` <- renderTable(read_qc_dir())
    output$`phosphosite_table` <- renderTable(read_phosphosite_dir(datatype = input$datatype))
    output$`kinase_table` <- renderTable(read_kinase_dir())
  })

  observe({
    if (is.null(input$params)) return()

    params_list <- read_yaml(input$params$datapath)
    load_saved_params(params_list)
    showNotification("Parameters loaded!", type = "message")

  })

  observeEvent(input$save, {
    save_params()
    showNotification("Parameters saved.", type = "message")
    output$downloadParams <- renderUI({
      downloadButton("download_params", "Download Report Parameters")
    })
    enable("knit")
  })

  output$download_params <- downloadHandler(
        filename <- function() { "params.yml" },
        content <- function(file) {
          file.copy("params.yml", file)
        }
    )

  observeEvent(input$knit, {
    withProgress(message = "Knitting Report...", {
      source("main.R")
    })
    showNotification("Knit complete!", type = "message")
    zip_name <- make_report_zip()

    output$download <- renderUI({
      downloadButton("download_report", "Download Report.zip")
    })

  })

  output$download_report <- downloadHandler(

        filename = function() {
          paste0(format(input$date, "%y%m%d"), "_Report.zip")
        },

        content <- function(file) {
          file.copy(paste0("output/", format(input$date, "%y%m%d"), "_Report.zip"), file)
        }
    )
}

# Run the application
shinyApp(ui = ui, server = server)