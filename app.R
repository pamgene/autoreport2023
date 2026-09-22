APP_VERSION <- "v1.14"

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
                div(
                  style = "display: flex; align-items: center;",
                  tags$label(style = "margin-right: 10px;", "Report date:"),
                  dateInput("date", label = NULL, width = "160px")
                ),
                textAreaInput("aim", "Project Aim", rows = 2),
                textAreaInput("comparisons", "Experiment Comparisons", rows = 2),
                numericInput("fscore_thr", "Final Score threshold", 1.3, 0, 10, step = 0.1),
                numericInput("spec_thr", "Specificity Score threshold", 0.7, 0, 10, step = 0.1),
                helpText("Score thresholds affect: main report top kinase table, coral tree dotsize, text."),
                numericInput("psite_p_thr", "Phosphosite significance p-value threshold", 0.05, min = 0, max = 1, step = 0.001),
                helpText("Affects: main report phosphosite analysis table, Supplement peptide volcano/heatmap."),
                checkboxGroupInput("kinase_analysis", "Kinase Analysis",
                                   choices = c("Coral Tree" = "tree")),
                helpText("The below outputs are deprecated and should be only used when necessary - not as default!"),
                checkboxGroupInput("kinase_analysis_old", "Deprecated Kinase outputs",
                                   choices = c("Score Table" = "table", "Score Plot - family" = "splotf", "Score Plot - specificity" = "splots")),
                helpText("Coral KS thresholds are derived from the 0.1 and 0.9 percentile of the data for each comparison."),
                radioButtons(
                  "coral_ks_thrs", "Coral Kinase Statistics thresholds", choices = c("Automatic" = "coral_auto", "Manual" = "coral_man")),
                helpText("The same manual thresholds apply for all comparisons."),
                conditionalPanel("input.coral_ks_thrs == 'coral_man'",
                                 numericInput("coral_min", "Coral KS min", -5, -30, 30),
                                 numericInput("coral_max", "Coral KS max", 5, -30, 30)
                ),
                radioButtons("xax_scale", "Same X axis for all score plots",
                             choices = c("No" = "no", "Yes" = "yes")),
                radioButtons(
                  "stk_qc_method", "STK QC method", choices = c("LOD" = "LOD", "Nominal CV" = "nom_cv")),
                checkboxGroupInput("normalizations", "Normalizations-BioNav input", choices = c("VSN" = "vsn", "ComBat Correction" = "combat")),
                helpText("Only used for BioNavigator QC files. For Tercen QC files, normalization is detected automatically from the uploaded files."),
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
        sprintf("PamDx Automated Report Version %s", APP_VERSION),
        br(),
        br()
    )
)

server <- function(input, output, session) {
  # csUKA is a property of the uploaded UKA file (which column-naming convention it
  # uses), not a user preference - detect_csUKA() sets this when kinase files are
  # processed, and save_params() persists it for the Rmds instead of a UI toggle.
  detected_csUKA <- reactiveVal(FALSE)
  # Significant Peptide Heatmap creation likewise follows from whether phosphosite
  # analysis files (Limma/MTvC/TT) were uploaded, not a separate user toggle.
  has_phosphosite_files <- reactiveVal(FALSE)
  # ...and the Overall Signal Heatmap text follows from whether QC files were uploaded.
  has_qc_files <- reactiveVal(FALSE)

  make_data_folders <- function() {
    folders <- c("01_Basic Processing", "02_Phosphosite Analysis", "03_Kinase Analysis", "01_REPORTS", "02_DATA", "03_FIGURES", "unzipped", "output")
    for (folder in folders) {
      if (!dir.exists(folder)) {
        dir.create(folder)
      }
    }
  }

  remove_data_folders <- function() {
    folders <- c("01_Basic Processing", "02_Phosphosite Analysis", "03_Kinase Analysis", "01_REPORTS", "02_DATA", "03_FIGURES", "unzipped", "output")
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
    updateRadioButtons(session, 'stk_qc_method', selected = params_list$`stk_qc_method`)
    updateCheckboxGroupInput(session, "kinase_analysis", selected = params_list$`kinase_analysis`)
    updateCheckboxGroupInput(session, "kinase_analysis_old", selected = params_list$`kinase_analysis_old`)
    updateRadioButtons(session, 'coral_ks_thrs', selected = params_list$`coral_ks_thrs`)
    updateNumericInput(session, "fscore_thr", selected = params_list$`fscore_thr`)
    updateNumericInput(session, "spec_thr", selected = params_list$`spec_thr`)
    psite_p_thr_val <- params_list$`psite_p_thr`
    if (is.null(psite_p_thr_val)) psite_p_thr_val <- 0.05
    updateNumericInput(session, "psite_p_thr", value = psite_p_thr_val)
    updateRadioButtons(session, "datatype", selected = params_list$datatype)
  }

  save_params <- function() {
    # Keep the phosphosite p-value threshold inside a valid p-value range, at most 3 decimals.
    psite_p_thr <- suppressWarnings(as.numeric(input$`psite_p_thr`))
    if (is.na(psite_p_thr)) psite_p_thr <- 0.05
    psite_p_thr <- round(min(max(psite_p_thr, 0), 1), 3)
    updateNumericInput(session, "psite_p_thr", value = psite_p_thr)

    params <- list(
            "author" = input$author,
            "date" = format(input$date, "%B %d, %Y"),
            "aim" = input$aim,
            "comparisons" = input$comparisons,
            "normalizations" = input$normalizations,
            "stk_qc_method" = input$`stk_qc_method`,
            "signal_heatmap" = if (has_qc_files()) "yes" else "no",
            "phosphosite_heatmap" = if (has_phosphosite_files()) "heatmap" else character(0),
            "kinase_analysis" = input$`kinase_analysis`,
            "kinase_analysis_old" = input$`kinase_analysis_old`,
            "csUKA" = detected_csUKA(),
            "fscore_thr" = input$`fscore_thr`,
            "spec_thr" = input$`spec_thr`,
            "psite_p_thr" = psite_p_thr,
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
      } else if (grepl("^MTvC|TT|Limma", stripped_name)) {
        file.copy(from = aFile$datapath, to = file.path("02_Phosphosite Analysis", stripped_name))
      } else if (grepl("^UKA", stripped_name)) {
        file.copy(from = aFile$datapath, to = file.path("03_Kinase Analysis", stripped_name))
      } else {
        print(paste("File not recognized:", aFile$name))
      }
    }
  }

  make_report_zip <- function() {
    date_str <- format(input$date, "%y%m%d")
    file_name <- paste0(date_str, "_Report.zip")
    zip_contents <- c(
      "01_REPORTS/",
      "02_DATA/", 
      "03_FIGURES/",
      "params.yml"
    )
    zip(paste0("output/", file_name), zip_contents)
    return(file_name)
  }

  observe({
    if (is.null(input$reportFiles)) return()
    remove_data_folders()
    make_data_folders()

    filetypes <- input$reportFiles %>% distinct(type) %>% pull()

    # updateRadioButtons() only queues a client message - input$datatype won't reflect
    # it until a round-trip with the browser completes, so code later in this same
    # observer tick must use current_datatype, not input$datatype, to see the corrected
    # value immediately.
    current_datatype <- input$datatype
    if (filetypes == "text/csv") {
      showNotification("Tercen output detected!", type = "message")
      updateRadioButtons(session, "datatype", selected = "tercen")
      current_datatype <- "tercen"
    } else if (filetypes == "text/plain") {
      showNotification("BioNavigator output detected!", type = "message")
      updateRadioButtons(session, "datatype", selected = "bionav")
      current_datatype <- "bionav"
    }

    if (length(input$reportFiles$name) == 1 && grepl(".zip$", input$reportFiles$name)) {
      files <- unzip(input$reportFiles$datapath, exdir = "unzipped")
      filenames <- basename(files)
      file_df <- data.frame(name = filenames, datapath = files)
      move_uploaded_files(file_df)
    } else {
      move_uploaded_files(input$reportFiles)
    }

    output$download <- renderUI({ NULL })
    qc_table <- read_qc_dir()
    has_qc_files(nrow(qc_table) > 0)
    output$`qc_table` <- renderTable(qc_table)
    phosphosite_table <- read_phosphosite_dir(datatype = current_datatype)
    has_phosphosite_files(nrow(phosphosite_table) > 0)
    output$`phosphosite_table` <- renderTable(phosphosite_table)

    # Make kinase table and save it
    if (!dir.exists("temp")){
      dir.create("temp")
    }
    
    kinase_table <- read_kinase_dir()
    if (!is.null(attr(kinase_table, "csUKA"))) {
      detected_csUKA(attr(kinase_table, "csUKA"))
    }

    if (nrow(kinase_table) > 0) {
      write_csv(kinase_table, file = "temp/kinase_files.csv")
      output$`kinase_table` <- renderTable(kinase_table)
    } else {
      write_csv(data.frame(), file = "temp/kinase_files.csv")
    }
    
    
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