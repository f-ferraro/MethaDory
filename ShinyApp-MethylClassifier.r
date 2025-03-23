options(shiny.maxRequestSize=250*1024^2)
# app.R
library(shiny)
library(shinydashboard)
library(DT)
library(ggplotify)
library(umap)
library(patchwork)
library(PCAtools)
library(data.table)
library(tidyverse)
library(caret)
library(BiocParallel)
library(ComplexHeatmap)
library(circlize)
library(Rtsne)
library(shinyFiles)
library(fs)

paste0(getwd(),"/ChAMPdata_2.38.0.tar.gz")
install.packages(paste0(getwd(),"/ChAMPdata_2.38.0.tar.gz"),
                 repos = NULL, type="source")


library(methyLImp2)
source("ShinyApp-MethylClassifier-helper_functions.R")

ui <- dashboardPage(
  dashboardHeader(title = "MethaDory"),
  dashboardSidebar(
    
    shinyDirButton("modelDir", "Select Model Folder", "Please select folder containing SVM models"),
    verbatimTextOutput("modelDirPath"),
    hr(),
    fileInput("dataFile", "Upload Test Data (.tsv)",
              accept = c(".tsv")),
    actionButton("loadData", "Load Data", class = "btn-primary"),
    hr(),
    selectInput("proband", "Select Proband for plotting",
                choices = NULL),
    hr(),

    fluidRow(
      column(6, actionButton("selectAll", "Select All")),
      column(6, actionButton("deselectAll", "Deselect All"))
    ),
    checkboxGroupInput("signatures", "Select Signatures for plotting",
                       choices = NULL),
    hr(),
    downloadButton("downloadResults", "Download Results"),
    downloadButton("downloadPlots", "Download Plots")
  ),
  dashboardBody(
    fluidRow(
      tabBox(
        id = "tabset1", height = "600px", width = 30, 
        tabPanel("Help",  includeMarkdown("help.md") ),
       
        tabPanel("Prediction results plots", 
                 div(style = "overflow-x: auto; width: 100%;",
                     div(style = "min-width: 1000px;", 
                         plotOutput("predictionPlot", height = "550px", width = "3000px") 
                     )
                 )
        ),
        tabPanel("Prediction table", DTOutput("predictionTable")),
        tabPanel("Cell proportion deconvolution", 
                 plotOutput("cellpropPlot", height = "550px")),
        tabPanel("Console Output", verbatimTextOutput("console"))
        
      )
    ),
    fluidRow(
      box(
        title = "Dimension Reduction Plots",
        width = 12,
        # Pagination controls
        fluidRow(
          column(4, numericInput("plotsPerPage", "Plots per page:", 2, min = 1, max = 10)),
          column(4, uiOutput("pageControls"))
        ),
        hr(),
        uiOutput("dimReductionPlots")
      )
    )
  )
)

server <- function(input, output, session) {
  # Initialize shinyDirChoose with current working directory
  volumes <- c(Home = fs::path_home(),
               "Working Directory" = getwd(),
               "R Installation" = R.home(),
               getVolumes()())
  
  shinyDirChoose(input, "modelDir",
                 roots = volumes,
                 defaultRoot = "Working Directory",
                 defaultPath = getwd(),
                 session = session)
  
  # Reactive values
  values <- reactiveValues(
    model_dir = NULL,
    data_list = NULL,
    background_data = NULL,
    background_data_cells = NULL,
    imputed_data = NULL,
    inference_data = NULL,
    results = NULL,
    beta_sig_data = NULL,
    real_cases_data = NULL,
    plot_data = NULL,
    plot_metadata = NULL,
    current_page = 1,
    console_messages = character(0)
  )
  
  # Display selected model directory
  output$modelDirPath <- renderText({
    if (!is.null(input$modelDir)) {
      values$model_dir <- parseDirPath(volumes, input$modelDir)
      paste("Selected model directory:", values$model_dir)
    }
  })
  # Console logging function
  log_message <- function(message) {
    values$console_messages <- c(values$console_messages, paste(Sys.time(), "-", message))
  }
  
  # Progress bar
  output$progressBox <- renderUI({
    if (!is.null(input$loadData) && input$loadData > 0) {
      div(
        class = "progress",
        div(class = "progress-bar", 
            style = "width: 0%",
            id = "progressBar",
            "0%"
        )
      )
    }
  })
  
  # Console output
  output$console <- renderPrint({
    cat(paste(rev(values$console_messages), collapse = "\n"))
  })
  
  # Load data when button is clicked
  observeEvent(input$loadData, {
    req(input$dataFile, values$model_dir)
    
    # Validate model directory
    if (!dir.exists(values$model_dir)) {
      log_message("Error: Invalid model directory")
      return()
    }
    
    
    # Create progress object
    withProgress(message = 'Loading data', value = 0, {
      
      # Load test data
      log_message("Loading test data...")
      incProgress(0.1)
      values$data_list <- load_test_data(input$dataFile$datapath)
      # Load background data
      log_message("Loading background data...")
      incProgress(0.2)
      values$background_data <- load_background_data(values$model_dir)
      
      # Prepare and impute data
      log_message("Preparing imputation data...")
      incProgress(0.3)
      imputation_data <- prepare_imputation_data(values$data_list$test_data,
                                                 values$background_data$imputation_background,
                                                 values$background_data$svm)
      
      log_message("Performing imputation...")
      incProgress(0.4)
      
      values$imputed_data <- perform_imputation(imputation_data$test_data,
                                                values$data_list$test_data_ids)
      
      
      # Load additional data for plotting
      log_message("Loading beta signatures...")
      incProgress(0.5)
      values$beta_sig_data <- load_beta_signatures()
      
      # Prepare inference data and make predictions
      log_message("Preparing inference data...")
      incProgress(0.6)
      
      
      values$inference_data <- prepare_inference_data(values$imputed_data, 
                                                      values$background_data$svm)
      
      log_message("Making predictions...")
      incProgress(0.7)
      values$results <- make_predictions(values$inference_data,
                                         values$background_data$svm,
                                         values$data_list$test_data_ids)
      
      
      log_message("Loading real cases...")
      incProgress(0.8)
      values$real_cases_data <- load_real_cases(values$beta_sig_data$signatures)
      
      # Prepare plot data
      log_message("Preparing plot data...")
      incProgress(0.9)
      
      values$plot_data <- prepare_plot_data(values$beta_sig_data$signatures,
                                            values$beta_sig_data$insilico_beta,
                                            values$imputed_data,
                                            values$real_cases_data$real_cases_beta)
      
      values$plot_metadata <- prepare_plot_metadata(values$beta_sig_data$signatures,
                                                    values$plot_data,
                                                    values$imputed_data,
                                                    values$real_cases_data$real_cases_meta)
      
      incProgress(1.0)
      log_message("Data loading complete!")
    })
    
    # Update UI elements
    updateSelectInput(session, "proband",
                      choices = values$data_list$test_data_ids)
    
    updateCheckboxGroupInput(session, "signatures",
                             choices = names(values$beta_sig_data$signatures))
  })
  
  # Pagination controls
  output$pageControls <- renderUI({
    req(input$signatures)
    total_pages <- ceiling(length(input$signatures) / input$plotsPerPage)
    
    if (total_pages > 1) {
      div(
        actionButton("prevPage", "Previous"),
        span(paste("Page", values$current_page, "of", total_pages)),
        actionButton("nextPage", "Next")
      )
    }
  })
  
  # Handle pagination
  observeEvent(input$prevPage, {
    values$current_page <- max(1, values$current_page - 1)
  })
  
  observeEvent(input$nextPage, {
    total_pages <- ceiling(length(input$signatures) / input$plotsPerPage)
    values$current_page <- min(total_pages, values$current_page + 1)
  })
  
  # Handle Select All button
  observeEvent(input$selectAll, {
    req(values$beta_sig_data)
    updateCheckboxGroupInput(session, "signatures",
                             choices = names(values$beta_sig_data$signatures),
                             selected = names(values$beta_sig_data$signatures))
  })
  
  # Handle Deselect All button
  observeEvent(input$deselectAll, {
    req(values$beta_sig_data)
    updateCheckboxGroupInput(session, "signatures",
                             choices = names(values$beta_sig_data$signatures),
                             selected = character(0))
  })
  
  # Create prediction plot
  output$predictionTable <- renderDT({datatable(values$results)})
  
  # Create prediction plot
  output$predictionPlot <- renderPlot({
    req(values$results, input$proband, input$signatures)
    
    create_prediction_plot(values$results, input$proband)
  })
  
  # Modify the cell proportion plot to fit the tab properly
  output$cellpropPlot <- renderPlot({
    req(values$imputed_data)
    
    cell_prop_input = values$data_list$test_data 
    rownames(cell_prop_input) = cell_prop_input$IlmnID
    cell_prop_input$IlmnID = NULL
    BloodFrac.m <- epidish(beta.m = cell_prop_input,
                           ref.m = centDHSbloodDMC.m,
                           method = "RPC")$estF
    
    bf = as.data.frame(t(BloodFrac.m))
    bf$CellType = rownames(bf)
    
    bf = pivot_longer(bf, -"CellType",
                      names_to="Sample",
                      values_to="CellProp")
    
    
    # Modified plot with better sizing and theme adjustments
    ggplot(bf, aes(CellType, CellProp, color=Sample)) + 
      geom_point(size = 3) +
      theme_classic() +
      theme(
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.margin = margin(20, 20, 20, 20)
      ) +
      ylab("Cell proportion") +
      xlab("") +
      ylim(0, NA) # Ensure y-axis starts at 0
    
  })
  
  
  # Create dimension reduction plots with pagination
  output$dimReductionPlots <- renderUI({
    req(values$plot_data, input$proband, input$signatures, input$plotsPerPage)
    
    # Calculate which plots to show on current page
    start_idx <- (values$current_page - 1) * input$plotsPerPage + 1
    end_idx <- min(values$current_page * input$plotsPerPage, length(input$signatures))
    current_signatures <- input$signatures[start_idx:end_idx]
    
    plot_output_list <- lapply(current_signatures, function(s) {
      plotname <- paste("plot", s, sep = "")
      plotOutput(plotname, height = 800)
    })
    
    do.call(tagList, plot_output_list)
  })
  
  # Render individual dimension reduction plots
  observe({
    req(values$plot_data, input$proband, input$signatures, input$plotsPerPage)
    
    # Calculate which plots to show on current page
    start_idx <- (values$current_page - 1) * input$plotsPerPage + 1
    end_idx <- min(values$current_page * input$plotsPerPage, length(input$signatures))
    current_signatures <- input$signatures[start_idx:end_idx]

    
    for(s in current_signatures) {
      local({
        s_local <- s

        output[[paste("plot", s_local, sep = "")]] <- renderPlot({
          create_dimension_reduction_plots(values$plot_data[[s_local]],
                                           values$plot_metadata[[s_local]],
                                           input$proband,
                                           s_local)
        })
      })
    }
  })
  
  # Download handlers remain the same
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste("results_", Sys.Date(), ".tsv", sep = "\tsv")
    },
    content = function(file) {
      write.table(values$results, file, quote = F, row.names = F)
    }
  )
  
  output$downloadPlots <- downloadHandler(
    filename = function() {
      paste("plots_", Sys.Date(), ".pdf", sep = "")
    },
    content = function(file) {
      pdf(file, width = 35, height = 8)
      
      # Print prediction plot
      print(create_prediction_plot(values$results, input$proband))
      
      # Print dimension reduction plots for selected signatures
      for(s in input$signatures) {
        

        dim_plots <- create_dimension_reduction_plots(values$plot_data[[s]],
                                                      values$plot_metadata[[s]],
                                                      input$proband,
                                                      s)
        print(dim_plots)
      }
      
      dev.off()
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))