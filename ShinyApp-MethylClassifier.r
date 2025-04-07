# app.R
# Max user input size
options(shiny.maxRequestSize=300*1024^2)

source("ShinyApp-MethylClassifier-helper_functions.R")

ui <- dashboardPage(
  dashboardHeader(title = "MethaDory"),
  dashboardSidebar(
    
    shinyDirButton("modelDir", "Select Model Folder", "Please select folder containing SVM models",
                   style="margin: 5px 5px 5px 45px;"),
    verbatimTextOutput("modelDirPath"),
    hr(),
    
    fileInput("dataFile", "Upload Test Data (.tsv)
              
              [Max 300MB]",
              accept = c(".tsv", ".txt" )),
    actionButton("loadData", "Load Data", class = "btn-primary",
                 style="margin: 0px 0px 0px 70px;"),
    
    checkboxGroupInput("proband", "Select Proband(s) for plotting",
                       choices = NULL),
    fluidRow(
      column(5, actionButton("selectAllprobands", "Select All")),
      column(4, actionButton("deselectAllprobands", "Deselect All"))
    ),
    
    hr(),
    
    fluidRow(
      column(5, actionButton("selectAllsignatures", "Select All")),
      column(4, actionButton("deselectAllsignatures", "Deselect All"))
    ),
    checkboxGroupInput("signatures", "Select Signature(s) for plotting",
                       choices = NULL),
    hr(),
    downloadButton("downloadResults", "Download Results",
                   style="margin: 5px 5px 5px 35px;"),
    downloadButton("downloadPlots", "Download Plots",
                   style="margin: 5px 5px 5px 42.5px;")
  ),
  dashboardBody(
    fluidRow(style='height:60vh',
      tabBox(
        id = "tabset1", height = "800px", width =  "600px", 
        tabPanel("Welcome",  includeMarkdown("help.md") ),
        
        tabPanel("Missing values in input", DTOutput("missingValTable")),
        
        tabPanel("Prediction results plots", 
                 div(style ='width:1400px; overflow-x: scroll; height:800px;',# width: 100%;',
                     div(style = "mix-width: 800px; overflow-x: scroll;",
                         plotlyOutput("predictionPlot") %>% 
                           layout(yaxis = list(
                             fixedrange=TRUE))
                 )
        )),
        
        tabPanel("SVM prediction table", DTOutput("predictionTable")),
        tabPanel("Cell proportion deconvolution", 
                 plotOutput("cellpropPlot", height = "500px")),
        
        tabPanel("Chromosomal sex prediction", 
                 DTOutput("chrSexTable", height = "500px"))#,
        # tabPanel("Methylation age prediction", 
        #          DTOutput("methAgeTable", height = "500px"))
        
      )
    ),
    fluidRow(
      box(
        title = "Dimension Reduction Plots",
        width = 12,
        # Pagination controls
        fluidRow(
          column(4, numericInput("plotsPerPage", "Plots per page:", 1, min = 1, max = 10)),
          column(4, uiOutput("pageControls"))
        ),
        hr(),
        uiOutput("dimReductionPlots")
      )
    )
  )
)

server <- function(input, output, session) {

  #Initialize shinyDirChoose with current working directory
  volumes <- c("Home" = fs::path_home(),
               "MethaDory Directory" = getwd(),
               "R Installation" = R.home()
  )
  
  shinyDirChoose(input, 
                 'modelDir',
                 defaultRoot = "Home",
                 defaultPath = "MethaDory Directory",
                 session = session,
                 roots = volumes
  )
  
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
    current_page = 1
  )
  
  # Display selected model directory
  output$modelDirPath <- renderText({
    if (!is.null(input$modelDir)) {
      values$model_dir <- parseDirPath(volumes, input$modelDir)
      paste(values$model_dir)
    }
  })

  
  
  # Load data when button is clicked
  observeEvent(input$loadData, {
    req(input$dataFile, values$model_dir)
    
    # Validate model directory
    if (!dir.exists(values$model_dir)) {
      log_message("Error: Invalid model directory")
      return()
    }
    
    
    # Initial proband selection update after data load
    observeEvent(input$loadData, {
      req(values$data_list$test_data_ids)
      
      # Update proband choices with all available test data IDs
      updateCheckboxGroupInput(session, "proband",
                               choices = values$data_list$test_data_ids,
                               selected = values$data_list$test_data_ids)  
    })
    show_modal_spinner(spin = "rotating-plane",
                       color =  "#279CED",
                       text = "Please wait...")
    # Load test data
    values$data_list <- load_test_data(input$dataFile$datapath)

    # Load background data
    values$background_data <- load_background_data(values$model_dir)

    # Prepare and impute data
    imputation_data <- prepare_imputation_data(values$data_list$test_data,
                                               values$background_data$imputation_background,
                                               values$background_data$svm)

    values$imputed_data <- perform_imputation(imputation_data$test_data,
                                              values$data_list$test_data_ids)
    

    # Load additional data for plotting
    values$beta_sig_data <- load_beta_signatures()

    # Prepare inference data and make predictions
    values$inference_data <- prepare_inference_data(values$imputed_data, 
                                                    values$background_data$svm)

    values$results <- make_predictions(values$inference_data,
                                       values$background_data$svm,
                                       values$data_list$test_data_ids)
    
    values$real_cases_data <- load_real_cases(values$beta_sig_data$signatures)
    
    # Prepare plot data
    values$plot_data <- prepare_plot_data(values$beta_sig_data$signatures,
                                          values$beta_sig_data$insilico_beta,
                                          values$imputed_data,
                                          values$real_cases_data$real_cases_beta)

    values$plot_metadata <- prepare_plot_metadata(values$beta_sig_data$signatures,
                                                  values$plot_data,
                                                  values$imputed_data,
                                                  values$real_cases_data$real_cases_meta)

    remove_modal_spinner()
    # Update UI elements
    updateCheckboxGroupInput(session, "proband",
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
  
  
  
  # Handle Select All button for probands 
  observeEvent(input$selectAllprobands, {
    req(values$data_list$test_data_ids)
    updateCheckboxGroupInput(session, "proband",
                             choices = values$data_list$test_data_ids,
                             selected = values$data_list$test_data_ids)
  })
  
  # Handle Deselect All button for probands 
  observeEvent(input$deselectAllprobands, {
    req(values$data_list$test_data_ids)
    updateCheckboxGroupInput(session, "proband",
                             choices = values$data_list$test_data_ids,
                             selected = character(0))
  })
  
  
  
  # Handle Select All button signatures 
  observeEvent(input$selectAllsignatures, {
    req(values$beta_sig_data)
    updateCheckboxGroupInput(session, "signatures",
                             choices = names(values$beta_sig_data$signatures),
                             selected = names(values$beta_sig_data$signatures))
  })
  
  # Handle Deselect All button signatures 
  observeEvent(input$deselectAllsignatures, {
    req(values$beta_sig_data)
    updateCheckboxGroupInput(session, "signatures",
                             choices = names(values$beta_sig_data$signatures),
                             selected = character(0))
  })
  
  
  
  show_modal_spinner()
  


  
  # Create dynamic missing value table
  output$missingValTable <- renderDT({
    req(values$data_list$test_data, 
        input$proband,
        values$beta_sig_data$signatures,
        input$signatures)
    
    datatable(
      count_missing_data(values$data_list$test_data, 
                         input$proband,
                         values$beta_sig_data$signatures,
                         input$signatures)
    )
  })
  
  
  
  
  
  # Create dynamic predicted age table
  output$chrSexTable <- renderDT({
    req(values$imputed_data, input$proband
    )
    
    datatable(
      predict_chr_sex(values$data_list$test_data, 
                      input$proband)
    )
  })
  
  
  
  # Create dynamic prediction table
  output$predictionTable <- renderDT({
    req(values$results, input$proband, input$signatures)
    
    datatable(
      values$results[values$results$SampleID %in% input$proband &
                       values$results$SVM %in% gsub("_", " ",input$signatures), ]
    )
  })
  
  
  
  # Create prediction plot
  output$predictionPlot <- renderPlotly({
    req(values$results, input$proband, input$signatures)
    
    # Calculate plot width based on number of signatures
    plot_width <- 100 + (length(input$signatures) * 100)
    
    # Create prediction plot with dynamic width
    plot <- create_prediction_plot(values$results, input$proband, input$signatures, plot_width)
    
  })
  
  # Cell proportion plot 
  output$cellpropPlot <- renderPlot({
    req(values$data_list$test_data, 
        values$background_data$cellprops,
        input$proband)
    
    create_cell_deconv_plot(values$data_list$test_data, 
                            values$background_data$cellprops,
                            input$proband)
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
      plotOutput(plotname, height = 1200)
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
  
  remove_modal_spinner() 
  output$downloadResults <- downloadHandler(
    filename = function() {
      paste0("MethaDory_results_", Sys.Date(), ".tsv")
    },
    content = function(file) {
      
      export_table = values$results
      
      export_table = pivot_wider(export_table, id_cols = "SampleID",  
                                 names_from = "SVM", 
                                 values_from = c("pSVM_average", "pSVM_sd"))
      
      
      names(export_table) = gsub("pSVM_", "pSVM", names(export_table))
      names(export_table) = gsub("_", " ", names(export_table))
      
      
      write.table(export_table, file, quote = F, row.names = F, sep = "\t")
    }
  )
  output$downloadPlots <- downloadHandler(
    filename = function() {
      paste0("MethaDory_plots_", Sys.Date(), ".pdf")
    },
    content = function(file) {
      pdf(file, width = 20, height = 8)
      plot_width <- min(1000, 100 + (length(input$signatures) * 100))
      print(create_prediction_plot2(values$results, input$proband, input$signatures, plot_width))
      
      # Print dimension reduction plots for selected signatures
      for(s in input$signatures) {
        
        s_local <- s
        
        
        
        print(create_dimension_reduction_plots(values$plot_data[[s_local]],
                                               values$plot_metadata[[s_local]],
                                               input$proband,
                                               s_local))
        
      }
      
      dev.off()
      
    })
}
  



# Run the app
shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))