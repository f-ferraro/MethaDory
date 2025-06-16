library(shiny)
library(shinydashboard)
library(shinyFiles)
library(DT)
library(plotly)

options(shiny.maxRequestSize = 300 * 1024^2)

source("ShinyApp-MethylClassifier-helper_functions.R")

ui <- dashboardPage(
  dashboardHeader(title = "MethaDory"),
  dashboardSidebar(
    shinyDirButton("modelDir", "Select Model Folder", "Please select folder containing SVM models",
                   style="margin: 5px 5px 5px 45px;"),
    verbatimTextOutput("modelDirPath"),
    hr(),
    fileInput("dataFile", "Upload Test Data (.tsv)\n\n[Max 300MB]",
              accept = c(".tsv", ".txt")),
    actionButton("loadData", "Load Data", class = "btn-primary", style = "margin: 0px 0px 0px 70px;"),
    checkboxGroupInput("proband", "Select Proband(s) for plotting", choices = NULL),
    fluidRow(
      column(5, actionButton("selectAllprobands", "Select All")),
      column(4, actionButton("deselectAllprobands", "Deselect All"))
    ),
    hr(),
    fluidRow(
      column(5, actionButton("selectAllsignatures", "Select All")),
      column(4, actionButton("deselectAllsignatures", "Deselect All"))
    ),
    checkboxGroupInput("signatures", "Select Signature(s) for plotting", choices = NULL),
    hr(),
    downloadButton("downloadResults", "Download Results", style="margin: 5px 5px 5px 35px;"),
    downloadButton("downloadPlots", "Download Plots", style="margin: 5px 5px 5px 42.5px;")
  ),
  dashboardBody(
    fluidRow(style='height:60vh',
             tabBox(
               id = "tabset1", height = "800px", width =  "600px", 
               tabPanel("Welcome", includeMarkdown("help.md")),
               
               #tabPanel("Missing values in input", DTOutput("missingValTable")),
               
               tabPanel("Prediction results plot", 
                        div(style = 'overflow-x: auto; white-space: nowrap;',
                            plotlyOutput("predictionPlot", height = "800px")  
                        )
               ),
               
               
               tabPanel("SVM prediction table", DTOutput("predictionTable")),
               tabPanel("Cell proportion deconvolution", plotOutput("cellpropPlot", height = "700px")),
               tabPanel("Chromosomal sex prediction", plotOutput("chrSexPlot", height = "800px")),
               tabPanel("Methylation age prediction", DTOutput("methAgeTable", height = "800px")),
               tabPanel("References", includeMarkdown("references.md"))
             )
    ),
    fluidRow(
      box(
        title = "Dimension Reduction Plots",
        width = 12,
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
  volumes <- c("Home" = fs::path_home(), "MethaDory Directory" = getwd(), "R Installation" = R.home())
  shinyDirChoose(input, 'modelDir', defaultRoot = "Home", defaultPath = "MethaDory Directory", session = session, roots = volumes)
  
  values <- reactiveValues(data = NULL, model_dir = NULL, current_page = 1)
  
  output$modelDirPath <- renderText({
    if (!is.null(input$modelDir)) {
      values$model_dir <- parseDirPath(volumes, input$modelDir)
      paste(values$model_dir)
    }
  })
  
  observeEvent(input$loadData, {
    req(input$dataFile, values$model_dir)
    show_modal_spinner(spin = "rotating-plane", color = "#279CED", text = "Please wait...")
    
    data_list <- load_test_data(input$dataFile$datapath)
    background_data <- load_background_data(values$model_dir)
    imputation_data <- prepare_imputation_data(data_list$test_data, background_data$imputation_background, background_data$svm)
    imputed_data <- perform_imputation(imputation_data$test_data, data_list$test_data_ids)
    
    beta_sig_data <- load_beta_signatures()
    inference_data <- prepare_inference_data(imputed_data, background_data$svm)
    results <- make_predictions(inference_data, background_data$svm, data_list$test_data_ids)
    real_cases_data <- load_real_cases(beta_sig_data$signatures)
    insilico_meta <- readRDS("data/syntheticcases/samples_for_insilico.meta.rds")
    insilico_meta = insilico_meta[,c("geo_accession", "Platform")]
    names(insilico_meta) = c("IDs", "Platform")
      
    plot_data <- prepare_plot_data(beta_sig_data$signatures, beta_sig_data$insilico_beta, imputed_data, real_cases_data$real_cases_beta)
    plot_metadata <- prepare_plot_metadata(beta_sig_data$signatures, plot_data, imputed_data, real_cases_data$real_cases_meta, insilico_meta)
    
    cell_props <- create_cell_deconv_table(data_list$test_data)
    chr_sex_table <- predict_chr_sex_table(data_list$test_data)
    age_table <- predict_age(data_list$test_data)
    
    values$data <- list(data_list=data_list, background_data=background_data, imputed_data=imputed_data,
                        beta_sig_data=beta_sig_data, inference_data=inference_data, results=results,
                        real_cases_data=real_cases_data, plot_data=plot_data, plot_metadata=plot_metadata,
                        cell_props=cell_props, chr_sex_table=chr_sex_table, age_table=age_table)
    
    updateCheckboxGroupInput(session, "proband", choices = data_list$test_data_ids, selected = data_list$test_data_ids)
    updateCheckboxGroupInput(session, "signatures", choices = names(beta_sig_data$signatures), selected = character(0))
    
    remove_modal_spinner()
  })
  
  observeEvent(input$selectAllprobands, {
    req(values$data)
    updateCheckboxGroupInput(session, "proband", choices = values$data$data_list$test_data_ids, selected = values$data$data_list$test_data_ids)
  })
  observeEvent(input$deselectAllprobands, {
    req(values$data)
    updateCheckboxGroupInput(session, "proband", choices = values$data$data_list$test_data_ids, selected = character(0))
  })
  observeEvent(input$selectAllsignatures, {
    req(values$data)
    updateCheckboxGroupInput(session, "signatures", choices = names(values$data$beta_sig_data$signatures), selected = names(values$data$beta_sig_data$signatures))
  })
  observeEvent(input$deselectAllsignatures, {
    req(values$data)
    updateCheckboxGroupInput(session, "signatures", choices = names(values$data$beta_sig_data$signatures), selected = character(0))
  })
  
  output$missingValTable <- renderDT({
    req(values$data, input$proband, input$signatures)
    datatable(count_missing_data(values$data$data_list$test_data, input$proband, values$data$beta_sig_data$signatures, input$signatures))
  })
  output$chrSexPlot <- renderPlot({
    req(values$data, input$proband)
    predict_chr_sex_plot(values$data$chr_sex_table, input$proband)
  })
  output$methAgeTable <- renderDT({
    req(values$data, input$proband)
    datatable(values$data$age_table[values$data$age_table$Proband %in% input$proband, ])
  })
  output$predictionTable <- renderDT({
    req(values$data, input$proband, input$signatures)
    datatable(values$data$results[values$data$results$SampleID %in% input$proband &
                                    values$data$results$SVM %in% gsub("_", " ", input$signatures), ])
  })
  output$predictionPlot <- renderPlotly({
    req(values$data, input$proband, input$signatures)
    plot_width <- 100 + (length(input$signatures) * 100)
    
    create_prediction_plot(values$data$results, input$proband, input$signatures, plot_width) 
  })
  output$cellpropPlot <- renderPlot({
    req(values$data, input$proband)
    create_cell_deconv_plot(values$data$cell_props, values$data$background_data$cellprops, input$proband) 
    
  })
  
  output$pageControls <- renderUI({
    req(input$signatures)
    total_pages <- ceiling(length(input$signatures) / input$plotsPerPage)
    if (total_pages > 1) {
      div(actionButton("prevPage", "Previous"), span(paste("Page", values$current_page, "of", total_pages)), actionButton("nextPage", "Next"))
    }
  })
  observeEvent(input$prevPage, {
    values$current_page <- max(1, values$current_page - 1)
  })
  observeEvent(input$nextPage, {
    total_pages <- ceiling(length(input$signatures) / input$plotsPerPage)
    values$current_page <- min(total_pages, values$current_page + 1)
  })
  
  output$dimReductionPlots <- renderUI({
    req(values$data, input$proband, input$signatures, input$plotsPerPage)
    start_idx <- (values$current_page - 1) * input$plotsPerPage + 1
    end_idx <- min(values$current_page * input$plotsPerPage, length(input$signatures))
    current_signatures <- input$signatures[start_idx:end_idx]
    plot_output_list <- lapply(current_signatures, function(s) plotOutput(paste0("plot", s), height = 1200))
    do.call(tagList, plot_output_list)
  })
  observe({
    req(values$data, input$proband, input$signatures, input$plotsPerPage)
    start_idx <- (values$current_page - 1) * input$plotsPerPage + 1
    end_idx <- min(values$current_page * input$plotsPerPage, length(input$signatures))
    current_signatures <- input$signatures[start_idx:end_idx]
    for (s in current_signatures) {
      local({
        s_local <- s
        output[[paste0("plot", s_local)]] <- renderPlot({
          create_dimension_reduction_plots(values$data$plot_data[[s_local]], values$data$plot_metadata[[s_local]], input$proband, s_local)
        })
      })
    }
  })
  
  output$downloadResults <- downloadHandler(
    filename = function() paste0("MethaDory_results_", Sys.Date(), ".xlsx"),
    content = function(file) {
      export_table <- values$data$results
      export_table <- pivot_wider(export_table, id_cols = "SampleID", names_from = "SVM", values_from = c("pSVM_average", "pSVM_sd"))
      names(export_table) <- gsub("pSVM_", "pSVM", names(export_table))
      names(export_table) <- gsub("_", " ", names(export_table))
      cell_props_exports <- pivot_wider(values$data$cell_props[values$data$cell_props$Proband %in% input$proband,],
                                        id_cols = "Proband", names_from = "CellType", values_from = "CellProp")
      openxlsx::write.xlsx(list("MethaDory Predictions" = export_table,
                                "Cell deconvolutions" = cell_props_exports,
                                "Predicted age" = values$data$age_table[values$data$age_table$Proband %in% input$proband,],
                                "Predicted chr sex" = values$data$chr_sex_table[values$data$age_table$Proband %in% input$proband,]),
                           file, quote = FALSE, row.names = FALSE, sep = "\t")
    })
  
  output$downloadPlots <- downloadHandler(
    filename = function() paste0("MethaDory_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      pdf(file, width = 20, height = 8)
      print(create_prediction_plot2(values$data$results, input$proband, input$signatures, 100 + (length(input$signatures) * 100)))
      for (s in input$signatures) {
        print(create_dimension_reduction_plots(values$data$plot_data[[s]], values$data$plot_metadata[[s]], input$proband, s))
      }
      dev.off()
    })
}

shinyApp(ui = ui, server = server, options = list(launch.browser = TRUE))
