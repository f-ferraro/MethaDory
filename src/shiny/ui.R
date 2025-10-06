ui <- dashboardPage(
  dashboardHeader(title = "MethaDory"),
  dashboardSidebar(
    shinyDirButton("modelDir", "Select Model Folder", "Please select folder containing SVM models (typically data/svm-models/)"),
    verbatimTextOutput("modelDirPath"),
    hr(),
    fileInput("dataFile", "Upload Test Data (.tsv)\n\n[Max 500MB]",
              accept = c(".tsv", ".txt")),
    numericInput("nImputationSamples", "Number of closest samples for imputation:",
                 value = 20, min = 5, max = 100, step = 5,
                 width = "200px"),
    actionButton("loadData", "Analyze", class = "btn-primary"),
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
    selectizeInput("signatures", "Select Signature(s) for plotting",
                   choices = NULL,
                   multiple = TRUE,
                   options = list(
                     maxItems = 500,
                     placeholder = 'Search and select signatures...',
                     searchField = 'label',
                     closeAfterSelect = FALSE,
                     highlight = TRUE,
                     maxOptions = 1000,
                     dropdownParent = 'body'
                   )),
    hr(),
    numericInput("minPSVMForPlots", "Minimum pSVM for dimension plots:",
                 value = 0.05, min = 0, max = 1, step = 0.01,
                 width = "200px"),
    numericInput("nSamplesPerGroup", "Number of additional samples to use for PCA and heatmap:",
                 value = 20, min = 5, max = 50, step = 5,
                 width = "200px"),
    hr(),
    h4("Export Options"),
    downloadButton("downloadResults", "Download Results"),
    downloadButton("downloadPlots", "Download Plots")
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .selectize-control.multi .selectize-input {
          max-height: 150px;
          overflow-y: auto;
          border: 1px solid #ddd;
        }

        .selectize-dropdown {
          max-height: 200px;
          overflow-y: auto;
        }

        .sidebar .selectize-control {
          font-size: 12px;
        }

        .selectize-input.items.not-full.has-options.has-items {
          min-height: 40px;
        }

        /* Center all sidebar text */
        .main-sidebar .sidebar {
          text-align: center;
        }

        /* Center labels and inputs */
        .main-sidebar label,
        .main-sidebar .control-label,
        .main-sidebar .checkbox,
        .main-sidebar h4,
        .main-sidebar hr {
          text-align: center;
        }

        /* Center form groups */
        .main-sidebar .form-group {
          text-align: center;
        }

        /* Center input fields */
        .main-sidebar input[type='number'],
        .main-sidebar .selectize-input {
          margin-left: auto;
          margin-right: auto;
        }

        /* Center buttons and file inputs */
        .main-sidebar .btn,
        .main-sidebar .btn-file,
        .main-sidebar .shiny-input-container,
        .main-sidebar .shinyDirectories {
          margin-left: auto;
          margin-right: auto;
          display: block;
        }

        /* Center checkbox groups */
        .main-sidebar .checkbox-group,
        .main-sidebar .shiny-options-group {
          text-align: center;
        }

        /* Center verbatim output */
        .main-sidebar pre {
          text-align: center;
        }

        /* Center fluidRow contents */
        .main-sidebar .row {
          text-align: center;
        }
      "))
    ),
    fluidRow(style='height:80vh',
             tabBox(
               id = "tabset1", height = "1000px", width =  "1200px",
               tabPanel("Welcome", includeMarkdown("html_imports/help.md")),

               tabPanel("Prediction results plot",
                        div(
                          numericInput("plotThreshold", "Minimum SVM Score Threshold:",
                                     value = 0.0, min = 0, max = 1, step = 0.05,
                                     width = "300px"),
                          div(style = 'overflow-x: auto; white-space: nowrap;',
                              plotlyOutput("predictionPlot", height = "800px")
                          )
                        )
               ),

               tabPanel("SVM prediction table",
                        numericInput("minPSVM", "Minimum pSVM_average:", value = 0.25, min = 0, max = 1, step = 0.01),
                        DTOutput("predictionTable")),
               tabPanel("Cell proportion deconvolution",
                        div(
                          h4("Cell Proportion Quality Control", style = "margin-bottom: 15px;"),
                          p("The table shows whether each sample's cell proportions are within the distributions observed in the training samples"),
                          div(
                            checkboxInput("showAllOutliers", "Show all results", value = FALSE),
                            style = "margin-bottom: 10px;"
                          ),
                          DTOutput("cellPropOutlierTable", height = "300px"),
                          br(),
                          plotOutput("cellpropPlot", height = "500px")
                        )),
               tabPanel("Chromosomal sex prediction", plotOutput("chrSexPlot", height = "800px")),
               tabPanel("Methylation age prediction", DTOutput("methAgeTable", height = "800px")),
               tabPanel("References", includeMarkdown("html_imports/references.md"))
             )
    ),
    fluidRow(
      box(
        width = 12,
        p(strong("Note:"), "Only signatures with pSVM ≥ threshold will be plotted."),
        fluidRow(
          column(4, numericInput("plotsPerPage", "Plots per page:", 1, min = 1, max = 10)),
          column(4, uiOutput("pageControls")),
          column(4, uiOutput("jumpToPlot"))
        ),
        hr(),
        uiOutput("dimReductionPlots")
      )
    )
  )
)