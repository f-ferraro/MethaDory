ui <- dashboardPage(
  dashboardHeader(title = "MethaDory"),
  dashboardSidebar(
    shinyDirButton("modelDir", "Select Model Folder", "Please select folder containing SVM/ and NNET/ subfolders"),
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
    numericInput("minPSVMForPlots", "Minimum mean(SVM,NNET) for dimension plots:",
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

               # All sample QC in one panel: missing values, methylation age, cell
               # proportions | chromosomal sex, PCA vs controls, beta-value distribution
               tabPanel("QC",
                        # Scrolls inside the fixed-height tabBox instead of spilling over it
                        div(style = "max-height: 920px; overflow-y: auto; padding-right: 10px;",
                          h4("Missing values before imputation"),
                          uiOutput("qcMissing"),

                          h4("Methylation age prediction", style = "margin-top: 25px;"),
                          DTOutput("methAgeTable"),

                          h4("Cell proportions and chromosomal sex", style = "margin-top: 25px;"),
                          p("The table shows whether each sample's cell proportions are within the distributions observed in the training samples"),
                          div(
                            checkboxInput("showAllOutliers", "Show all results", value = FALSE),
                            style = "margin-bottom: 10px;"
                          ),
                          DTOutput("cellPropOutlierTable"),
                          br(),
                          plotOutput("qcCellSexPlot", height = "700px"),

                          h4("PCA against controls", style = "margin-top: 25px;"),
                          p("Principal component analysis of each selected sample together with the control samples, computed on the beta values before imputation using the top 1% most variable autosomal CpGs measured in the sample. Controls are shown in grey and the sample in red. A sample lying far from the controls may be of poor quality (or come from a different tissue or platform) and its predictions should be interpreted with caution."),
                          plotOutput("qcPcaPlot", height = "auto"),

                          h4("Beta-value distribution against controls", style = "margin-top: 25px;"),
                          p("Distribution of the beta values before imputation over all autosomal CpGs measured in the sample. Controls are shown in grey and the sample in red. A distribution that departs from the two peaks (near 0 and 1) of the controls may indicate poor quality, or data from a different platform that is not on the array scale."),
                          plotOutput("qcDensityPlot", height = "auto")
                        )),

               tabPanel("Prediction results plot",
                        div(
                          numericInput("plotThreshold", "Minimum mean(SVM,NNET) Score Threshold:",
                                     value = 0.0, min = 0, max = 1, step = 0.05,
                                     width = "300px"),
                          div(style = 'overflow-x: auto; white-space: nowrap;',
                              plotlyOutput("predictionPlot", height = "800px")
                          )
                        )
               ),

               tabPanel("Prediction table",
                        numericInput("minPSVM", "Minimum mean(SVM,NNET):", value = 0.25, min = 0, max = 1, step = 0.01),
                        DTOutput("predictionTable")),
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