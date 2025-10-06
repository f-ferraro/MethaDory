source("../core/data_processing.R")
source("../core/svm_prediction.R")
source("../core/methylation_analysis.R")
source("../visualization/prediction_plots.R")
source("../visualization/dimension_plots.R")
source("../visualization/analysis_plots.R")
source("../export/html_export.R")
source("../export/table_export.R")

server <- function(input, output, session) {
  # Set up directory browser 
  methadory_root <- file.path("..", "..")
  volumes <- c("Home" = fs::path_home(),
               "MethaDory Directory" = normalizePath(methadory_root),
               "Data Directory" = normalizePath(file.path(methadory_root, "data")),
               "R Installation" = R.home())
  shinyDirChoose(input, 'modelDir', defaultRoot = "Data Directory", defaultPath = "", session = session, roots = volumes)

  values <- reactiveValues(data = NULL, model_dir = NULL, current_page = 1,
                           high_scoring_cache = NULL, plot_cache = list())

  output$modelDirPath <- renderText({
    if (!is.null(input$modelDir)) {
      values$model_dir <- parseDirPath(volumes, input$modelDir)
      paste(values$model_dir)
    }
  })

  observeEvent(input$loadData, {
    req(input$dataFile, values$model_dir, input$nImputationSamples)
    show_modal_spinner(spin = "rotating-plane", color = "#279CED", text = "Please wait...")

    data_list <- load_test_data(input$dataFile$datapath)
    background_data <- load_background_data(values$model_dir)
    cat("Using", input$nImputationSamples, "closest samples for imputation\n")
    imputation_data <- prepare_imputation_data(data_list$test_data, background_data$imputation_background,
                                               background_data$svm, n_closest = input$nImputationSamples)
    imputed_data <- perform_imputation(imputation_data$test_data, data_list$test_data_ids)

    beta_sig_data <- load_beta_signatures()
    inference_data <- prepare_inference_data(imputed_data, background_data$svm)
    results <- make_predictions(inference_data, background_data$svm, data_list$test_data_ids)
    real_cases_data <- load_real_cases(beta_sig_data$signatures)

    # Use insilico_meta from beta_sig_data and prepare it for plotting
    insilico_meta <- beta_sig_data$insilico_meta[,c("geo_accession", "platform_id", "Sex", "AgeGroup")]
    names(insilico_meta) <- c("IDs", "Platform", "Sex", "AgeGroup")

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
    updateSelectizeInput(session, "signatures", choices = names(beta_sig_data$signatures), selected = character(0))

    remove_modal_spinner()
  })

  # Cached high-scoring signatures with reactive invalidation
  high_scoring_signatures <- reactive({
    req(values$data, input$proband, input$signatures, input$minPSVMForPlots)

    # Create cache key including the threshold
    cache_key <- paste(c(paste(sort(input$proband), collapse = "_"),
                        paste(sort(input$signatures), collapse = "_"),
                        as.character(input$minPSVMForPlots)),
                      collapse = "_")

    # Check if cache is valid
    if (!is.null(values$high_scoring_cache) &&
        values$high_scoring_cache$key == cache_key) {
      cat("Using cached high-scoring signatures:", length(values$high_scoring_cache$signatures), "\n")
      return(values$high_scoring_cache$signatures)
    }

    # Compute and cache using user-defined threshold
    # cat("Computing high-scoring signatures with threshold:", input$minPSVMForPlots, "\n")
    # cat("Selected signatures:", length(input$signatures), "\n")
    # cat("Selected probands:", length(input$proband), "\n")

    signatures <- get_high_scoring_signatures_with_threshold(
      values$data$results,
      input$proband,
      input$signatures,
      input$minPSVMForPlots
    )

    # cat("Found", length(signatures), "high-scoring signatures\n")
    # if(length(signatures) > 0) {
    #   cat("High-scoring signatures:", paste(head(signatures, 5), collapse=", "), "...\n")
    # }

    values$high_scoring_cache <- list(key = cache_key, signatures = signatures)

    # Clear plot cache when signatures change
    values$plot_cache <- list()

    return(signatures)
  })

  # Debounced reactive for plot generation
  plot_trigger <- debounce(reactive({
    list(
      probands = input$proband,
      page = values$current_page,
      plots_per_page = input$plotsPerPage
    )
  }), 300) # 300ms delay

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
    updateSelectizeInput(session, "signatures", choices = names(values$data$beta_sig_data$signatures), selected = names(values$data$beta_sig_data$signatures))
  })
  observeEvent(input$deselectAllsignatures, {
    req(values$data)
    updateSelectizeInput(session, "signatures", choices = names(values$data$beta_sig_data$signatures), selected = character(0))
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
    req(values$data, input$proband, input$minPSVM)

    # Filter by proband and minimum pSVM threshold
    filtered_data <- values$data$results[values$data$results$SampleID %in% input$proband &
                                        values$data$results$pSVM_average >= input$minPSVM, ]

    # Additionally filter by signatures if any are selected
    if (!is.null(input$signatures) && length(input$signatures) > 0) {
      filtered_data <- filtered_data[filtered_data$SVM %in% gsub("_", " ", input$signatures), ]
    }

    datatable(filtered_data, options = list(pageLength = 25)) %>%
      formatRound(columns = c("pSVM_average", "pSVM_sd"), digits = 2) %>%
      formatStyle(
        columns = "pSVM_average",
        backgroundColor = styleInterval(c(0.25, 0.5), c("white", "#FAD302", "#9E1E05")),
        color = styleInterval(0.5, c("black", "white")),
        fontWeight = styleInterval(0.5, c("normal", "bold"))
      )
  })
  output$predictionPlot <- renderPlotly({
    req(values$data, input$proband, input$signatures, input$plotThreshold)

    # Calculate plot width based on filtered data
    if(input$plotThreshold > 0) {
      # Count how many signatures will remain after filtering
      filtered_data <- values$data$results[values$data$results$SampleID %in% input$proband &
                                         values$data$results$SVM %in% gsub("_", " ", input$signatures) &
                                         values$data$results$pSVM_average >= input$plotThreshold, ]
      n_signatures <- length(unique(filtered_data$SVM))
    } else {
      n_signatures <- length(input$signatures)
    }

    plot_width <- 100 + (n_signatures * 100)

    create_prediction_plot_interactive(values$data$results, input$proband, input$signatures, plot_width, input$plotThreshold)
  })


  output$cellPropOutlierTable <- renderDT({
    req(values$data, input$proband)

    outlier_table <- create_cell_prop_outlier_table(
      values$data$cell_props,
      values$data$background_data$cellprops,
      input$proband
    )

    if (nrow(outlier_table) == 0) {
      return(datatable(data.frame(Message = "No samples selected"),
                      options = list(dom = 't', ordering = FALSE)))
    }

    # Get the current state of showAllOutliers 
    show_all <- if(is.null(input$showAllOutliers)) FALSE else input$showAllOutliers

    # print(paste("Unique status values:", paste(unique(outlier_table$status), collapse = ", ")))
    # print(paste("Show all:", show_all))
    # print(paste("Total rows before filtering:", nrow(outlier_table)))

    # Filter to show only WARNING|FAIL rows if checkbox is unchecked
    if (!show_all) {
      outlier_table <- outlier_table[outlier_table$status == "WARNING"| outlier_table$status == "FAIL", ]
      print(paste("Rows after WARNING filtering:", nrow(outlier_table)))
    }

    # Check if any rows remain after filtering
    if (nrow(outlier_table) == 0) {
      return(datatable(data.frame(Message = "Cell proportion table not shown - all cell proportions are within the expected range"),
                      options = list(dom = 't', ordering = FALSE)))
    }

    # Sort data to put FAIL rows first when showing all results
    if (show_all) {
      outlier_table <- outlier_table[order(outlier_table$status, decreasing = TRUE), ]
    }

    datatable(outlier_table,
              options = list(
                pageLength = 15,
                dom = 'tp',
                ordering = TRUE,
                order = if(show_all) list(list(3, 'desc')) else list(),
                columnDefs = list(
                  list(targets = c(3), className = 'dt-center'),
                  list(targets = c(2, 4), className = 'dt-right')
                )
              ),
              colnames = c("Sample", "Cell Type", "Proportion", "Status", "SD from Mean"),
              rownames = FALSE) %>%
      formatStyle("status",
                  backgroundColor = styleEqual(c("WARNING", "FAIL"), c("#FFEBCD", "#DE4C35"))) 
  })

  output$cellpropPlot <- renderPlot({
    req(values$data, input$proband)
    create_cell_deconv_plot(values$data$cell_props, values$data$background_data$cellprops, input$proband)

  })

  output$pageControls <- renderUI({
    req(input$plotsPerPage)
    signatures <- high_scoring_signatures()
    total_pages <- ceiling(length(signatures) / input$plotsPerPage)
    if (total_pages > 1) {
      div(actionButton("prevPage", "Previous"), span(paste("Page", values$current_page, "of", total_pages)), actionButton("nextPage", "Next"))
    }
  })

  output$jumpToPlot <- renderUI({
    signatures <- high_scoring_signatures()
    if (length(signatures) > 1) {
      div(
        selectInput("selectedSignature", "Jump to plot:",
                   choices = c("Select signature..." = "", setNames(signatures, signatures)),
                   selected = ""),
        style = "margin-top: 0px;"
      )
    }
  })
  observeEvent(input$prevPage, {
    values$current_page <- max(1, values$current_page - 1)
  })
  observeEvent(input$nextPage, {
    signatures <- high_scoring_signatures()
    total_pages <- ceiling(length(signatures) / input$plotsPerPage)
    values$current_page <- min(total_pages, values$current_page + 1)
  })

  observeEvent(input$selectedSignature, {
    req(input$selectedSignature)
    if (input$selectedSignature != "") {
      signatures <- high_scoring_signatures()
      if (input$selectedSignature %in% signatures) {
        # Find which page contains the selected signature
        signature_index <- which(signatures == input$selectedSignature)
        target_page <- ceiling(signature_index / input$plotsPerPage)
        values$current_page <- target_page

        # Reset the dropdown
        updateSelectInput(session, "selectedSignature", selected = "")
      }
    }
  })

  output$dimReductionPlots <- renderUI({
    req(input$plotsPerPage)

    # Get cached high-scoring signatures
    signatures <- high_scoring_signatures()

    if(length(signatures) == 0) {
      threshold_msg <- paste("No signatures with pSVM ≥", input$minPSVMForPlots, "found.")
      return(div(class = "alert alert-warning",
                 threshold_msg))
    }

    start_idx <- (values$current_page - 1) * input$plotsPerPage + 1
    end_idx <- min(values$current_page * input$plotsPerPage, length(signatures))
    current_signatures <- signatures[start_idx:end_idx]
    plot_output_list <- lapply(current_signatures, function(s) plotOutput(paste0("plot", s), height = 1200))
    do.call(tagList, plot_output_list)
  })
  observe({
    req(values$data)

    # Use debounced trigger to avoid rapid re-renders
    trigger <- plot_trigger()
    req(trigger$plots_per_page)

    # Get cached high-scoring signatures
    signatures <- high_scoring_signatures()

    if(length(signatures) > 0) {
      start_idx <- (values$current_page - 1) * trigger$plots_per_page + 1
      end_idx <- min(values$current_page * trigger$plots_per_page, length(signatures))
      current_signatures <- signatures[start_idx:end_idx]

      for (s in current_signatures) {
        local({
          s_local <- s

          # Check plot cache first
          cache_key <- paste(s_local, paste(sort(trigger$probands), collapse = "_"), sep = "_")

          if (!is.null(values$plot_cache[[cache_key]])) {
            output[[paste0("plot", s_local)]] <- values$plot_cache[[cache_key]]
          } else {
            # Generate and cache plot with progress indication
            plot_render <- renderPlot({
              withProgress(message = paste("Generating plot for", s_local), {
                tryCatch({
                  create_dimension_reduction_plots(values$data$plot_data[[s_local]],
                                                 values$data$plot_metadata[[s_local]],
                                                 trigger$probands, s_local,
                                                 values$data$age_table,
                                                 values$data$chr_sex_table,
                                                 n_samples_per_group = input$nSamplesPerGroup)
                }, error = function(e) {
                  cat("ERROR generating plot for", s_local, ":", e$message, "\n")
                  print(str(values$data$plot_data[[s_local]]))
                  print(str(values$data$plot_metadata[[s_local]]))
                  # Return an error plot
                  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
                  text(1, 1, paste("Error generating plot:\n", e$message), cex = 1.2, col = "red")
                })
              })
            })

            values$plot_cache[[cache_key]] <- plot_render
            output[[paste0("plot", s_local)]] <- plot_render
          }
        })
      }
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
                                "Predicted chr sex" = values$data$chr_sex_table[values$data$chr_sex_table$Proband %in% input$proband,]),
                           file, quote = FALSE, row.names = FALSE, sep = "\t")
    })

  output$downloadPlots <- downloadHandler(
    filename = function() paste0("MethaDory_plots_", Sys.Date(), ".pdf"),
    content = function(file) {
      # Filter to only high-scoring signatures
      high_scoring <- values$data$results[values$data$results$pSVM_average >= input$minPSVMForPlots, ]
      high_scoring_sigs <- unique(high_scoring$SVM)
      high_scoring_sigs <- gsub(" ", "_", high_scoring_sigs)

      # Only include signatures that are selected AND high-scoring AND have plot data
      sigs_to_plot <- intersect(input$signatures, high_scoring_sigs)
      sigs_to_plot <- sigs_to_plot[sigs_to_plot %in% names(values$data$plot_data)]

      pdf(file, width = 25, height = 16)
      print(create_prediction_plot(values$data$results, input$proband, input$signatures, 100 + (length(input$signatures) * 100)))

      for (s in sigs_to_plot) {
        tryCatch({
          print(create_dimension_reduction_plots(values$data$plot_data[[s]], values$data$plot_metadata[[s]],
                                                input$proband, s, values$data$age_table, values$data$chr_sex_table,
                                                n_samples_per_group = input$nSamplesPerGroup))
        }, error = function(e) {
          cat("ERROR generating plot for", s, "in PDF export:", e$message, "\n")
          # Create an error page in the PDF
          plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
          text(1, 1, paste("Error generating plot for", s, ":\n", e$message), cex = 1.2, col = "red")
        })
      }
      dev.off()
    })
}