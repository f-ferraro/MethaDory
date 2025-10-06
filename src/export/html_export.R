#' Create HTML export of MethaDory results
#'
#' @param values Reactive values from Shiny app
#' @param input Input values from Shiny app
#' @param file_path Output file path
#' @return Creates self-contained HTML file
create_html_export <- function(values, input, file_path) {

  # Get export options
  export_options <- input$htmlExportOptions
  include_dim_plots <- "dim_plots" %in% export_options
  include_cell_plot <- "cell_plot" %in% export_options
  include_chr_sex_plot <- "chr_sex_plot" %in% export_options
  # Set default image compression level
  compress_level <- 0.9

  # Read markdown files
  welcome_content <- ""
  references_content <- ""

  # Load markdown files 
  tryCatch({
    if (file.exists("apps/shiny/html_imports/help.md")) {
      welcome_content <- paste(readLines("apps/shiny/html_imports/help.md"), collapse = "\n")
    }
  }, error = function(e) {
    welcome_content <- "<h3>Welcome to MethaDory</h3><p>Methylation analysis dashboard</p>"
  })

  tryCatch({
    if (file.exists("apps/shiny/html_imports/references.md")) {
      references_content <- paste(readLines("apps/shiny/html_imports/references.md"), collapse = "\n")
    }
  }, error = function(e) {
    references_content <- "<h3>References</h3><p>References not available</p>"
  })

  # Convert markdown to HTML
  welcome_html <- markdown::markdownToHTML(text = welcome_content, fragment.only = TRUE)
  references_html <- markdown::markdownToHTML(text = references_content, fragment.only = TRUE)

  # Filter data based on current selections
  filtered_results <- values$data$results[
    values$data$results$SampleID %in% input$proband &
    values$data$results$SVM %in% gsub("_", " ", input$signatures), ]

  # Create static prediction plot with filtering
  prediction_plot <- create_prediction_plot_static_filtered(values$data$results, input$proband, input$signatures)

  # Conditionally create static plots and outlier table
  cell_prop_plot <- if(include_cell_plot) {
    create_cell_deconv_plot(values$data$cell_props, values$data$background_data$cellprops, input$proband)
  } else NULL

  cell_prop_outlier_table <- if(include_cell_plot) {
    create_cell_prop_outlier_table(values$data$cell_props, values$data$background_data$cellprops, input$proband)
  } else NULL

  chr_sex_plot <- if(include_chr_sex_plot) {
    predict_chr_sex_plot(values$data$chr_sex_table, input$proband)
  } else NULL

  # Filter signatures with pSVM_average >= threshold for dimension plots
  high_scoring_signatures <- character(0)
  min_psvm_threshold <- if (!is.null(input$minPSVMForPlots)) input$minPSVMForPlots else 0.05

  if(include_dim_plots) {
    # Get signatures with pSVM_average >= threshold
    high_scoring_results <- filtered_results[filtered_results$pSVM_average >= min_psvm_threshold, ]
    high_scoring_signatures <- unique(high_scoring_results$SVM)
    high_scoring_signatures <- gsub(" ", "_", high_scoring_signatures)
    high_scoring_signatures <- high_scoring_signatures[high_scoring_signatures %in% input$signatures]
  }

  # Create dimension reduction plots using disk caching for memory efficiency
  dim_plot_files <- list()
  temp_files_to_cleanup <- c()

  if(include_dim_plots && length(high_scoring_signatures) > 0) {
    for (s in high_scoring_signatures) {
      # Generate plot
      plot_obj <- create_dimension_reduction_plots(
        values$data$plot_data[[s]],
        values$data$plot_metadata[[s]],
        input$proband,
        s,
        values$data$age_table,
        values$data$chr_sex_table,
        n_samples_per_group = input$nSamplesPerGroup
      )

      # Save plot to temporary file immediately
      temp_file <- tempfile(pattern = paste0("dimplot_", s, "_"), fileext = ".jpg")
      ggsave(temp_file, plot_obj, width = 16, height = 16, dpi = 100,
             device = "jpeg", quality = 90)

      # Store file path instead of plot object
      dim_plot_files[[s]] <- temp_file
      temp_files_to_cleanup <- c(temp_files_to_cleanup, temp_file)

      # Remove plot object from memory immediately
      rm(plot_obj)
      gc(verbose = FALSE)  # Force garbage collection
    }
  }

  # Convert plot to base64
  prediction_plot_base64 <- plot_to_base64(prediction_plot, width = 12, height = 8)

  # Convert ggplot to base64 images
  cell_prop_base64 <- if(!is.null(cell_prop_plot)) {
    plot_to_base64(cell_prop_plot, width = 10, height = 7)
  } else ""

  chr_sex_base64 <- if(!is.null(chr_sex_plot)) {
    plot_to_base64(chr_sex_plot, width = 10, height = 8)
  } else ""

  # Convert cached dimension plots to base64 by reading from disk
  dim_plots_base64 <- list()
  for (s in names(dim_plot_files)) {
    tryCatch({
      # Read image file and convert to base64 directly
      img_data <- readBin(dim_plot_files[[s]], "raw", file.info(dim_plot_files[[s]])$size)
      base64_string <- base64enc::base64encode(img_data)
      dim_plots_base64[[s]] <- paste0("data:image/jpeg;base64,", base64_string)
    }, error = function(e) {
      warning(paste("Failed to read cached plot for", s, ":", e$message))
      # Skip this plot if it can't be read
    })
  }

  # Create filtered tables
  age_table_filtered <- values$data$age_table[values$data$age_table$Proband %in% input$proband, ]

  # Clean column names for DataTables compatibility 
  names(age_table_filtered) <- gsub("\\.", "_", names(age_table_filtered))

  # Generate HTML
  html_content <- generate_html_template(
    welcome_html = welcome_html,
    references_html = references_html,
    prediction_plot_base64 = prediction_plot_base64,
    prediction_table_data = filtered_results,
    cell_prop_base64 = cell_prop_base64,
    chr_sex_base64 = chr_sex_base64,
    age_table_data = age_table_filtered,
    dim_plots_base64 = dim_plots_base64,
    signatures = input$signatures,
    include_cell_plot = include_cell_plot,
    include_chr_sex_plot = include_chr_sex_plot,
    include_dim_plots = include_dim_plots
  )

  writeLines(html_content, file_path)

  # Clean up temporary files
  for (temp_file in temp_files_to_cleanup) {
    if (file.exists(temp_file)) {
      tryCatch({
        unlink(temp_file)
      }, error = function(e) {
        warning(paste("Failed to clean up temporary file:", temp_file, "-", e$message))
      })
    }
  }
}

#' Convert ggplot to base64 image with optimization
#'
#' @param plot ggplot object
#' @param width Width in inches
#' @param height Height in inches
#' @return Base64 encoded image string
plot_to_base64 <- function(plot, width = 10, height = 8) {
  temp_file <- tempfile(fileext = ".jpg")

  # Use JPEG with compression for smaller file size
  ggsave(temp_file, plot, width = width, height = height, dpi = 100,
         device = "jpeg", quality = 90)

  img_data <- readBin(temp_file, "raw", file.info(temp_file)$size)
  base64_string <- base64enc::base64encode(img_data)

  unlink(temp_file)
  paste0("data:image/jpeg;base64,", base64_string)
}