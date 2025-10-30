#!/usr/bin/env Rscript

#' MethaDory Command Line HTML Report Generator
#'
#' This script generates HTML reports for methylation analysis using MethaDory
#' It processes all SVM classifiers in a given folder and generates a comprehensive report
#'
#' Usage: Rscript generate_methadory_report.R <model_folder> <sample_file> <output_path> [options]
#'
#' Arguments:
#'   model_folder: Path to folder containing SVM model .rds files
#'   sample_file:  Path to .tsv file containing sample data
#'   output_path:  Path where HTML report should be saved
#'
#' Options:
#'   --include-dim-plots     Include dimension reduction plots (default: TRUE)
#'   --include-cell-plots    Include cell deconvolution plots (default: TRUE)
#'   --include-chr-sex       Include chromosomal sex prediction plots (default: TRUE)
#'   --min-psvm              Minimum pSVM threshold for plots (default: 0.05)
#'   --n-imputation-samples  Number of closest samples for imputation (default: 20)
#'   --n-samples-plots       Number of additional samples for visualization (default: 20)
#'   --help                  Show this help message


options(timeout = 2000)
Sys.setenv(R_DEFAULT_INTERNET_TIMEOUT = "2000")

# Add packages that are broken in pixi
packages <- c("FDb.InfiniumMethylation.hg19", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
              "ChAMPdata",  "GenomeInfoDb")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

# Load required libraries
suppressPackageStartupMessages({
  library(BiocParallel)
  library(caret)
  library(circlize)
  library(ComplexHeatmap)
  library(data.table)
  library(DT)
  library(butcher)
  library(EpiDISH)
  library(wateRmelon)
  library(fs)
  library(ggplotify)
  library(methyLImp2)
  library(patchwork)
  library(PCAtools)
  library(plotly)
  library(shiny)
  library(shinybusy)
  library(shinydashboard)
  library(shinyFiles)
  library(tidyverse)
  library(scales)
  library(ggrepel)
  library(htmlwidgets)
  library(base64enc)
  library(jsonlite)
  library(markdown)
})

# Source the modular functions (from MethaDory root)
source("../../src/core/data_processing.R")
source("../../src/core/svm_prediction.R")
source("../../src/core/methylation_analysis.R")
source("../../src/visualization/prediction_plots.R")
source("../../src/visualization/dimension_plots.R")
source("../../src/visualization/analysis_plots.R")
source("../../src/export/html_export.R")
source("../../src/export/table_export.R")

#' Parse command line arguments
parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0 || "--help" %in% args) {
    cat("MethaDory Command Line HTML Report Generator\n\n")
    cat("Usage: Rscript generate_methadory_report.R <model_folder> <sample_file> <output_path> [options]\n\n")
    cat("Arguments:\n")
    cat("  model_folder: Path to folder containing SVM model .rds files\n")
    cat("  sample_file:  Path to .tsv file containing sample data\n")
    cat("  output_path:  Path where HTML report should be saved\n\n")
    cat("Options:\n")
    cat("  --include-dim-plots      Include dimension reduction plots (default: TRUE)\n")
    cat("  --include-cell-plots     Include cell deconvolution plots (default: TRUE)\n")
    cat("  --include-chr-sex        Include chromosomal sex prediction plots (default: TRUE)\n")
    cat("  --min-psvm               Minimum pSVM threshold for plots (default: 0.05)\n")
    cat("  --n-imputation-samples   Number of closest samples for imputation (default: 20)\n")
    cat("  --n-samples-plots        Number of additional samples for visualization (default: 20)\n")
    cat("  --help                   Show this help message\n\n")
    cat("Example:\n")
    cat("  Rscript generate_methadory_report.R ./models sample.tsv report.html\n")
    quit(status = 0)
  }
  
  if (length(args) < 3) {
    stop("ERROR: Missing required arguments. Use --help for usage information.")
  }
  
  # Parse positional arguments
  parsed <- list(
    model_folder = args[1],
    sample_file = args[2],
    output_path = args[3],
    include_dim_plots = TRUE,
    include_cell_plots = TRUE,
    include_chr_sex = TRUE,
    min_psvm = 0.05,
    n_imputation_samples = 20,
    n_samples_plots = 20
  )
  
  # Parse optional arguments
  optional_args <- args[4:length(args)]
  
  if ("--include-dim-plots" %in% optional_args) {
    idx <- which(optional_args == "--include-dim-plots")
    if (idx < length(optional_args)) {
      parsed$include_dim_plots <- as.logical(optional_args[idx + 1])
    }
  }
  
  if ("--include-cell-plots" %in% optional_args) {
    idx <- which(optional_args == "--include-cell-plots")
    if (idx < length(optional_args)) {
      parsed$include_cell_plots <- as.logical(optional_args[idx + 1])
    }
  }
  
  if ("--include-chr-sex" %in% optional_args) {
    idx <- which(optional_args == "--include-chr-sex")
    if (idx < length(optional_args)) {
      parsed$include_chr_sex <- as.logical(optional_args[idx + 1])
    }
  }
  
  
  if ("--min-psvm" %in% optional_args) {
    idx <- which(optional_args == "--min-psvm")
    if (idx < length(optional_args)) {
      parsed$min_psvm <- as.numeric(optional_args[idx + 1])
    }
  }

  if ("--n-imputation-samples" %in% optional_args) {
    idx <- which(optional_args == "--n-imputation-samples")
    if (idx < length(optional_args)) {
      parsed$n_imputation_samples <- as.integer(optional_args[idx + 1])
    }
  }

  if ("--n-samples-plots" %in% optional_args) {
    idx <- which(optional_args == "--n-samples-plots")
    if (idx < length(optional_args)) {
      parsed$n_samples_plots <- as.integer(optional_args[idx + 1])
    }
  }

  return(parsed)
}

#' Validate input arguments
validate_arguments <- function(args) {
  # Check model folder exists
  if (!dir.exists(args$model_folder)) {
    stop(paste("ERROR: Model folder does not exist:", args$model_folder))
  }
  
  # Check sample file exists
  if (!file.exists(args$sample_file)) {
    stop(paste("ERROR: Sample file does not exist:", args$sample_file))
  }
  
  # Check output directory exists
  output_dir <- dirname(args$output_path)
  if (!dir.exists(output_dir)) {
    stop(paste("ERROR: Output directory does not exist:", output_dir))
  }
  
  # Check file extensions
  if (!grepl("\\.(tsv|txt)$", args$sample_file, ignore.case = TRUE)) {
    warning("Sample file should be a .tsv or .txt file")
  }
  
  if (!grepl("\\.html?$", args$output_path, ignore.case = TRUE)) {
    warning("Output file should have .html extension")
  }
  
  # Validate min_psvm range
  if (args$min_psvm < 0 || args$min_psvm > 1) {
    stop("ERROR: min_psvm must be between 0 and 1")
  }
  
  cat(" Arguments validated successfully\n")
}

#' Create CLI-specific HTML export function
create_cli_html_export <- function(data_list, model_dir, output_path, options) {
  
  cat("Generating HTML report...\n")
  
  # Load and encode MethaDory logo
  logo_base64 <- ""
  if (file.exists("../../src/shiny/html_imports/methadory.png")) {
    tryCatch({
      img_data <- readBin("../../src/shiny/html_imports/methadory.png", "raw", file.info("../../src/shiny/html_imports/methadory.png")$size)
      logo_base64 <- base64enc::base64encode(img_data)
    }, error = function(e) {
      cat("Warning: Could not load methadory.png logo\n")
    })
  }
  
  # Create welcome content with logo
  logo_html <- if (nchar(logo_base64) > 0) {
    paste0('<div style="text-align: center; margin-bottom: 30px;">',
           '<img src="data:image/png;base64,', logo_base64, '" ',
           'style="max-width: 400px; height: auto;" alt="MethaDory Logo">',
           '</div>')
  } else ""
  
  # Create welcome content
  # Convert to HTML and add logo
  welcome_html_content <- markdown::markdownToHTML(file = "../../src/shiny/html_imports/welcome_cli.md", fragment.only = TRUE)
  welcome_html <- paste0(logo_html, welcome_html_content)
  references_html <- markdown::markdownToHTML(file = "../../src/shiny/html_imports/references.md", fragment.only = TRUE)
  
  # Filter results for all samples and signatures
  filtered_results <- data_list$results
  all_signatures <- unique(gsub(" ", "_", filtered_results$SVM))
  
  cat("Creating plots...\n")
  
  # Create static prediction plot with all data
  prediction_plot <- create_prediction_plot_static_filtered(
    data_list$results,
    data_list$test_data_ids,
    all_signatures
  )
  
  # Create other plots based on options
  cell_prop_plot <- if(options$include_cell_plots) {
    create_cell_deconv_plot(
      data_list$cell_props,
      data_list$background_data$cellprops,
      data_list$test_data_ids
    )
  } else NULL
  
  chr_sex_plot <- if(options$include_chr_sex) {
    predict_chr_sex_plot(data_list$chr_sex_table, data_list$test_data_ids)
  } else NULL
  
  # Filter signatures with pSVM >= min_psvm for dimension plots
  high_scoring_signatures <- character(0)
  if(options$include_dim_plots) {
    high_scoring_results <- filtered_results[filtered_results$pSVM_average >= options$min_psvm, ]
    high_scoring_signatures <- unique(high_scoring_results$SVM)
    high_scoring_signatures <- gsub(" ", "_", high_scoring_signatures)
    high_scoring_signatures <- high_scoring_signatures[high_scoring_signatures %in% all_signatures]
    
    cat(paste("Found", length(high_scoring_signatures), "signatures with pSVM >=", options$min_psvm, "\n"))
  }
  
  # Create dimension reduction plots using disk caching
  dim_plot_files <- list()
  temp_files_to_cleanup <- c()
  
  if(options$include_dim_plots && length(high_scoring_signatures) > 0) {
    cat("Generating dimension reduction plots (using disk caching)...\n")
    
    for (i in seq_along(high_scoring_signatures)) {
      s <- high_scoring_signatures[i]
      cat(paste("  Processing plot", i, "of", length(high_scoring_signatures), ":", s, "\n"))
      
      # Generate plot with error handling
      tryCatch({
        # Check if required data exists
        if (is.null(data_list$plot_data[[s]]) || is.null(data_list$plot_metadata[[s]])) {
          stop(paste("Missing plot data or metadata for signature:", s))
        }
        
        # Suppress non-critical warnings during plot generation
        plot_obj <- suppressWarnings({
          create_dimension_reduction_plots(
            data_list$plot_data[[s]],
            data_list$plot_metadata[[s]],
            data_list$test_data_ids,
            s,
            data_list$age_table,
            data_list$chr_sex_table,
            n_samples_per_group = options$n_samples_plots
          )
        })
        
        # Validate plot object before saving
        if (is.null(plot_obj) || !inherits(plot_obj, c("gg", "ggplot", "patchwork"))) {
          stop("Generated plot object is invalid or NULL")
        }
        
        # Save to temporary file immediately
        temp_file <- tempfile(pattern = paste0("dimplot_", s, "_"), fileext = ".jpg")
        ggsave(filename = temp_file, plot = plot_obj, width = 16, height = 16, dpi = 100,
               device = "jpeg", quality = 90)
        
        # Store file path
        dim_plot_files[[s]] <- temp_file
        temp_files_to_cleanup <- c(temp_files_to_cleanup, temp_file)
        
        # Clean up memory
        rm(plot_obj)
        gc(verbose = FALSE)
        
        cat("    Successfully generated and cached\n")
        
      }, error = function(e) {
        cat(paste("    Error generating plot for", s, ":", e$message, "\n"))
        cat("    Skipping this signature and continuing...\n")
      })
    }
    
    # Summary of plot generation
    successful_plots <- length(dim_plot_files)
    cat(paste(" Successfully generated", successful_plots, "of", length(high_scoring_signatures), "dimension reduction plots\n"))
  }
  
  cat("Converting plots to base64...\n")
  
  # Convert plots to base64
  prediction_plot_base64 <- plot_to_base64(prediction_plot, width = 12, height = 8)
  
  cell_prop_base64 <- if(!is.null(cell_prop_plot)) {
    plot_to_base64(cell_prop_plot, width = 10, height = 7)
  } else ""
  
  chr_sex_base64 <- if(!is.null(chr_sex_plot)) {
    plot_to_base64(chr_sex_plot, width = 10, height = 8)
  } else ""
  
  # Convert cached dimension plots to base64
  dim_plots_base64 <- list()
  for (s in names(dim_plot_files)) {
    tryCatch({
      img_data <- readBin(dim_plot_files[[s]], "raw", file.info(dim_plot_files[[s]])$size)
      base64_string <- base64enc::base64encode(img_data)
      dim_plots_base64[[s]] <- paste0("data:image/jpeg;base64,", base64_string)
    }, error = function(e) {
      warning(paste("Failed to read cached plot for", s, ":", e$message))
    })
  }
  
  cat("Creating data tables...\n")
  
  # Create filtered tables
  age_table_filtered <- data_list$age_table[data_list$age_table$Proband %in% data_list$test_data_ids, ]
  names(age_table_filtered) <- gsub("\\.", "_", names(age_table_filtered))

  cat("Generating HTML content...\n")
  
  # Generate HTML using existing template function
  html_content <- generate_html_template(
    welcome_html = welcome_html,
    references_html = references_html,
    prediction_plot_base64 = prediction_plot_base64,
    prediction_table_data = filtered_results,
    cell_prop_base64 = cell_prop_base64,
    chr_sex_base64 = chr_sex_base64,
    age_table_data = age_table_filtered,
    dim_plots_base64 = dim_plots_base64,
    signatures = all_signatures,
    include_cell_plot = options$include_cell_plots,
    include_chr_sex_plot = options$include_chr_sex,
    include_dim_plots = options$include_dim_plots
  )
  
  # Write HTML file
  writeLines(html_content, output_path)
  
  # Clean up temporary files
  cat("Cleaning up temporary files...\n")
  for (temp_file in temp_files_to_cleanup) {
    if (file.exists(temp_file)) {
      tryCatch({
        unlink(temp_file)
      }, error = function(e) {
        warning(paste("Failed to clean up temporary file:", temp_file))
      })
    }
  }
  
  cat(paste(" HTML report successfully generated:", output_path, "\n"))
}

main <- function() {
  cat("=== MethaDory Command Line HTML Report Generator ===\n\n")
  
  # Parse and validate arguments
  args <- parse_arguments()
  validate_arguments(args)
  
  cat("Starting analysis...\n")
  
  # Load test data
  cat("Loading sample data...\n")
  data_list <- load_test_data(args$sample_file)
  data_list$sample_file <- args$sample_file  # Store for reference
  
  # Load background data and models
  cat("Loading background data and SVM models...\n")
  background_data <- load_background_data(args$model_folder)
  
  # Prepare and perform imputation
  cat("Performing data imputation...\n")
  cat("Using", args$n_imputation_samples, "closest samples for imputation\n")
  imputation_data <- prepare_imputation_data(
    data_list$test_data,
    background_data$imputation_background,
    background_data$svm,
    n_closest = args$n_imputation_samples
  )
  imputed_data <- perform_imputation(imputation_data$test_data, data_list$test_data_ids)
  
  # Load additional data for plotting
  cat("Loading signature data...\n")
  beta_sig_data <- load_beta_signatures()
  
  # Prepare inference and make predictions
  cat("Making SVM predictions...\n")
  inference_data <- prepare_inference_data(imputed_data, background_data$svm)
  results <- make_predictions(inference_data, background_data$svm, data_list$test_data_ids)
  
  # Load real cases data for plotting
  real_cases_data <- load_real_cases(beta_sig_data$signatures)
  
  # Use insilico metadata from beta_sig_data and prepare it for plotting
  insilico_meta <- beta_sig_data$insilico_meta[, c("geo_accession", "platform_id", "Sex", "AgeGroup")]
  names(insilico_meta) <- c("IDs", "Platform", "Sex", "AgeGroup")
  
  # Prepare plot data
  cat("Preparing plot data...\n")
  plot_data <- prepare_plot_data(
    beta_sig_data$signatures,
    beta_sig_data$insilico_beta,
    imputed_data,
    real_cases_data$real_cases_beta
  )
  plot_metadata <- prepare_plot_metadata(
    beta_sig_data$signatures,
    plot_data,
    imputed_data,
    real_cases_data$real_cases_meta,
    insilico_meta
  )
  
  # Create additional analyses
  cat("Performing additional analyses...\n")
  cell_props <- create_cell_deconv_table(data_list$test_data)
  chr_sex_table <- predict_chr_sex_table(data_list$test_data)
  age_table <- predict_age(data_list$test_data)
  
  # Combine all data
  full_data <- list(
    data_list = data_list,
    background_data = background_data,
    imputed_data = imputed_data,
    beta_sig_data = beta_sig_data,
    inference_data = inference_data,
    results = results,
    real_cases_data = real_cases_data,
    plot_data = plot_data,
    plot_metadata = plot_metadata,
    cell_props = cell_props,
    chr_sex_table = chr_sex_table,
    age_table = age_table,
    test_data_ids = data_list$test_data_ids,
    sample_file = args$sample_file
  )
  
  # Generate HTML report
  create_cli_html_export(full_data, args$model_folder, args$output_path, args)
  
  cat("\n=== Analysis Complete ===\n")
}

# Execute main function if script is run directly
if (!interactive()) {
  tryCatch({
    main()
  }, error = function(e) {
    cat(paste("ERROR:", e$message, "\n"))
    quit(status = 1)
  })
}
