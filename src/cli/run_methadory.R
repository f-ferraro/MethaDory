#!/usr/bin/env Rscript

#' MethaDory Command Line Analysis Tool
#'
#' This script runs the complete MethaDory methylation analysis pipeline
#' It processes SVM classifiers and generates Excel tables, PDF plots, and optionally imputed data
#'
#' Usage: Rscript run_methadory.R <model_folder> <sample_file> <output_prefix> [options]
#'
#' Arguments:
#'   model_folder: Path to folder containing SVM model .rds files
#'   sample_file:  Path to .tsv file containing sample data
#'   output_prefix: Prefix for output files (will create .xlsx and .pdf files)
#'
#' Options:
#'   --include-dim-plots     Include dimension reduction plots in PDF (default: TRUE)
#'   --include-cell-plots    Include cell deconvolution plots in PDF (default: TRUE)
#'   --include-chr-sex       Include chromosomal sex prediction plots in PDF (default: TRUE)
#'   --min-psvm              Minimum pSVM threshold for dimension plots (default: 0.05)
#'   --n-imputation-samples  Number of closest samples for imputation (default: 20)
#'   --n-samples-plots       Number of additional samples for visualization (default: 20)
#'   --export-imputed        Export imputed methylation data (default: FALSE)
#'   --help                  Show this help message

# Load packages required for MethaDory
# Add packages that are broken in pixi
packages <- c("ChAMPdata", "FDb.InfiniumMethylation.hg19", "GenomeInfoDb", "IlluminaHumanMethylation450kanno.ilmn12.hg19")

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
  library(openxlsx)
})

# Source the modular functions
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
    cat("MethaDory Command Line Analysis Tool\n\n")
    cat("Usage: Rscript run_methadory.R <model_folder> <sample_file> <output_prefix> [options]\n\n")
    cat("Arguments:\n")
    cat("  model_folder:   Path to folder containing SVM model .rds files\n")
    cat("  sample_file:    Path to .tsv file containing sample data\n")
    cat("  output_prefix:  Prefix for output files\n\n")
    cat("Options:\n")
    cat("  --include-dim-plots      Include dimension reduction plots in PDF (default: TRUE)\n")
    cat("  --include-cell-plots     Include cell deconvolution plots in PDF (default: TRUE)\n")
    cat("  --include-chr-sex        Include chromosomal sex prediction plots in PDF (default: TRUE)\n")
    cat("  --min-psvm               Minimum pSVM threshold for dimension plots (default: 0.05)\n")
    cat("  --n-imputation-samples   Number of closest samples for imputation (default: 20)\n")
    cat("  --n-samples-plots        Number of additional samples for visualization (default: 20)\n")
    cat("  --export-imputed         Export imputed methylation data including controls and cases (default: FALSE)\n")
    cat("  --help                   Show this help message\n\n")
    cat("Output files:\n")
    cat("  <output_prefix>.xlsx            Excel workbook with all results\n")
    cat("  <output_prefix>.pdf             PDF with all plots\n")
    cat("  <output_prefix>_imputed.tsv     Imputed data: user samples + controls + real cases (if --export-imputed TRUE)\n\n")
    cat("Example:\n")
    cat("  Rscript run_methadory.R ./models sample.tsv results --min-psvm 0.1 --export-imputed TRUE\n")
    quit(status = 0)
  }

  if (length(args) < 3) {
    stop("ERROR: Missing required arguments. Use --help for usage information.")
  }

  # Parse positional arguments
  parsed <- list(
    model_folder = args[1],
    sample_file = args[2],
    output_prefix = args[3],
    include_dim_plots = TRUE,
    include_cell_plots = TRUE,
    include_chr_sex = TRUE,
    min_psvm = 0.05,
    n_imputation_samples = 20,
    n_samples_plots = 20,
    export_imputed = FALSE
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

  if ("--export-imputed" %in% optional_args) {
    idx <- which(optional_args == "--export-imputed")
    if (idx < length(optional_args)) {
      parsed$export_imputed <- as.logical(optional_args[idx + 1])
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
  output_dir <- dirname(args$output_prefix)
  if (output_dir != "." && !dir.exists(output_dir)) {
    stop(paste("ERROR: Output directory does not exist:", output_dir))
  }

  # Check file extensions
  if (!grepl("\\.(tsv|txt)$", args$sample_file, ignore.case = TRUE)) {
    warning("Sample file should be a .tsv or .txt file")
  }

  # Validate min_psvm range
  if (args$min_psvm < 0 || args$min_psvm > 1) {
    stop("ERROR: min_psvm must be between 0 and 1")
  }

  cat(" Arguments validated successfully\n")
}

#' Create Excel tables export
create_excel_export <- function(data_list, output_path) {
  cat("Generating Excel tables...\n")

  # Prepare prediction results table
  prediction_results <- data_list$results

  # Create wide format for easier reading
  prediction_wide <- prediction_results %>%
    select(SampleID, SVM, pSVM_average, pSVM_sd) %>%
    pivot_wider(
      id_cols = SampleID,
      names_from = SVM,
      values_from = c(pSVM_average, pSVM_sd),
      names_sep = "_"
    )

  # Clean column names
  names(prediction_wide) <- gsub("pSVM_average_", "pSVM_", names(prediction_wide))
  names(prediction_wide) <- gsub("pSVM_sd_", "pSVM_SD_", names(prediction_wide))

  # Prepare cell proportions table
  cell_props_wide <- data_list$cell_props %>%
    pivot_wider(
      id_cols = Proband,
      names_from = CellType,
      values_from = CellProp
    )

  # Prepare age predictions
  age_table_clean <- data_list$age_table
  names(age_table_clean) <- gsub("\\.", "_", names(age_table_clean))

  # Prepare chromosomal sex predictions
  chr_sex_clean <- data_list$chr_sex_table

  # Create signatures summary
  signatures_summary <- prediction_results %>%
    group_by(SVM) %>%
    summarise(
      Max_pSVM = round(max(pSVM_average, na.rm = TRUE), 3),
      Min_pSVM = round(min(pSVM_average, na.rm = TRUE), 3),
      Avg_pSVM = round(mean(pSVM_average, na.rm = TRUE), 3),
      N_Samples = n(),
      High_Confidence = sum(pSVM_average >= 0.5),
      Medium_Confidence = sum(pSVM_average >= 0.25 & pSVM_average < 0.5),
      Low_Confidence = sum(pSVM_average >= 0.05 & pSVM_average < 0.25),
      Very_Low = sum(pSVM_average < 0.05),
      .groups = "drop"
    ) %>%
    arrange(desc(Avg_pSVM))

  # Create workbook with separate sheets
  wb <- createWorkbook()

  # Create separate sheets
  addWorksheet(wb, "SVM_Predictions_Wide")
  addWorksheet(wb, "SVM_Predictions_Long")
  addWorksheet(wb, "Signatures_Summary")
  addWorksheet(wb, "Cell_Proportions")
  addWorksheet(wb, "Methylation_Age")
  addWorksheet(wb, "Chromosomal_Sex")

  # Write data to sheets
  writeData(wb, "SVM_Predictions_Wide", prediction_wide, rowNames = FALSE)
  writeData(wb, "SVM_Predictions_Long", prediction_results, rowNames = FALSE)
  writeData(wb, "Signatures_Summary", signatures_summary, rowNames = FALSE)
  writeData(wb, "Cell_Proportions", cell_props_wide, rowNames = FALSE)
  writeData(wb, "Methylation_Age", age_table_clean, rowNames = FALSE)
  writeData(wb, "Chromosomal_Sex", chr_sex_clean, rowNames = FALSE)

  # Add conditional formatting for pSVM values in summary
  conditionalFormatting(wb, "Signatures_Summary", cols = 2:4, rows = 2:(nrow(signatures_summary) + 1),
                       rule = ">=0.5", style = createStyle(bgFill = "#9E1E05", fontColour = "white"))
  conditionalFormatting(wb, "Signatures_Summary", cols = 2:4, rows = 2:(nrow(signatures_summary) + 1),
                       rule = ">=0.25", style = createStyle(bgFill = "#FAD302"))

  # Save workbook
  saveWorkbook(wb, output_path, overwrite = TRUE)
  cat(paste("Excel file saved:", output_path, "\n"))
}

#' Create PDF plots export
create_pdf_export <- function(data_list, output_path, options) {
  cat("Generating PDF plots...\n")

  # Open PDF device
  pdf(output_path, width = 16, height = 12)

  # Create prediction plot
  all_signatures <- unique(gsub(" ", "_", data_list$results$SVM))
  prediction_plot <- create_prediction_plot(
    data_list$results,
    data_list$test_data_ids,
    all_signatures,
    100 + (length(all_signatures) * 100)
  )

  print(prediction_plot + ggtitle("MethaDory SVM Prediction Results"))

  # Add cell deconvolution plot if requested
  if (options$include_cell_plots) {
    cell_plot <- create_cell_deconv_plot(
      data_list$cell_props,
      data_list$background_data$cellprops,
      data_list$test_data_ids
    )
    print(cell_plot)
  }

  # Add chromosomal sex plot if requested
  if (options$include_chr_sex) {
    chr_sex_plot <- predict_chr_sex_plot(data_list$chr_sex_table, data_list$test_data_ids)
    print(chr_sex_plot)
  }

  # Add dimension reduction plots if requested
  if (options$include_dim_plots) {
    # Filter signatures with pSVM >= min_psvm
    high_scoring_results <- data_list$results[data_list$results$pSVM_average >= options$min_psvm, ]
    high_scoring_signatures <- unique(high_scoring_results$SVM)
    high_scoring_signatures <- gsub(" ", "_", high_scoring_signatures)
    high_scoring_signatures <- high_scoring_signatures[high_scoring_signatures %in% all_signatures]

    cat(paste("  Including", length(high_scoring_signatures), "dimension reduction plots (pSVM >=", options$min_psvm, ")\n"))

    if (length(high_scoring_signatures) > 0) {
      for (i in seq_along(high_scoring_signatures)) {
        s <- high_scoring_signatures[i]
        cat(paste("    Processing plot", i, "of", length(high_scoring_signatures), ":", s, "\n"))

        tryCatch({
          # Check if required data exists
          if (is.null(data_list$plot_data[[s]]) || is.null(data_list$plot_metadata[[s]])) {
            cat(paste("      Skipping", s, "- missing plot data\n"))
            next
          }

          # Generate dimension reduction plot
          dim_plot <- suppressWarnings({
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

          if (!is.null(dim_plot)) {
            print(dim_plot)
            cat("      Successfully added to PDF\n")
          } else {
            cat(paste("      Skipping", s, "- plot generation returned NULL\n"))
          }

        }, error = function(e) {
          cat(paste("      Error generating plot for", s, ":", e$message, "\n"))
        })
      }
    } else {
      # Add a message page if no high-scoring signatures
      grid.newpage()
      grid.text(paste("No signatures with pSVM >=", options$min_psvm, "found.\nDimension reduction plots not generated."),
               x = 0.5, y = 0.5, just = "center", gp = gpar(fontsize = 16))
    }
  }

  dev.off()
  cat(paste(" PDF file saved:", output_path, "\n"))
}


main <- function() {
  cat("=== MethaDory Command Line Analysis Tool ===\n\n")

  # Parse and validate arguments
  args <- parse_arguments()
  validate_arguments(args)

  cat("Starting analysis...\n")

  # Load test data
  cat("Loading sample data...\n")
  data_list <- load_test_data(args$sample_file)
  data_list$sample_file <- args$sample_file

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

  # Use insilico metadata from beta_sig_data
  # Include Sex and AgeGroup from the metadata
  insilico_meta <- beta_sig_data$insilico_meta[, c("geo_accession", "platform_id", "Sex", "AgeGroup")]
  names(insilico_meta) <- c("IDs", "Platform", "Sex", "AgeGroup")

  # Prepare plot data 
  if (args$include_dim_plots) {
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
  } else {
    plot_data <- NULL
    plot_metadata <- NULL
  }

  # Create additional analyses
  cat("Performing additional analyses...\n")
  cell_props <- create_cell_deconv_table(data_list$test_data)
  chr_sex_table <- predict_chr_sex_table(data_list$test_data)
  age_table <- predict_age(data_list$test_data)

  # Combine all data
  full_data <- list(
    data_list = data_list,
    background_data = background_data,
    results = results,
    plot_data = plot_data,
    plot_metadata = plot_metadata,
    cell_props = cell_props,
    chr_sex_table = chr_sex_table,
    age_table = age_table,
    test_data_ids = data_list$test_data_ids
  )

  # Generate outputs
  excel_path <- paste0(args$output_prefix, ".xlsx")
  pdf_path <- paste0(args$output_prefix, ".pdf")

  # Create Excel export
  create_excel_export(full_data, excel_path)

  # Create PDF export
  create_pdf_export(full_data, pdf_path, args)

  # Export imputed data if requested 
  if (args$export_imputed) {
    imputed_path <- paste0(args$output_prefix, "_imputed.tsv")
    cat("\nExporting imputed methylation data (user samples + controls + cases)...\n")

    # Combine user imputed data with background controls and real cases
    combined_imputed <- imputed_data

    # Add control samples from insilico_beta 
    if (!is.null(beta_sig_data$insilico_beta)) {
      # Remove IlmnID column if present
      controls_beta <- beta_sig_data$insilico_beta
      if ("IlmnID" %in% names(controls_beta)) {
        ilmn_col <- controls_beta$IlmnID
        controls_beta$IlmnID <- NULL
      } else {
        ilmn_col <- rownames(controls_beta)
      }
      controls_beta$IlmnID <- ilmn_col

      # Merge with user data
      combined_imputed <- merge(combined_imputed, controls_beta, by = "IlmnID", all = TRUE)
      cat("  Added", ncol(controls_beta) - 1, "control samples\n")
    }

    # Add real cases beta values 
    if (!is.null(real_cases_data$real_cases_beta) && length(real_cases_data$real_cases_beta) > 0) {
      # Real cases are stored per signature, need to get unique samples across all signatures
      # Extract all unique samples from the first signature's data
      first_sig <- real_cases_data$real_cases_beta[[1]]
      if (!is.null(first_sig) && nrow(first_sig) > 0) {
        # Load the full real cases beta matrix 
        tryCatch({
          full_real_cases_beta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.beta.rds")
          full_real_cases_meta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.meta.rds")

          # Filter to only real cases
          real_case_ids <- full_real_cases_meta[full_real_cases_meta$RealLabel != "control", ]$geo_accession
          real_cases_only <- full_real_cases_beta[, colnames(full_real_cases_beta) %in% real_case_ids, drop = FALSE]

          real_cases_only$IlmnID <- rownames(real_cases_only)

          # Merge with combined data
          combined_imputed <- merge(combined_imputed, real_cases_only, by = "IlmnID", all = TRUE)
          cat("  Added", ncol(real_cases_only) - 1, "real case samples\n")
        }, error = function(e) {
          cat("  Warning: Could not load real cases data:", e$message, "\n")
        })
      }
    }

    # Write combined data
    write.table(combined_imputed, file = imputed_path, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(paste(" Imputed data saved:", imputed_path, "\n"))
    cat(paste("  Total samples:", ncol(combined_imputed) - 1, "\n"))
    cat(paste("  Total CpGs:", nrow(combined_imputed), "\n"))
  }

  cat("\n=== Analysis Complete ===\n")
  cat(paste("Files created:\n"))
  cat(paste("  Excel tables:", excel_path, "\n"))
  cat(paste("  PDF plots:  ", pdf_path, "\n"))
  if (args$export_imputed) {
    cat(paste("  Imputed data:", paste0(args$output_prefix, "_imputed.tsv"), "\n"))
  }
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