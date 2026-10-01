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
#'   output_path:  Where the report(s) should be saved. With a single sample in the
#'                 input this is the .html file to write. With several samples it is
#'                 the output DIRECTORY: the pipeline is run once and each sample is
#'                 then rendered on its own into <output_dir>/<SampleName>.MethaDory-output.html
#'
#' Next to the report(s) the result tables are written as one .xlsx workbook
#' covering every sample of the run: <report>.xlsx for a single sample,
#' <output_dir>/<input file name>.MethaDory-output.xlsx (all samples together)
#' for a multi-sample input. Disable with --export-xlsx FALSE.
#' With --export-imputed TRUE the imputed beta matrix is written there too, as
#' <same name>_imputed.tsv.
#'
#' Options:
#'   --include-dim-plots     Include dimension reduction plots (default: TRUE)
#'   --include-cell-plots    Include cell deconvolution plots (default: TRUE)
#'   --include-chr-sex       Include chromosomal sex prediction plots (default: TRUE)
#'   --min-p                 Keep a signature for the per-signature dimension plots when its combined
#'                           score (pCombined, the mean of the SVM and NNET scores) is at or above
#'                           this value, between 0 and 1 (default: 0.20). Tables and the prediction
#'                           plot always show every signature.
#'   --n-imputation-samples  Number of closest samples for imputation (default: 20)
#'   --n-samples-plots       Number of additional samples for visualization (default: 20)
#'   --export-xlsx           Also write the result tables as an .xlsx workbook (default: TRUE)
#'   --export-imputed        Also write the imputed methylation data: user samples + controls +
#'                           real cases, as <report>_imputed.tsv (default: FALSE)
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
  # tidyverse components actually used (avoid pulling the meta-package).
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(tibble)
  library(scales)
  library(ggrepel)
  library(htmlwidgets)
  library(base64enc)
  library(jsonlite)
  library(markdown)
  library(reticulate)
})

# Source the modular functions (from MethaDory root)
source("../../src/core/data_processing.R")
source("../../src/core/svm_prediction.R")
source("../../src/core/nnet_prediction.R")
source("../../src/core/methylation_analysis.R")
source("../../src/core/sample_qc.R")
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
    cat("  output_path:  Single-sample input: the .html file to write.\n")
    cat("                Multi-sample input: the output DIRECTORY. One report per sample is\n")
    cat("                written as <output_dir>/<SampleName>.MethaDory-output.html\n")
    cat("                The result tables are also written as one .xlsx workbook next to the\n")
    cat("                report(s): <report>.xlsx, or for a multi-sample input\n")
    cat("                <output_dir>/<input file name>.MethaDory-output.xlsx\n\n")
    cat("Options:\n")
    cat("  --include-dim-plots      Include dimension reduction plots (default: TRUE)\n")
    cat("  --include-cell-plots     Include cell deconvolution plots (default: TRUE)\n")
    cat("  --include-chr-sex        Include chromosomal sex prediction plots (default: TRUE)\n")
    cat("  --min-p                  Keep a signature for the per-signature dimension plots when its combined\n")
    cat("                           score (pCombined, the mean of the SVM and NNET scores) is at or above\n")
    cat("                           this value, between 0 and 1 (default: 0.20). Tables and the prediction\n")
    cat("                           plot always show every signature.\n")
    cat("  --n-imputation-samples   Number of closest samples for imputation (default: 20)\n")
    cat("  --n-samples-plots        Number of additional samples for visualization (default: 20)\n")
    cat("  --export-xlsx            Also write the result tables as an .xlsx workbook (default: TRUE)\n")
    cat("  --export-imputed         Also write the imputed methylation data (user samples + controls +\n")
    cat("                           real cases) as <report>_imputed.tsv (default: FALSE)\n")
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
    min_p = 0.20,
    n_imputation_samples = 20,
    n_samples_plots = 20,
    export_xlsx = TRUE,
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
  
  
  # Renamed from --min-psvm: the threshold is on the combined mean(SVM, NNET)
  # score, not on the SVM alone. The old spelling is refused rather than
  # ignored, so a script still using it cannot silently run at the default.
  if ("--min-psvm" %in% optional_args) {
    stop("ERROR: --min-psvm was renamed to --min-p")
  }

  if ("--min-p" %in% optional_args) {
    idx <- which(optional_args == "--min-p")
    if (idx < length(optional_args)) {
      parsed$min_p <- as.numeric(optional_args[idx + 1])
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

  if ("--export-xlsx" %in% optional_args) {
    idx <- which(optional_args == "--export-xlsx")
    if (idx < length(optional_args)) {
      parsed$export_xlsx <- as.logical(optional_args[idx + 1])
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
  
  # The output argument is either an .html file (single sample) or a directory
  # (multi-sample). Which one applies is only known once the input is loaded, so
  # here we just check that the enclosing directory is reachable.
  output_dir <- if (grepl("\\.html?$", args$output_path, ignore.case = TRUE)) {
    dirname(args$output_path)
  } else {
    args$output_path
  }
  if (!dir.exists(output_dir)) {
    stop(paste("ERROR: Output directory does not exist:", output_dir))
  }
  
  # Check file extensions
  if (!grepl("\\.(tsv|txt)$", args$sample_file, ignore.case = TRUE)) {
    warning("Sample file should be a .tsv or .txt file")
  }

  if (is.na(args$export_xlsx)) {
    stop("ERROR: --export-xlsx must be TRUE or FALSE")
  }
  if (is.na(args$export_imputed)) {
    stop("ERROR: --export-imputed must be TRUE or FALSE")
  }
  
  # Validate min_p range
  if (args$min_p < 0 || args$min_p > 1) {
    stop("ERROR: --min-p must be between 0 and 1")
  }
  
  cat(" Arguments validated successfully\n")
}

#' Resolve the output file for each sample.
#'
#' A single-sample input keeps the historical behaviour: `output_path` is the
#' .html file to write. As soon as the input holds more than one sample the
#' reports have to be kept apart, so `output_path` is taken as the output
#' directory and each sample gets `<SampleName>.MethaDory-output.html`. Passing
#' an .html path together with a multi-sample input is tolerated - its directory
#' is used and the per-sample naming applies.
#'
#' @param output_path the third positional argument.
#' @param sample_ids  sample IDs from the input file.
#' @return character vector of file paths, named by sample ID.
resolve_output_paths <- function(output_path, sample_ids) {
  looks_like_html <- grepl("\\.html?$", output_path, ignore.case = TRUE)

  if (length(sample_ids) == 1 && looks_like_html) {
    return(setNames(output_path, sample_ids))
  }

  output_dir <- if (looks_like_html) dirname(output_path) else output_path
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  if (looks_like_html) {
    cat("Multi-sample input: '", basename(output_path), "' is ignored, one report ",
        "per sample is written to ", output_dir, "\n", sep = "")
  }

  # Sample IDs come from a user-supplied header, so keep only characters that
  # are safe in a filename on every platform.
  safe_ids <- gsub("[^A-Za-z0-9._-]+", "_", sample_ids)
  setNames(file.path(output_dir, paste0(safe_ids, ".MethaDory-output.html")),
           sample_ids)
}

#' Resolve the .xlsx workbook written next to the report(s).
#'
#' The workbook holds every sample of the run, so there is one per run rather
#' than one per report. A single-sample run names
#' it after its report (report.html -> report.xlsx). A multi-sample run has no
#' single report to be named after, so it takes the input file's name and the
#' same suffix the per-sample reports carry.
#'
#' @param output_paths per-sample report paths, from resolve_output_paths().
#' @param sample_file  the input .tsv, used to name a multi-sample workbook.
#' @return path of the .xlsx file.
resolve_xlsx_path <- function(output_paths, sample_file) {
  if (length(output_paths) == 1) {
    return(sub("\\.html?$", ".xlsx", output_paths[[1]], ignore.case = TRUE))
  }
  input_name <- tools::file_path_sans_ext(basename(sample_file))
  file.path(dirname(output_paths[[1]]),
            paste0(gsub("[^A-Za-z0-9._-]+", "_", input_name),
                   ".MethaDory-output.xlsx"))
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
  
  # Restrict the report to the sample(s) it is being written for. With a
  # multi-sample input this function is called once per sample, so everything
  # downstream (tables, prediction plot, dimension plots) covers that sample
  # alone; the other plotting helpers filter on test_data_ids themselves.
  filtered_results <- data_list$results[
    data_list$results$SampleID %in% data_list$test_data_ids, , drop = FALSE]
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
  
  qc_pca_plot <- create_qc_pca_plot(data_list$qc_pca, data_list$test_data_ids)
  qc_density_plot <- create_qc_density_plot(data_list$qc_pca, data_list$test_data_ids)

  # Filter signatures with combined mean(SVM, NNET) score >= min_p for dimension plots
  high_scoring_signatures <- character(0)
  if(options$include_dim_plots) {
    score_col <- if ("mean_case" %in% names(filtered_results)) "mean_case" else "pSVM_average"
    high_scoring_results <- filtered_results[!is.na(filtered_results[[score_col]]) &
                                                filtered_results[[score_col]] >= options$min_p, ]
    high_scoring_signatures <- unique(high_scoring_results$SVM)
    high_scoring_signatures <- gsub(" ", "_", high_scoring_signatures)
    high_scoring_signatures <- high_scoring_signatures[high_scoring_signatures %in% all_signatures]

    cat(paste("Found", length(high_scoring_signatures), "signatures with", score_col, ">=", options$min_p, "\n"))
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
        ggsave(filename = temp_file, plot = plot_obj,
               width = 16,
               height = 16 * DIMENSION_PLOT_HEIGHT_SCALE, dpi = 100,
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
  
  # QC tab: cell proportions with the chromosomal sex prediction to their right
  qc_cell_sex_plot <- create_qc_cell_sex_plot(cell_prop_plot, chr_sex_plot)
  qc_cell_sex_base64 <- if(!is.null(qc_cell_sex_plot)) {
    plot_to_base64(qc_cell_sex_plot,
                   width = if(!is.null(cell_prop_plot) && !is.null(chr_sex_plot)) 20 else 10,
                   height = 8)
  } else ""
  
  qc_pca_base64 <- if(!is.null(qc_pca_plot)) {
    plot_to_base64(qc_pca_plot, width = 14, height = 7)
  } else ""
  qc_density_base64 <- if(!is.null(qc_density_plot)) {
    plot_to_base64(qc_density_plot, width = 14, height = 6)
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
    age_table_data = age_table_filtered,
    dim_plots_base64 = dim_plots_base64,
    signatures = all_signatures,
    include_dim_plots = options$include_dim_plots,
    signature_version = signature_version(),
    qc_missing_html = create_qc_missing_html(data_list$qc_missing, data_list$test_data_ids),
    qc_cell_sex_base64 = qc_cell_sex_base64,
    qc_pca_base64 = qc_pca_base64,
    qc_density_base64 = qc_density_base64
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
  
  # Sample QC: PCA of each proband with the controls, on pre-imputation betas.
  # A failure here costs the QC tab only, not the run.
  cat("Computing sample QC PCA...\n")
  qc_pca <- tryCatch(
    compute_qc_pca(data_list$test_data, background_data$imputation_background,
                   data_list$test_data_ids),
    error = function(e) {
      cat(paste("ERROR computing the sample QC PCA:", e$message, "\n"))
      list()
    })

  # Prepare and perform imputation
  cat("Performing data imputation...\n")
  cat("Using", args$n_imputation_samples, "closest samples for imputation\n")
  extra_nnet_cpgs <- get_nnet_required_cpgs(background_data$nnet_files)
  all_signature_cpgs <- unique(read.delim(
    "../../data/support_files/merged_signatures_90DMRs.tsv",
    header = TRUE)$ProbeID)
  imputation_data <- prepare_imputation_data(
    data_list$test_data,
    background_data$imputation_background,
    background_data$svm,
    n_closest = args$n_imputation_samples,
    extra_cpgs = c(extra_nnet_cpgs, all_signature_cpgs)
  )
  # QC: share of the model CpGs missing before imputation (PASS/WARNING/FAIL)
  qc_missing <- compute_pre_imputation_missing(imputation_data$test_data, data_list$test_data_ids)
  imputed_data <- perform_imputation(imputation_data$test_data, data_list$test_data_ids)

  # Load additional data for plotting
  cat("Loading signature data...\n")
  beta_sig_data <- load_beta_signatures()

  # QC: pre-imputation NA fraction per (sample x signature).
  na_pct_pre <- compute_pre_imputation_na_pct(imputation_data$test_data,
                                              beta_sig_data$signatures,
                                              data_list$test_data_ids)

  # SVM + NNET predictions, combined into the metapredictor summary
  cat("Making SVM predictions...\n")
  inference_data <- prepare_inference_data(imputed_data, background_data$svm)
  svm_results <- make_predictions(inference_data, background_data$svm, data_list$test_data_ids)
  cat("Making NNET predictions on", length(background_data$nnet_files), "checkpoints...\n")
  nnet_results <- make_nnet_predictions(imputed_data, background_data$nnet_files,
                                        data_list$test_data_ids)
  results <- combine_svm_nnet_results(svm_results, nnet_results, na_pct = na_pct_pre)
  
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
    qc_pca = qc_pca,
    qc_missing = qc_missing,
    test_data_ids = data_list$test_data_ids,
    sample_file = args$sample_file
  )
  
  # Generate the HTML report(s). The pipeline above ran once over the whole
  # input; the rendering is then done one sample at a time so that each report
  # is self-contained and covers a single proband.
  sample_ids <- data_list$test_data_ids
  output_paths <- resolve_output_paths(args$output_path, sample_ids)

  # The result tables, once for the whole run and before the rendering: they
  # need nothing the reports produce, and a report failing further down should
  # not cost the numbers. A failure here is reported but does not stop the
  # reports either.
  xlsx_path <- NULL
  if (isTRUE(args$export_xlsx)) {
    xlsx_path <- resolve_xlsx_path(output_paths, args$sample_file)
    tryCatch({
      create_excel_export(full_data, xlsx_path)
    }, error = function(e) {
      cat(paste("ERROR writing the Excel tables:", e$message, "\n"))
      xlsx_path <<- NULL
    })
  }

  # Optional: the imputed beta matrix, named after the workbook.
  imputed_path <- NULL
  if (isTRUE(args$export_imputed)) {
    imputed_path <- sub("\\.xlsx$", "_imputed.tsv",
                        resolve_xlsx_path(output_paths, args$sample_file))
    tryCatch({
      export_imputed_data(imputed_data, beta_sig_data$insilico_beta, imputed_path)
    }, error = function(e) {
      cat(paste("ERROR writing the imputed data:", e$message, "\n"))
      imputed_path <<- NULL
    })
  }

  cat("\nRendering", length(sample_ids), "report(s)...\n")

  failed <- character(0)
  for (i in seq_along(sample_ids)) {
    id <- sample_ids[i]
    cat(sprintf("\n--- Sample %d of %d: %s ---\n", i, length(sample_ids), id))

    sample_data <- full_data
    sample_data$test_data_ids <- id

    tryCatch({
      create_cli_html_export(sample_data, args$model_folder, output_paths[[id]], args)
    }, error = function(e) {
      cat(paste("ERROR generating report for", id, ":", e$message, "\n"))
      failed <<- c(failed, id)
    })

    # Each report holds its own base64-encoded figures; drop them before moving on.
    rm(sample_data)
    gc(verbose = FALSE)
  }

  cat("\n=== Analysis Complete ===\n")
  cat(paste(length(sample_ids) - length(failed), "of", length(sample_ids),
            "report(s) written\n"))
  if (!is.null(xlsx_path)) cat(paste("Excel tables:", xlsx_path, "\n"))
  if (!is.null(imputed_path)) cat(paste("Imputed data:", imputed_path, "\n"))
  if (length(failed) > 0) {
    cat(paste("Failed:", paste(failed, collapse = ", "), "\n"))
    quit(status = 1)
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
