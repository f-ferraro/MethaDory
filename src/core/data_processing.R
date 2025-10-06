#' Load and prepare test data
#'
#' @param file_path Path to the test data file
#' @return List containing the test data and test data IDs
load_test_data <- function(file_path) {
  test_data_user <- read.delim(file_path, header = TRUE)

  test_data <- test_data_user %>% relocate(IlmnID)

  return(list(
    test_data = test_data,
    test_data_ids = setdiff(names(test_data), "IlmnID")
  ))
}

#' Load background data and SVM models
#'
#' @param model_dir Directory containing the SVM models
#' @return List containing imputation background and SVM models
load_background_data <- function(model_dir) {
  tryCatch({
    # Load imputation background
    imputation_background = readRDS("../../data/imputationsamples/samples.beta.rds")
    # cat("Loaded imputation background:", dim(imputation_background), "\n")

    # Load models from the selected directory
    models <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)

    if (length(models) == 0) {
      stop("No .rds model files found in the selected directory")
    }

    # cat("Loading", length(models), "SVM models...\n")

    # Load SVM models
    svm <- lapply(models, function(model_path) {
      tryCatch({
        butcher(readRDS(model_path))
      }, error = function(e) {
        warning(paste("Error loading model:", basename(model_path), "-", e$message))
        NULL
      })
    })

    # Remove any NULL entries from failed loads
    svm <- Filter(Negate(is.null), svm)

    if (length(svm) == 0) {
      stop("No valid SVM models could be loaded from the selected directory")
    }

    # Name the models based on their filenames
    valid_model_names <- gsub("\\.rds$", "", basename(models[1:length(svm)]))
    names(svm) <- valid_model_names

    cat("Successfully loaded", length(svm), "SVM models\n")

    # Read cell props from the background
    cellprops = readRDS("../../data/support_files/background_training.cellprops.rds")

    return(list(
      imputation_background = imputation_background,
      svm = svm,
      cellprops = cellprops
    ))

  }, error = function(e) {
    stop(paste("Error loading background data:", e$message))
  })
}

#' Prepare data for imputation
#'
#' @param test_data Test data frame
#' @param imputation_background Background data for imputation
#' @param svm SVM models
#' @param n_closest Number of closest samples to select for imputation
#' @return List containing prepared test data and merged signatures
prepare_imputation_data <- function(test_data, imputation_background, svm, n_closest = 20) {
  merged_signatures <- lapply(svm, predictors)

  # Add CpGs for cell deconvolution
  cpgs=unique(unlist(merged_signatures))

  test_data_cpgs <- test_data[test_data$IlmnID %in% cpgs,]
  imputation_background <- imputation_background[imputation_background$IlmnID %in% cpgs,]

  # Select closest samples from background for imputation
  # cat("Selecting", n_closest, "closest samples from background for imputation...\n")

  # Get test sample IDs
  test_sample_ids <- setdiff(names(test_data_cpgs), "IlmnID")
  background_sample_ids <- setdiff(names(imputation_background), "IlmnID")

  cat("Total background samples available:", length(background_sample_ids), "\n")

  if(length(background_sample_ids) > n_closest && length(test_sample_ids) > 0) {
    # Calculate distances between test samples and all background samples using all non empty CpGs
    test_data_matrix <- test_data_cpgs[, test_sample_ids, drop = FALSE]
    rownames(test_data_matrix) <- test_data_cpgs$IlmnID

    background_matrix <- imputation_background[, background_sample_ids, drop = FALSE]
    rownames(background_matrix) <- imputation_background$IlmnID

    # Find common CpGs between test and background
    common_cpgs <- intersect(rownames(test_data_matrix), rownames(background_matrix))

    if(length(common_cpgs) > 0) {
      # cat("Using", length(common_cpgs), "CpGs for distance calculation\n")

      # Calculate pairwise distances between each test sample and each background sample
      # Then aggregate by taking the mean distance across test samples
      all_distances <- matrix(NA, nrow = length(test_sample_ids), ncol = length(background_sample_ids))
      rownames(all_distances) <- test_sample_ids
      colnames(all_distances) <- background_sample_ids

      for(i in seq_along(test_sample_ids)) {
        test_sample <- test_sample_ids[i]
        test_values <- test_data_matrix[common_cpgs, test_sample]

        for(j in seq_along(background_sample_ids)) {
          bg_sample <- background_sample_ids[j]
          bg_values <- background_matrix[common_cpgs, bg_sample]

          # Find positions where both test and background have non-missing values
          valid_positions <- !is.na(test_values) & !is.na(bg_values)

          if(sum(valid_positions) > 0) {
            # Calculate Euclidean distance using only non-missing CpGs
            all_distances[i, j] <- sqrt(sum((test_values[valid_positions] - bg_values[valid_positions])^2))
          }
        }
      }

      # Take mean distance across all test samples for each background sample
      mean_distances <- colMeans(all_distances, na.rm = TRUE)

      # cat("Calculated distances for", sum(!is.na(mean_distances)), "background samples\n")

      # Select closest n_closest samples
      closest_sample_ids <- names(sort(mean_distances))[1:min(n_closest, sum(!is.na(mean_distances)))]

      # cat("Selected", length(closest_sample_ids), "closest background samples\n")

      # Subset background to closest samples
      imputation_background <- imputation_background[, c("IlmnID", closest_sample_ids)]
    } else {
      cat("Warning: No common CpGs found, using first", n_closest, "background samples\n")
      imputation_background <- imputation_background[, c("IlmnID", head(background_sample_ids, n_closest))]
    }
  } else {
    # cat("Using all", length(background_sample_ids), "background samples (less than threshold)\n")
  }

  test_data <- merge(test_data_cpgs,
                     imputation_background,
                     by = "IlmnID",
                     all = TRUE)

  rownames(test_data) <- test_data$IlmnID
  test_data$IlmnID <- NULL

  return(list(
    test_data = test_data,
    merged_signatures = merged_signatures
  ))
}

#' Perform imputation on test data
#'
#' @param test_data Prepared test data
#' @param test_sample_ids IDs of test samples
#' @return Imputed data frame
perform_imputation <- function(test_data, test_sample_ids) {
  # Reload manifest
  manifest <- readRDS("../../data/support_files/manifest.qc_filtered.rds")
  manifest$MAPINFO <- NULL
  manifest <- as.data.frame(manifest)
  names(manifest) <- c('cpg', 'chr')

  test_data = test_data[rownames(test_data) %in% manifest$cpg,]
  
  # Impute missing data
  beta_SE_imputed <- methyLImp2(input = t(test_data),
                                type = "user",
                                annotation = manifest,
                                BPPARAM = SnowParam(exportglobals = FALSE,
                                                    workers = 1))
  df <- as.data.frame(t(beta_SE_imputed))
  
  # cat("Imputed dataset head")
  # cat(head(df))

  df$IlmnID <- rownames(df)

  df <- df[, names(df) %in% c("IlmnID", test_sample_ids)]

  return(df)
}

#' Prepare data for inference
#'
#' @param test_data Test data frame
#' @param svm SVM models
#' @return List of data frames for inference
prepare_inference_data <- function(test_data, svm) {
  lapply(svm, function(x) {
    z <- test_data[match(predictors(x), test_data$IlmnID), ,drop=F]
    z$IlmnID <- NULL
    return(z)
  })
}

#' Load beta values and signatures
#'
#' @return List containing insilico beta, meta, and signatures
load_beta_signatures <- function() {
  
  message("Loading controls and signatures")

  insilico_beta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.beta.rds")
  insilico_meta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.meta.rds")
  # Subset controls
  insilico_meta <- insilico_meta[insilico_meta$RealLabel == "control",]
  insilico_beta <- insilico_beta[, which(names(insilico_beta) %in% insilico_meta$geo_accession)]

  insilico_beta$IlmnID <- rownames(insilico_beta)
  signatures <- read.delim('../../data/support_files/merged_signatures_90DMRs.tsv', header = TRUE)
  signatures <- split(signatures, f = as.factor(paste(signatures$Label)))

  return(list(
    insilico_beta = insilico_beta,
    insilico_meta = insilico_meta,
    signatures = signatures
  ))
}

#' Load and prepare patient samples data
#'
#' @param signatures List of signatures
#' @return List containing patient samples beta and meta data
load_real_cases <- function(signatures) {
  
  message("Loading patient samples")

  tryCatch({
    real_cases_beta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.beta.rds")
    real_cases_meta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.meta.rds")
    real_cases_meta <- real_cases_meta[real_cases_meta$RealLabel != "control",]
    
    real_cases_beta <- real_cases_beta[, which(names(real_cases_beta) %in% real_cases_meta$geo_accession)]
    cat("Loaded patient samples data: beta matrix", dim(real_cases_beta), ", meta table", dim(real_cases_meta), "\n")

    # Check if data is empty
    if (nrow(real_cases_beta) == 0 || nrow(real_cases_meta) == 0) {
      warning("patient samples data is empty, creating empty placeholder")
      return(create_empty_real_cases(signatures))
    }

    # Filter metadata to samples that exist in beta matrix
    real_cases_meta <- real_cases_meta[real_cases_meta$geo_accession %in% colnames(real_cases_beta),]

    # Filter beta matrix to samples that exist in metadata
    real_cases_beta <- real_cases_beta[, colnames(real_cases_beta) %in% real_cases_meta$geo_accession, drop = FALSE]

    # cat("After filtering: beta matrix", dim(real_cases_beta), ", meta table", dim(real_cases_meta), "\n")

    # Add IlmnID column
    real_cases_beta$IlmnID <- rownames(real_cases_beta)

    # Rename columns and include Sex and AgeGroup for heatmap annotations
    # Keep all metadata columns, just rename the key ones
    real_cases_meta_processed <- data.frame(
      IDs = real_cases_meta$geo_accession,
      Status = real_cases_meta$RealLabel,
      Platform = real_cases_meta$platform_id,
      Source = "literature",
      stringsAsFactors = FALSE
    )

    # Add Sex and AgeGroup if they exist in the metadata
    if("Sex" %in% names(real_cases_meta)) {
      real_cases_meta_processed$Sex <- real_cases_meta$Sex
    }
    if("AgeGroup" %in% names(real_cases_meta)) {
      real_cases_meta_processed$AgeGroup <- real_cases_meta$AgeGroup
    }

    real_cases_meta <- real_cases_meta_processed

    # Process beta data for each signature
    real_cases_beta <- lapply(signatures, function(x) {
      y <- real_cases_beta[real_cases_beta$IlmnID %in% x$ProbeID, , drop = FALSE]
      return(as.data.frame(y))
    })

    return(list(
      real_cases_beta = real_cases_beta,
      real_cases_meta = real_cases_meta
    ))

  }, error = function(e) {
    warning(paste("Error loading patient samples data:", e$message, "- using empty placeholder"))
    return(create_empty_real_cases(signatures))
  })
}

#' Create empty patient samples data structure when real data is unavailable
#'
#' @param signatures List of signatures
#' @return Empty patient samples data structure
create_empty_real_cases <- function(signatures) {
  # Create empty beta data for each signature
  empty_beta <- lapply(signatures, function(x) {
    empty_df <- data.frame(IlmnID = character(0))
    return(empty_df)
  })

  # Create empty metadata
  empty_meta <- data.frame(
    IDs = character(0),
    Status = character(0),
    Platform = character(0),
    Source = character(0),
    stringsAsFactors = FALSE
  )

  return(list(
    real_cases_beta = empty_beta,
    real_cases_meta = empty_meta
  ))
}