# Path to the signature-definition file shipped with the app. Single source of
# truth so the signature *version* (a content hash, below) and the loaders agree.
SIGNATURE_FILE <- "../../data/support_files/merged_signatures_90DMRs.tsv"

#' Signature version: a short content hash (md5) of the signature file shipped
#' with the app. Recorded in the HTML report and the output tables so results
#' can be traced to the exact signature set / classifier version used.
#'
#' @param path Path to the signature file (defaults to SIGNATURE_FILE)
#' @param n    Number of leading hex chars to keep (default 10)
#' @return Character scalar (e.g. "634d32f155"), or NA if the file is missing
signature_version <- function(path = SIGNATURE_FILE, n = 10) {
  if (!file.exists(path)) return(NA_character_)
  substr(unname(tools::md5sum(path)), 1, n)
}

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

#' Load background data, SVM models, and NNET checkpoints
#'
#' The model directory is expected to contain two subfolders:
#'   <model_dir>/SVM/   -> caret SVM models as .rds
#'   <model_dir>/NNET/  -> PyTorch FlexibleNNet checkpoints as .pth
#' For backward compatibility, if no SVM/ subfolder is present but .rds files
#' live directly under <model_dir>, those are used as the SVM models and NNET
#' inference is disabled.
#'
#' @param model_dir Directory containing SVM/ and NNET/ subfolders
#' @return List with imputation background, svm models, nnet checkpoint paths, cellprops
load_background_data <- function(model_dir) {
  tryCatch({
    imputation_background <- readRDS("../../data/imputationsamples/samples.beta.rds")

    svm_dir  <- file.path(model_dir, "SVM")
    nnet_dir <- file.path(model_dir, "NNET")

    if (!dir.exists(svm_dir)) {
      # Backward-compatible fallback: flat .rds layout
      warning("No 'SVM' subfolder found in ", model_dir,
              " - falling back to flat .rds layout (NNET inference disabled)")
      svm_dir <- model_dir
      nnet_dir <- NULL
    }

    svm_files <- list.files(svm_dir, pattern = "\\.rds$", full.names = TRUE)
    if (length(svm_files) == 0) {
      stop("No .rds SVM model files found in ", svm_dir)
    }

    svm <- lapply(svm_files, function(p) {
      tryCatch(butcher(readRDS(p)),
               error = function(e) {
                 warning(paste("Error loading SVM:", basename(p), "-", e$message))
                 NULL
               })
    })
    svm <- Filter(Negate(is.null), svm)
    if (length(svm) == 0) stop("No valid SVM models could be loaded from ", svm_dir)
    names(svm) <- gsub("\\.rds$", "", basename(svm_files[seq_along(svm)]))
    cat("Successfully loaded", length(svm), "SVM models\n")

    nnet_files <- character(0)
    if (!is.null(nnet_dir) && dir.exists(nnet_dir)) {
      nnet_files <- list.files(nnet_dir, pattern = "\\.pth$", full.names = TRUE)
      cat("Found", length(nnet_files), "NNET checkpoints\n")
    } else if (!is.null(nnet_dir)) {
      warning("NNET folder '", nnet_dir, "' does not exist - NNET inference disabled")
    }

    cellprops <- readRDS("../../data/support_files/background_training.cellprops.rds")

    list(
      imputation_background = imputation_background,
      svm = svm,
      nnet_files = nnet_files,
      cellprops = cellprops
    )

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
prepare_imputation_data <- function(test_data, imputation_background, svm, n_closest = 20,
                                    extra_cpgs = character(0)) {
  merged_signatures <- lapply(svm, predictors)

  # Union of SVM predictors and any extra CpGs required by NNET checkpoints
  cpgs <- unique(c(unlist(merged_signatures), extra_cpgs))

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

#' Compute the percentage of NA probes per (sample x signature) BEFORE
#' imputation, so it can be surfaced as a QC column in the prediction table.
#'
#' @param test_data_pre  Pre-imputation beta matrix (probes x samples, rownames = IlmnID)
#' @param signatures     Named list keyed by signature label, each element a
#'                       data.frame with a `ProbeID` column (output of
#'                       `load_beta_signatures()$signatures`).
#' @param sample_ids     Test sample IDs (columns of `test_data_pre`).
#' @return data.frame with columns SampleID, .sig_key, pct_na_pre.
compute_pre_imputation_na_pct <- function(test_data_pre, signatures, sample_ids) {
  if ("IlmnID" %in% names(test_data_pre)) {
    rownames(test_data_pre) <- test_data_pre$IlmnID
    test_data_pre$IlmnID <- NULL
  }
  sample_ids <- intersect(sample_ids, names(test_data_pre))
  if (length(sample_ids) == 0) {
    return(data.frame(SampleID = character(0),
                      .sig_key = character(0),
                      pct_na_pre = numeric(0)))
  }

  rows <- list()
  for (sig_name in names(signatures)) {
    probes <- unique(signatures[[sig_name]]$ProbeID)
    n_total <- length(probes)
    if (n_total == 0) next
    present <- intersect(probes, rownames(test_data_pre))
    for (s in sample_ids) {
      n_na_present <- if (length(present) > 0)
                        sum(is.na(test_data_pre[present, s]))
                      else 0L
      n_missing_probes <- n_total - length(present)   # probes not on the array count as NA
      pct <- 100 * (n_na_present + n_missing_probes) / n_total
      rows[[length(rows) + 1]] <- data.frame(
        SampleID   = s,
        .sig_key   = sig_name,
        pct_na_pre = round(pct, 1),
        stringsAsFactors = FALSE
      )
    }
  }
  if (length(rows) == 0) {
    return(data.frame(SampleID = character(0),
                      .sig_key = character(0),
                      pct_na_pre = numeric(0)))
  }
  do.call(rbind, rows)
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

#' Map raw GPL platform IDs to the friendly labels expected by the heatmap
#' annotation color palette (EpicV1 / EpicV2 / 450k). Values that don't match
#' a known GPL ID are returned unchanged.
gpl_to_friendly_platform <- function(x) {
  m <- c(
    "GPL13534"  = "450k",
    "GPL16304"  = "450k",
    "GPL21145"  = "EpicV1",
    "GPL23976"  = "EpicV1",
    "GPL33022"  = "EpicV2",
    "EPIC"      = "EpicV1",
    "EPICv2"    = "EpicV2",
    "EpicV1"    = "EpicV1",
    "EpicV2"    = "EpicV2",
    "450k"      = "450k"
  )
  out <- m[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}

#' Load beta values and signatures
#'
#' @return List containing insilico beta, meta, and signatures
load_beta_signatures <- function() {
  
  message("Loading controls and signatures")

  insilico_beta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.beta.rds")
  insilico_meta <- readRDS("../../data/affectedindividuals/affectedindividuals_methadory.meta.rds")
  # Subset controls. Files generated by 12_methadory_app_data.R store CpG IDs
  # in an `IlmnID` column (rownames are integer indices), so preserve that
  # column when subsetting to control samples; the older code assumed
  # CpG IDs lived in rownames(insilico_beta) and silently produced "1","2"...
  # IlmnIDs for the new files.
  insilico_meta <- insilico_meta[insilico_meta$RealLabel == "control",]
  keep_cols <- c("IlmnID",
                 intersect(names(insilico_beta), insilico_meta$geo_accession))
  insilico_beta <- insilico_beta[, keep_cols]

  # If the file stored CpG IDs in rownames instead (legacy layout), promote
  # them; otherwise the IlmnID column is already correct.
  if (all(insilico_beta$IlmnID == seq_len(nrow(insilico_beta)))) {
    insilico_beta$IlmnID <- rownames(insilico_beta)
  }

  # Normalise platform_id to the friendly labels used by the heatmap palette.
  if ("platform_id" %in% names(insilico_meta)) {
    insilico_meta$platform_id <- gpl_to_friendly_platform(insilico_meta$platform_id)
  }
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
    
    # Preserve the IlmnID column (CpG IDs stored as a column, not as rownames,
    # by 12_methadory_app_data.R) when subsetting samples.
    keep_cols <- c("IlmnID",
                   intersect(names(real_cases_beta), real_cases_meta$geo_accession))
    real_cases_beta <- real_cases_beta[, keep_cols]
    cat("Loaded patient samples data: beta matrix", dim(real_cases_beta), ", meta table", dim(real_cases_meta), "\n")

    # Check if data is empty
    if (nrow(real_cases_beta) == 0 || nrow(real_cases_meta) == 0) {
      warning("patient samples data is empty, creating empty placeholder")
      return(create_empty_real_cases(signatures))
    }

    # Filter metadata to samples that exist in beta matrix
    real_cases_meta <- real_cases_meta[real_cases_meta$geo_accession %in% colnames(real_cases_beta),]

    # Filter beta matrix to samples that exist in metadata (keep IlmnID col)
    real_cases_beta <- real_cases_beta[, c("IlmnID",
                                           intersect(colnames(real_cases_beta), real_cases_meta$geo_accession)),
                                       drop = FALSE]

    # Backfill IlmnID from rownames only if it's a placeholder integer index
    # (legacy file layout where CpG IDs lived in rownames).
    if (!is.null(real_cases_beta$IlmnID) &&
        all(real_cases_beta$IlmnID == seq_len(nrow(real_cases_beta)))) {
      real_cases_beta$IlmnID <- rownames(real_cases_beta)
    }

    # Rename columns and include Sex and AgeGroup for heatmap annotations
    # Keep all metadata columns, just rename the key ones
    real_cases_meta_processed <- data.frame(
      IDs = real_cases_meta$geo_accession,
      Status = real_cases_meta$RealLabel,
      Platform = gpl_to_friendly_platform(real_cases_meta$platform_id),
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

    # Process beta data for each signature. Mirror the controls_plot_beta
    # layout so the downstream full_join in prepare_plot_data lines up on
    # CpG IDs: store CpG IDs in rownames and DROP the IlmnID column so the
    # subsequent `x$IlmnID = rownames(x)` re-injects the correct CpG IDs.
    # (Previously the IlmnID column survived with the original ProbeIDs but
    # was overwritten by integer row indices, causing zero CpG overlap
    # between real cases and other cohorts and silently dropping all real
    # cases from the dimension-reduction plots.)
    real_cases_beta <- lapply(signatures, function(x) {
      y <- real_cases_beta[real_cases_beta$IlmnID %in% x$ProbeID, , drop = FALSE]
      rownames(y) <- y$IlmnID
      y$IlmnID    <- NULL
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