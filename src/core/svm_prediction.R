#' Make predictions using SVM models
#'
#' @param imputed_data Prepared inference data
#' @param svm SVM models
#' @param test_data_ids IDs of test samples
#' @return Data frame of prediction results
make_predictions <- function(imputed_data, svm, test_data_ids) {
  
  results <- list()
  
  message("Processing samples for DNAm testing...")
  for (i in names(imputed_data)) {
    
    for (j in names(svm)) {
      if (grepl(i, j)) {
        
        tryCatch({
          
          # Get prediction for each model matching by name
          results[[paste(i, j)]] <- predict(svm[[j]],
                                            newdata = as.data.frame(t(imputed_data[[i]])),
                                            type = "prob")
          
          results[[paste(i, j)]]$SampleID <- names(imputed_data[[i]])
        }, error=function(err){
          print(paste("Failure processing", err))
        })
      }
    }
  }
  
  process_results(results, test_data_ids)
}

#' Process prediction results
#'
#' @param results Raw prediction results
#' @param test_data_ids IDs of test samples
#' @return Processed prediction results
process_results <- function(results, test_data_ids) {
  results <- lapply(results, as.data.frame)
  results <- bind_rows(results, .id = "Model")
  
  results <- pivot_longer(results,
                          -c("SampleID", "Model"),
                          names_to = "Signature",
                          values_to = "SVM_score")
  
  results <- results[results$Signature != "control",]
  results <- results[results$SampleID %in% test_data_ids,]
  results$SVM <- paste(str_split(results$Model, "_", simplify = TRUE)[,3], 
                       str_split(results$Model, "_", simplify = TRUE)[,4])
  
  # Calculate prediction scores after excluding the highest and lowest scores.
  # `n_svm` is captured *before* the trim filter so the QC column reflects how
  # many classifiers were actually loaded for this signature, not how many
  # contributed to the mean (post-trim is always n-2 when n >= 3).
  results %>%
    group_by(SampleID, SVM) %>%
    mutate(n_svm = dplyr::n(),
           rank  = rank(SVM_score, ties.method = "first")) %>%
    filter(n_svm < 3 | (rank != min(rank) & rank != max(rank))) %>%
    summarise(pSVM_average = mean(SVM_score),
              pSVM_sd      = sd(SVM_score),
              n_svm        = dplyr::first(n_svm),
              .groups      = "drop")
}

#' Normalised signature key for joining SVM and NNET results.
#' SVMs name signatures like "Kabuki KMT2D"; NNET checkpoints store them as
#' "Kabuki_KMT2D" (or sometimes lowercase). Strip separators and case so both
#' sides land on the same key.
normalize_signature_key <- function(x) {
  tolower(gsub("[_[:space:]]+", "", as.character(x)))
}

#' Combine SVM + NNET per-(sample, signature) summaries into the metapredictor
#' table used by plots, tables and exports.
#'
#' @param svm_results  Output of make_predictions() (SampleID, SVM, pSVM_average, pSVM_sd)
#' @param nnet_results Output of make_nnet_predictions() (SampleID, Signature, pNNET_average, pNNET_sd)
#' @return Tibble with SVM column (display label) plus pSVM_average/sd,
#'         pNNET_average/sd, mean_case, whisker_low/high, abs_diff.
combine_svm_nnet_results <- function(svm_results, nnet_results, na_pct = NULL) {
  svm_results$.key  <- normalize_signature_key(svm_results$SVM)
  nnet_results$.key <- normalize_signature_key(nnet_results$Signature)
  nnet_results$Signature <- NULL

  merged <- dplyr::full_join(svm_results, nnet_results,
                             by = c("SampleID", ".key"))

  # If NNET was disabled, the join still needs the columns to exist downstream.
  if (!"pNNET_average" %in% names(merged)) merged$pNNET_average <- NA_real_
  if (!"pNNET_sd"      %in% names(merged)) merged$pNNET_sd      <- NA_real_
  if (!"n_nnet"        %in% names(merged)) merged$n_nnet        <- 0L
  if (!"n_svm"         %in% names(merged)) merged$n_svm         <- 0L
  # NA counts (from a full_join row that only had NNET, no SVM, or vice versa)
  # are zero loaded classifiers on that side.
  merged$n_svm[is.na(merged$n_svm)]   <- 0L
  merged$n_nnet[is.na(merged$n_nnet)] <- 0L

  # If a NNET-only signature comes through (no SVM label), fall back to its key.
  merged$SVM[is.na(merged$SVM)] <- merged$.key[is.na(merged$SVM)]
  merged$.key <- NULL

  merged$mean_case    <- rowMeans(merged[, c("pSVM_average", "pNNET_average")], na.rm = TRUE)
  merged$whisker_low  <- pmin(merged$pSVM_average, merged$pNNET_average, na.rm = TRUE)
  merged$whisker_high <- pmax(merged$pSVM_average, merged$pNNET_average, na.rm = TRUE)
  merged$abs_diff     <- abs(merged$pSVM_average - merged$pNNET_average)

  # NaN from rowMeans of all-NA rows -> NA
  merged$mean_case[is.nan(merged$mean_case)] <- NA_real_

  # Optional: attach pre-imputation NA% per (sample x signature). `na_pct`
  # is keyed by signature label (e.g. "ADNPc_ArefEshghi2020"); results$SVM
  # uses a space-separated display form ("ADNPc ArefEshghi2020"). Normalise
  # both sides with the same key as the SVM/NNET merge above.
  if (!is.null(na_pct) && nrow(na_pct) > 0) {
    na_pct$.key <- normalize_signature_key(na_pct$.sig_key)
    na_pct$.sig_key <- NULL
    merged$.key <- normalize_signature_key(merged$SVM)
    merged <- dplyr::left_join(merged, na_pct, by = c("SampleID", ".key"))
    merged$.key <- NULL
  } else {
    merged$pct_na_pre <- NA_real_
  }

  tibble::as_tibble(merged)
}

#' Get signatures with pSVM >= threshold for dimension plot filtering
#'
#' @param results Prediction results data frame
#' @param probands Selected probands
#' @param signatures Selected signatures
#' @param threshold Minimum pSVM threshold
#' @return Vector of signature names with pSVM >= threshold
get_high_scoring_signatures_with_threshold <- function(results, probands, signatures, threshold = 0.05) {
  # Filter results for selected probands and signatures
  filtered_results <- results[results$SampleID %in% probands &
                                results$SVM %in% gsub("_", " ", signatures), ]

  # Score column used for thresholding is the metapredictor mean (SVM+NNET);
  # fall back to pSVM_average if mean_case isn't present (NNET-disabled run).
  score <- if ("mean_case" %in% names(filtered_results))
             filtered_results$mean_case
           else filtered_results$pSVM_average
  high_scoring_results <- filtered_results[!is.na(score) & score >= threshold, ]
  high_scoring_signatures <- unique(high_scoring_results$SVM)
  
  # Convert back to signature format 
  high_scoring_signatures <- gsub(" ", "_", high_scoring_signatures)
  
  # Return only signatures that are in the original selection
  high_scoring_signatures[high_scoring_signatures %in% signatures]
}
