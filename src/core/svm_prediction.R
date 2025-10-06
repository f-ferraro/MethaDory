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
  
  # Calculate prediction scores after excluding the highest and lowest scores
  results %>%
    group_by(SampleID, SVM) %>%
    mutate(rank = rank(SVM_score, ties.method = "first")) %>%
    filter(rank != min(rank) & rank != max(rank)) %>%
    summarise(pSVM_average = mean(SVM_score),
              pSVM_sd = sd(SVM_score),
              .groups = "drop")
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
  
  # Get signatures with pSVM_average >= threshold
  high_scoring_results <- filtered_results[filtered_results$pSVM_average >= threshold, ]
  high_scoring_signatures <- unique(high_scoring_results$SVM)
  
  # Convert back to signature format 
  high_scoring_signatures <- gsub(" ", "_", high_scoring_signatures)
  
  # Return only signatures that are in the original selection
  high_scoring_signatures[high_scoring_signatures %in% signatures]
}
