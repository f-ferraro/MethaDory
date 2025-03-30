# Helper functions for MethaDory app
# Load required libraries
library(BiocParallel)
library(caret)
paste0(getwd(),"/ChAMPdata_2.38.0.tar.gz")
install.packages(paste0(getwd(),"/ChAMPdata_2.38.0.tar.gz"),
                 repos = NULL, type="source")

library(circlize)
library(ComplexHeatmap)
library(data.table)
library(DT)
library(EpiDISH)
library(fs)
library(ggplotify)
library(methyLImp2)
library(patchwork)
library(PCAtools)
library(plotly)
library(Rtsne)
library(shiny)
library(shinybusy)
library(shinydashboard)
library(shinyFiles)
library(tidyverse)
library(umap)
library(scales)
#' Load and prepare test data
#'
#' @param file_path Path to the test data file
#' @return List containing the test data and test data IDs
load_test_data <- function(file_path) {
  test_data_user <- read.table(file_path, header = TRUE)
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
  # Load imputation background
  imputation_background <- lapply(list.files("data/imputation", "rds", full.names = TRUE), readRDS)
  imputation_background <- bind_rows(imputation_background)
  imputation_background$IlmnID <- rownames(imputation_background)
  
  # Load models from the selected directory
  models <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(models) == 0) {
    stop("No .rds model files found in the selected directory")
  }
  
  # Load SVM models
  svm <- lapply(models, function(model_path) {
    tryCatch({
      readRDS(model_path)
    }, error = function(e) {
      warning(paste("Error loading model:", model_path))
      NULL
    })
  })
  
  # Remove any NULL entries from failed loads
  svm <- Filter(Negate(is.null), svm)
  
  if (length(svm) == 0) {
    stop("No valid SVM models could be loaded from the selected directory")
  }
  
  # Name the models based on their filenames
  names(svm) <- gsub("\\.rds$", "", basename(models))
  
  
  # Read cell props from the background 
  cellprops = readRDS("data/background_training.cellprops.rds")
  
  
  return(list(
    imputation_background = imputation_background,
    svm = svm,
    cellprops = cellprops
  ))
}

#' Prepare data for imputation
#'
#' @param test_data Test data frame
#' @param imputation_background Background data for imputation
#' @param svm SVM models
#' @return List containing prepared test data and merged signatures
prepare_imputation_data <- function(test_data, imputation_background, svm) {
  merged_signatures <- lapply(svm, predictors)
  test_data <- test_data[test_data$IlmnID %in% unlist(merged_signatures),]
  imputation_background <- imputation_background[imputation_background$IlmnID %in% unlist(merged_signatures),]
  
  test_data <- merge(test_data,
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
  manifest <- readRDS("data/manifest.qc_filtered.rds")
  manifest$MAPINFO <- NULL
  names(manifest) <- c('cpg', 'chr')
  
  beta_SE_imputed <- methyLImp2(input = t(test_data),
                                type = "user",
                                annotation = manifest,
                                BPPARAM = SnowParam(exportglobals = FALSE,
                                                    workers = 1))
  df <- as.data.frame(t(beta_SE_imputed))
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
    z <- test_data[match(predictors(x), test_data$IlmnID), ]
    z$IlmnID <- NULL
    return(z)
  })
}

#' Make predictions using SVM models
#'
#' @param test_for_inference Prepared inference data
#' @param svm SVM models
#' @param test_data_ids IDs of test samples
#' @return Data frame of prediction results
make_predictions <- function(test_for_inference, svm, test_data_ids) {
  results <- list()
  
  for (i in names(test_for_inference)) {
    for (j in names(svm)) {
      if (grepl(i, j)) {
        results[[paste(i, j)]] <- predict(svm[[j]],
                                          newdata = t(test_for_inference[[i]]),
                                          type = "prob")
        results[[paste(i, j)]]$SampleID <- names(test_for_inference[[i]])
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
  
  results %>%
    group_by(SampleID, SVM) %>%
    mutate(rank = rank(SVM_score, ties.method = "first")) %>%
    filter(rank != min(rank) & rank != max(rank)) %>%
    summarise(pSVM_average = mean(SVM_score),
              pSVM_sd = sd(SVM_score),
              .groups = "drop")
}

#' Load beta values and signatures
#'
#' @return List containing insilico beta, meta, and signatures
load_beta_signatures <- function() {
  insilico_beta <- readRDS("data/syntheticcases/samples_for_insilico.beta.rds")
  insilico_meta <- readRDS("data/syntheticcases/samples_for_insilico.meta.rds")
  signatures <- read.table('data/merged_signatures_90DMRs.tsv', header = TRUE)
  signatures <- split(signatures, f = as.factor(paste(signatures$Label)))
  
  return(list(
    insilico_beta = insilico_beta,
    insilico_meta = insilico_meta,
    signatures = signatures
  ))
}

#' Load and prepare real cases data
#'
#' @param signatures List of signatures
#' @return List containing real cases beta and meta data
load_real_cases <- function(signatures) {
  real_cases_beta <- readRDS("data/affectedindividuals/affectedindividuals_methadory.beta.rds")
  real_cases_meta <- readRDS("data/affectedindividuals/affectedindividuals_methadory.meta.rds")
  real_cases_meta <- real_cases_meta[, c("geo_accession", "RealLabel")]
  
  real_cases_meta <- real_cases_meta[real_cases_meta$geo_accession %in% colnames(real_cases_beta),]
  real_cases_beta <- real_cases_beta[, colnames(real_cases_beta) %in% real_cases_meta$geo_accession]
  real_cases_beta$IlmnID <- rownames(real_cases_beta)
  names(real_cases_meta) <- c('IDs', "Status")
  real_cases_meta$Source <- "literature"
  
  real_cases_beta <- lapply(signatures, function(x) {
    y <- real_cases_beta[real_cases_beta$IlmnID %in% x$ProbeID,]
    return(as.data.frame(y))
  })
  
  return(list(
    real_cases_beta = real_cases_beta,
    real_cases_meta = real_cases_meta
  ))
}

#' Prepare in-silico cases data
#'
#' @param signatures List of signatures
#' @param insilico_beta In-silico beta values
#' @return List of prepared in-silico cases data
prepare_insilico_cases <- function(signatures, insilico_beta) {
  lapply(signatures, function(x) {
    # Subset beta values for relevant probes
    y <- insilico_beta[insilico_beta$IlmnID %in% x$ProbeID, ]
    rownames(y) <- y$IlmnID
    y$IlmnID <- NULL
    
    # Match the order of probes to the signature
    y <- y[match(x$ProbeID, rownames(y)), ]
    
    # Apply delta beta and clamp values between 0 and 1
    y <- y + x$deltaBeta
    y[y < 0] <- 0
    y[y > 1] <- 1
    
    # Add suffix to column names to identify as in-silico cases
    names(y) <- paste0(names(y), '_isc')
    
    # Add back IlmnID column
    y$IlmnID <- rownames(y)
    
    return(as.data.frame(y))
  })
}

#' Prepare user data for plots
#'
#' @param signatures List of signatures
#' @param test_data_user User test data
#' @return List of prepared user data for plots
prepare_user_plots <- function(signatures, test_data_user) {
  lapply(signatures, function(x) {
    # Subset user data for relevant probes
    y <- test_data_user[test_data_user$IlmnID %in% x$ProbeID, ]
    rownames(y) <- y$IlmnID
    y <- y[match(x$ProbeID, rownames(y)), ]
    
    return(as.data.frame(y))
  })
}

#' Prepare data for plots
#'
#' @param signatures List of signatures
#' @param insilico_beta In-silico beta values
#' @param test_data_user_imputed Imputed user test data
#' @param real_cases_beta Real cases beta values
#' @return List of prepared data for plots
prepare_plot_data <- function(signatures, insilico_beta, test_data_user_imputed, real_cases_beta) {
  insilicocases_plot_beta <- prepare_insilico_cases(signatures, insilico_beta)
  user_plot_beta <- prepare_user_plots(signatures, test_data_user_imputed)
  
  data_for_plots_beta <- lapply(names(signatures), function(s) {
    lst <- list(insilicocases_plot_beta[[s]],
                user_plot_beta[[s]],
                real_cases_beta[[s]]) %>%
      purrr::reduce(full_join, by = "IlmnID")
    
    lst$IlmnID <- NULL
    return(lst)
  })
  
  names(data_for_plots_beta) <- names(signatures)
  return(data_for_plots_beta)
}

#' Prepare metadata for plots
#'
#' @param signatures List of signatures
#' @param insilicocases_plot_beta In-silico cases plot beta values
#' @param user_plot_beta User plot beta values
#' @param real_cases_meta Real cases metadata
#' @return List of prepared metadata for plots
prepare_plot_metadata <- function(signatures, insilicocases_plot_beta, user_plot_beta, real_cases_meta) {
  data_for_plots_meta <- lapply(names(signatures), function(s) {
    mt_ids_iscases <- setdiff(colnames(insilicocases_plot_beta[[s]]), "IlmnID")
    mt_ids_usertests <- setdiff(colnames(user_plot_beta), "IlmnID")

    lsd <- data.frame(
      'IDs' = c(mt_ids_iscases, mt_ids_usertests),
      'Status' = rep("in_silico_case", length(c(mt_ids_iscases, mt_ids_usertests))),
      'Source' = rep("literature", length(c(mt_ids_iscases, mt_ids_usertests)))
    )

    lsd <- lsd[!lsd$IDs %in% real_cases_meta$IDs, ]
    lsd <- rbind(lsd, real_cases_meta)
    lsd$Status <- ifelse(lsd$IDs %in% mt_ids_usertests, "proband", lsd$Status)
    lsd$Source <- ifelse(lsd$IDs %in% mt_ids_usertests, "user_sample", lsd$Source)

    lsd <- lsd[!duplicated(lsd),]
    rownames(lsd) <- lsd$IDs
    return(lsd)
  })

  names(data_for_plots_meta) <- names(signatures)
  return(data_for_plots_meta)
}

#' Create prediction plot
#'
#' @param results Prediction results
#' @param proband Proband ID
#' @return ggplot object of prediction plot
create_prediction_plot <- function(results, proband, signatures, plot_width) {
  
  plot_data = results[results$SampleID %in% proband & 
                        results$SVM %in% gsub("_", " ",signatures),]
  plot_data$pSVMmin = ifelse(plot_data$pSVM_average - plot_data$pSVM_sd < 0, 0, plot_data$pSVM_average - plot_data$pSVM_sd)
  plot_data$pSVMmax = ifelse(plot_data$pSVM_average + plot_data$pSVM_sd > 1, 1, plot_data$pSVM_average + plot_data$pSVM_sd)
  
  p = ggplot(plot_data,
         aes(SVM, pSVM_average, color = SampleID)) +
    geom_hline(yintercept = c(0, 0.5, 1), color = "darkgrey") +
    geom_pointrange(aes(ymin = pSVMmin,
                        ymax = pSVMmax),
                    position = position_jitter(width = 0.20, height = 0)) +
    theme_minimal() +
    scale_x_discrete(labels = scales::label_wrap(20)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none") +
    ylab("SVM score") +
    xlab("")

  p
}



#' Create prediction plot for pdf
#'
#' @param results Prediction results
#' @param proband Proband ID
#' @return ggplot object of prediction plot
create_prediction_plot2 <- function(results, proband, signatures, plot_width) {
  
  plot_data = results[results$SampleID %in% proband & 
                        results$SVM %in% gsub("_", " ",signatures),]
  plot_data$pSVMmin = ifelse(plot_data$pSVM_average - plot_data$pSVM_sd < 0, 0, plot_data$pSVM_average - plot_data$pSVM_sd)
  plot_data$pSVMmax = ifelse(plot_data$pSVM_average + plot_data$pSVM_sd > 1, 1, plot_data$pSVM_average + plot_data$pSVM_sd)
  
  ggplot(plot_data,
             aes(SVM, pSVM_average, color = SampleID)) +
    geom_hline(yintercept = c(0, 0.5, 1), color = "darkgrey") +
    geom_pointrange(aes(ymin = pSVMmin,
                        ymax = pSVMmax),
                    position = position_jitter(width = 0.20, height = 0)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    ylab("SVM score") +
    xlab("Predictive Model")

}


#' Create color scheme for plots
#'
#' @param test_id Test ID
#' @return List of annotation colors
create_annotation_colors <- function(test_id) {
  ann_colors <- list(
    Status = c('#FEB24C',
               "#FCDE9C",
               '#089099',
               "#DC3977"), 
    Source = c("#A599CA",
               "#DC3977")
  )
  
  names(ann_colors$Status) <- c(ifelse(length(test_id) == 0, 'NA', test_id),
                                "in_silico_case",
                                'control',
                                'proband')
  
  names(ann_colors$Source) <- c("literature",
                                'user_sample')
  
  return(ann_colors)
}

#' Create PCA plot
#'
#' @param pca_object PCA object
#' @return ggplot object of PCA plot
create_pca_plot <- function(pca_object) {
  # Get test_id for color scheme
  test_id <- setdiff(unique(pca_object$metadata$Status), 
                     c("in_silico_case", "control", "proband"))
  
  # Create color scheme
  ann_colors <- create_annotation_colors(test_id)
  
  # Create data frame for plotting
  pca_plot_data <- data.frame(
    PC1 = pca_object$rotated$PC1,
    PC2 = pca_object$rotated$PC2,
    Status = pca_object$metadata$Status
  )
  
  # Calculate variance explained percentages
  variance_pct <- round(pca_object$variance / sum(pca_object$variance) * 100)
  
  # Create the plot
  ggplot() +
    geom_point(pca_plot_data[pca_plot_data$Status != "proband",], 
               mapping=aes(PC1, PC2, 
                            color = Status,
                            shape = (Status == "proband"), 
                            size = 1,
                            alpha = 0.85)) +
    geom_point(pca_plot_data[pca_plot_data$Status == "proband",], 
               mapping=aes(PC1, PC2, 
                            color = Status,
                            shape = (Status == "proband"), 
                            size = 1,
                            alpha = 0.85)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlab(paste0("PC1 (", variance_pct[1], "%)")) +
    ylab(paste0("PC2 (", variance_pct[2], "%)")) +
    ggtitle("PCA") +
    scale_colour_manual(values = ann_colors$Status) +
    guides(size = "none", shape = "none", alpha = "none")
}

#' Create t-SNE plot
#'
#' @param data_beta Beta values
#' @param data_meta Metadata
#' @param perplexity Perplexity parameter for t-SNE
#' @param max_iter Maximum iterations for t-SNE
#' @return ggplot object of t-SNE plot
create_tsne_plot <- function(data_beta, data_meta, perplexity = 4, max_iter = 1500) {
  # Get test_id for color scheme
  test_id <- setdiff(unique(data_meta$Status), 
                     c("in_silico_case", "control", "proband"))
  
  # Create color scheme
  ann_colors <- create_annotation_colors(test_id)
  
  # Run t-SNE
  tsne_results <- Rtsne(t(data_beta),
                        dims = 2,
                        perplexity = perplexity,
                        verbose = TRUE,
                        max_iter = max_iter)
  
  # Create data frame for plotting
  tsne_df <- data.frame(
    X = tsne_results$Y[, 1],
    Y = tsne_results$Y[, 2],
    Status = data_meta$Status
  )
  
  # Create the plot
  ggplot(tsne_df) +
    geom_point(tsne_df[tsne_df$Status != "proband",], 
               mapping=aes(X, Y,
                           color = Status,
                           shape = Status == "proband", 
                           size = 1,
                           alpha = 0.85)) +
    geom_point(tsne_df[tsne_df$Status == "proband",], 
               mapping=aes(X, Y,
                           color = Status,
                           shape = Status == "proband", 
                           size = 1,
                           alpha = 0.85)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    ggtitle("t-SNE") +
    scale_colour_manual(values = ann_colors$Status) +
    guides(size = "none", shape = "none", alpha = "none")
}

#' Create UMAP plot
#'
#' @param data_beta Beta values
#' @param data_meta Metadata
#' @return ggplot object of UMAP plot
create_umap_plot <- function(data_beta, data_meta) {
  # Get test_id for color scheme
  test_id <- setdiff(unique(data_meta$Status), 
                     c("in_silico_case", "control", "proband"))
  
  # Create color scheme
  ann_colors <- create_annotation_colors(test_id)
  
  # Run UMAP
  umap_results <- umap(t(data_beta))
  
  # Create data frame for plotting
  umap_df <- data.frame(
    IDs = rownames(umap_results$layout),
    UMAP1 = umap_results$layout[,1],
    UMAP2 = umap_results$layout[,2]
  )
  
  # Merge with metadata
  umap_df <- merge(umap_df,
                   data_meta,
                   by = "IDs")
  
  # Create the plot
  ggplot(umap_df) +
    geom_point(umap_df[umap_df$Status != "proband",], 
               mapping=aes(UMAP1, UMAP2,
                           color = Status,
                           shape = Status == "proband", 
                           size = 1,
                           alpha = 0.85)) +
    geom_point(umap_df[umap_df$Status == "proband",], 
               mapping=aes(UMAP1, UMAP2,
                           color = Status,
                           shape = Status == "proband", 
                           size = 1,
                           alpha = 0.85)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("UMAP") +
    scale_colour_manual(values = ann_colors$Status) +
    guides(size = "none", shape = "none", alpha = "none")
}

#' Create heatmap
#'
#' @param data_beta Beta values
#' @param data_meta Metadata
#' @return Heatmap grob
create_heatmap <- function(data_beta, data_meta) {
  # Get test_id for color scheme
  test_id <- setdiff(unique(data_meta$Status), 
                     c("in_silico_case", "control", "proband"))
  
  # Create color scheme
  ann_colors <- create_annotation_colors(test_id)
  
  # Normalize beta values
  norm_beta <- as.data.frame(t(apply(data_beta, 1, function(x) {
    (x - mean(x)) / sd(x)
  })))
  
  # Create source/status annotation
  ta <- HeatmapAnnotation(
    df = data_meta[, c("Source", "Status")],
    col = ann_colors
  )
  
  #Add sample names at the bottom of the heatmap
  ba = columnAnnotation(foo = anno_mark(at = which(data_meta$Status == "proband"), 
                                        labels = data_meta[ which(data_meta$Status == "proband"),]$IDs,
                                        side="bottom"))
  
  
  # Create heatmap
  htm <- Heatmap(
    norm_beta,
    top_annotation = ta,
    bottom_annotation = ba,
    col = colorRamp2(
      seq(-2, 2, length = 3),
      c("#58b0ff", "black", "#ffc000")
    ),
    show_column_names = FALSE,
    show_row_names = F,
    name = " ",
    heatmap_legend_param = list(direction = "horizontal")
  ) 
  # Convert to grob for compatibility with ggplot2
  grid.grabExpr(draw(htm,
                     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
                     ))
}

#' Create dimension reduction plots
#'
#' @param data_beta Beta values
#' @param data_meta Metadata
#' @param proband Proband ID
#' @param signature_name Signature name
#' @return Combined plot object
create_dimension_reduction_plots <- function(data_beta, data_meta, proband, signature_name) {
  goi <- str_split(signature_name, "_", simplify = TRUE)[,1]

  data_meta <- data_meta[data_meta$Status == goi |
                         data_meta$Status %in% c("in_silico_case", "control") |
                         data_meta$IDs %in% proband,]

  data_beta <- data_beta[, colnames(data_beta) %in% rownames(data_meta)]
  data_beta <- as.data.frame(t(data_beta))
  data_beta <- data_beta[!duplicated(data_beta), ]
  data_beta <- as.data.frame(t(data_beta))
  data_meta <- data_meta[rownames(data_meta) %in% colnames(data_beta),]
  data_beta <- data_beta[, match(rownames(data_meta), colnames(data_beta))]
  # PCA plot
  
  p <- pca(data_beta, data_meta)
  pca_plot <- create_pca_plot(p)

  # t-SNE plot
  tsne_plot <- create_tsne_plot(data_beta, data_meta)

  # UMAP plot
  umap_plot <- create_umap_plot(data_beta, data_meta)

  # Heatmap
  heatmap_plot <- create_heatmap(data_beta, data_meta)

  # Combine plots
  
  design <- "AABBCC
             AABBCC
             DDDDDD
             DDDDDD
             DDDDDD"
  
  # combined_plot <- (pca_plot + tsne_plot + umap_plot) / heatmap_plot +
  combined_plot <- pca_plot + tsne_plot + umap_plot + heatmap_plot +
    plot_layout(design = design)+
    plot_annotation(title = signature_name,
                   subtitle = paste("Proband:", proband))

  return(combined_plot)
}




#' Create UMAP plot
#'
#' @param cell_prop_input Input to calculate cell proportions, ie. beta values of sample(s) to test
#' @param cellproportions_background Cell proportions fo the samples used for the tests 
#' @return ggplot for cell proportions
create_cell_deconv_plot <- function(samples_of_interest_beta, cellproportions_background, probands) {
  
  rownames(samples_of_interest_beta) = samples_of_interest_beta$IlmnID
  samples_of_interest_beta$IlmnID = NULL
  BloodFrac.m <- epidish(beta.m = samples_of_interest_beta,
                         ref.m = centDHSbloodDMC.m,
                         method = "RPC")$estF
  
  bf = as.data.frame(t(BloodFrac.m))
  bf$CellType = rownames(bf)
  
  bf = pivot_longer(bf, -"CellType",
                    names_to="Sample",
                    values_to="CellProp")
  
  ggplot() + 
    geom_jitter(cellproportions_background, mapping=aes(CellType, CellProp),
                size = 1, width = 0.25, height = 0, alpha=0.20) +
    geom_jitter(bf[bf$Sample %in% probands,], mapping=aes(CellType, CellProp, color=Sample),
                size = 5,  width = 0.25, height = 0) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    ggtitle("Cell proportions of tested samples compared to that of samples in SVM training set")+
    ylab("Deconvoluted cell proportions") +
    xlab("") +
    ylim(0, NA) # Ensure y-axis starts at 0
  
  
}
