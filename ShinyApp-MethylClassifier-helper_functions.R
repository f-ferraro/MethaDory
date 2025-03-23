# Load required libraries
library(ggplotify)
library(umap)
library(patchwork)
library(PCAtools)
library(data.table)
library(tidyverse)
library(caret)
library(methyLImp2)
library(BiocParallel)
library(ComplexHeatmap)
library(circlize)
library(Rtsne)
library(EpiDISH)

# Function to load and prepare test data
load_test_data <- function(file_path) {
  test_data_user <- read.table(file_path, header = TRUE)
  test_data <- test_data_user %>% relocate(IlmnID)
  return(list(
    test_data = test_data,
    test_data_ids = setdiff(names(test_data), "IlmnID")
  ))
}

# Modified load_background_data function in helper_functions.R
load_background_data <- function(model_dir) {
  # Load imputation background
  imputation_background <- lapply(list.files("imputation", "qs",full.names = T),qs::qread)
  imputation_background <- bind_rows(imputation_background)
  imputation_background$IlmnID = rownames(imputation_background)
  
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
  
  return(list(
    imputation_background = imputation_background,
    svm = svm
  ))
}


# Function to merge and prepare data for imputation
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

# Function to perform imputation
perform_imputation <- function(test_data, test_sample_ids) {
  manifest <- qs::qread("manifest.qc_filtered.qs")
  manifest$MAPINFO <- NULL
  names(manifest) <- c('cpg', 'chr')
  
  beta_SE_imputed <- methyLImp2(input = t(test_data),
                                type = "user",
                                annotation = manifest,
                                BPPARAM = SnowParam(exportglobals = FALSE,
                                                    workers = 1))
  df =as.data.frame(t(beta_SE_imputed))
  df$IlmnID = rownames(df)
  
  df = df[, names(df) %in% c("IlmnID", test_sample_ids)]
 
  return(df)
}

# Function to prepare data for inference
prepare_inference_data <- function(test_data, svm) {
  lapply(svm, function(x) {
    z = test_data[match(predictors(x), test_data$IlmnID), ]
    z$IlmnID = NULL
    return(z)
  })
}

# Function to make predictions
make_predictions <- function(test_for_inference, svm, test_data_ids) {
  results <- list()
  
  for (i in names(test_for_inference)) {

    for (j in names(svm)) {
      if (grepl(i, j)) {
        # print("Processing",i)
        results[[paste(i, j)]] <- predict(svm[[j]],
                                          newdata = t(test_for_inference[[i]]),
                                          type = "prob")
        results[[paste(i, j)]]$SampleID <- names(test_for_inference[[i]])
      }
    }
  }

  process_results(results,test_data_ids)
}

# Function to process prediction results
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
    summarise(pCase = mean(SVM_score),
              pCase_sd = sd(SVM_score))
}

# Function to load beta values and signatures
load_beta_signatures <- function() {
  insilico_beta <- qs::qread("samples_for_insilico.beta.qs")
  insilico_meta <- qs::qread("samples_for_insilico.meta.qs")
  signatures <- read.table('../0b.Signatures_Processing/age70signatures_signature500/merged_signatures_90DMRs.tsv', header = TRUE)
  signatures <- split(signatures, f = as.factor(paste(signatures$Label)))
  
  return(list(
    insilico_beta = insilico_beta,
    insilico_meta = insilico_meta,
    signatures = signatures
  ))
}

# Function to load and prepare real cases
load_real_cases <- function(signatures) {
  real_cases_beta <- qs::qread("realcases/subset.realcases.plots.beta.qs")
  real_cases_meta <- qs::qread("realcases/subset.realcases.plots.meta.qs")
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

# Function to prepare insilico cases data
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

# Helper function to prepare user plots
prepare_user_plots <- function(signatures, test_data_user) {
  lapply(signatures, function(x) {
    # Subset user data for relevant probes
    y = test_data_user[test_data_user$IlmnID %in% x$ProbeID, ]
    rownames(y) = y$IlmnID
    y <- y[match(x$ProbeID, rownames(y)), ]
    
    
    return(as.data.frame(y))
  })
}
# Function to prepare plot data
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

# Function to prepare plot metadata
prepare_plot_metadata <- function(signatures, insilicocases_plot_beta, user_plot_beta, real_cases_meta) {
  
  data_for_plots_meta <- lapply(names(signatures), function(s) {
    
    mt_ids_iscases <- setdiff(colnames(insilicocases_plot_beta[[s]]), "IlmnID")
    
    mt_ids_usertests <- setdiff(colnames(user_plot_beta), "IlmnID")
    
    lsd <- data.frame(
      'IDs' = c(mt_ids_iscases, mt_ids_usertests),
      'Status' = c(rep("in_silico_case", length(c(mt_ids_iscases, mt_ids_usertests)))),
      'Source' = c(rep("literature", length(c(mt_ids_iscases, mt_ids_usertests))))        
    )           
    
    lsd <- lsd[!lsd$IDs %in% real_cases_meta$IDs, ]
    lsd <- rbind(lsd, real_cases_meta)
    lsd$Status = ifelse(lsd$IDs %in% mt_ids_usertests, "proband", lsd$Status)
    lsd$Source = ifelse(lsd$IDs %in% mt_ids_usertests, "user_sample", lsd$Source)
    
    lsd <- lsd[!duplicated(lsd),]
    rownames(lsd) <- lsd$IDs
    return(lsd)
  })
  
  names(data_for_plots_meta) <- names(signatures)
  return(data_for_plots_meta)
}

# Function to create prediction plots
create_prediction_plot <- function(results, proband) {
  # & results$SVM %in% signature
  ggplot(results[results$SampleID == proband,],
         aes(SVM, pCase, color = pCase)) +
    geom_hline(yintercept = c(0, 0.5, 1), color = "darkgrey") +
    geom_pointrange(aes(ymin = ifelse(pCase - pCase_sd < 0, 0,pCase - pCase_sd),
                        ymax = ifelse(pCase + pCase_sd > 1, 1, pCase + pCase_sd)),
                    position = position_jitter(width = 0, height = 0),
                    show.legend = FALSE) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    ggtitle(proband) +
    ylab("Probability") +
    xlab("Predictive Model") +
    scale_colour_gradient(low = "#132B43", high = "darkred")
}

# Function to create color scheme for plots
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

# Function to create PCA plot
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
  ggplot(pca_plot_data) +
    geom_point(aes(PC1, PC2, 
                   color = Status,
                   shape = (Status == "proband"), 
                   size=1,
                   alpha=0.85)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlab(paste0("PC1 (", variance_pct[1], "%)")) +
    ylab(paste0("PC2 (", variance_pct[2], "%)")) +
    ggtitle("PCA") +
    scale_colour_manual(values = ann_colors$Status)+
    guides(size = "none",shape="none",alpha="none")
}

# Function to create t-SNE plot
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
    geom_point(aes(X, Y,
                   color = Status,
                   shape = Status == "proband", 
                   size=1,
                   alpha=0.85)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlab("t-SNE 1") +
    ylab("t-SNE 2") +
    ggtitle("t-SNE") +
    scale_colour_manual(values = ann_colors$Status)+
    guides(size = "none",shape="none",alpha="none")
}

# Function to create UMAP plot
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
    geom_point(aes(UMAP1, UMAP2,
                   color = Status,
                   shape = Status == "proband", 
                   size=1,
                   alpha=0.85)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin()
    ) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    ggtitle("UMAP") +
    scale_colour_manual(values = ann_colors$Status)+
    guides(size = "none",shape="none",alpha="none")
}

# Function to create heatmap
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
  
  # Create annotation
  ha <- HeatmapAnnotation(
    df = data_meta[, c("Source", "Status")],
    col = ann_colors
  )
  
  # Create heatmap
  htm <- Heatmap(
    norm_beta,
    top_annotation = ha,
    col = colorRamp2(
      seq(-2, 2, length = 3),
      c("#58b0ff", "black", "#ffc000")
    ),
    show_column_names = FALSE,
    name = " "
  )
  
  # Convert to grob for compatibility with ggplot2
  grid.grabExpr(draw(htm))
}

# Function to create dimension reduction plots
create_dimension_reduction_plots <- function(data_beta, data_meta, proband, signature_name) {
  
  goi <- str_split(signature_name, "_", simplify = T)[,1]
  
  data_meta <- data_meta[data_meta$Status == goi |
                           data_meta$Status %in%  c("in_silico_case","control") |
                           data_meta$IDs == proband,]
  
  data_beta <- data_beta[, colnames(data_beta) %in% rownames(data_meta)] 
  data_beta <- data_beta[, match(rownames(data_meta), colnames(data_beta) )]
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
  combined_plot <- (pca_plot + heatmap_plot) / (tsne_plot + umap_plot) +
    plot_annotation(title = signature_name,
                    subtitle = paste("Proband:", proband))
  
  return(combined_plot)
}