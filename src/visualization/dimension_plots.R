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

  # Prepare control samples
  controls_plot_beta <- lapply(signatures, function(x) {
    y <- insilico_beta[insilico_beta$IlmnID %in% x$ProbeID, ]
    rownames(y) <- y$IlmnID
    y <- y[match(x$ProbeID, rownames(y)), ]
    return(as.data.frame(y))
  })

  data_for_plots_beta <- lapply(names(signatures), function(s) {
    # Clean column names before joining
    clean_df <- function(df) {
      if(is.null(df) || nrow(df) == 0) return(df)
      # Remove any completely empty columns 
      df <- df[, !sapply(df, function(x) {
        tryCatch({
          all(is.na(x) | x == "")
        }, error = function(e) {
          FALSE
        })
      }), drop = FALSE]
      # Ensure column names are not empty
      colnames(df)[colnames(df) == ""] <- paste0("col_", seq_along(colnames(df)[colnames(df) == ""]))
      return(df)
    }

    df_list <- list(
      clean_df(controls_plot_beta[[s]]),
      clean_df(insilicocases_plot_beta[[s]]),
      clean_df(user_plot_beta[[s]]),
      clean_df(real_cases_beta[[s]])
    )

    # Remove NULL or empty data frames
    df_list <- df_list[sapply(df_list, function(x) !is.null(x) && nrow(x) > 0)]

    df_list <- lapply(df_list, function(x) {x$IlmnID = rownames(x); return(x)})

    if(length(df_list) > 0) {
      lst <- purrr::reduce(df_list, full_join, by = "IlmnID")
      lst$IlmnID <- NULL
      return(lst)
    } else {
      return(data.frame())
    }
  })

  names(data_for_plots_beta) <- names(signatures)
  return(data_for_plots_beta)
}

#' Prepare metadata for plots
#'
#' @param signatures List of signatures
#' @param plot_data Combined plot data
#' @param user_plot_beta User plot beta values
#' @param real_cases_meta Real cases metadata
#' @param insl In-silico metadata
#' @return List of prepared metadata for plots
prepare_plot_metadata <- function(signatures, plot_data, user_plot_beta, real_cases_meta, insl) {
  data_for_plots_meta <- lapply(names(signatures), function(s) {
    # Get all sample IDs from the combined plot data for this signature
    all_sample_ids <- colnames(plot_data[[s]])

    # cat("Processing metadata for signature:", s, "\n")
    # cat("Total sample IDs in plot_data:", length(all_sample_ids), "\n")

    # User test IDs
    mt_ids_usertests <- setdiff(colnames(user_plot_beta), "IlmnID")
    # cat("User test IDs:", length(mt_ids_usertests), "\n")

    # In-silico cases have '_isc' suffix
    mt_ids_iscases <- all_sample_ids[grepl("_isc$", all_sample_ids)]
    # cat("In-silico case IDs:", length(mt_ids_iscases), "\n")

    # Control IDs are those in insl that are in plot_data but NOT in-silico cases or user tests
    mt_ids_controls <- all_sample_ids[all_sample_ids %in% insl$IDs &
                                      !grepl("_isc$", all_sample_ids) &
                                      !all_sample_ids %in% mt_ids_usertests]
    # cat("Control IDs:", length(mt_ids_controls), "\n")

    # Real cases are any remaining samples 
    mt_ids_realcases <- all_sample_ids[all_sample_ids %in% real_cases_meta$IDs]
    # cat("Real case IDs:", length(mt_ids_realcases), "\n")

    # Create metadata for all samples found in plot_data
    lsd <- data.frame(
      'IDs' = c(mt_ids_controls, mt_ids_iscases, mt_ids_usertests, mt_ids_realcases),
      'Status' = c(rep("control", length(mt_ids_controls)),
                   rep("in_silico_case", length(mt_ids_iscases)),
                   rep("in_silico_case", length(mt_ids_usertests)),
                   rep("case", length(mt_ids_realcases))),
      'Source' = c(rep("literature", length(mt_ids_controls)),
                   rep("literature", length(mt_ids_iscases)),
                   rep("literature", length(mt_ids_usertests)),
                   rep("literature", length(mt_ids_realcases)))
    )

    # For in-silico cases, create a mapping to their source control IDs so to inherit Sex and AgeGroup from the source sample
    lsd$source_ID <- gsub("_isc$", "", lsd$IDs)

    cols_to_merge <- c("IDs", "Platform")
    if("Sex" %in% names(insl)) cols_to_merge <- c(cols_to_merge, "Sex")
    if("AgeGroup" %in% names(insl)) cols_to_merge <- c(cols_to_merge, "AgeGroup")

    insl_for_merge <- insl[, cols_to_merge]
    names(insl_for_merge)[names(insl_for_merge) == "IDs"] <- "source_ID"

    # Merge to get Platform, Sex, and AgeGroup from source controls
    lsd <- merge(lsd, insl_for_merge, by = "source_ID", all.x = TRUE)

    # Remove the temporary source_ID column
    lsd$source_ID <- NULL

    # cat("After merge, lsd has", nrow(lsd), "rows with IDs:", paste(head(lsd$IDs, 10), collapse=", "), "...\n")

    lsd <- lsd[!lsd$IDs %in% real_cases_meta$IDs, ]
    lsd <- bind_rows(lsd, real_cases_meta)

    lsd$Status <- ifelse(lsd$IDs %in% mt_ids_usertests, "proband", lsd$Status)
    lsd$Source <- ifelse(lsd$IDs %in% mt_ids_usertests, "user_sample", lsd$Source)

    lsd <- lsd[!duplicated(lsd),]
    rownames(lsd) <- lsd$IDs
    return(lsd)
  })

  names(data_for_plots_meta) <- names(signatures)
  return(data_for_plots_meta)
}

#' Create color scheme for plots
#'
#' @param test_id Test ID
#' @return List of annotation colors
create_annotation_colors <- function(test_id) {
  base_test_ids <- unique(test_id[!test_id %in% c("in_silico_case", "control", "proband")])

  status_colors <- c('#FB8C00',
                     "#FCDE9C", 
                     '#089099',
                     "#DC3977")

  status_names <- c(ifelse(length(base_test_ids) > 0, base_test_ids[1], "cases"),
                    "in_silico_case",
                    'control',
                    'proband')

  if(length(base_test_ids) > 1) {
    # In case more labels matching the syndrome
    additional_colors <- rep('#FB8C00', length(base_test_ids) - 1)
    status_colors <- c(status_colors[1], additional_colors, status_colors[-1])
    status_names <- c(base_test_ids, "in_silico_case", 'control', 'proband')
  }

  ann_colors <- list(
    Status = status_colors[1:length(status_names)],
    Platform = c("#089099", "#045275", "#7C1D6F", "#DC3977")
  )

  names(ann_colors$Status) <- status_names
  names(ann_colors$Platform) <- c('EpicV2', 'EpicV1', '450k', 'proband')

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

  # Ensure all status values in metadata have colors
  unique_statuses <- unique(pca_object$metadata$Status)
  missing_statuses <- setdiff(unique_statuses, names(ann_colors$Status))

  if(length(missing_statuses) > 0) {
    # Add default colors for missing status values
    default_colors <- rainbow(length(missing_statuses))
    names(default_colors) <- missing_statuses
    ann_colors$Status <- c(ann_colors$Status, default_colors)
  }

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
                            alpha = 1)) +
    geom_point(pca_plot_data[pca_plot_data$Status == "proband",],
               mapping=aes(PC1, PC2,
                            color = Status,
                            shape = (Status == "proband"),
                            size = 1,
                            alpha = 1)) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(),
      panel.border = element_rect(color = "black", fill = NA, size = 1)
    ) +
    xlab(paste0("PC1 (", variance_pct[1], "%)")) +
    ylab(paste0("PC2 (", variance_pct[2], "%)")) +
    ggtitle("PCA") +
    scale_colour_manual(values = ann_colors$Status) +
    guides(size = "none", shape = "none", alpha = "none")
}

#' Create heatmap
#'
#' @param data_beta Beta values
#' @param data_meta Metadata
#' @return Heatmap grob
create_heatmap <- function(data_beta, data_meta, age_table = NULL, chr_sex_table = NULL) {
  # Get test_id for color scheme
  test_id <- setdiff(unique(data_meta$Status),
                     c("in_silico_case", "control", "proband"))

  # Create color scheme
  ann_colors <- create_annotation_colors(test_id)

  # Ensure all status values in data_meta have colors
  unique_statuses <- unique(data_meta$Status)
  missing_statuses <- setdiff(unique_statuses, names(ann_colors$Status))

  if(length(missing_statuses) > 0) {
    # Add default colors for missing status values
    default_colors <- rainbow(length(missing_statuses))
    names(default_colors) <- missing_statuses
    ann_colors$Status <- c(ann_colors$Status, default_colors)
  }

  # Add colors for Sex and AgeGroup if present in metadata
  if("Sex" %in% names(data_meta)) {
    # Sex values: "Male" or "Female"
    ann_colors$Sex <- c("Male" = "#4169E1", "Female" = "#FF69B4")
  }

  if("AgeGroup" %in% names(data_meta)) {
    ann_colors$AgeGroup <- c(
      "Infant new born" = "#FFF4E6",
      "Infant" = "#FFE0B2",
      "Preschool child" = "#FFCC80",
      "Child" = "#FFB74D",
      "Adolescent" = "#FFA726",
      "Adult" = "#FF9800",
      "MiddleAged" = "#FB8C00",
      "Aged65plus" = "#E65100"
    )
    
    age_levels <- c("Infant new born", "Infant", "Preschool child", "Child",
                    "Adolescent", "Adult", "MiddleAged", "Aged65plus")
    data_meta$AgeGroup <- factor(data_meta$AgeGroup, levels = age_levels, ordered = TRUE)
  }

  # Normalize beta values
  norm_beta <- as.data.frame(t(apply(data_beta, 1, function(x) {
    (x - mean(x)) / sd(x)
  })))

  # annotation  for plot
  annot_cols <- c("Status", "Platform")
  if("Sex" %in% names(data_meta)) annot_cols <- c(annot_cols, "Sex")
  if("AgeGroup" %in% names(data_meta)) annot_cols <- c(annot_cols, "AgeGroup")

  # Create annotation with Status at top and double height
  annotation_heights <- rep(1, length(annot_cols))
  annotation_heights[1] <- (3)

  ta <- HeatmapAnnotation(
    df = data_meta[, annot_cols, drop = FALSE],
    col = ann_colors[names(ann_colors) %in% annot_cols],
    annotation_height = annotation_heights
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

#' Calculate pairwise distances between test samples and candidate samples
#'
#' @param test_beta Beta matrix for test samples
#' @param candidate_beta Beta matrix for candidate samples
#' @return Named vector of mean distances for each candidate sample
calculate_sample_distances <- function(test_beta, candidate_beta) {
  # Calculate distances between each test sample and each candidate sample
  n_test <- ncol(test_beta)
  n_candidates <- ncol(candidate_beta)

  all_distances <- matrix(NA, nrow = n_test, ncol = n_candidates)
  rownames(all_distances) <- colnames(test_beta)
  colnames(all_distances) <- colnames(candidate_beta)

  for(i in 1:n_test) {
    test_values <- test_beta[, i]

    for(j in 1:n_candidates) {
      candidate_values <- candidate_beta[, j]

      # Find positions where both have non-missing values
      valid_positions <- !is.na(test_values) & !is.na(candidate_values)

      if(sum(valid_positions) > 0) {
        # Calculate Euclidean distance using only non-missing CpGs
        all_distances[i, j] <- sqrt(sum((test_values[valid_positions] - candidate_values[valid_positions])^2))
      }
    }
  }

  # Return mean distance across all test samples for each candidate
  colMeans(all_distances, na.rm = TRUE)
}

#' Create dimension reduction plots
#'
#' @param data_beta Beta values
#' @param data_meta Metadata
#' @param proband Proband ID
#' @param signature_name Signature name
#' @param age_table Optional age prediction table
#' @param chr_sex_table Optional chromosomal sex prediction table
#' @param n_samples_per_group Number of samples per group 
#' @return Combined plot object
create_dimension_reduction_plots <- function(data_beta, data_meta, proband, signature_name,
                                            age_table = NULL, chr_sex_table = NULL,
                                            n_samples_per_group = 20) {

  # Add debugging information
  # cat("Creating dimension reduction plot for:", signature_name, "\n")
  # cat("Data beta dimensions:", dim(data_beta), "\n")
  # cat("Data meta dimensions:", dim(data_meta), "\n")
  # cat("Proband(s):", proband, "\n")

  # Check for required data
  if(is.null(data_beta) || ncol(data_beta) == 0 || nrow(data_beta) == 0) {
    stop("data_beta is NULL or empty for signature: ", signature_name)
  }

  if(is.null(data_meta) || nrow(data_meta) == 0) {
    stop("data_meta is NULL or empty for signature: ", signature_name)
  }

  goi <- str_split(signature_name, "_", simplify = TRUE)[,1]

  # Get proband data for distance calculation
  proband_samples <- data_meta[data_meta$IDs %in% proband, ]

  if(nrow(proband_samples) == 0) {
    stop("No proband samples found in metadata for: ", paste(proband, collapse=", "))
  }

  # cat("Proband IDs from metadata:", paste(proband_samples$IDs, collapse=", "), "\n")
  # cat("Available columns in data_beta:", paste(head(colnames(data_beta), 20), collapse=", "), "...\n")

  proband_beta <- data_beta[, colnames(data_beta) %in% proband_samples$IDs, drop = FALSE]

  if(ncol(proband_beta) == 0) {
    cat("ERROR: Proband(s) not found in beta data\n")
    cat("Requested proband IDs:", paste(proband, collapse=", "), "\n")
    cat("IDs in metadata:", paste(data_meta$IDs, collapse=", "), "\n")
    cat("Columns in beta data:", paste(colnames(data_beta), collapse=", "), "\n")
    stop("No proband beta values found for: ", paste(proband, collapse=", "))
  }

  # Separate different types of samples
  original_insilico_cases <- data_meta[data_meta$Status == "in_silico_case", ]

  signature_name_base <- str_split(signature_name, "_", simplify = TRUE)[,1]

  real_cases_mask <- sapply(data_meta$Status, function(status) {
    status_escaped <- gsub("([.\\-\\+\\*\\?\\[\\]\\(\\)\\{\\}\\^\\$\\|\\\\])", "\\\\\\1", status)
    grepl(paste0("^", status_escaped, "($|[._])"), signature_name, ignore.case = TRUE)
  }) & !data_meta$IDs %in% proband &
       data_meta$Status != "in_silico_case" &
       data_meta$Status != "control"

  real_cases_available <- data_meta[real_cases_mask, ]

  # cat("Signature name:", signature_name, "\n")
  # cat("Found", nrow(real_cases_available), "real cases matching signature\n")
  if(nrow(real_cases_available) > 0) {
    cat("Real case statuses:", paste(unique(real_cases_available$Status), collapse=", "), "\n")
  }

  controls_available <- data_meta[data_meta$Status == "control", ]

  # Calculate distances for controls using all CpGs and only non-missing values
  if(nrow(controls_available) > 0) {
    controls_beta <- data_beta[, colnames(data_beta) %in% controls_available$IDs, drop = FALSE]

    # cat("Calculating distances to", ncol(controls_beta), "control samples\n")
    controls_distances <- calculate_sample_distances(proband_beta, controls_beta)

    # Select closest N controls
    closest_controls_ids <- names(sort(controls_distances))[1:min(n_samples_per_group, length(controls_distances))]
    controls_meta <- controls_available[controls_available$IDs %in% closest_controls_ids, ]
    cat("Selected", nrow(controls_meta), "closest controls\n")
  } else {
    controls_meta <- data.frame()
  }

  # Calculate distances for real cases using all CpGs and only non-missing values
  if(nrow(real_cases_available) > 0) {
    real_cases_beta <- data_beta[, colnames(data_beta) %in% real_cases_available$IDs, drop = FALSE]

    # cat("Calculating distances to", ncol(real_cases_beta), "real case samples\n")
    real_cases_distances <- calculate_sample_distances(proband_beta, real_cases_beta)

    closest_real_cases_ids <- names(sort(real_cases_distances))[1:min(n_samples_per_group, length(real_cases_distances))]
    real_cases_meta <- real_cases_available[real_cases_available$IDs %in% closest_real_cases_ids, ]
    n_real_cases <- nrow(real_cases_meta)
    # cat("Selected", n_real_cases, "closest real cases\n")
  } else {
    real_cases_meta <- data.frame()
    n_real_cases <- 0
  }

  # Include probands
  proband_meta <- data_meta[data_meta$IDs %in% proband, ]

  # Only include in-silico cases if we need more cases to reach N total
  if(n_real_cases < n_samples_per_group) {
    # Calculate how many in-silico cases we need
    n_insilico_needed <- n_samples_per_group - n_real_cases

    if(nrow(original_insilico_cases) > 0) {
      # Select in-silico cases based on distance of their source controls to proband
      # In-silico case IDs have "_isc" suffix, so we need to strip it to get the control ID
      insilico_source_ids <- gsub("_isc$", "", original_insilico_cases$IDs)

      # Check if we have control distances calculated
      if(exists("controls_distances") && length(controls_distances) > 0) {
        # Match in-silico cases to their source control distances
        insilico_distances <- controls_distances[insilico_source_ids]
        # Remove NAs (in-silico cases whose source control wasn't in controls_distances)
        valid_insilico <- !is.na(insilico_distances)
        insilico_distances <- insilico_distances[valid_insilico]
        valid_insilico_meta <- original_insilico_cases[valid_insilico, ]

        # Select closest in-silico cases based on their source control distances
        if(length(insilico_distances) > 0) {
          closest_insilico_indices <- order(insilico_distances)[1:min(n_insilico_needed, length(insilico_distances))]
          insilico_cases_selected <- valid_insilico_meta[closest_insilico_indices, ]
          cat("Selected", nrow(insilico_cases_selected), "closest in-silico cases based on source control distances:\n")
          cat("  In-silico case IDs:", paste(insilico_cases_selected$IDs, collapse=", "), "\n")
          cat("  Source control IDs:", paste(gsub("_isc$", "", insilico_cases_selected$IDs), collapse=", "), "\n")
        } else {
          insilico_cases_selected <- data.frame()
        }
      } else {
        # Fallback: if no control distances, take first N in-silico cases
        cat("Warning: No control distances available, using first", n_insilico_needed, "in-silico cases\n")
        insilico_cases_selected <- original_insilico_cases[1:min(n_insilico_needed, nrow(original_insilico_cases)), ]
      }
    } else {
      insilico_cases_selected <- data.frame()
    }
  } else {
    # If we have N or more real cases, take only the first N and no in-silico cases
    real_cases_meta <- real_cases_meta[1:n_samples_per_group, ]
    insilico_cases_selected <- data.frame()
  }

  # Combine metadata with consistent columns
  # Ensure all data frames have the same columns in the same order
  all_meta_list <- list(proband_meta, controls_meta, real_cases_meta, insilico_cases_selected)

  # Remove empty data frames
  all_meta_list <- all_meta_list[sapply(all_meta_list, nrow) > 0]

  if(length(all_meta_list) > 0) {
    # Get common column names
    all_cols <- unique(unlist(lapply(all_meta_list, colnames)))

    # Ensure all data frames have the same columns
    all_meta_list <- lapply(all_meta_list, function(df) {
      missing_cols <- setdiff(all_cols, colnames(df))
      if(length(missing_cols) > 0) {
        df[missing_cols] <- NA
      }
      return(df[, all_cols, drop = FALSE])
    })

    data_meta <- do.call(rbind, all_meta_list)
  } else {
    data_meta <- data.frame()
  }

  data_beta <- data_beta[, colnames(data_beta) %in% rownames(data_meta)]
  data_beta <- as.data.frame(t(data_beta))
  data_beta <- data_beta[!duplicated(data_beta), ]
  data_beta <- as.data.frame(t(data_beta))
  data_meta <- data_meta[rownames(data_meta) %in% colnames(data_beta),]
  data_beta <- data_beta[, match(rownames(data_meta), colnames(data_beta))]
  
  data_beta = na.omit(data_beta)

  # PCA
  p <- pca(data_beta, data_meta, rank = 2)
  pca_plot <- create_pca_plot(p)

  # Heatmap
  heatmap_plot <- create_heatmap(data_beta, data_meta, age_table, chr_sex_table)

  # Combine plots
  design <- "AA
             BB"

  combined_plot <- pca_plot + heatmap_plot +
    plot_layout(design = design)+
    plot_annotation(title = signature_name)

  return(combined_plot)
}