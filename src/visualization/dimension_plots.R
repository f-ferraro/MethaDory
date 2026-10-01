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
    
    # Beta distribution for noise in effect size
    a <- 10
    b <- 1
    beta_mode <- (a - 1) / (a + b - 2)
    set.seed(42);beta_values <- rbeta(n = ncol(y), shape1 = a, shape2 = b) + (1 - beta_mode)
    
    # Apply delta beta and clamp values between 0 and 1
    y <- y + outer(x$deltaBeta, beta_values)
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

#' Case labels a signature accepts as its own real cases.
#'
#' A signature is named <genes>_<Study>, and <genes> may cover several genes
#' joined by "-" (KMT2D-KDM6A, NIPBL-RAD21-SMC3-SMC1A) and carry a subtype after
#' a "." (SRCAP.FLHS, ARID1A-ARID1B.c.6200). The real cases are labelled with
#' ONE gene or subtype each (KMT2D, KDM6A, SRCAP.FLHS), so the accepted labels
#' are: the whole <genes> part, each "-"-separated gene in it, and every
#' "."-prefix of those (SMARCA2.BISS also accepts SMARCA2).
#'
#' The whole part is kept next to its pieces because "-" is also part of some
#' gene names: RNU2-2 must accept the label RNU2-2, and its pieces "RNU2" and
#' "2" are harmless since no case carries them.
#'
#' The previous rule took a label only when it was a prefix of the signature
#' name followed by the end, "." or "_". After KMT2D in KMT2D-KDM6A comes "-",
#' so every multi-gene signature lost all of its real cases.
#'
#' @param signature_name "<genes>_<Study>", possibly with a suffix after that.
#' @return character vector of accepted labels, lower case.
signature_case_labels <- function(signature_name) {
  genes <- sub("_.*$", "", signature_name)
  parts <- unique(c(genes, strsplit(genes, "-", fixed = TRUE)[[1]]))
  dot_prefixes <- function(x) {
    pieces <- strsplit(x, ".", fixed = TRUE)[[1]]
    vapply(seq_along(pieces), function(i) paste(pieces[seq_len(i)], collapse = "."),
           character(1))
  }
  labels <- unique(unlist(lapply(parts, dot_prefixes)))
  tolower(labels[nzchar(labels)])
}

#' Which Status values are real cases of this signature.
status_matches_signature <- function(status, signature_name) {
  tolower(as.character(status)) %in% signature_case_labels(signature_name)
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
    guides(size = "none", shape = "none", alpha = "none") +
    # Force a square panel so the plot's x-axis matches its y-axis. The
    # surrounding patchwork row keeps its allotted height; the empty space
    # ends up flanking the (now-narrower) PCA panel.
    theme(aspect.ratio = 1)
}

#' Create similarity-to-median-profile scatter (NSD1-paper style).
#'
#' For a signature's probe set, builds a median CASE profile (real cases +
#' in-silico synthetic cases - the same samples shown in the PCA/heatmap) and a
#' median CONTROL profile, then plots every sample by its Pearson similarity to
#' the control profile (x) vs the case profile (y). Uses the same Status colour
#' scheme as the PCA/heatmap; proband(s) are drawn as a top layer so their
#' landing point is unambiguous.
#'
#' @param data_beta probes x samples matrix/data.frame (already aligned to
#'   data_meta).
#' @param data_meta metadata aligned to data_beta columns (needs Status).
#' @return ggplot object, or NULL if there aren't enough cases/controls/probes.
create_similarity_plot <- function(data_beta, data_meta) {
  mat <- as.matrix(data_beta)
  status <- data_meta[colnames(mat), "Status"]

  # Case = real cases (any diagnosis label) + synthetic in-silico cases.
  test_id   <- setdiff(unique(status), c("in_silico_case", "control", "proband"))
  case_cols <- colnames(mat)[status %in% c(test_id, "in_silico_case")]
  ctrl_cols <- colnames(mat)[status == "control"]

  # Need a stable median on each side and enough probes for a correlation.
  if (length(case_cols) < 2 || length(ctrl_cols) < 2 || nrow(mat) < 3) {
    return(NULL)
  }

  med_case <- apply(mat[, case_cols, drop = FALSE], 1, median, na.rm = TRUE)
  med_ctrl <- apply(mat[, ctrl_cols, drop = FALSE], 1, median, na.rm = TRUE)

  sim_case <- apply(mat, 2, function(x) suppressWarnings(cor(x, med_case, use = "complete.obs")))
  sim_ctrl <- apply(mat, 2, function(x) suppressWarnings(cor(x, med_ctrl, use = "complete.obs")))

  df <- data.frame(
    Sample   = colnames(mat),
    sim_ctrl = sim_ctrl,
    sim_case = sim_case,
    Status   = status,
    stringsAsFactors = FALSE
  )
  df <- df[is.finite(df$sim_ctrl) & is.finite(df$sim_case), , drop = FALSE]
  if (nrow(df) < 2) return(NULL)

  # Same colour scheme as the PCA / heatmap Status annotation.
  ann_colors <- create_annotation_colors(test_id)
  missing_statuses <- setdiff(unique(df$Status), names(ann_colors$Status))
  if (length(missing_statuses) > 0) {
    default_colors <- rainbow(length(missing_statuses))
    names(default_colors) <- missing_statuses
    ann_colors$Status <- c(ann_colors$Status, default_colors)
  }

  # Shared square axis limits (x == y) so the diagonal is meaningful.
  lim <- range(c(df$sim_ctrl, df$sim_case), na.rm = TRUE)
  pad <- diff(lim) * 0.03
  lim <- c(lim[1] - pad, lim[2] + pad)

  ggplot() +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    # Base layer: everything except the proband.
    geom_point(df[df$Status != "proband", ],
               mapping = aes(sim_ctrl, sim_case,
                             color = Status,
                             shape = (Status == "proband")),
               size = 6, alpha = 0.75) +
    # Top layer: proband(s), larger so their landing point is unambiguous.
    geom_point(df[df$Status == "proband", ],
               mapping = aes(sim_ctrl, sim_case,
                             color = Status,
                             shape = (Status == "proband")),
               size = 12, alpha = 1) +
    scale_colour_manual(values = ann_colors$Status) +
    coord_cartesian(xlim = lim, ylim = lim) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      aspect.ratio = 1
    ) +
    xlab("Similarity to control DNAm profile") +
    ylab("Similarity to case DNAm profile") +
    ggtitle("Similarity to median profiles") +
    guides(size = "none", shape = "none", alpha = "none", color = "none")
}

#' Status colour scheme shared by the PCA / heatmap, padded for unknown labels.
#'
#' @param status character vector of Status values present in the plot.
#' @return named colour vector covering every value in `status`.
get_status_colors <- function(status) {
  test_id <- setdiff(unique(status), c("in_silico_case", "control", "proband"))
  status_colors <- create_annotation_colors(test_id)$Status
  missing_statuses <- setdiff(unique(status), names(status_colors))
  if (length(missing_statuses) > 0) {
    default_colors <- rainbow(length(missing_statuses))
    names(default_colors) <- missing_statuses
    status_colors <- c(status_colors, default_colors)
  }
  status_colors
}

#' Collapse Status into control / case / proband.
#'
#' Case = real cases (any diagnosis label) + synthetic in-silico cases, the same
#' definition create_similarity_plot() uses.
simplify_status <- function(status) {
  ifelse(status %in% c("control", "proband"), status, "case")
}

#' Empty panel with a message, so the combined layout keeps its shape when a
#' view can't be drawn (e.g. too few cases for a median profile).
placeholder_panel <- function(title, msg) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 5, colour = "gray50") +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    ggtitle(title)
}

#' Centre each probe on the median control beta (row means if < 2 controls).
#'
#' Raw betas are dominated by each probe's baseline level, so every sample
#' correlates highly with every other; centring leaves only the deviation from
#' the control profile, which is what the signature is about.
center_on_controls <- function(mat, status) {
  ctrl_cols <- which(status == "control")
  ref <- if (length(ctrl_cols) >= 2) {
    apply(mat[, ctrl_cols, drop = FALSE], 1, median, na.rm = TRUE)
  } else {
    rowMeans(mat, na.rm = TRUE)
  }
  mat - ref
}

#' PCA fitted on the reference samples only, with the proband(s) projected in.
#'
#' In create_pca_plot() the proband takes part in the fit, so a noisy or
#' platform-shifted proband can define a component by itself. Here the axes
#' come from case-vs-control variation alone and the proband is placed onto
#' them afterwards.
#'
#' @param data_beta probes x samples, complete (no NA), aligned to data_meta.
#' @param data_meta metadata aligned to data_beta columns (needs Status).
#' @param orient_like optional PC1/PC2 scores of the joint PCA (rownames = sample
#'   IDs). Component signs are arbitrary, so each axis is flipped to agree with
#'   it and the two PCA panels can be compared side by side.
#' @return ggplot object, or NULL if fewer than 3 reference samples.
create_pca_projected_plot <- function(data_beta, data_meta, orient_like = NULL) {
  mat <- as.matrix(data_beta)
  status <- data_meta[colnames(mat), "Status"]
  is_proband <- status == "proband"
  if (sum(!is_proband) < 3 || !any(is_proband)) return(NULL)

  # Same preprocessing as PCAtools::pca(): centre probes, no scaling.
  fit <- prcomp(t(mat[, !is_proband, drop = FALSE]), center = TRUE, scale. = FALSE, rank. = 2)
  projected <- predict(fit, t(mat[, is_proband, drop = FALSE]))

  if (!is.null(orient_like)) {
    for (k in 1:2) {
      r <- suppressWarnings(cor(fit$x[, k], orient_like[rownames(fit$x), k]))
      if (is.finite(r) && r < 0) {
        fit$x[, k] <- -fit$x[, k]
        projected[, k] <- -projected[, k]
      }
    }
  }

  pca_plot_data <- data.frame(
    PC1 = c(fit$x[, 1], projected[, 1]),
    PC2 = c(fit$x[, 2], projected[, 2]),
    Status = c(status[!is_proband], status[is_proband])
  )
  variance_pct <- round(fit$sdev^2 / sum(fit$sdev^2) * 100)

  status_colors <- get_status_colors(status)

  ggplot() +
    geom_point(pca_plot_data[pca_plot_data$Status != "proband",],
               mapping = aes(PC1, PC2, color = Status, shape = (Status == "proband")),
               size = 3) +
    geom_point(pca_plot_data[pca_plot_data$Status == "proband",],
               mapping = aes(PC1, PC2, color = Status, shape = (Status == "proband")),
               size = 5) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.margin = margin(),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      aspect.ratio = 1
    ) +
    xlab(paste0("PC1 (", variance_pct[1], "%)")) +
    ylab(paste0("PC2 (", variance_pct[2], "%)")) +
    ggtitle("PCA on cases + controls, proband projected") +
    scale_colour_manual(values = status_colors) +
    guides(size = "none", shape = "none", alpha = "none")
}

#' MDS (classical, via cmdscale) on pairwise-complete 1 - Spearman distances.
#'
#' Works on the matrix BEFORE na.omit(), so probes the proband is missing still
#' inform the distances among the reference samples. Betas are centred on the
#' control median first (see center_on_controls()).
#'
#' @param data_beta_full probes x samples, may contain NA, aligned to data_meta.
#' @param data_meta metadata aligned to data_beta_full columns (needs Status).
#' @return ggplot object, or NULL if there are too few usable samples.
create_pcoa_plot <- function(data_beta_full, data_meta) {
  mat <- as.matrix(data_beta_full)
  status <- data_meta[colnames(mat), "Status"]

  cm <- suppressWarnings(cor(center_on_controls(mat, status), use = "pairwise.complete.obs",
                             method = "spearman"))
  # A sample with no overlapping probes yields NA distances; cmdscale can't take them.
  keep <- rowSums(is.na(cm)) == 0
  if (sum(keep) < 3) return(NULL)
  cm <- cm[keep, keep]
  status <- status[keep]

  mds <- cmdscale(as.dist(1 - cm), k = 2, eig = TRUE)
  pos_eig <- mds$eig[mds$eig > 0]
  variance_pct <- round(pos_eig[1:2] / sum(pos_eig) * 100)

  pcoa_plot_data <- data.frame(PCo1 = mds$points[, 1], PCo2 = mds$points[, 2], Status = status)

  ggplot() +
    geom_point(pcoa_plot_data[pcoa_plot_data$Status != "proband",],
               mapping = aes(PCo1, PCo2, color = Status, shape = (Status == "proband")),
               size = 3) +
    geom_point(pcoa_plot_data[pcoa_plot_data$Status == "proband",],
               mapping = aes(PCo1, PCo2, color = Status, shape = (Status == "proband")),
               size = 5) +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      aspect.ratio = 1
    ) +
    xlab(paste0("PCo1 (", variance_pct[1], "%)")) +
    ylab(paste0("PCo2 (", variance_pct[2], "%)")) +
    ggtitle("MDS") +
    scale_colour_manual(values = get_status_colors(status))
}

#' Delta-concordance scatter: one point per CpG.
#'
#' x = signature effect (median case - median control), y = the proband's
#' deviation from the median control. A typical case lies on the diagonal
#' (slope 1), a control on the horizontal (slope 0). The fitted slope is the
#' fraction of the signature the proband carries; the intercept is a global
#' offset (typically platform) that is independent of the signature.
#'
#' @param data_beta_full probes x samples, may contain NA, aligned to data_meta.
#' @param data_meta metadata aligned to data_beta_full columns (needs Status).
#' @return ggplot object, or NULL if there aren't enough cases/controls/probes.
create_delta_concordance_plot <- function(data_beta_full, data_meta) {
  mat <- as.matrix(data_beta_full)
  status <- data_meta[colnames(mat), "Status"]
  group <- simplify_status(status)

  if (sum(group == "case") < 2 || sum(group == "control") < 2) return(NULL)

  med_case <- apply(mat[, group == "case", drop = FALSE], 1, median, na.rm = TRUE)
  med_ctrl <- apply(mat[, group == "control", drop = FALSE], 1, median, na.rm = TRUE)

  proband_ids <- colnames(mat)[group == "proband"]
  df <- do.call(rbind, lapply(proband_ids, function(id) {
    data.frame(Proband = id, effect = med_case - med_ctrl, deviation = mat[, id] - med_ctrl)
  }))
  df <- df[is.finite(df$effect) & is.finite(df$deviation), , drop = FALSE]
  if (nrow(df) < 3) return(NULL)

  fit_label <- sapply(split(df, df$Proband), function(d) {
    if (nrow(d) < 3) return(NA_character_)
    fit <- coef(lm(deviation ~ effect, d))
    sprintf("%s: slope %.2f, offset %+.3f, r %.2f",
            d$Proband[1], fit[2], fit[1], cor(d$effect, d$deviation))
  })
  fit_label <- paste(na.omit(fit_label), collapse = "\n")

  status_colors <- get_status_colors(status)
  many <- length(proband_ids) > 1

  p <- ggplot(df, aes(effect, deviation)) +
    geom_hline(yintercept = 0, color = status_colors[["control"]]) +
    geom_abline(slope = 1, intercept = 0, color = status_colors[[1]])
  p <- if (many) {
    # Proband hues stay clear of the control / case colours used for the guides.
    proband_colors <- colorRampPalette(c(status_colors[["proband"]], "#7C1D6F", "#2D2D2D"))(length(proband_ids))
    p + geom_point(aes(color = Proband), size = 2, alpha = 0.6) +
      geom_smooth(aes(color = Proband), method = "lm", formula = y ~ x, se = FALSE) +
      scale_colour_manual(values = setNames(proband_colors, proband_ids))
  } else {
    p + geom_point(size = 2, alpha = 0.6, color = "gray25") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  color = status_colors[["proband"]], fill = status_colors[["proband"]], alpha = 0.15)
  }
  p +
    theme_minimal() +
    theme(
      legend.position = if (many) "bottom" else "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      aspect.ratio = 1
    ) +
    xlab("Signature effect (median case - median control)") +
    ylab("Proband - median control") +
    labs(title = "Delta concordance (one point per CpG)", subtitle = fit_label)
}

#' Sample-sample correlation heatmap.
#'
#' Spearman (rank) correlation of control-centred betas (see
#' center_on_controls()), pairwise-complete, clustered on 1 - rho. Proband(s)
#' are marked the same way as in create_heatmap().
#'
#' Centred, because raw betas share each probe's baseline level: every pair of
#' samples then correlates near 1 and the group structure is squeezed into a
#' narrow band. On the deviations from the control profile the scale has a
#' meaningful zero - no shared deviation - so the colours run over the fixed
#' -1..1 range and are comparable between figures.
#'
#' Ranks rather than Pearson, so a few probes with a large deviation cannot
#' carry the correlation on their own.
#'
#' @param data_beta_full probes x samples, may contain NA, aligned to data_meta.
#' @param data_meta metadata aligned to data_beta_full columns (needs Status).
#' @return Heatmap grob, or NULL if there are too few usable samples.
create_sample_correlation_heatmap <- function(data_beta_full, data_meta) {
  mat <- as.matrix(data_beta_full)
  status <- data_meta[colnames(mat), "Status"]

  cm <- suppressWarnings(cor(center_on_controls(mat, status),
                             use = "pairwise.complete.obs", method = "spearman"))
  keep <- rowSums(is.na(cm)) == 0
  if (sum(keep) < 3) return(NULL)
  cm <- cm[keep, keep]
  status <- status[keep]

  hc <- hclust(as.dist(1 - cm), method = "average")
  proband_at <- which(status == "proband")


  ta <- HeatmapAnnotation(
    Status = status,
    col = list(Status = get_status_colors(status)),
    show_legend = FALSE
  )
  ra <- rowAnnotation(foo = anno_mark(at = proband_at, labels = colnames(cm)[proband_at]))

  htm <- Heatmap(
    cm,
    top_annotation = ta,
    right_annotation = ra,
    cluster_rows = hc,
    cluster_columns = hc,
    show_row_dend = FALSE,
    col = colorRamp2(c(-1, 0, 1), c("#2C5F9E", "white", "#B5361C")),
    show_column_names = FALSE,
    show_row_names = FALSE,
    name = "Spearman rho",
    column_title = "Sample-sample correlation (control-centred, Spearman)",
    column_title_gp = gpar(fontsize = 13),
    heatmap_legend_param = list(direction = "horizontal")
  )
  grid.grabExpr(draw(htm, heatmap_legend_side = "bottom"))
}

#' Ranked-neighbour strip.
#'
#' Every displayed reference sample ordered by its distance to the proband
#' (nearest on the left): a colour strip of their Status on top, the actual
#' distances below. Unlike the heatmap dendrogram this is centred on the
#' proband, has no arbitrary leaf order, and shows magnitude - a gap between the
#' last case and the first control is the feature to look for. Distance is the
#' mean absolute beta difference over the probes both samples have, so missing
#' probes don't shrink it.
#'
#' The proband itself is drawn at rank 0, distance 0, in the proband colour and
#' labelled: it is the origin every other point is measured from, so the first
#' neighbour's height reads directly as "how far is the nearest sample", and the
#' strip starts with the sample the panel is about.
#'
#' @param data_beta_full probes x samples, may contain NA, aligned to data_meta.
#' @param data_meta metadata aligned to data_beta_full columns (needs Status).
#' @return ggplot object (one facet per proband), or NULL if too few samples.
create_ranked_neighbour_plot <- function(data_beta_full, data_meta) {
  mat <- as.matrix(data_beta_full)
  status <- data_meta[colnames(mat), "Status"]
  is_proband <- status == "proband"
  if (!any(is_proband) || sum(!is_proband) < 3) return(NULL)

  ref <- mat[, !is_proband, drop = FALSE]
  ref_status <- status[!is_proband]
  n_case <- sum(simplify_status(ref_status) == "case")

  df <- do.call(rbind, lapply(colnames(mat)[is_proband], function(id) {
    d <- colMeans(abs(ref - mat[, id]), na.rm = TRUE)
    out <- data.frame(Proband = id, Status = ref_status, d = d)
    out <- out[is.finite(out$d), , drop = FALSE]
    out <- out[order(out$d), , drop = FALSE]
    out$rank <- seq_len(nrow(out))
    # Strip sits above the points, in this facet's own y range.
    span <- diff(range(out$d))
    if (span == 0) span <- max(out$d, 1e-6)
    out$strip_y <- max(out$d) + 0.22 * span
    out$strip_h <- 0.16 * span
    # Of the n_case nearest neighbours, how many are cases (all of them for a typical case).
    k <- min(n_case, nrow(out))
    n_hit <- sum(simplify_status(out$Status[seq_len(k)]) == "case")
    out$Facet <- if (k > 0) sprintf("%s: %d of the %d nearest are cases", id, n_hit, k) else id
    # The proband itself, at rank 0. Added after the ranking and the strip
    # geometry, so neither the neighbour ranks nor the hit count include it.
    self <- out[1, , drop = FALSE]
    self$Status <- "proband"
    self$d <- 0
    self$rank <- 0L
    out$is_self <- FALSE
    self$is_self <- TRUE
    rbind(self, out)
  }))
  if (is.null(df) || nrow(df) < 3) return(NULL)

  status_colors <- get_status_colors(status)

  ggplot(df) +
    geom_tile(aes(rank, strip_y, height = strip_h, fill = Status), width = 0.9) +
    geom_point(data = df[!df$is_self, , drop = FALSE],
               aes(rank, d, color = Status), size = 2.5) +
    # The proband: a larger diamond at the origin, named, so it cannot be read
    # as one more reference sample.
    geom_point(data = df[df$is_self, , drop = FALSE],
               aes(rank, d, color = Status), shape = 18, size = 5) +
    geom_text(data = df[df$is_self, , drop = FALSE],
              aes(rank, d, label = Proband, color = Status),
              hjust = -0.15, vjust = -0.9, size = 3.5, show.legend = FALSE) +
    facet_wrap(~ Facet, nrow = 1, scales = "free") +
    scale_fill_manual(values = status_colors) +
    scale_colour_manual(values = status_colors) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    xlab("Proband (rank 0), then reference samples ranked by distance to it (nearest on the left)") +
    ylab("Mean |beta difference|") +
    ggtitle("Ranked neighbours")
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

    # The source metadata codes Sex inconsistently - most samples use
    # "Male"/"Female", but some carry numeric codes (e.g. 0/1). Map any
    # unmapped non-NA levels to grey so the heatmap still renders instead of
    # erroring with "cannot map colors to some of the levels". (NA is handled
    # by ComplexHeatmap's default na_col.)
    observed_sex <- unique(as.character(data_meta$Sex))
    observed_sex <- observed_sex[!is.na(observed_sex)]
    missing_sex <- setdiff(observed_sex, names(ann_colors$Sex))
    if(length(missing_sex) > 0) {
      default_sex <- rep("#BDBDBD", length(missing_sex))
      names(default_sex) <- missing_sex
      ann_colors$Sex <- c(ann_colors$Sex, default_sex)
    }
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

  #Add sample names at the bottom of the heatmap, rotated 30 degrees so
  # they don't overlap each other for multi-proband runs.
  ba = columnAnnotation(foo = anno_mark(at = which(data_meta$Status == "proband"),
                                        labels = data_meta[ which(data_meta$Status == "proband"),]$IDs,
                                        side = "bottom",
                                        labels_rot = 30))

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

# Relative row heights of the combined per-signature figure: three rows of
# paired square panels, the ranked-neighbour strip, then the probe heatmap.
DIMENSION_PLOT_ROW_HEIGHTS <- c(1, 1, 1, 0.7, 2)
DIMENSION_PLOT_UNITS <- sum(DIMENSION_PLOT_ROW_HEIGHTS)
# The figure used to be 3 units tall (PCA row + 2x heatmap). Front-ends scale
# their old device height by this so every row keeps its previous size.
DIMENSION_PLOT_HEIGHT_SCALE <- DIMENSION_PLOT_UNITS / 3

#' Keep a single QC plot at its old size on a page sized for the per-signature
#' figure (one PDF device has one page size): the plot takes the top 3 units.
fit_to_dimension_page <- function(p) {
  p / patchwork::plot_spacer() +
    plot_layout(heights = c(3, DIMENSION_PLOT_UNITS - 3))
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

  real_cases_mask <- status_matches_signature(data_meta$Status, signature_name) &
       !data_meta$IDs %in% proband &
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
  }
  # The displayed controls are chosen further down, once the in-silico cases
  # are known, so no individual is shown both as a control and as a case.

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

  # Select the closest N controls from the remainder: an in-silico case is its
  # source control plus the signature effect, so that control is not also shown
  # as a control. Source controls are reused (closest first) only when the
  # remainder can't fill N.
  if(nrow(controls_available) > 0) {
    sorted_control_ids <- names(sort(controls_distances))
    used_as_source <- sorted_control_ids %in% gsub("_isc$", "", insilico_cases_selected$IDs)
    closest_controls_ids <- head(c(sorted_control_ids[!used_as_source],
                                   sorted_control_ids[used_as_source]),
                                 n_samples_per_group)
    controls_meta <- controls_available[controls_available$IDs %in% closest_controls_ids, ]
    n_reused <- sum(closest_controls_ids %in% sorted_control_ids[used_as_source])
    cat("Selected", nrow(controls_meta), "closest controls not used for in-silico cases",
        if (n_reused > 0) paste0("(", n_reused, " reused: not enough remaining controls)") else "", "\n")
  } else {
    controls_meta <- data.frame()
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

  # Drop duplicate samples (identical beta vectors). Transposing through a
  # matrix is safe for the dedup, but we re-anchor sample IDs from the
  # matrix dimnames instead of trusting check.names to leave hyphens alone
  # (signature labels like "KMT2D-KDM6A" otherwise get mangled to "KMT2D.KDM6A").
  mat <- as.matrix(data_beta[, colnames(data_beta) %in% rownames(data_meta), drop = FALSE])
  tmat <- t(mat)
  tmat <- tmat[!duplicated(tmat), , drop = FALSE]
  mat  <- t(tmat)

  data_beta <- as.data.frame(mat, check.names = FALSE, stringsAsFactors = FALSE)

  # If the proband has no methylation data at this signature's CpGs (common
  # for newly added signatures whose probes aren't in the user's array set,
  # or for ONT/PacBio samples missing coverage), na.omit below would drop
  # every row. Bail out with a clear message instead of letting PCA fail
  # with a cryptic rank-deficient error downstream.
  proband_cols <- intersect(colnames(data_beta),
                            data_meta$IDs[data_meta$IDs %in% proband])
  if (length(proband_cols) > 0) {
    proband_nonNA <- sum(rowSums(!is.na(data_beta[, proband_cols, drop = FALSE])) > 0)
    if (proband_nonNA == 0) {
      stop("Proband has no methylation values at any of the ",
           nrow(data_beta), " ", signature_name,
           " signature CpGs - skipping plot.")
    }
  }

  # Keep the pre-na.omit matrix for the pairwise-complete views (MDS, delta
  # concordance, correlation heatmap, ranked neighbours): they can use probes the
  # proband is missing, which na.omit drops for the PCA / heatmap.
  data_beta_full <- data_beta
  data_beta <- na.omit(data_beta)

  # Enforce strict alignment that PCAtools::pca() requires:
  # colnames(data_beta) must equal rownames(data_meta) exactly, same order.
  common <- intersect(colnames(data_beta), rownames(data_meta))
  if (length(common) == 0) {
    stop("No overlapping samples between beta matrix and metadata for signature: ",
         signature_name)
  }
  data_beta <- data_beta[, common, drop = FALSE]
  data_beta_full <- data_beta_full[, common, drop = FALSE]
  data_meta <- data_meta[common, , drop = FALSE]

  # Need >= 3 samples for a meaningful 2-component PCA. With fewer samples
  # the SVD goes rank-deficient and the downstream heatmap merge fails with
  # "differing number of rows".
  if (ncol(data_beta) < 3) {
    stop("Only ", ncol(data_beta), " sample(s) survived alignment for ",
         signature_name, " - need >= 3 for PCA. Likely too many NA probes ",
         "from the proband's input data.")
  }

  # PCA
  p <- pca(data_beta, data_meta, rank = 2)
  pca_plot <- create_pca_plot(p)

  # Heatmap
  heatmap_plot <- create_heatmap(data_beta, data_meta, age_table, chr_sex_table)

  # Similarity-to-median-profile scatter (may be NULL if too few cases/controls).
  similarity_plot <- create_similarity_plot(data_beta, data_meta)

  # Views that may be NULL when there are too few cases/controls get a
  # placeholder so the layout keeps its shape.
  too_few <- "Not enough cases / controls"
  or_placeholder <- function(p, title) if (is.null(p)) placeholder_panel(title, too_few) else p

  pca_projected_plot <- or_placeholder(create_pca_projected_plot(data_beta, data_meta,
                                                                 orient_like = p$rotated[, 1:2]),
                                       "PCA on cases + controls, proband projected")
  pcoa_plot          <- or_placeholder(create_pcoa_plot(data_beta_full, data_meta), "MDS")
  delta_plot         <- or_placeholder(create_delta_concordance_plot(data_beta_full, data_meta),
                                       "Delta concordance")
  similarity_plot    <- or_placeholder(similarity_plot, "Similarity to median profiles")
  correlation_plot   <- or_placeholder(create_sample_correlation_heatmap(data_beta_full, data_meta),
                                       "Sample-sample correlation")
  neighbour_plot     <- or_placeholder(create_ranked_neighbour_plot(data_beta_full, data_meta),
                                       "Ranked neighbours")

  # Row 1: current PCA (A)        | PCA with proband projected (B)
  # Row 2: MDS (C)                | delta concordance (D)
  # Row 3: similarity scatter (E) | sample-sample correlation heatmap (F)
  # Row 4: ranked-neighbour strip (G), full width
  # Row 5: probe heatmap (H), full width, 2x the height of a scatter row
  design <- "AB
             CD
             EF
             GG
             HH"
  combined_plot <- pca_plot + pca_projected_plot +
    pcoa_plot + delta_plot +
    similarity_plot + patchwork::wrap_elements(full = correlation_plot) +
    neighbour_plot +
    patchwork::wrap_elements(full = heatmap_plot) +
    plot_layout(design = design, heights = DIMENSION_PLOT_ROW_HEIGHTS) +
    plot_annotation(title = signature_name)

  return(combined_plot)
}
