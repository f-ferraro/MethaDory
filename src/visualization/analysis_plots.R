#' Create cell proportion outlier table
#'
#' @param sample_cell_props Sample cell proportions table
#' @param background_cell_props Background cell proportions data
#' @param probands Selected probands
#' @return Data frame showing outlier status for each sample and cell type
create_cell_prop_outlier_table <- function(sample_cell_props, background_cell_props, probands) {

  # Filter for selected probands
  sample_data <- sample_cell_props[sample_cell_props$Proband %in% probands, ]

  if (nrow(sample_data) == 0) {
    return(data.frame())
  }

  # Calculate mean and SD for each cell type from background
  bg_stats <- background_cell_props %>%
    group_by(CellType) %>%
    summarise(
      min_prop = min(CellProp, na.rm=TRUE),
      max_prop = max(CellProp, na.rm=TRUE),
      mean_prop = mean(CellProp, na.rm = TRUE),
      sd_prop = sd(CellProp, na.rm = TRUE),
      .groups = "drop"
    )

  # Create outlier check table
  outlier_results <- sample_data %>%
    left_join(bg_stats, by = "CellType") %>%
    mutate(
      lower_bound = mean_prop - 3 * sd_prop,
      upper_bound = mean_prop + 3 * sd_prop,
      within_3sd = CellProp >= lower_bound & CellProp <= upper_bound,
      out_of_range = CellProp < min_prop | CellProp > max_prop,
      status = ifelse(within_3sd, "PASS", "WARNING"),
      status = ifelse(out_of_range, "FAIL", status),
      deviation = abs((CellProp - mean_prop) / sd_prop)
    ) %>%
    select(Proband, CellType, CellProp, mean_prop, sd_prop, within_3sd, status, deviation)

  # Create summary table for display
  summary_table <- outlier_results %>%
    select(Proband, CellType, CellProp, status, deviation) %>%
    mutate(
      CellProp = round(CellProp, 3),
      deviation = round(deviation, 2)
    )

  return(summary_table)
}

#' Create cell deconvolution plot
#'
#' @param create_cell_deconv_table Input cell deconvolution table
#' @param cellproportions_background Cell proportions for the samples used for the tests
#' @param probands Selected probands
#' @return ggplot for cell proportions
create_cell_deconv_plot <- function(create_cell_deconv_table, cellproportions_background, probands) {

 ggplot() +
    # Static background points
    geom_jitter(data =  cellproportions_background,
                mapping = aes(CellType, CellProp),
                size = 1, width = 0.25, height = 0, alpha = 0.20) +
    # Interactive foreground points
    geom_jitter(data = create_cell_deconv_table[create_cell_deconv_table$Proband %in% probands,],
                mapping = aes(CellType, CellProp, color = Proband),
                size = 5, width = 0.25, height = 0) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    ylab("Deconvoluted cell proportions") +
    xlab("") +
    ylim(0, NA)
}

#' Predict chromosomal sex and make plot
#'
#' @param chr_sex_table Chromosomal sex table
#' @param probands Probands of interest
#' @return ggplot for chromosomal sex prediction
predict_chr_sex_plot <- function(chr_sex_table, probands) {
  print("Plotting predicted chromosomal sex")

  x = chr_sex_table[chr_sex_table$Proband %in% probands,]

  boundary = max(x$X, abs(x$X), x$Y, abs(x$Y))

  ggplot(x, aes(X, Y, color = Proband, label = Proband)) +
    geom_vline(xintercept = 0) +
    geom_hline(yintercept = 0) +
    annotate("text", label = expression(bold("47, XXY")), x = boundary + 5, y = boundary + 5, size = 8, colour = "darkred") +
    annotate("text", label = expression(bold("46 XX")), x = boundary + 5, y = -(boundary + 5), size = 8, colour = "black") +
    annotate("text", label = expression(bold("45, X0")), x = -(boundary + 5), y = -(boundary + 5), size = 8, colour = "darkred") +
    annotate("text", label = expression(bold("46, XY")), x = -(boundary + 5), y = boundary + 5, size = 8, colour = "black") +
    xlim(-boundary - 7, boundary + 7) +
    ylim(-boundary - 7, boundary + 7) +
    geom_point(size = 6) +
    geom_text_repel(size = 7, max.overlaps = Inf) +
    theme_minimal() +
    theme(legend.position = "none") +
    ggtitle("Chromosomal sex prediction")
}
#' Sample QC PCA plot: proband against the control background
#'
#' One row per proband, showing PC1 vs PC2 and PC3 vs PC4. Controls are drawn in
#' grey and the proband is drawn last, on top, so its position is always visible.
#'
#' @param qc_pca Output of compute_qc_pca()
#' @param probands Selected probands
#' @return patchwork of ggplots, or NULL if no selected proband has a QC PCA
create_qc_pca_plot <- function(qc_pca, probands) {
  probands <- intersect(probands, names(qc_pca))
  if (length(probands) == 0) return(NULL)

  pc_panel <- function(res, pcs) {
    scores <- res$scores
    axis_label <- function(pc) sprintf("%s (%.1f%%)", pc, res$var_explained[as.integer(sub("PC", "", pc))])
    ggplot(mapping = aes(.data[[pcs[1]]], .data[[pcs[2]]])) +
      geom_point(data = scores[scores$Group == "Control", ],
                 colour = "grey70", size = 3, alpha = 0.8) +
      # Proband layer last so it sits on top of the controls
      geom_point(data = scores[scores$Group == "Proband", ],
                 colour = "#DE4C35", size = 6) +
      geom_text_repel(data = scores[scores$Group == "Proband", ],
                      mapping = aes(label = SampleID),
                      colour = "#DE4C35", size = 5, fontface = "bold",
                      max.overlaps = Inf, box.padding = 0.8) +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")) +
      xlab(axis_label(pcs[1])) +
      ylab(axis_label(pcs[2]))
  }

  rows <- lapply(probands, function(id) {
    res <- qc_pca[[id]]
    (pc_panel(res, c("PC1", "PC2")) | pc_panel(res, c("PC3", "PC4"))) +
      plot_annotation(
        title = paste0(id, " (red) vs ", sum(res$scores$Group == "Control"), " controls (grey)"),
        subtitle = paste0("Pre-imputation betas, top ", format(res$n_cpgs_used, big.mark = ","),
                          " most variable of ", format(res$n_cpgs_available, big.mark = ","),
                          " measured autosomal CpGs"),
        theme = theme(plot.title = element_text(size = 16, face = "bold"),
                      plot.subtitle = element_text(size = 12)))
  })

  wrap_plots(lapply(rows, wrap_elements), ncol = 1)
}

#' Sample QC beta-value distribution plot: proband against the controls
#'
#' One row per proband. Each control is a grey density curve and the proband is
#' drawn last, on top.
#'
#' @param qc_pca Output of compute_qc_pca()
#' @param probands Selected probands
#' @return patchwork of ggplots, or NULL if no selected proband has QC data
create_qc_density_plot <- function(qc_pca, probands) {
  probands <- intersect(probands, names(qc_pca))
  if (length(probands) == 0) return(NULL)

  rows <- lapply(probands, function(id) {
    res <- qc_pca[[id]]
    dens <- res$density
    ggplot(mapping = aes(beta, density, group = SampleID)) +
      geom_line(data = dens[dens$Group == "Control", ],
                colour = "grey70", linewidth = 0.4, alpha = 0.6) +
      # Proband layer last so it sits on top of the controls
      geom_line(data = dens[dens$Group == "Proband", ],
                colour = "#DE4C35", linewidth = 1.4) +
      theme_classic() +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")) +
      xlab("Beta value") +
      ylab("Density") +
      ggtitle(paste0(id, " (red) vs ", sum(res$scores$Group == "Control"), " controls (grey)"),
              subtitle = paste0("Pre-imputation betas, ", format(res$n_cpgs_available, big.mark = ","),
                                " measured autosomal CpGs")) +
      theme(plot.title = element_text(size = 16, face = "bold"),
            plot.subtitle = element_text(size = 12))
  })

  wrap_plots(rows, ncol = 1)
}

#' QC status line(s): missing values before imputation
#'
#' Plain HTML with inline styles so the same markup works in the Shiny app and
#' in the exported report.
#'
#' @param qc_missing Output of compute_pre_imputation_missing()
#' @param probands Selected probands
#' @return HTML string, one line per proband ("" if none)
create_qc_missing_html <- function(qc_missing, probands) {
  if (is.null(qc_missing)) return("")
  qc_missing <- qc_missing[qc_missing$SampleID %in% probands, , drop = FALSE]
  if (nrow(qc_missing) == 0) return("")

  badge_style <- c(PASS = "background-color:#2E9E5B;color:white;",
                   WARNING = "background-color:#F0C419;color:black;",
                   FAIL = "background-color:#DE4C35;color:white;")
  lines <- sprintf(
    paste0('<div style="font-size:16px;margin-bottom:10px;display:flex;align-items:center;flex-wrap:wrap;">',
           '<span style="%sfont-size:32px;line-height:1.2;font-weight:bold;padding:4px 20px;border-radius:6px;margin-right:14px;">%s</span>',
           '<span><strong>%s</strong>: %.1f%% missing values before imputation ',
           '(%s of %s CpGs used by the models)</span></div>'),
    badge_style[qc_missing$status], qc_missing$status, qc_missing$SampleID,
    qc_missing$pct_missing,
    format(qc_missing$n_missing, big.mark = ","), format(qc_missing$n_cpgs, big.mark = ","))

  paste0(paste(lines, collapse = "\n"),
         '<div style="color:#666;margin-bottom:15px;">PASS &lt; ', QC_MISSING_WARN_PCT,
         '%, WARNING ', QC_MISSING_WARN_PCT, '-', QC_MISSING_FAIL_PCT,
         '%, FAIL &gt; ', QC_MISSING_FAIL_PCT, '%</div>')
}

#' QC figure: cell proportions with the chromosomal sex prediction to its right
#'
#' @param cell_plot Output of create_cell_deconv_plot(), or NULL
#' @param chr_sex_plot Output of predict_chr_sex_plot(), or NULL
#' @return patchwork / ggplot, or NULL if both are NULL
create_qc_cell_sex_plot <- function(cell_plot, chr_sex_plot) {
  if (is.null(cell_plot)) return(chr_sex_plot)
  if (is.null(chr_sex_plot)) return(cell_plot)
  cell_plot | chr_sex_plot
}
