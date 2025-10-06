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
    ggtitle("Cell proportions of tested samples compared to that of samples in SVM training set") +
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