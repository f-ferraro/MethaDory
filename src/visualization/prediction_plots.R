#' Create prediction static plot
#'
#' @param results Prediction results
#' @param proband Proband ID
#' @param signatures Selected signatures
#' @param plot_width Plot width
#' @param threshold Minimum SVM score threshold for filtering 
#' @return ggplot object of prediction plot
create_prediction_plot <- function(results, proband, signatures, plot_width, threshold = 0) {

  plot_data <- results[results$SampleID %in% proband &
                         results$SVM %in% gsub("_", " ", signatures), ]

  # Filter by threshold if specified
  if (threshold > 0) {
    plot_data <- plot_data[plot_data$pSVM_average >= threshold, ]
  }

  # Handle empty data after filtering
  if (nrow(plot_data) == 0) {
    # Create empty plot with message
    return(ggplot() +
           annotate("text", x = 0.5, y = 0.5,
                   label = paste("No signatures above threshold", threshold),
                   size = 6, color = "gray50") +
           theme_minimal() +
           theme(axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank()) +
           ylab("SVM score") + xlab(""))
  }

  plot_data$pSVMmin <- pmax(0, plot_data$pSVM_average - plot_data$pSVM_sd)
  plot_data$pSVMmax <- pmin(1, plot_data$pSVM_average + plot_data$pSVM_sd)

  # Drop unused factor levels to fix x-axis scaling
  plot_data$SVM <- droplevels(as.factor(plot_data$SVM))

  # Create plot
  p <- ggplot(plot_data, aes(SVM, pSVM_average, color = SampleID)) +
    geom_hline(yintercept = c(0, 0.5, 1), color = "darkgrey") +
    geom_pointrange(aes(ymin = pSVMmin, ymax = pSVMmax),
                    position = position_jitter(width = 0.20, height = 0)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    scale_x_discrete(labels = scales::label_wrap(20), drop = TRUE) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "none") +
    ylab("SVM score") +
    xlab("")

  # Add threshold line if threshold > 0
  if (threshold > 0) {
    p <- p + geom_hline(yintercept = threshold, color = "red", linetype = "dashed", alpha = 0.7)
  }

  return(p)
}

#' Create interactive prediction plot for Shiny GUI
#'
#' @param results Prediction results
#' @param proband Proband ID
#' @param signatures Selected signatures
#' @param plot_width Plot width
#' @param threshold Minimum SVM score threshold for filtering
#' @return plotly object for interactive display
create_prediction_plot_interactive <- function(results, proband, signatures, plot_width, threshold = 0) {

  # Get the base ggplot with threshold filtering
  p <- create_prediction_plot(results, proband, signatures, plot_width, threshold)

  # Convert to plotly with custom layout
  ggplotly(p) %>%
    layout(yaxis = list(fixedrange = TRUE),
           margin = list(b = 120))
}

#' Create static prediction plot with score filtering for HTML export
#'
#' @param results Prediction results
#' @param proband Proband ID
#' @param signatures Selected signatures
#' @return ggplot object with filtered data
create_prediction_plot_static_filtered <- function(results, proband, signatures) {

  plot_data <- results[results$SampleID %in% proband &
                         results$SVM %in% gsub("_", " ", signatures), ]

  # Filter out signatures with pSVM_average < 0.05
  plot_data <- plot_data[plot_data$pSVM_average >= 0.05, ]

  # If no data meets criteria, return a message plot
  if(nrow(plot_data) == 0) {
    return(ggplot() +
           annotate("text", x = 0.5, y = 0.5, label = "No signatures with pSVM ≥ 0.05", size = 6) +
           theme_void() +
           xlim(0, 1) + ylim(0, 1))
  }

  plot_data$pSVMmin <- pmax(0, plot_data$pSVM_average - plot_data$pSVM_sd)
  plot_data$pSVMmax <- pmin(1, plot_data$pSVM_average + plot_data$pSVM_sd)

  ggplot(plot_data, aes(SVM, pSVM_average, color = SampleID)) +
    geom_hline(yintercept = c(0, 0.5, 1), color = "darkgrey") +
    geom_pointrange(aes(ymin = pSVMmin, ymax = pSVMmax),
                    position = position_jitter(width = 0.20, height = 0)) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal() +
    scale_x_discrete(labels = scales::label_wrap(20)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5)) +
    ylab("SVM score") +
    xlab("") +
    ggtitle("Prediction Results (pSVM ≥ 0.05)")
}