#' @importFrom stats density
#' @importFrom ggplot2 ggplot aes geom_line geom_vline labs theme_minimal geom_point
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid grid.draw
#' @importFrom magrittr %>%
NULL

#' Compute Density and Derivatives for Ranking
#'
#' @description 
#' This function calculates the density of aggregated SHAP values and computes
#' slopes and window-based differences to help identify a cut-off (Knee/Elbow method).
#'
#' @param data Dataframe. The gene summary output from aggregation (must contain the score column).
#' @param column Character. The column name to analyze (default "aggregated_max").
#' @param window_size Numeric. The continuous window size for difference calculation.
#'
#' @return A list containing 'density_df' (x, y, slope), 'window_results', and the original 'data'.
#'
#' @export
compute_density_analysis <- function(data, column = "aggregated_max", window_size = 0.001) {
  
  vals <- data[[column]]
  if (is.null(vals)) stop("Column not found in data.")
  # 1. Compute density
  dens <- stats::density(vals)
  x_vals <- dens$x
  y_vals <- dens$y
  # 2. Compute first derivative (slope)
  # We pad the slope to match the length of x_vals
  slope <- diff(y_vals) / diff(x_vals)
  slope <- c(slope, slope[length(slope)]) 
  # 3. Window-based analysis
  results_window <- data.frame()
  start_index <- 1
  while (start_index <= length(x_vals)) {
    indices <- which(x_vals >= x_vals[start_index] & x_vals < (x_vals[start_index] + window_size))
    if (length(indices) > 1) {
      y_diff <- max(y_vals[indices]) - min(y_vals[indices])
      slope_diff <- max(slope[indices]) - min(slope[indices])
    } else {
      y_diff <- NA
      slope_diff <- NA
    }
    results_window <- rbind(results_window, data.frame(
      x_start = x_vals[start_index],
      x_end = x_vals[start_index] + window_size,
      y_diff = y_diff,
      slope_diff = slope_diff
    ))
    start_index <- max(indices) + 1
  }

  return(list(
    density_df = data.frame(x = x_vals, y = y_vals, slope = slope),
    window_results = results_window,
    original_data = data,
    column_used = column
  ))
}

#' Plot Density Diagnostic for Threshold Selection
#' 
#' @description 
#' This function takes the output from `compute_density_analysis` and generates two plots:
#' a density curve of the specified column and a line plot of the slope differences across windows
#'
#' @param analysis_list List. The output from `compute_density_analysis`.
#' @param suggested_threshold Numeric. Optional value to draw a vertical line.
#' @param output_file Character. Optional path to save the plot (e.g., "plot.pdf"). 
#' If NULL (default), the plot is displayed in the active device.
#'
#' @return A grob object (invisibly) containing the arranged plots.
#' 
#' @export
plot_density_diagnostic <- function(analysis_list, suggested_threshold = NULL, output_file = NULL) {
  
  dens_df <- analysis_list$density_df
  win_df <- analysis_list$window_results
  # Plot 1: Density Curve
  p1 <- ggplot2::ggplot(dens_df, ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_line(color = "blue") +
    ggplot2::labs(title = "Density Plot", x = analysis_list$column_used, y = "Density") +
    ggplot2::theme_minimal()
  # Plot 2: Slope Differences
  p2 <- ggplot2::ggplot(win_df, ggplot2::aes(x = .data$x_start, y = .data$slope_diff)) +
    ggplot2::geom_line(color = "red") +
    ggplot2::geom_point(color = "darkred", size = 0.5) +
    ggplot2::labs(title = "Window-Based Slope Differences", x = "X Value", y = "Slope Difference") +
    ggplot2::theme_minimal()
  if (!is.null(suggested_threshold)) {
    p1 <- p1 + ggplot2::geom_vline(xintercept = suggested_threshold, linetype = "dashed")
    p2 <- p2 + ggplot2::geom_vline(xintercept = suggested_threshold, linetype = "dashed")
  }
  # Combine plots into a "grob" (graphical object)
  combined_plot <- gridExtra::arrangeGrob(p1, p2, ncol = 1)
  # If the user provided a path, save the file
  if (!is.null(output_file)) {
    # Create the directory if it doesn't exist
    if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
    ggplot2::ggsave(filename = output_file, plot = combined_plot, width = 8, height = 10)
    cat("Plot saved to:", output_file, "\n")
  }
  # Always display the plot on the screen
  grid::grid.draw(combined_plot)

  # Return the object invisibly (so it doesn't clutter the console)
  return(invisible(combined_plot))
}

#' Select Candidate Genes by Threshold
#' 
#' @description
#' This function filters the original data based on a user-defined threshold for the specified column.
#' It returns a sorted dataframe of selected genes and prints the number and percentage of genes selected.
#'
#' @param analysis_list List. The output from `compute_density_analysis`.
#' @param threshold Numeric. The threshold decided by the user.
#'
#' @return A sorted dataframe of selected genes.
#' 
#' @export
get_candidate_genes <- function(analysis_list, threshold) {
  
  # Filter original data based on the threshold and sort by the specified column
  data <- analysis_list$original_data
  col <- analysis_list$column_used
  sgenes <- data[data[[col]] > threshold, ]
  sgenes <- sgenes[order(-sgenes[[col]]), ]
  perc <- (nrow(sgenes) / nrow(data)) * 100
  cat(sprintf("Selected %d genes (%.2f%% of total genes).\n", nrow(sgenes), perc))
  
  return(sgenes)
}