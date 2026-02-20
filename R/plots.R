#' @importFrom dplyr filter select mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_line labs theme_minimal scale_color_manual ggsave scale_color_gradient theme_bw theme element_text
#' @importFrom fmsb radarchart
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices pdf dev.off rainbow
#' @importFrom stats as.formula reorder
#' @importFrom stringr str_extract
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation draw
#' @importFrom circlize colorRamp2
#' @importFrom grid unit gpar
NULL

## --------------------------------------------------------------------------------------------------------------------
# RADAR CHART FUNCTIONS

#' Prepare data for Radar Chart
#'
#' @description 
#' Filters metrics and constructs the specific matrix format required by fmsb,
#' identifying the best performing input list for each classifier.
#'
#' @param results Dataframe containing processed metrics.
#' @param metric Character. The metric to plot (e.g., "normMCC").
#' @param estimate_column Character. The column name containing the values (e.g., "estimate" or "mean").
#' @param min Numeric. Minimum value for the radar scale.
#' @param max Numeric. Maximum value for the radar scale.
#' 
#' @return A list containing the radarchart dataframe, best selection names, and best values.
#' 
#' @export
prepare_data_for_radarchart <- function(results, metric, estimate_column, min = 0, max = 1) {

  # 1. Filter metrics by the desired metric
  results_metric <- results[results$metric == metric, ]
  if (nrow(results_metric) == 0) {
    stop("No data found for metric: ", metric)
  }
  # 2. Ensure values are numeric (avoids 0.0 errors)
  results_metric[[estimate_column]] <- as.numeric(results_metric[[estimate_column]])
  # 3. Identify Classifiers and Gene Lists
  classifiers <- unique(results_metric$classifier)
  mlinput <- unique(results_metric$input_list)
  # 4. Find the best selection for each classifier (Original logic)
  best_selection <- rep("", length(classifiers))
  best_value <- numeric(length(classifiers))
  names(best_selection) <- names(best_value) <- classifiers
  for (clf in classifiers) {
    # Data for this classifier
    sub_clf <- results_metric[results_metric$classifier == clf, ]
    if (nrow(sub_clf) > 0) {
      # Find the maximum value
      max_idx <- which.max(sub_clf[[estimate_column]])
      best_selection[clf] <- as.character(sub_clf$input_list[max_idx])
      best_value[clf] <- sub_clf[[estimate_column]][max_idx]
    }
  }
  # 5. Construct the matrix for fmsb (Manually pivoted to ensure order)
  # Create an empty matrix filled with NAs or zeros
  mat_data <- matrix(0, nrow = length(mlinput), ncol = length(classifiers))
  rownames(mat_data) <- mlinput
  colnames(mat_data) <- classifiers
  for (inp in mlinput) {
    for (clf in classifiers) {
      val <- results_metric[[estimate_column]][results_metric$input_list == inp & results_metric$classifier == clf]
      if (length(val) > 0) {
        mat_data[inp, clf] <- val[1] # Take the first value for safety
      }
    }
  }
  # 6. Add MAX and MIN rows (Requirement of fmsb::radarchart)
  values_radarchart <- rbind(
    max = rep(max, ncol(mat_data)),
    min = rep(min, ncol(mat_data)),
    mat_data
  )

  return(list(
    values_radarchart = as.data.frame(values_radarchart),
    best_selection = unname(best_selection),
    best_value = unname(best_value)
  ))
}

#' Plot Radar Chart
#'
#' @description
#' Generates a PDF radar chart with specific color logic for known and unknown input lists.
#'
#' @param data List. Output from `prepare_data_for_radarchart`.
#' @param metric Character. Metric name for the title.
#' @param plot_name Character. Name for the output file.
#' @param output_dir Character. Directory to save the PDF.
#' @param min Numeric. Minimum axis value.
#' @param max Numeric. Maximum axis value.
#' @param seq Numeric. Step for axis labels.
#'
#' @return Invisibly returns the plot object.
#'
#' @export
plot_radarchart <- function(data, metric, plot_name, output_dir = NULL, min, max, seq) {

  # Save current graphics parameters and ensure they are restored on exit
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par))
  # Extract values and best selections
  values_radarchart <- data$values_radarchart
  best_selection <- data$best_selection
  best_value <- data$best_value
  # Logic to manage colors (Custom vector from utils.R)
  row_names_to_plot <- rownames(values_radarchart)[3:nrow(values_radarchart)]
  # Reference the package color object
  known_names <- intersect(row_names_to_plot, names(genes_list_colors))
  unknown_names <- setdiff(row_names_to_plot, names(genes_list_colors))
  # Generate unique colors for unknowns
  fallback_colors <- character(0)
  if (length(unknown_names) > 0) {
    palette_size <- max(length(unknown_names), 3)
    if (length(unknown_names) <= 8) {
      fallback_colors <- RColorBrewer::brewer.pal(palette_size, "Set2")[1:length(unknown_names)]
    } else {
      fallback_colors <- grDevices::rainbow(length(unknown_names))
    }
    names(fallback_colors) <- unknown_names
    warning(sprintf("Assigning fallback colors for: %s", paste(unknown_names, collapse = ", ")))
  }
  # Merge known and fallback colors in correct order
  all_colors <- c(genes_list_colors[known_names], fallback_colors)
  radar_colors <- all_colors[row_names_to_plot]
  # Ensure the output directory exists
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    file_path <- file.path(output_dir, paste0(plot_name, "_", metric, ".pdf"))
    grDevices::pdf(file = file_path, width = 15, height = 10)
    # Asegurar que el dispositivo se cierra al terminar, incluso si hay error
    on.exit(grDevices::dev.off(), add = TRUE)
  } else {
    warning("No output directory specified. Plot will be displayed in the active device.")
  }
  # Create Plot
  # Plot radar chart with custom colors and settings
  fmsb::radarchart(
    values_radarchart,
    axistype = 5,
    # Line style
    pcol = radar_colors, plwd = 3, plty = 1,
    # Grid and axis
    cglcol = "grey",     cglwd = 2,
    axislabcol = radar_colors[best_selection],
    paxislabels = round(best_value, 3),
    palcex = 1.4, # Tamaño números mejor valor (14pt aprox)
    caxislabels = seq(min, max, seq),
    calcex = 1.3, # Tamaño de etiquetas de eje central
    vlcex = 1.5 # Etiquetas de las esquinas
  )
  # Position legend
  legend(
    x = 1.3, y = 1.3,
    legend = names(radar_colors),
    bty = "n", pch = 20,
    col = radar_colors,
    text.col = "black",
    cex = 1, pt.cex = 1.5
  )
  if (!is.null(output_dir)) {
      print_message(sprintf("Radar chart saved to: %s\n", file.path(output_dir, paste0(plot_name, "_", metric, ".pdf"))))
  }

  return(invisible(NULL))
}

## --------------------------------------------------------------------------------------------------------------------
# TIME ANALYSIS PLOT

#' Plot Computation Time vs Gene List Size
#' 
#' @description
#' This function creates a plot of computation time against the number of genes in the input list for each classifier.
#'
#' @param time_df Dataframe containing columns 'time', 'classifier', and 'input_list'.
#' @param gene_list_counts Named vector with the number of genes for each input list.
#' @param output_dir Character. Directory to save the plot.
#' @param colors Named vector. Optional colors for each classifier.
#' 
#' @return Invisibly returns the plot object.
#'
#' @export
plot_time_vs_genes <- function(time_df, gene_list_counts, output_dir, colors = NULL) {
  
  # Merge time data with gene counts
  plot_df <- time_df %>%
    dplyr::mutate(num_genes = gene_list_counts[.data$input_list])
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$time, y = .data$num_genes, color = .data$classifier)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_line(ggplot2::aes(group = .data$classifier), size = 1) +
    ggplot2::labs(title = "Efficiency: Time vs Number of Genes",
                  x = "Time",
                  y = "Number of Genes") +
    ggplot2::theme_minimal()
  
  if (!is.null(colors)) {
    p <- p + ggplot2::scale_color_manual(values = colors)
  }
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, "time_analysis.pdf"), plot = p, width = 8, height = 6)
}

## --------------------------------------------------------------------------------------------------------------------
# SUMMARY BARPLOT

#' Prepare data for Summary Barplot
#'
#' @description
#' Filters the performance data by a specific metric and calculates
#' mean, SD, and median per classifier.
#'
#' @param data Dataframe. Processed metrics (long format).
#' @param target_metric Character. The metric to summarize (e.g., "normMCC").
#' @param value_col Character. The name of the column containing the values.
#'
#' @return A summarized dataframe ready for plotting.
#'
#' @export
prepare_summary_barplot_data <- function(data, target_metric, value_col) {
  
  # 1. Make sure the values are numeric and filter by the target metric
  summary_df <- data %>%
    dplyr::filter(.data$metric == !!target_metric) %>%
    dplyr::mutate(value_num = as.numeric(!!rlang::sym(value_col))) %>%
    # 2. Group and summarize
    dplyr::group_by(.data$classifier) %>%
    dplyr::summarize(
      mean_value = mean(.data$value_num, na.rm = TRUE),
      sd = stats::sd(.data$value_num, na.rm = TRUE),
      median_value = stats::median(.data$value_num, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # 3. Reorder factor levels by mean (descending)
    dplyr::mutate(
      classifier = stats::reorder(.data$classifier, -.data$mean_value)
    )
  
  return(summary_df)
}

#' Plot Summary Barplot of Classifier Performance
#'
#' @description
#' Generates a barplot with error bars (SD), median points, and SD annotations.
#'
#' @param summary_data Dataframe. Output from `prepare_summary_barplot_data`.
#' @param metric_name Character. Metric name for titles and axes.
#' @param colors Named vector. Colors for each classifier.
#' @param base_size Numeric. Base font size (default 14).
#'
#' @return A ggplot object.
#'
#' @export
plot_summary_barplot <- function(summary_data, metric_name, colors = ml_models_colors, base_size = 14) {
  
  # Adjust size for internal labels (14pt approx = 5mm)
  number_size <- 9 
  p <- ggplot2::ggplot(summary_data, ggplot2::aes(x = .data$classifier, y = .data$mean_value)) +
    # Mean bars
    ggplot2::geom_col(ggplot2::aes(fill = .data$classifier), position = "dodge", show.legend = FALSE) +
    # Error bars (SD)
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = .data$mean_value - .data$sd, ymax = .data$mean_value + .data$sd),
      width = 0.2, color = "grey30"
    ) +
    # Median point
    ggplot2::geom_point(ggplot2::aes(y = .data$median_value), color = "#003049", size = 3) +
    # SD value annotation above the bar
    ggplot2::geom_text(
      ggplot2::aes(y = .data$mean_value + .data$sd, label = round(.data$sd, 3)),
      vjust = -0.8, size = number_size, fontface = "italic"
    ) +
    # Aesthetics and Colors
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = paste("Comparison by", metric_name),
      y = paste("Mean", metric_name),
      x = "Classifier"
    ) +
    # Theme for the manuscript (Title 18, Labels 14)
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 4, hjust = 0.5),
      axis.title = ggplot2::element_text(),
      axis.text.y = ggplot2::element_text(size = base_size + 4),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 0, size = base_size +1),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

## --------------------------------------------------------------------------------------------------------------------
# BARPLOT INPUT LISTS

#' Plot Barplot of input lists
#'
#' @description
#' This function generates a barplot showing the number of genes in each input list used for feature selection. It uses the 
#' package-wide color scheme for consistency.
#'
#' @param feature_selection_data Dataframe with columns 'input_list' and 'intersection_count'.
#' @param output_file Character. Path to save the output PDF.
#'
#' @return Invisibly returns the plot object.
#'
#' @export
barplot_feature_selection <- function(feature_selection_data, output_file = NULL) {
  
  # Validate data
  if (!all(c("input_list", "intersection_count") %in% colnames(feature_selection_data))) {
    stop("El dataframe debe contener las columnas 'input_list' y 'intersection_count'.")
  }
  # 2. Create the barplot
  p <- ggplot2::ggplot(feature_selection_data) +
    ggplot2::aes(
      x = input_list,
      fill = input_list,
      weight = intersection_count
    ) +
    ggplot2::labs(x = "Input selection")+ggplot2::labs(y = "")+
    ggplot2::geom_bar(position = "dodge") +
    ggplot2::geom_text(ggplot2::aes(label = intersection_count, y = intersection_count), 
              position = ggplot2::position_dodge(width = 0.9), 
              vjust = -0.5, size = 6) +  # Add labels on top of bars
    ggplot2::scale_fill_manual(values = genes_list_colors) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      fill = "Feature Selection Approach",  # Lengend title
      title = "Feature Selection Methods Number of Genes" # Title
    ) +
    ggplot2::theme(
      # Eje y
      axis.text.y = ggplot2::element_text(size = 15),
      # Eje x
      axis.text.x = ggplot2::element_blank(),  # Remove x-axis labels
      axis.title.x = ggplot2::element_blank(), # Remove x-axis title
      # Leyend
      legend.title = ggplot2::element_text(size = 18),
      legend.text = ggplot2::element_text(size = 17),
      # Title
      plot.title = ggplot2::element_text(face = "bold", 
                                         size = 20, 
                                         hjust = 0.5),
      )
  # 3. Save the plot
  if (!is.null(output_file)) {
    if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
    ggplot2::ggsave(
      filename = output_file, 
      plot = p, 
      device = "pdf", 
      width = 12, 
      height = 6
    )
    print_message(sprintf("Barplot saved to: %s\n", output_file))
  }

  return(invisible(p))
}

## --------------------------------------------------------------------------------------------------------------------
# FUNCTIONAL ENRICHMENT

#' Plot Functional Enrichment Dotplot
#'
#' @description
#' This function creates a dotplot for functional enrichment results, showing gene ratio, count, and adjusted p-value for the top terms.
#'
#' @param enrichment_df Dataframe from one database of Enrichr results.
#' @param title Character. Title for the plot.
#' @param n_input_genes Integer. Total number of genes used in the enrichment call.
#' @param top_n Integer. Number of top terms to show.
#' @param wrap_width Integer. Width for wrapping long term names.
#' @param color_limits Numeric vector of length 2. Limits for the color scale.
#' @param size_limits Numeric vector of length 2. Limits for the size scale.
#' @param aspect_ratio Numeric. Aspect ratio for the plot.
#' @param output_file Character. Optional path to save the plot (e.g., "enrichment_dotplot.pdf"). If NULL (default), the plot is displayed in the active device.
#'
#' @return A ggplot object representing the enrichment dotplot.
#'
#' @export
plot_enrichment_dotplot <- function(enrichment_df, title, n_input_genes, top_n = 15, wrap_width = 30, 
                                        color_limits = c(1, 5), size_limits = c(1, 30), aspect_ratio = 1, output_file = NULL) {
  
  if (nrow(enrichment_df) == 0) {
    message("No data to plot for: ", title)
    return(NULL)
  }
  # Prepare data for plotting
  plot_data <- enrichment_df %>%
    dplyr::arrange(.data$Adjusted.P.value) %>%
    head(top_n) %>%
    dplyr::mutate(
      # Extraer el numerador de "5/100" en la columna Overlap
      Count = as.numeric(stringr::str_extract(.data$Overlap, "^[0-9]+")),
      GeneRatio = .data$Count / n_input_genes,
      minusLog10FDR = -log10(.data$Adjusted.P.value),
      Term = stringr::str_wrap(.data$Term, width = wrap_width)
    )
  # Create the dotplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$GeneRatio, y = stats::reorder(.data$Term, .data$GeneRatio))) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Count, color = .data$minusLog10FDR)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    # ggplot2::scale_size_continuous(limits = size_limits, range = c(1, 8)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x = "Gene Ratio",
      y = NULL,
      color = "-log10(FDR)",
      size = "Count"
    ) +
    ggplot2::theme(
      # Titulo
      plot.title = ggplot2::element_text(face = "bold", 
                                         size = 18, 
                                         hjust = 0.5),
      # Eje y
      axis.text.y = ggplot2::element_text(size = 14),
      # Eje x
      axis.text.x = ggplot2::element_text(size = 14),
      axis.title.x = ggplot2::element_text(size = 15),
      # Leyend
      legend.title = ggplot2::element_text(size = 12),
      legend.text = ggplot2::element_text(size = 11),
      # Proporción del grid
      aspect.ratio = aspect_ratio)

  # Save the plot if output_file is provided
  if (!is.null(output_file)) {
    if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
    ggplot2::ggsave(filename = output_file, plot=p, device = "pdf", width = 10, height = 8)
    print_message("Enrichment dotplot saved to:", output_file)
  }

  return(invisible(p))
}

## --------------------------------------------------------------------------------------------------------------------
# GTEx HEATMAP

#' Plot GTEx Expression Heatmap
#'
#' @description
#' This function generates a ComplexHeatmap based on median TPM values with GTEx-style scaling.
#' It includes annotations for tissue macro-categories and uses a color scheme inspired by GTEx visualizations.
#'
#' @param tpm_matrix Matrix. Genes as rows, GTEx tissue IDs as columns.
#' @param tissue_macro_map Named vector mapping tissue IDs to macro-categories.
#' @param tissue_fine_map Named vector mapping tissue IDs to pretty names.
#' @param output_path Character. Path to save the PDF.
#' @param title Character. Title for the heatmap.
#'
#' @return A ComplexHeatmap object.
#'
#' @export
plot_gtex_heatmap <- function(tpm_matrix, tissue_macro_map, tissue_fine_map, 
                              title = "GTEx Expression") {

  # 1. Map tissue IDs to pretty names and filter out any tissues without mapping
  colnames(tpm_matrix) <- tissue_fine_map[colnames(tpm_matrix)]
  heatmap_matrix <- tpm_matrix[, !is.na(colnames(tpm_matrix)), drop = FALSE]
  stopifnot(!any(is.na(colnames(heatmap_matrix))))
  # 2. Map tissue IDs to macro-categories and order
  macro_vec <- tissue_macro_map[colnames(heatmap_matrix)]
  ord <- order(macro_vec)
  heatmap_matrix <- heatmap_matrix[, ord, drop = FALSE]
  macro_vec <- macro_vec[ord]
  # 3. Color scale (sqrt to highlight low values)
  # Cut at 2500 as per GTEx standard
  heatmap_matrix[heatmap_matrix > 2500] <- 2500
  col_fun <- circlize::colorRamp2(
    sqrt(c(0, 3, 15, 60, 240, 2500)),
    c("#f7fcb9", "#addd8e", "#41ab5d", "#1d91c0", "#225ea8", "#0c2c84")
  )
  # 4. Top annotation
  top_ha <- ComplexHeatmap::HeatmapAnnotation(
    Tissue = macro_vec,
    col = list(Tissue = gtex_group_colors),
    annotation_name_gp = grid::gpar(fontsize = 8),
    simple_anno_size = grid::unit(4, "mm")
  )
  tpm_breaks <- c(0, 3, 15, 60, 240, 1000, 2500)
  # 5. Create Heatmap
  p <- ComplexHeatmap::Heatmap(
    heatmap_matrix,
    name = "TPM",
    col = function(x) col_fun(sqrt(x)),
    rect_gp = grid::gpar(col = "#ffffff", lwd = 0.4),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_split = macro_vec,
    column_title_rot = 45,
    column_title_gp = grid::gpar(fontsize = 7),
    top_annotation = top_ha,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    heatmap_legend_param = list(
      title = "TPM",
      at = sqrt(tpm_breaks),
      labels = tpm_breaks,
      title_gp = grid::gpar(fontsize = 8)
    )
  )

  return(p)
}

## --------------------------------------------------------------------------------------------------------------------
# DEA CANDIDATE GENES PLOTS

#' Plot DEA Candidate Genes
#'
#' @description 
#' Generates either an MA plot or a Volcano plot for DEA candidate genes, 
#' with points colored by gene status and threshold lines for log fold change and adjusted p-value.
#'
#' @param deg_data Dataframe containing DEA results with columns 'AveExpr', 'logFC', 'adj.P.Val', and 'status'.
#' @param type Character. Type of plot to generate: "MA" or "Vol
#' @param lfc_thresh Numeric. Log fold change threshold for significance (default log2(1.5)).
#' @param p_val_thresh Numeric. Adjusted p-value threshold for significance (default 0.01).
#'
#' @return A ggplot object representing the MA or Volcano plot.
#'
#' @export
plot_dea_candidates <- function(deg_data, type = "MA", lfc_thresh = log2(1.5), p_val_thresh = 0.01) {

  if (type == "MA") {
    p <- ggplot2::ggplot(deg_data, ggplot2::aes(x = .data$AveExpr, y = .data$logFC)) +
      ggplot2::geom_hline(yintercept = c(lfc_thresh, -lfc_thresh), linetype = "dashed")
  } else {
    p <- ggplot2::ggplot(deg_data, ggplot2::aes(x = .data$logFC, y = -log10(.data$adj.P.Val))) +
      ggplot2::geom_vline(xintercept = c(lfc_thresh, -lfc_thresh), linetype = "dashed") +
      ggplot2::geom_hline(yintercept = -log10(p_val_thresh), linetype = "dashed")
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(color = .data$status), alpha = 0.5, size = 2) +
    # ggrepel::geom_text_repel(
    #   data = deg_data %>% dplyr::filter(.data$status == "Extended"),
    #   ggplot2::aes(label = .data$variable_name),
    #   size = 2.5, max.overlaps = Inf
    # ) +
    ggplot2::scale_color_manual(values = gene_status_colors) +
    ggplot2::theme_minimal(base_size = 20) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"))+
    ggplot2::labs(color = "Gene Status", title = paste(type, "Plot: Candidate Genes in DEA"))

  return(p)
}

## --------------------------------------------------------------------------------------------------------------------
# SHAP OCCURRENCE PLOT

#' Plot SHAP Magnitude vs Occurrence
#'
#' @description
#' This function creates a scatter plot of SHAP magnitude (aggregated_max) against occurrence
#' count for each gene, colored by gene status. An optional horizontal line can be added to indicate a threshold for SHAP magnitude.
#' Marginal density plots are included on the axes to show the distribution of occurrence counts and SHAP magnitudes.
#'
#' @param data Dataframe containing columns 'occurrence_count', 'aggregated_max', and 'status'.
#' @param y_threshold Numeric. Optional threshold for SHAP magnitude to be indicated with a horizontal line.
#' @param title Character. Title for the plot.
#'
#' @return A ggplot object with marginal density plots.
#'
#' @export
plot_shap_occurrence <- function(data, y_threshold = NULL, title = "SHAP Magnitude vs Occurrence") {

  # 1. Create the base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$occurrence_count,
                                          y = .data$aggregated_max,
                                          color = .data$status)) +
    ggplot2::geom_point(size = 2, alpha = 0.7) +
    ggplot2::scale_color_manual(values = gene_status_colors) +
    ggplot2::labs(
      x = "Occurrence Count",
      y = "Max normSHAP Value",
      color = "Gene Status",
      title = title
    ) +
    ggplot2::theme_minimal(base_size = 20) +
    ggplot2::theme(
      # # Leyend
      # legend.title = ggplot2::element_text(size = 12),
      # legend.text = ggplot2::element_text(size = 11),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
      )

  # 2. Add threshold line if provided
  if (!is.null(y_threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = y_threshold,
                                 linetype = "solid",
                                 color = "gray70",
                                 linewidth = 0.6)
  }

  # 3. Add marginal density plots
  # Note: ggMarginal returns an object that can be printed or saved
  p_final <- ggExtra::ggMarginal(p, type = "density", fill = "gray70", col = "gray70", alpha = 0.6)

  return(p_final)
}