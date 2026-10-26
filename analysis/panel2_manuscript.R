## ====================================================================================================================
# PANEL FIGURE — manuscrito
#
# Idea: cada figura se redibuja con EXACTAMENTE el mismo diseño que tenía (mismas proporciones entre texto,
# puntos, líneas y lienzo), pero reducida con un factor de escala elegido para que el texto de los ejes
# mida `target_pt` en TODAS las figuras. Así, colocando cada archivo en el editor al 100 % (sin reescalar),
# todas las figuras tienen tamaños de texto coherentes y se ven igual que las originales.
#
#   scale_factor = target_pt / tamaño original del texto de ejes de esa figura
#   lienzo final = lienzo original * scale_factor
#
# Ejecutar desde la raíz del repo de EBEx.
## ====================================================================================================================

devtools::load_all(".")   # o library(EBEx)
library(ggplot2)

## --------------------------------------------------------------------------------------------------------------------
# 0. PARÁMETROS COMUNES

# Datos de origen de la figura (ver test/manuscript_data/README.md)
data_dir   <- "test/manuscript_data"
out_dir    <- file.path("test", "panels_manuscript", "panel2")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

target_pt  <- 6      # escala del diseño (lienzo, puntos, líneas): texto de ejes a 6 pt
text_scale <- 7 / 6  # aumento SOLO del texto: ejes 6 -> 7 pt (el resto del plot no cambia)
dpi        <- 500

# Tamaños de texto de referencia para TODO el panel = los del dotplot lado a lado (figura D final).
# Definidos aquí arriba para que cada sección (A, B, E...) se pueda ejecutar por separado.
comb_text_scale <- 1                              # texto de ejes a 6 pt
s_ref  <- (target_pt / 14) * comb_text_scale
fs_ref <- list(
  title        = 18 * s_ref,   # ~7.7 pt bold
  axis_text    = 14 * s_ref,   # 6 pt
  axis_title   = 15 * s_ref,   # ~6.4 pt
  legend_title = 12 * s_ref,   # ~5.1 pt
  legend_text  = 11 * s_ref,   # ~4.7 pt
  key_lines    = 1.2 * s_ref   # alto de las cajitas de leyenda (en "lines")
)

## --------------------------------------------------------------------------------------------------------------------
# 1. HELPERS

#' Guarda un ggplot en PDF vectorial + PNG, al tamaño final (pulgadas)
save_gg <- function(p, name, width_in, height_in, dpi = 500) {
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p,
         width = width_in, height = height_in, units = "in", device = cairo_pdf)
  ggsave(file.path(out_dir, paste0(name, ".png")), p,
         width = width_in, height = height_in, units = "in", dpi = dpi, device = "png", bg = "white")
  message(sprintf("Saved: %-28s %5.1f x %5.1f mm", name, width_in * 25.4, height_in * 25.4))
}

#' Guarda un dibujo de grid (función sin argumentos) en PDF + PNG, tamaño en mm
save_grid <- function(draw_fun, name, size_mm, dpi = 500) {
  cairo_pdf(file.path(out_dir, paste0(name, ".pdf")),
            width = size_mm[1] / 25.4, height = size_mm[2] / 25.4)
  draw_fun()
  dev.off()
  png(file.path(out_dir, paste0(name, ".png")),
      width = size_mm[1], height = size_mm[2], units = "mm", res = dpi, type = "cairo", bg = "white")
  draw_fun()
  dev.off()
  message(sprintf("Saved: %-28s %5.1f x %5.1f mm", name, size_mm[1], size_mm[2]))
}

## --------------------------------------------------------------------------------------------------------------------
# 2. FUNCIÓN: tu plot_enrichment_dotplot original + `scale_factor` + `drop_go_id`
#    Con scale_factor = 1 y drop_go_id = FALSE el resultado es idéntico al original.
#    (Cuando funcione, puedes sustituir con ella la versión de EBEx.)

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
#' @param scale_factor Numeric. Uniform scaling of text, points, lines and margins (1 = original look). Use together with a canvas of 10 x 8 in * scale_factor to keep the original proportions.
#' @param drop_go_id Logical. Remove the "(GO:0000000)" suffix from term names.
#'
#' @return A ggplot object representing the enrichment dotplot.
#'
#' @export
plot_enrichment_dotplot <- function(enrichment_df, title, n_input_genes, top_n = 15, wrap_width = 30,
                                    color_limits = NULL, size_limits = NULL, aspect_ratio = 1, output_file = NULL,
                                    scale_factor = 1, text_scale = 1, drop_go_id = FALSE,
                                    label_lineheight = 0.9,     # 0.9 = valor por defecto de ggplot
                                    x_n_breaks = NULL) {        # nº aprox. de marcas en el eje x (NULL = automático)

  if (nrow(enrichment_df) == 0) {
    message("No data to plot for: ", title)
    return(NULL)
  }
  s  <- scale_factor
  ts <- scale_factor * text_scale   # solo para el texto

  # Prepare data for plotting
  plot_data <- enrichment_df %>%
    dplyr::arrange(.data$Adjusted.P.value) %>%
    head(top_n) %>%
    dplyr::mutate(
      # Extraer el numerador de "5/100" en la columna Overlap
      Count = as.numeric(stringr::str_extract(.data$Overlap, "^[0-9]+")),
      GeneRatio = .data$Count / n_input_genes,
      minusLog10FDR = -log10(.data$Adjusted.P.value),
      Term = if (drop_go_id) stringr::str_remove(.data$Term, "\\s*\\(GO:[0-9]+\\)") else .data$Term,
      Term = stringr::str_wrap(.data$Term, width = wrap_width)
    )
  # Create the dotplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$GeneRatio, y = stats::reorder(.data$Term, .data$GeneRatio))) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Count, color = .data$minusLog10FDR)) +
    # Límites opcionales: si se pasan, la escala es común entre plots (misma leyenda)
    ggplot2::scale_color_gradient(low = "blue", high = "red", limits = color_limits) +
    # Rango de tamaños por defecto de ggplot (1-6 mm), escalado
    ggplot2::scale_size_continuous(range = c(1, 6) * s, limits = size_limits) +
    (if (!is.null(x_n_breaks)) ggplot2::scale_x_continuous(n.breaks = x_n_breaks) else NULL) +
    # theme_bw() por defecto (base_size 11, líneas 0.5), escalado: márgenes, ticks y bordes se reducen igual
    ggplot2::theme_bw(base_size = 11 * s, base_line_size = 0.5 * s, base_rect_size = 0.5 * s) +
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
                                         size = 18 * ts,
                                         hjust = 0.5),
      # Eje y
      axis.text.y = ggplot2::element_text(size = 14 * ts, lineheight = label_lineheight),
      # Eje x
      axis.text.x = ggplot2::element_text(size = 14 * ts),
      axis.title.x = ggplot2::element_text(size = 15 * ts),
      # Leyend
      legend.title = ggplot2::element_text(size = 12 * ts),
      legend.text = ggplot2::element_text(size = 11 * ts),
      legend.key.size = ggplot2::unit(1.2 * s, "lines"),
      # Proporción del grid
      aspect.ratio = aspect_ratio)

  # Save the plot if output_file is provided
  if (!is.null(output_file)) {
    if (!dir.exists(dirname(output_file))) dir.create(dirname(output_file), recursive = TRUE)
    ggplot2::ggsave(filename = output_file, plot = p, device = "pdf", width = 10 * s, height = 8 * s)
    print_message("Enrichment dotplot saved to:", output_file)
  }

  return(invisible(p))
}

## --------------------------------------------------------------------------------------------------------------------
# D. FUNCTIONAL ENRICHMENT DOTPLOTS

D1 <- readRDS(file.path(data_dir, "fig3_D1_reactome.rds"))
D2 <- readRDS(file.path(data_dir, "fig3_D2_kegg.rds"))
D3 <- readRDS(file.path(data_dir, "fig3_D3_go.rds"))

# Diseño original: texto de ejes 14 pt; aspect_ratio y lienzo (in) de cada plot tal como los usabas
s_D <- target_pt / 14

D_specs <- list(
  D1_enrichment_reactome = list(obj = D1, title = "REACTOME: Extended Candidate Genes",
                                aspect = 3,   canvas = c(10, 9.7), drop_go_id = FALSE, wrap = 30),
  D2_enrichment_kegg     = list(obj = D2, title = "KEGG: Extended Candidate Genes",
                                aspect = 2.8, canvas = c(11, 9.7), drop_go_id = FALSE, wrap = 30),
  D3_enrichment_go       = list(obj = D3, title = "GO: Extended Candidate Genes",
                                aspect = 4,   canvas = c(10, 9.7), drop_go_id = TRUE,  wrap = 30)
)

top_n_D    <- 15
D_height_mm <- 105   # más alto para que las etiquetas de 3 líneas (GO) no se solapen con 15 términos

# Escalas comunes: mismos límites de color (-log10 FDR) y tamaño (Count) en los tres plots,
# para que la leyenda sea idéntica y baste con dejar UNA en el editor
D_top <- lapply(D_specs, function(d) {
  d$obj$data %>% dplyr::arrange(.data$Adjusted.P.value) %>% head(top_n_D)
})
all_fdr   <- unlist(lapply(D_top, function(x) -log10(x$Adjusted.P.value)))
all_count <- unlist(lapply(D_top, function(x) as.numeric(stringr::str_extract(x$Overlap, "^[0-9]+"))))

D_color_limits <- c(floor(min(all_fdr) * 2) / 2, ceiling(max(all_fdr) * 2) / 2)   # redondeado a 0.5
D_size_limits  <- range(all_count)
message("Escala común D -> -log10(FDR): ", paste(D_color_limits, collapse = "-"),
        " | Count: ", paste(D_size_limits, collapse = "-"))

for (nm in names(D_specs)) {
  d <- D_specs[[nm]]
  p <- plot_enrichment_dotplot(d$obj$data, d$title, d$obj$n_input_genes,
                               aspect_ratio = d$aspect, scale_factor = s_D, text_scale = text_scale,
                               drop_go_id = d$drop_go_id, label_lineheight = 0.75,
                               wrap_width = d$wrap, x_n_breaks = 3, top_n = top_n_D,
                               color_limits = D_color_limits, size_limits = D_size_limits)
  # Alto fijo = alto disponible en el panel; el panel se ajusta manteniendo aspect_ratio
  if (!is.null(p)) save_gg(p, nm, d$canvas[1] * s_D, D_height_mm / 25.4, dpi)
}

## --------------------------------------------------------------------------------------------------------------------
# D (ALTERNATIVA). UN SOLO DOTPLOT CON LAS TRES BASES DE DATOS
# Mismo estilo que plot_enrichment_dotplot (theme_bw, colores, tamaños, escalado), pero con un facet por base de
# datos apilado en vertical: un único eje x, una única leyenda y escalas comunes de forma natural.

#' Plot combined enrichment dotplot (one facet per database)
#'
#' @param enr_list Named list of Enrichr data.frames; names are used as facet labels (in this order).
#' @param n_input_genes Integer. Total number of genes used in the enrichment call.
#' @param top_n Integer. Number of top terms per database.
#' @param wrap_width Integer. Width for wrapping long term names.
#' @param scale_factor,text_scale,label_lineheight,x_n_breaks,drop_go_id As in plot_enrichment_dotplot().
#' @param title Character. Optional overall title.
#'
#' @return A ggplot object.
plot_enrichment_dotplot_combined <- function(enr_list, n_input_genes, top_n = 15, wrap_width = 40,
                                             scale_factor = 1, text_scale = 1, label_lineheight = 0.9,
                                             x_n_breaks = NULL, x_breaks = NULL, drop_go_id = TRUE,
                                             title = NULL, show_legend = TRUE,
                                             layout = c("rows", "cols"),     # rows = apilados; cols = lado a lado
                                             strip_box = TRUE,               # FALSE = títulos de panel sin recuadro gris
                                             title_hjust = 0.5) {            # 0.5 centrado, 1 a la derecha
  layout <- match.arg(layout)
  s  <- scale_factor
  ts <- scale_factor * text_scale

  plot_data <- dplyr::bind_rows(lapply(names(enr_list), function(db) {
    enr_list[[db]] %>%
      dplyr::arrange(.data$Adjusted.P.value) %>%
      head(top_n) %>%
      dplyr::mutate(Database = db)
  })) %>%
    dplyr::mutate(
      Count         = as.numeric(stringr::str_extract(.data$Overlap, "^[0-9]+")),
      GeneRatio     = .data$Count / n_input_genes,
      minusLog10FDR = -log10(.data$Adjusted.P.value),
      Term          = if (drop_go_id) stringr::str_remove(.data$Term, "\\s*\\(GO:[0-9]+\\)") else .data$Term,
      Term          = stringr::str_wrap(.data$Term, width = wrap_width),
      Database      = factor(.data$Database, levels = names(enr_list)),
      # clave única término+base (un mismo término podría aparecer en dos bases)
      Key           = paste(.data$Term, .data$Database, sep = "___")
    ) %>%
    dplyr::arrange(.data$Database, .data$GeneRatio) %>%
    dplyr::mutate(Key = factor(.data$Key, levels = unique(.data$Key)))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$GeneRatio, y = .data$Key)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Count, color = .data$minusLog10FDR)) +
    ggplot2::scale_color_gradient(low = "blue", high = "red") +
    ggplot2::scale_size_continuous(range = c(1, 6) * s) +
    ggplot2::scale_y_discrete(labels = function(x) sub("___.*$", "", x)) +
    (if (!is.null(x_breaks)) ggplot2::scale_x_continuous(breaks = x_breaks)
     else if (!is.null(x_n_breaks)) ggplot2::scale_x_continuous(n.breaks = x_n_breaks) else NULL) +
    (if (layout == "rows")
      ggplot2::facet_grid(rows = ggplot2::vars(.data$Database), scales = "free_y", space = "free_y")
     else  # lado a lado: cada panel con sus etiquetas, eje x y leyenda comunes
       ggplot2::facet_wrap(ggplot2::vars(.data$Database), nrow = 1, scales = "free_y")) +
    ggplot2::theme_bw(base_size = 11 * s, base_line_size = 0.5 * s, base_rect_size = 0.5 * s) +
    ggplot2::labs(title = title, x = "Gene Ratio", y = NULL,
                  color = "-log10(FDR)", size = "Count") +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 18 * ts, hjust = title_hjust),
      plot.title.position = "plot",   # alinea el título respecto a toda la figura, no solo a los paneles
      strip.background    = if (strip_box) NULL else ggplot2::element_blank(),   # NULL = recuadro de theme_bw
      axis.text.y     = ggplot2::element_text(size = 14 * ts, lineheight = label_lineheight),
      axis.text.x     = ggplot2::element_text(size = 14 * ts),
      axis.title.x    = ggplot2::element_text(size = 15 * ts),
      strip.text      = ggplot2::element_text(size = 14 * ts, face = "bold"),
      legend.title    = ggplot2::element_text(size = 12 * ts),
      legend.text     = ggplot2::element_text(size = 11 * ts),
      legend.key.size = ggplot2::unit(1.2 * s, "lines"),
      legend.position = if (show_legend) "right" else "none")
}

#' Extrae la leyenda de un ggplot como grob (para guardarla aparte)
extract_legend <- function(p) {
  g     <- ggplot2::ggplotGrob(p)
  boxes <- g$grobs[grepl("guide-box", g$layout$name)]
  boxes[!vapply(boxes, inherits, logical(1), "zeroGrob")][[1]]
}

comb_width_mm   <- 90     # tamaño máximo en el panel (sin leyenda: va en un archivo aparte)
comb_height_mm  <- 175

comb_args <- list(
  enr_list         = list("Reactome" = D1$data, "KEGG" = D2$data, "GO BP" = D3$data),
  n_input_genes    = D1$n_input_genes,
  top_n            = 15,
  wrap_width       = 48,      # más ancho: menos términos en 2 líneas (ya no hay leyenda robando ancho)
  scale_factor     = s_D,
  text_scale       = comb_text_scale,
  label_lineheight = 0.75,
  x_breaks         = seq(0, 1, by = 0.1),   # marcas cada 0.1: sin números solapados
  title            = NULL                   # el título va en el editor / pie de figura: más alto para los términos
)

# Plot SIN leyenda: todo el ancho para etiquetas y panel
p_D_comb <- do.call(plot_enrichment_dotplot_combined, c(comb_args, show_legend = FALSE))
save_gg(p_D_comb, "D_enrichment_combined", comb_width_mm / 25.4, comb_height_mm / 25.4, dpi)

# Leyenda en archivo aparte (mismas fuentes): colócala en el editor en un hueco libre
lgd_D <- extract_legend(do.call(plot_enrichment_dotplot_combined, c(comb_args, show_legend = TRUE)))
lgd_size_mm <- c(grid::convertWidth(sum(lgd_D$widths), "mm", valueOnly = TRUE),
                 grid::convertHeight(sum(lgd_D$heights), "mm", valueOnly = TRUE)) + 2
save_grid(function() { grid::grid.newpage(); grid::grid.draw(lgd_D) },
          "D_enrichment_combined_legend", lgd_size_mm)

## --------------------------------------------------------------------------------------------------------------------
# D (ALTERNATIVA 2). TRES PANELES LADO A LADO, ESCALAS COMUNES Y UNA ÚNICA LEYENDA — 180 x 80 mm

side_width_mm  <- 180
side_height_mm <- 85

p_D_side <- plot_enrichment_dotplot_combined(
  enr_list         = list("Reactome" = D1$data, "KEGG" = D2$data, "GO BP" = D3$data),
  n_input_genes    = D1$n_input_genes,
  top_n            = 15,
  wrap_width       = 30,                    # como en tus dotplots originales
  scale_factor     = s_D,
  text_scale       = comb_text_scale,       # texto de ejes a 6 pt, igual que la combinada vertical
  label_lineheight = 0.75,
  x_breaks         = seq(0, 1, by = 0.1),
  title            = "Functional Enrichment Extended genes",
  title_hjust      = 0,                     # alineado a la izquierda
  strip_box        = FALSE,                 # títulos Reactome / KEGG / GO BP sin recuadro gris
  show_legend      = TRUE,                  # una sola leyenda, a la derecha (muévela en el editor si quieres)
  layout           = "cols"
)

# El título se AÑADE encima: se mide su alto y se suma al lienzo para que los paneles mantengan el tamaño
g_side         <- ggplot2::ggplotGrob(p_D_side)
title_row      <- g_side$layout$t[g_side$layout$name == "title"]
title_height_mm <- grid::convertHeight(g_side$heights[title_row], "mm", valueOnly = TRUE)

save_gg(p_D_side, "D_enrichment_side_by_side",
        side_width_mm / 25.4, (side_height_mm + title_height_mm) / 25.4, dpi)

## ====================================================================================================================
# TAMAÑOS DE TEXTO DE REFERENCIA (definidos en la sección 0, fs_ref)

message(sprintf("Texto de referencia (pt): título %.1f | ejes %.1f | título ejes %.1f | leyenda %.1f / %.1f",
                fs_ref$title, fs_ref$axis_text, fs_ref$axis_title, fs_ref$legend_title, fs_ref$legend_text))

## --------------------------------------------------------------------------------------------------------------------
# A. SHAP MAGNITUDE VS OCCURRENCE
# Tu función original + 3 argumentos opcionales (con los valores por defecto el plot es idéntico al original):
#   - colors:       paleta (por defecto gene_status_colors de EBEx)
#   - scale_factor: escala puntos, líneas y márgenes (1 = original)
#   - font_sizes:   tamaños de texto explícitos (NULL = los de theme_minimal(base_size = 20))

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
#' @param colors Named character vector of colors for gene status.
#' @param scale_factor Numeric. Uniform scaling of points, lines and margins (1 = original look).
#' @param font_sizes Optional list (title, axis_text, axis_title, legend_title, legend_text, key_lines) of text sizes in pt.
#'
#' @return A ggplot object with marginal density plots.
#'
#' @export
plot_shap_occurrence <- function(data, y_threshold = NULL, title = "SHAP Magnitude vs Occurrence",
                                 colors = gene_status_colors, scale_factor = 1, font_sizes = NULL) {
  s <- scale_factor

  # 1. Create the base plot
  p <- ggplot2::ggplot(data, ggplot2::aes(x = .data$occurrence_count,
                                          y = .data$aggregated_max,
                                          color = .data$status)) +
    ggplot2::geom_point(size = 2 * s, alpha = 0.7) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::labs(
      x = "Occurrence Count",
      y = "Max normSHAP Value",
      color = "Gene Status",
      title = title
    ) +
    # base_size escala también márgenes y grosor de líneas (base_line_size = base_size / 22)
    ggplot2::theme_minimal(base_size = 20 * s) +
    ggplot2::theme(
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )

  # Tamaños de texto explícitos (mismos que el resto del panel)
  if (!is.null(font_sizes)) {
    p <- p + ggplot2::theme(
      plot.title          = ggplot2::element_text(face = "bold", size = font_sizes$title, hjust = 0),
      plot.title.position = "plot",
      axis.text           = ggplot2::element_text(size = font_sizes$axis_text),
      axis.title          = ggplot2::element_text(size = font_sizes$axis_title),
      legend.title        = ggplot2::element_text(size = font_sizes$legend_title),
      legend.text         = ggplot2::element_text(size = font_sizes$legend_text),
      legend.key.size     = ggplot2::unit(font_sizes$key_lines, "lines")
    )
  }

  # 2. Add threshold line if provided
  if (!is.null(y_threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = y_threshold,
                                 linetype = "solid",
                                 color = "gray70",
                                 linewidth = 0.6 * s)
  }

  # 3. Add marginal density plots
  # Note: ggMarginal returns an object that can be printed or saved
  p_final <- ggExtra::ggMarginal(p, type = "density", fill = "gray70", col = "gray70", alpha = 0.6)

  return(p_final)
}

A <- readRDS(file.path(data_dir, "fig3_A_shap_occurrence.rds"))

A_size_mm <- c(95, 85)                      # hueco de A en el panel
s_A       <- A_size_mm[1] / (15 * 25.4)     # original: 15 x 12 in -> puntos y líneas en la misma proporción

p_A <- plot_shap_occurrence(
  data         = A$data,
  y_threshold  = A$y_threshold,
  title        = "Candidate Genes: Importance vs Robustness",
  colors       = A$colors,
  scale_factor = s_A,
  font_sizes   = fs_ref
)
save_gg(p_A, "A_shap_occurrence", A_size_mm[1] / 25.4, A_size_mm[2] / 25.4, dpi)

## --------------------------------------------------------------------------------------------------------------------
# B. UPSET (ORIGINS)
# Se reutilizan los argumentos originales guardados en el .rds (upset_args) y solo se cambian los tamaños:
#   - text.scale: UpSetR multiplica unos tamaños base fijos (títulos de ejes 8.3 pt, números de ejes y nombres de
#     sets 7 pt, números encima de las barras 2.9 mm). Se calcula el multiplicador para que coincidan con fs_ref.
#   - point.size / line.size: se reducen en la misma proporción que el texto respecto al original (text.scale = 2).

library(UpSetR)

B <- readRDS(file.path(data_dir, "fig3_B_upset_origins.rds"))
message("---- B contiene: ", paste(names(B), collapse = ", "))

# Matriz de presencia/ausencia (ya renombrada si el .rds la trae; si no, se construye desde la lista de genes)
upset_data_B <- if (!is.null(B$upset_data)) B$upset_data else
  if (!is.null(B$data))       B$data       else fromList(B$gene_list)
upset_data_B <- as.data.frame(upset_data_B)

# Argumentos originales (sin los datos); si no se guardaron, los de tu llamada original
args_B <- if (!is.null(B$upset_args)) B$upset_args else list(
  sets        = colnames(upset_data_B),
  nintersects = NA,
  keep.order  = TRUE,
  order.by    = "freq",
  mb.ratio    = c(0.5, 0.5)
)
args_B$data <- NULL
if (is.null(args_B$sets)) args_B$sets <- if (!is.null(B$set_names)) B$set_names else colnames(upset_data_B)

# Tamaños de texto = fs_ref
pt_to_mm <- 1 / ggplot2::.pt   # 1 pt en mm
text_scale_B <- c(
  fs_ref$axis_title / 8.3,                 # título eje "Intersection Size"
  fs_ref$axis_text  / 7,                   # números del eje de intersecciones
  fs_ref$axis_title / 8.3,                 # título eje "Set Size"
  fs_ref$axis_text  / 7,                   # números del eje de set size
  fs_ref$axis_text  / 7,                   # nombres de los sets
  (fs_ref$axis_text * pt_to_mm) / 2.9      # números encima de las barras
)

# Puntos y líneas: misma reducción que el texto respecto a tu original (text.scale = 2 -> números a 14 pt)
s_B <- fs_ref$axis_text / (7 * 2)
args_B$point.size <- (if (!is.null(args_B$point.size)) args_B$point.size else 2.2) * s_B
args_B$line.size  <- (if (!is.null(args_B$line.size))  args_B$line.size  else 0.7) * s_B
args_B$text.scale <- text_scale_B

B_size_mm <- c(90, 60)   # alto máximo 60 mm
# Reparto del alto entre barras de intersección (1º) y matriz de sets (2º)
args_B$mb.ratio <- c(0.62, 0.38)   # filas de la matriz más juntas: más alto para las barras

draw_upset_B <- function() {
  u <- do.call(upset, c(list(data = upset_data_B), args_B))
  # UpSetR dibuja las barras con un margen inferior NEGATIVO (-0.65 líneas) para pegarlas a la matriz.
  # En una figura pequeña ese solape tapa la parte baja de las barras (las de 1, 2, 5 desaparecen).
  # Se pone a 0: las barras arrancan justo encima de la matriz y se ven todas.
  # upset() guarda las barras ya convertidas a gtable: el margen inferior es la última fila de alturas.
  if (inherits(u$Main_bar, "gtable")) {
    n <- length(u$Main_bar$heights)
    u$Main_bar$heights[n] <- grid::unit(0, "pt")
  } else if (inherits(u$Main_bar, "ggplot")) {
    u$Main_bar <- u$Main_bar + ggplot2::theme(plot.margin = grid::unit(c(0.5, 0.5, 0, 0.5), "lines"))
  }
  print(u)
}
save_grid(draw_upset_B, "B_upset_origins", B_size_mm)

## --------------------------------------------------------------------------------------------------------------------
# E. HEATMAP DE SUBGRUPOS — 180 x 44 mm (+ leyendas en archivo aparte)
# Tu Heatmap + HeatmapAnnotation originales, con:
#   - textos = fs_ref (nombres de anotaciones 6 pt, leyendas 5.1 / 4.7 pt, números de cluster 6.4 pt)
#   - pistas de anotación algo más gruesas (anno_track_mm)
#   - leyendas FUERA del heatmap, en un archivo aparte: en 44 mm no caben junto con dendrograma + 8 pistas
#     (se montaban encima de las anotaciones y empujaban el dendrograma fuera del lienzo)
#   - sin nombres ni dendrograma de filas: el alto disponible es para el dendrograma de columnas y las pistas

library(ComplexHeatmap)
library(circlize)

E <- readRDS(file.path(data_dir, "fig3_E_heatmap.rds"))
message("---- E contiene: ", paste(names(E), collapse = ", "))
message("---- E$colors contiene: ", paste(names(E$colors), collapse = ", "))

# Paleta por nombre, con error claro si el .rds usa otro nombre
get_col <- function(nm) {
  hit <- intersect(nm, names(E$colors))
  if (length(hit) == 0) stop("No encuentro ", paste(nm, collapse = " / "), " en E$colors. Disponibles: ",
                             paste(names(E$colors), collapse = ", "), call. = FALSE)
  E$colors[[hit[1]]]
}

E_size_mm      <- c(180, 44)
E_dend_mm      <- 6        # alto del dendrograma de columnas
anno_track_mm  <- 2.8      # alto de cada pista de anotación
anno_na        <- "grey80"
grid_mm        <- 2.2      # cajitas de leyenda: mismo tamaño que en las leyendas de los ggplot (D)

ann <- E$annotation_data

# Colores de cada anotación (nombre que verás en la figura -> paleta)
E_anno_cols <- list(
  'Disease condition' = get_col("dis_condition"),
  'Sex'               = get_col("sex"),
  'Smoker'            = get_col("smoker"),
  'Age'               = get_col("col_age"),
  'Emphysema'         = get_col("col_emphysema"),
  'predicted DLCO'    = get_col("col_pred_dlco"),
  'predicted FEV1'    = get_col(c("col_pred_fev_postbd", "col_pred_fev1")),
  'predicted FVC'     = get_col(c("col_pred_fvc_postbd", "col_pred_fvc"))
)

ha_E <- HeatmapAnnotation(
  'Disease condition' = ann$dis_condition,
  'Sex'               = ann$sex,
  'Smoker'            = ann$smoker,
  'Age'               = ann$age,
  'Emphysema'         = ann$emphysema,
  'predicted DLCO'    = ann$pred_dlco,
  'predicted FEV1'    = ann$pred_fev_postbd,
  'predicted FVC'     = ann$pred_fvc_postbd,
  col                  = E_anno_cols,
  show_legend          = FALSE,                     # leyendas en archivo aparte
  annotation_name_side = "left",
  simple_anno_size     = unit(anno_track_mm, "mm"),
  gap                  = unit(0.4, "mm"),
  na_col               = anno_na,
  annotation_name_gp   = gpar(fontsize = fs_ref$axis_text)
)

p_E <- Heatmap(
  E$heatmap_matrix,
  name                 = paste0("Scaled ", E$heatmap_data, " data"),
  na_col               = "white",
  col                  = get_col("default_colors"),
  cluster_columns      = E$hclust,
  column_split         = E$column_split,
  row_split            = E$row_split,
  column_gap           = unit(0.5, "mm"),
  row_gap              = unit(0.3, "mm"),
  show_row_names       = FALSE,
  show_column_names    = FALSE,
  show_row_dend        = FALSE,
  row_title            = NULL,
  column_title_gp      = gpar(fontsize = fs_ref$axis_title),   # números de cluster 1-5
  column_dend_height   = unit(E_dend_mm, "mm"),
  top_annotation       = ha_E,
  show_heatmap_legend  = FALSE                     # leyenda en archivo aparte
)

draw_heatmap_E <- function() {
  draw(p_E, padding = unit(c(1, 1, 1, 1), "mm"))
}
save_grid(draw_heatmap_E, paste0("E_heatmap_", E$heatmap_data), E_size_mm)

## ---- Leyendas de E (archivo aparte, mismas fuentes; recolócalas en el editor)

lgd_title_gp  <- gpar(fontsize = fs_ref$legend_title, fontface = "bold")
lgd_labels_gp <- gpar(fontsize = fs_ref$legend_text)

make_lgd <- function(title, col) {
  if (is.function(col)) {
    # continua: barra horizontal
    Legend(col_fun = col, title = title, direction = "horizontal",
           legend_width = unit(15, "mm"), grid_height = unit(grid_mm, "mm"),
           title_position = "topleft", title_gp = lgd_title_gp, labels_gp = lgd_labels_gp)
  } else {
    # discreta: niveles en una fila
    Legend(labels = names(col), legend_gp = gpar(fill = unname(col)), title = title, nrow = 1,
           grid_height = unit(grid_mm, "mm"), grid_width = unit(grid_mm, "mm"),
           title_gp = lgd_title_gp, labels_gp = lgd_labels_gp)
  }
}

E_legends <- c(
  mapply(make_lgd, names(E_anno_cols), E_anno_cols, SIMPLIFY = FALSE),
  list(make_lgd(paste0("Scaled ", E$heatmap_data, " data"), get_col("default_colors")),
       Legend(labels = "Missing", title = "NA", legend_gp = gpar(fill = anno_na),
              grid_height = unit(grid_mm, "mm"), grid_width = unit(grid_mm, "mm"),
              title_gp = lgd_title_gp, labels_gp = lgd_labels_gp))
)

E_lgd_pack <- packLegend(list = unname(E_legends), direction = "horizontal",
                         max_width = unit(E_size_mm[1] - 4, "mm"),   # si no caben, pasan a una 2ª fila
                         column_gap = unit(4, "mm"), row_gap = unit(2, "mm"))

E_lgd_size_mm <- c(E_size_mm[1],
                   grid::convertHeight(ComplexHeatmap:::height(E_lgd_pack), "mm", valueOnly = TRUE) + 2)

save_grid(function() {
  grid::grid.newpage()
  draw(E_lgd_pack, x = unit(2, "mm"), just = "left")
}, paste0("E_heatmap_", E$heatmap_data, "_legend"), E_lgd_size_mm)

## --------------------------------------------------------------------------------------------------------------------
# C. TABLA DE GENES EXTENDIDOS ENRIQUECIDOS POR VARIABLE CLÍNICA
# Se dibuja en R (gridExtra::tableGrob) con las mismas fuentes que el resto del panel y se guarda en PDF + SVG
# al tamaño justo de la tabla: en Inkscape se importa al 100 % y el texto sigue siendo editable.
# La tabla (variable, nº de genes, genes; "\n" = salto de línea en la celda) se carga de test/manuscript_data.

library(gridExtra)

tab_C <- readRDS(file.path(data_dir, "fig3_C_table.rds"))

lwd_C <- 0.5   # grosor de las líneas de la tabla (pt), como los bordes de los paneles de D

theme_C <- ttheme_default(
  core = list(
    fg_params = list(fontsize = fs_ref$axis_text,
                     # Variable y nº en redonda, genes en cursiva
                     fontface = matrix(rep(c("plain", "plain", "italic"), each = nrow(tab_C)), ncol = 3),
                     lineheight = 0.9),
    bg_params = list(fill = "white", col = "black", lwd = lwd_C),
    padding   = unit(c(3, 1.6), "mm")
  ),
  colhead = list(
    fg_params = list(fontsize = fs_ref$axis_title, fontface = "bold", lineheight = 0.9),
    bg_params = list(fill = "white", col = "black", lwd = lwd_C),
    padding   = unit(c(3, 1.6), "mm")
  )
)

g_C <- tableGrob(tab_C, rows = NULL, theme = theme_C)

# Tamaño fijo 85 x 30 mm: anchos de columna y altos de fila repartidos para llenar exactamente ese hueco
C_size_mm  <- c(85, 30)
C_col_mm   <- c(22, 16, 47)                           # Variable | nº genes | genes (suman 85)
# Alto de cada fila proporcional a sus líneas de texto + 1 (relleno):
# cabecera (2 líneas: "# enriched / genes"), emphysema, DLCO (2 líneas), FEV1
C_row_lines <- c(2, 1, 2, 1)
C_row_mm    <- C_size_mm[2] * (C_row_lines + 1) / sum(C_row_lines + 1)

# La tabla se dibuja 1 mm más pequeña que el lienzo y centrada: así los bordes exteriores
# (izquierdo, derecho, arriba y abajo) no quedan cortados por el borde del archivo
C_inset_mm  <- 0.5
g_C$widths  <- unit(C_col_mm * (C_size_mm[1] - 2 * C_inset_mm) / C_size_mm[1], "mm")
g_C$heights <- unit(C_row_mm * (C_size_mm[2] - 2 * C_inset_mm) / C_size_mm[2], "mm")

# Marco exterior explícito alrededor de toda la tabla
g_C <- gtable::gtable_add_grob(g_C, grid::rectGrob(gp = grid::gpar(fill = NA, col = "black", lwd = lwd_C)),
                               t = 1, b = nrow(g_C), l = 1, r = ncol(g_C), name = "outer_border", z = Inf)
message(sprintf("Tabla C: %.0f x %.0f mm | filas: %s mm", C_size_mm[1], C_size_mm[2],
                paste(round(C_row_mm, 1), collapse = ", ")))

draw_table_C <- function() { grid::grid.newpage(); grid::grid.draw(g_C) }
save_grid(draw_table_C, "C_table_enriched_genes", C_size_mm)

# SVG además del PDF: en Inkscape mantiene el texto como texto editable
svg(file.path(out_dir, "C_table_enriched_genes.svg"),
    width = C_size_mm[1] / 25.4, height = C_size_mm[2] / 25.4)
draw_table_C()
dev.off()
