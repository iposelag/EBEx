#!/usr/bin/env Rscript
######################################
## Batch correction
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(devtools)
devtools::load_all()
devtools::document()
library(EBEx)
library(Biobase)
library(sva); library(ggplot2);library(patchwork)
# Global parameters
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
# Relative paths
raw_data_dir <- "test/raw_data/"
output_dir <- "test/analysis/batch_correction/"

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. DATA PREPARATION (project-specific)
# Load ExpressionSet
cat("Loading ExpressionSet...\n")
eset_path <- file.path(raw_data_dir, "eset.Rda")
load(eset_path)  # Load the ExpressionSet object
# Relevel condition factor and extract expression, metadata, and gene names
eset$dis_condition <- relevel(factor(eset$dis_condition), ref = "CTRL")
expression <- Biobase::exprs(eset)
# Phenotype metadata
phenotypic_ctics <- Biobase::pData(eset)

# Filter and clean phenotypes
pheno_data <- phenotypic_ctics %>%
  dplyr::select(dis_condition, sex, smoker, age, emphysema, pred_dlco, pred_fev_postbd, pred_fvc_postbd, platform_id) %>%
#   dplyr::filter(dis_condition == "COPD") %>%
  dplyr::mutate(
    across(c(sex, smoker, platform_id), as.factor),
    across(c(age, emphysema, pred_dlco, pred_fev_postbd, pred_fvc_postbd), as.numeric),
  ) %>%
  dplyr::mutate(sample_id = rownames(.))

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. PCA and Visualization before batch correction
pca <- perform_pca(t(expression), num_components = 30, # samples x genes
                       save_eig_plot = TRUE, eig_ncp = 12,
                       output_dir = output_dir, plot_filename="pca_eigenvalues_pre_combat.png")
pca_pheno_data <- join_pca_with_metadata(pca$pca_data, pheno_data)

tema_labels <- theme(
  axis.title   = element_text(size = 20),  # títulos de ejes (PC1, PC2...)
  axis.text    = element_text(size = 16),  # números de los ejes
  legend.title = element_text(size = 18),
  legend.text  = element_text(size = 16)
)

pdf(file.path(output_dir, "pca_pre_combat.pdf"), width = 7, height = 6)

for(pc in c("PC2","PC3")){
  for(feature in c("dis_condition","platform_id", "sex", "smoker", "age",
                  "emphysema", "pred_dlco", "pred_fev_postbd", "pred_fvc_postbd")){
    if(is.numeric(pca_pheno_data[[feature]])){
      print(ggplot(
        pca_pheno_data,
        aes(
          x = PC1,
          y = !!sym(pc),
          shape = platform_id,
          color = !!sym(feature) # label points with different colors for each `subgroup`
        )
      ) +
        geom_point() + # Plot individual points to make a scatterplot
        theme_classic(
        ) + # Format as a classic-looking plot with no gridlines
        scale_fill_continuous(type = "gradient",na.value="white")+
        tema_labels)
    } else {
      palette <- get(feature)  # paleta con nombre desde utils.R (dis_condition/platform_id/sex/smoker)
      print(ggplot(
        pca_pheno_data,
        aes(
          x = PC1,
          y = !!sym(pc),
          shape = dis_condition,
          color = !!sym(feature) # label points with different colors for each `subgroup`
        )
      ) +
        geom_point() + # Plot individual points to make a scatterplot
        theme_classic(

        ) + # Format as a classic-looking plot with no gridlines
        ggplot2::scale_color_manual(values = palette)+
        tema_labels)
    }
  }
}

dev.off()

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. BATCH effect correction (platform): ComBAT
batch = pheno_data$platform_id
#' model w effect of interest
modcombat = model.matrix(~ dis_condition, data=pheno_data)
combat_edata = ComBat(dat=expression, batch=batch, mod=modcombat,
                      par.prior=TRUE, prior.plots=FALSE)

# PCA after batch correction
pca_combat <- perform_pca(t(combat_edata), num_components = 30, # samples x genes
                       save_eig_plot = TRUE, eig_ncp = 12,
                       output_dir = output_dir, plot_filename="pca_eigenvalues_post_combat.png")
pca_combat_pheno_data <- join_pca_with_metadata(pca_combat$pca_data, pheno_data)
pca_summary <- summary(pca_combat$pca_summary)

## Plot PCA after batch correction
pdf(file.path(output_dir, "pca_post_combat.pdf"), width = 7, height = 6)

for(pc in c("PC2","PC3")){
  for(feature in c("dis_condition","platform_id", "sex", "smoker", "age",
                  "emphysema", "pred_dlco", "pred_fev_postbd", "pred_fvc_postbd")){
    if(is.numeric(pca_combat_pheno_data[[feature]])){
      print(ggplot(
        pca_combat_pheno_data,
        aes(
          x = PC1,
          y = !!sym(pc),
          shape = platform_id,
          color = !!sym(feature) # label points with different colors for each `subgroup`
        )
      ) +
        geom_point() + # Plot individual points to make a scatterplot
        theme_classic() + # Format as a classic-looking plot with no gridlines
        scale_fill_continuous(type = "gradient",na.value="white")+
        tema_labels)
    } else {
      palette <- get(feature)  # paleta con nombre desde utils.R (dis_condition/platform_id/sex/smoker)
      print(ggplot(
        pca_combat_pheno_data,
        aes(
          x = PC1,
          y = !!sym(pc),
          shape = dis_condition,
          color = !!sym(feature) # label points with different colors for each `subgroup`
        )
      ) +
        geom_point() + # Plot individual points to make a scatterplot
        theme_classic() + # Format as a classic-looking plot with no gridlines
        ggplot2::scale_color_manual(values = palette)+
        tema_labels)
    }
  }
}

dev.off()


## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. Plot for manuscript
vars <- c("dis_condition", "platform_id", "sex", "emphysema")

tema_sup <- theme_classic(base_size = 12) +
  theme(legend.key.size = unit(0.35, "cm"),
        plot.title = element_text(hjust = 0.5))

nombres <- c(dis_condition = "disease condition",
             platform_id   = "platform ID",
             sex           = "sex",
             emphysema     = "emphysema")

plot_pc_stage <- function(data, feature, y, titulo = NULL) {
  df <- data[order(!is.na(data[[feature]])), ]          # NA debajo
  if (is.numeric(df[[feature]])) {
    shape_var <- "platform_id"
    p <- ggplot(df, aes(PC1, .data[[y]], shape = platform_id, color = .data[[feature]])) +
      scale_color_continuous(type = "gradient", na.value = "grey85")
  } else {
    shape_var <- "dis_condition"
    p <- ggplot(df, aes(PC1, .data[[y]], shape = dis_condition, color = .data[[feature]])) +
      scale_color_manual(values = get(feature)) +
      guides(color = guide_legend(override.aes = list(size = 2.5)))
  }
  p + geom_point(size = 1.1) +
    scale_shape_manual(values = c(16, 17)) +
    guides(shape = guide_legend(override.aes = list(size = 2.5))) +
    labs(y = y, title = titulo,
         color = nombres[[feature]],
         shape = nombres[[shape_var]]) +
    tema_sup
}

filas <- lapply(seq_along(vars), function(i) {
  f <- vars[i]
  tit <- function(t) if (i == 1) t else NULL            # títulos solo en la primera fila
  wrap_plots(
    plot_pc_stage(pca_pheno_data,        f, "PC2", tit("Pre-ComBat")),
    plot_pc_stage(pca_pheno_data,        f, "PC3", tit("Pre-ComBat")),
    plot_pc_stage(pca_combat_pheno_data, f, "PC2", tit("Post-ComBat")),
    plot_pc_stage(pca_combat_pheno_data, f, "PC3", tit("Post-ComBat")),
    ncol = 4
  ) + plot_layout(guides = "collect")
})

figura_pca <- wrap_plots(filas, ncol = 1)

ggsave(file.path(output_dir, "FigS_PCA_pre_post_ComBat.svg"), figura_pca,
       width = 11.69, height = 8.27, units = "in")

ggsave(file.path(output_dir, "FigS_PCA_pre_post_ComBat.jpeg"), figura_pca,
       width = 11.69, height = 8.27, units = "in",
       dpi = 600, quality = 100, bg = "white")
