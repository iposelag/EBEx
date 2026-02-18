#!/usr/bin/env Rscript
######################################
## Heatmap of candidate genes in GTEx expression data
## Iria Pose
## 19 nov 2025
######################################

## --------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(devtools)
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
library(data.table)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
output_dir <- "test/plots"
raw_data_dir <- "test/raw_data/"
candidate_results_dir <- "test/candidate_genes"

## --------------------------------------------------------------------------------------------------------------------
# 1. LOAD DATA
print_message("Loading data...")
# GTEx gene mapping
expr <- data.table::fread(file.path(raw_data_dir, "GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"),
              skip = 2)
# Candidate genes
candidate_genes <- read.csv(file.path(candidate_results_dir, "candidate_genes_density.csv"))$variable_name
canonical_genes <- scan(file.path(candidate_results_dir, "canonical_genes.txt"), what = "", sep = ",")
extended_genes <- scan(file.path(candidate_results_dir, "extended_genes.txt"), what = "", sep = ",")
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. PREPARE DATA
print_message("Preparing data for heatmap...")
# Map gene names with GTEx names
expr$Description <- ifelse(expr$Description %in% names(gtex_gene_map), 
                           gtex_gene_map[expr$Description], 
                           expr$Description)
print_message("________________________________")
## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. CREATE HEATMAPS
# Canonical genes
print_message("Creating heatmap for canonical genes...")
expression_sub <- expr[expr$Description %in% canonical_genes, ]
heatmap_data <- as.data.frame(expression_sub[, -c(1,2)])
rownames(heatmap_data) <- expression_sub$Description
# Convert to matrix
heatmap_data <- as.matrix(heatmap_data)
# Plot heatmap
h_gtex_known <- plot_gtex_heatmap(
  tpm_matrix = heatmap_data,
  tissue_macro_map = tissue_macro_map, # Este objeto debe estar cargado o en utils
  tissue_fine_map = tissue_fine_map    # Este objeto debe estar cargado o en utils
)
# Save output
pdf(file.path(output_dir, "heatmap_GTEx_canonical.pdf"), width = 17, height = 14)
ComplexHeatmap::draw(h_gtex_known)
dev.off()

# Extended genes
print_message("Creating heatmap for extended genes...")
expression_sub <- expr[expr$Description %in% extended_genes, ]
heatmap_data <- as.data.frame(expression_sub[, -c(1,2)])
rownames(heatmap_data) <- expression_sub$Description
# Convert to matrix
heatmap_data <- as.matrix(heatmap_data)
# Plot heatmap
h_gtex_new <- plot_gtex_heatmap(
    tpm_matrix = heatmap_data,
    tissue_macro_map = tissue_macro_map, # Este objeto debe estar cargado o en utils
    tissue_fine_map = tissue_fine_map    # Este objeto debe estar cargado o en utils
)
# Save output
pdf(file.path(output_dir, "heatmap_GTEx_extended.pdf"), width = 17, height = 14)
ComplexHeatmap::draw(h_gtex_new)
dev.off()

## --------------------------------------------------------------------------------------------------------------------
# FINISH
print_message(paste0("Heatmaps created and saved in ", output_dir))