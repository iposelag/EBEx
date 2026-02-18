#!/usr/bin/env Rscript
######################################
## Analysis of origins of candidate genes: generate canonical and extended lists
## Plots: UpSet plots of gene list overlaps
## Iria Pose
## 2025
######################################

## --------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
library(tidyverse)
library(UpSetR)
library(grid)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
files_to_load <- c(
  "disease-related curated"          = "disease_related",
  "disease-related entire list"      = "disease_related_entire_list",
  "disease-related GUILDify" = "guildify_disease_related",
  "disease-related OmniPath" = "omnipath_disease_related",
  "mRMR"                     = "mrmr",
  "DEA"                      = "dea",
  "data-driven"              = "data_driven",
  "data-driven GUILDify"     = "guildify_data_driven",
  "data-driven OmniPath"     = "omnipath_data_driven"
)
candidate_results_dir <- "test/candidate_genes"
features_lists_dir <- "test/feature_selection"
plots_dir <- "test/plots"

## --------------------------------------------------------------------------------------------------------------------
# 2. LOAD DATA
# 2.1. Load candidate genes selected by density
cat("Loading candidate genes from density ranking...\n")
candidate_genes_df <- read.csv(file.path(candidate_results_dir, "candidate_genes_density.csv"))
top_genes <- candidate_genes_df$variable_name
# 2.2. Load all origin lists automatically into an R list
cat("Loading gene lists for all origins...\n")
gene_lists <- lapply(files_to_load, load_gene_list, base_path = features_lists_dir)
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 3. PREPARE DATA
cat("Preparing data for UpSet plot...\n")
# All unique genes across all lists
all_unique_genes <- unique(unlist(gene_lists))
# Create binary dataframe (1 if in list, 0 if not)
feature_lists <- data.frame(Gene = all_unique_genes)
for (list_name in names(gene_lists)) {
  feature_lists[[list_name]] <- as.integer(all_unique_genes %in% gene_lists[[list_name]])
}
# Filter the matrix to keep only candidate genes (those in the ranking)
upset_data_full <- feature_lists %>%
  filter(Gene %in% top_genes)
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 4. IDENTIFY CANONICAL AND EXTENDED GENES
# 4.1. Define canonical genes as those already present in DEA or the curated DisGeNET list
cat("Identifying canonical and extended genes...\n")
canonical_genes <- upset_data_full %>% 
  filter(DEA == 1 | `disease-related curated` == 1) %>% 
  pull(Gene)
extended_genes <- setdiff(top_genes, canonical_genes)
# Save resulting lists using package helpers
write(canonical_genes, file = file.path(candidate_results_dir, "canonical_genes.txt"), ncolumns = length(canonical_genes), sep = ",")
write(extended_genes,   file = file.path(candidate_results_dir, "extended_genes.txt"),   ncolumns = length(extended_genes),   sep = ",")
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 5. PREPARE DATA FOR PLOTTING
# Group expansions into two categories: "Expansions disease-related" and "Expansions data-driven"
upset_data_grouped <- upset_data_full %>%
  mutate(
    'Expansions disease-related' = as.integer(`disease-related entire list` | `disease-related GUILDify` | `disease-related OmniPath`),
    'Expansions data-driven'    = as.integer(`data-driven GUILDify` | `data-driven OmniPath`)
  ) %>%
  select(Gene, DEA, `disease-related curated`, mRMR, `Expansions disease-related`, `Expansions data-driven`)

## --------------------------------------------------------------------------------------------------------------------
# 6. PLOTTING
cat("Generating UpSet plots...\n")
pdf(file.path(plots_dir, "upset_origins_genes.pdf"), width = 10, height = 7)
# Create a new page and divide it into two viewports (top and bottom)
# Top Panel: Detailed
  print(upset(upset_data_full[, -1], sets = colnames(upset_data_full)[-1],
              nintersects = NA, order.by = "freq", keep.order = TRUE,
              text.scale = c(2.5, 2.5, 2.5, 2.5, 2, 2.5)))
dev.off()
pdf(file.path(plots_dir, "upset_origins_genes_ced_grouped.pdf"), width = 10, height = 7)
# Create a new page and divide it into two viewports (top and bottom)
# Top Panel: Detailed
  print(upset(upset_data_grouped[, -1], sets = colnames(upset_data_grouped)[-1],
              nintersects = NA, order.by = "freq", keep.order = TRUE,
              text.scale = c(2.5, 2.5, 2.5, 2.5, 1.8, 2.5)))
dev.off()
cat("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 7. FINISH
cat("UpSet plots saved in:", plots_dir, "\n")