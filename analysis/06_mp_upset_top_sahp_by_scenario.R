#!/usr/bin/env Rscript
######################################
## Upset plot of top SHAP values by scenario
## Iria Pose
## 19 nov 2025
######################################

## --------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(EBEx)
library(dplyr)
library(UpSetR)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
candidate_results_dir <- "test/candidate_genes"
output_dir <- "test/plots"
classifiers <- c("rf", "svm_r", "svm_p", "glm", "knn", "xgb")
input <- "disease_related"
if(input != "disease_related"){
  input_lists <- c("data_driven",
                   "dea",
                   "mrmr",
                   "omnipath_data_driven",
                   "guildify_data_driven")
}else{
  input_lists <- c("disease_related",
                      "disease_related_entire_list",
                      "omnipath_disease_related",
                       "guildify_disease_related")}

## ----------------------------------------------------------------------------------------------------------------------------------------
# LOAD DATA
print_message("Loading data...")
# 1. Load normalized shap data
normalized_shap_list <- readRDS(file.path(candidate_results_dir, "normalized_shap_list.Rds"))
# 2. Filter input_lists
normalized_shap_list_filtered <- normalized_shap_list[
  sapply(names(normalized_shap_list), function(x) {
    any(sapply(input_lists, function(pat) grepl(pat, x)))
  })
]
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# PREPARE DATA
# 3. Aggregated shap by gene for each combination
print_message("Aggregating SHAP values by gene for each scenario...")
aggregated_shap <- map(normalized_shap_list_filtered, function(res){
  res %>%
    group_by(variable_name) %>%
    summarize(mean_norm_shap = mean(abs(norm_shap)), .groups = "drop") %>%
    arrange(desc(mean_norm_shap)) %>% filter(mean_norm_shap !=0)
})

# 4. Extract top 30 genes
print_message("Extracting top 30 genes for each scenario...")
top_genes <- map(aggregated_shap, function(res){
    res %>%
    slice_max(mean_norm_shap, n=30) %>% select(variable_name)
})

# 4.1. Combine all `variable_name` columns into a single vector
all_genes <- unlist(lapply(top_genes, function(df) df$variable_name))
# 4.2. Count occurrences of each gene
gene_counts <- sort(table(all_genes))
print(gene_counts)
# saveRDS(gene_counts, paste0("../COPD/gene_counts_",input,".Rds"))

# 4.3. Extract the list of genes from the tibbles
gene_list <- lapply(top_genes, function(tibble) tibble$variable_name)

## ----------------------------------------------------------------------------------------------------------------------------------------
# PLOT
# 5. Plot upset plot
print_message("Creating UpSet plot for top SHAP ",input," genes...")
# 5.1. Convert the list of genes into a presence/absence matrix
upset_data <- fromList(gene_list)
colnames(upset_data) <- rename_input_classifier(colnames(upset_data))

# Open a PDF device
pdf(file.path(output_dir, paste0("upset_plot_top_shap_",input,"_genes.pdf")), width = 22, height = 12)
upset(
  upset_data,
  sets = colnames(upset_data), # All columns except the "gene" column
  nintersects = NA,            # Display all intersections
  keep.order = TRUE,            # Keep the order of sets as in the data
  order.by = "freq",           # Order intersections by frequency
  text.scale = c(2, 2, 2, 2, 2, 2), # Scale text size (axes, labels, etc.)
  mb.ratio = c(0.5, 0.5)          # Adjust the ratio of the main bar plot to the set size plot (default is c(0.5, 0.5))
)
dev.off()
