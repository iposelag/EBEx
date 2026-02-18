#!/usr/bin/env Rscript
######################################
## Functional Enrichment of candidate genes
## Iria Pose
## 19 nov 2025
######################################

## --------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(devtools)
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
library(ggplot2)
library(dplyr)

## --------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
functional_enrichment_results_dir <- "test/analysis/functional_enrichment"
output_dir <- "test/plots"
candidate_results_dir <- "test/candidate_genes"

## --------------------------------------------------------------------------------------------------------------------
# 1. LOAD DATA
print_message("Loading candidate (Canonical and Extended) genes...")
candidate_genes <- read.csv(file.path(candidate_results_dir, "candidate_genes_density.csv"))$variable_name
canonical_genes <- load_gene_list("canonical_genes", file.path(candidate_results_dir))
extended_genes <- load_gene_list("extended_genes", file.path(candidate_results_dir))
print_message("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 2. Run enrichment
# 2.1 Load available databases
print_message("Running functional enrichment analysis...\n")
dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_Pathways_2024")
print_message("Databases loaded: ", paste(dbs, collapse = ", "), "\n")
# 2.2. Run enrichment for extended genes and canonical genes
print_message("Enriching Extended Genes...\n")
enr_extended <- run_functional_enrichment(extended_genes, dbs = dbs)
# Save results into csv files
for (db_name in names(enr_extended)) {
    df <- enr_extended[[db_name]]
    # Only save if there are results (rows > 0)
    if (nrow(df) > 0) {
      #   Create a clean file name (e.g., extended_genes_GO_Biological_Process_2023.csv)
      clean_db_name <- gsub(" ", "_", db_name)
      filename <- paste0("extended_genes_", clean_db_name, ".csv")
      save_helper_csv(
        data = df,
        filename = filename,
        directory = functional_enrichment_results_dir
      )
    }
}
print_message("Enriching Canonical Genes...\n")
enr_canonical <- run_functional_enrichment(canonical_genes, dbs = dbs)
# Save results into csv files
for (db_name in names(enr_canonical)) {
    df <- enr_canonical[[db_name]]
    # Only save if there are results (rows > 0)
    if (nrow(df) > 0) {
      #   Create a clean file name (e.g., canonical_genes_GO_Biological_Process_2023.csv)
      clean_db_name <- gsub(" ", "_", db_name)
      filename <- paste0("canonical_genes_", clean_db_name, ".csv")
      save_helper_csv(
        data = df,
        filename = filename,
        directory = output_directory
      )
    }
}
print_message("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 3. Visualization
print_message("Creating enrichment dotplots...\n")
# 3.1. Dotplot para GO (extended Genes)
p_go_extended <- plot_enrichment_dotplot(extended_go,
                                      title = "GO: Extended Candidate Genes",
                                      n_input_genes = length(extended_genes),
                                      aspect_ratio = 4)
if (!is.null(p_go_extended)) {
  ggplot2::ggsave(file.path(output_plots_dir, "dotplot_GO_Extended.pdf"), p_go_extended, width = 10, height = 9.7)
}
# 3.2. Dotplot para KEGG (extended Genes)
p_kegg_extended <- plot_enrichment_dotplot(extended_kegg,
                                          title = "KEGG: Extended Candidate Genes",
                                          n_input_genes = length(extended_genes),
                                          aspect_ratio = 2.8)
if (!is.null(p_kegg_extended)) {
  ggplot2::ggsave(file.path(output_plots_dir, "dotplot_KEGG_Extended.pdf"), p_kegg_extended, width = 11, height = 9.7)
}
# 3.3. Dotplot para Reactome (extended Genes)
p_reactome_extended <- plot_enrichment_dotplot(extended_reactome,
                                        title = "REACTOME: Extended Candidate Genes",
                                        n_input_genes = length(extended_genes),
                                        aspect_ratio = 3)
if (!is.null(p_reactome_extended)) {
  ggplot2::ggsave(file.path(output_plots_dir, "dotplot_Reactome_Extended.pdf"), p_reactome_extended, width = 10, height = 9.7)
}
# 3.4. Dotplot para GO (canonical Genes)
p_go_canonical <- plot_enrichment_dotplot(canonical_go,
                                    title = "GO: Canonical Candidate Genes",
                                    n_input_genes = length(canonical_genes),
                                    aspect_ratio = 2)
if (!is.null(p_go_canonical)) {
  ggsave(filename = file.path(output_plots_dir, "dotplot_go_canonical.pdf"), plot = p_go_canonical, device = "pdf", width = 10, height = 9.7)
}
# 3.5. Dotplot para KEGG (canonical Genes)

p_kegg_canonical <- plot_enrichment_dotplot(canonical_kegg,
                                      title = "KEGG: Canonical Candidate Genes",
                                      n_input_genes = length(canonical_genes),
                                      aspect_ratio = 3)
if (!is.null(p_kegg_canonical)) {
  ggsave(filename = file.path(output_plots_dir, "dotplot_kegg_canonical.pdf"), plot = p_kegg_canonical, device = "pdf", width = 10, height = 9.7)
}
# 3.6. Dotplot para Reactome (canonical Genes)
p_reactome_canonical <- plot_enrichment_dotplot(knwon_reactome,
                                          title = "REACTOME: Canonical Candidate Genes",
                                          n_input_genes = length(canonical_genes),
                                          aspect_ratio = 3)
if (!is.null(p_reactome_canonical)) {                                          
  ggsave(filename = file.path(output_plots_dir, "dotplot_Reactome_canonical.pdf"), plot = p_reactome_canonical, device = "pdf", width = 10, height = 9.7)
}
print_message("________________________________\n")

## --------------------------------------------------------------------------------------------------------------------
# 4. FINISH
print_message("Enrichment analysis completed. Results saved in:", functional_enrichment_results_dir)
print_message("Plots saved in:", output_dir)
print_message("________________________________\n")

# p_reactome_extended <- plot_enrichment_dotplot(extended_reactome,
#                                         title = "REACTOME: Extended Candidate Genes",
#                                         n_input_genes = length(extended_genes),
#                                         aspect_ratio = 3)
# p_reactome_extended
# if (!is.null(p_reactome_extended)) {
#   ggplot2::ggsave(file.path(output_plots_dir, "dotplot_Reactome_Extended.pdf"), p_reactome_extended, width = 10, height = 9.7)
# }

# p_kegg_extended <- plot_enrichment_dotplot(extended_kegg,
#                                           title = "KEGG: Extended Candidate Genes",
#                                           n_input_genes = length(extended_genes),
#                                           aspect_ratio = 2.8)
# p_kegg_extended
# if (!is.null(p_kegg_extended)) {
#   ggplot2::ggsave(file.path(output_plots_dir, "dotplot_KEGG_Extended.pdf"), p_kegg_extended, width = 11, height = 9.7)
# }

# p_go_extended <- plot_enrichment_dotplot(extended_go,
#                                       title = "GO: Extended Candidate Genes",
#                                       n_input_genes = length(extended_genes),
#                                       aspect_ratio = 4)
# p_go_extended
# if (!is.null(p_go_extended)) {
#   ggplot2::ggsave(file.path(output_plots_dir, "dotplot_GO_Extended.pdf"), p_go_extended, width = 10, height = 9.7)
# }

# p_reactome_canonical <- plot_enrichment_dotplot(knwon_reactome,
#                                           title = "REACTOME: Canonical Candidate Genes",
#                                           n_input_genes = length(canonical_genes),
#                                           aspect_ratio = 3)
# ggsave(filename = file.path(output_plots_dir, "dotplot_Reactome_canonical.pdf"), plot = p_reactome_canonical, device = "pdf", width = 10, height = 9.7)

# p_kegg_canonical <- plot_enrichment_dotplot(canonical_kegg,
#                                       title = "KEGG: Canonical Candidate Genes",
#                                       n_input_genes = length(canonical_genes),
#                                       aspect_ratio = 3)
# ggsave(filename = file.path(output_plots_dir, "dotplot_kegg_canonical.pdf"), plot = p_kegg_canonical, device = "pdf", width = 10, height = 9.7)
# p_go_canonical <- plot_enrichment_dotplot(canonical_go,
#                                     title = "GO: Canonical Candidate Genes",
#                                     n_input_genes = length(canonical_genes),
#                                     aspect_ratio = 2)
# ggsave(filename = file.path(output_plots_dir, "dotplot_go_canonical.pdf"), plot = p_go_canonical, device = "pdf", width = 10, height = 9.7)

