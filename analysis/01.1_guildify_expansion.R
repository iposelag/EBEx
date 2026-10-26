#!/usr/bin/env Rscript
######################################
## Feature selection: GUILDify expansion of data_driven and disease_related lists
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(EBEx)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
base_dir <- "test"
feature_selection_dir <- file.path(base_dir, "feature_selection")
preprocessing_dir <- file.path(base_dir, "analysis/preprocessing")
input_lists <- c("data_driven", "disease_related")
top_fraction <- 0.01

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. LOAD DATA
load(file.path(preprocessing_dir, "expression_train.Rda"))

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. EXPANSION + CHECKS
guildify_results <- list()
for (name in input_lists) {
  print_message("________________________________")
  print_message("GUILDify expansion:", name)

  original <- load_gene_list(name, feature_selection_dir)
  # NOTE: Gene.Symbol is not trimmed on purpose (reproduces original lists; e.g. "BAG6 " is lost in the filter)
  guild <- read.delim(file.path(feature_selection_dir, paste0("guild_scores_", name, ".txt")), check.names = TRUE)
  if (is.unsorted(rev(guild$GUILD.Score))) {
    print_message("WARNING:", name, ": file not sorted by GUILD.Score, reordering")
    guild <- guild[order(guild$GUILD.Score, decreasing = TRUE), ]
  }

  # Seeds: Seed column, cross-checked by membership (some rows hold several symbols, e.g. "CT45A10 ; CT45A5")
  row_symbols <- lapply(strsplit(guild$Gene.Symbol, ";"), trimws)
  is_member <- sapply(row_symbols, function(s) any(s %in% original))
  if (!identical(guild$Seed == 1, is_member)) {
    print_message("WARNING:", name, ": Seed column differs from membership in original list: ",
            paste(guild$Gene.Symbol[(guild$Seed == 1) != is_member], collapse = ", "))
  }
  seeds <- guild$Gene.Symbol[guild$Seed == 1]
  seeds_missing <- setdiff(original, unlist(row_symbols))

  candidates <- guild$Gene.Symbol[guild$Seed != 1]
  top_n <- round(top_fraction * length(candidates))
  # NOTE: top taken from the position of the last seed on purpose (reproduces original lists, one-position shift)
  top_genes <- guild$Gene.Symbol[length(seeds):(length(seeds) + top_n - 1)]

  seed_pos <- which(guild$Gene.Symbol %in% seeds)
  if (max(seed_pos) != length(seeds)) {
    print_message("WARNING:", name, ": seeds are not the first ", length(seeds), " positions (max seed position = ", max(seed_pos), ")")
  }

  final <- union(original, top_genes)
  n_expected <- length(original) + top_n - length(intersect(original, top_genes))
  stopifnot(length(final) == n_expected)

  final_filtered <- intersect(final, colnames(expression_train))
  lost <- setdiff(final, final_filtered)
  stopifnot(!anyDuplicated(final_filtered), !anyNA(final_filtered))

  cat(sprintf("  Genes in GUILDify file:        %d\n", nrow(guild)))
  cat(sprintf("  Original list:                 %d\n", length(original)))
  cat(sprintf("  Seeds found in GUILDify:       %d\n", length(seeds)))
  cat(sprintf("  Seeds not found:               %d  %s\n", length(seeds_missing), paste(seeds_missing, collapse = ", ")))
  cat(sprintf("  Non-seed candidates:           %d\n", length(candidates)))
  cat(sprintf("  top_n (%.0f%%):                  %d\n", 100 * top_fraction, top_n))
  cat(sprintf("  Union (expected / obtained):   %d / %d\n", n_expected, length(final)))
  cat(sprintf("  After expression_train filter: %d (lost %d)\n", length(final_filtered), length(lost)))
  cat(sprintf("    lost original genes: %s\n", paste(intersect(lost, original), collapse = ", ")))
  cat(sprintf("    lost top genes:      %s\n", paste(intersect(lost, top_genes), collapse = ", ")))

  guildify_results[[name]] <- list(guild = guild, top_genes = top_genes, final = final_filtered)
}
