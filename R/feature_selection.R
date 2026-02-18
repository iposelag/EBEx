#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom OmnipathR import_all_interactions
#' @importFrom stats ks.test qnorm sd
#' @importFrom utils read.csv read.table
NULL

#' Load a gene list from a file
#'
#' @description
#' This function loads a gene list from a specified text file. The file is expected to
#' contain gene names separated by commas.
#'
#' @param genes_list Character string of the filename.
#' @param base_path Base directory. Defaults to "../COPD/raw_data/".
#'
#' @return A character vector of gene names.
#'
#' @export
load_gene_list <- function(genes_list, base_path = "../COPD/feature_selection/") {

  file_path <- file.path(base_path, paste0(genes_list, ".txt"))
  genes <- scan(file_path, what = "character", sep = ",")

  return(genes)
}

#' Main Function for Feature Selection
#'
#' @description 
#' This function serves as the main entry point for performing feature selection based on 
#' various procedures. It takes in the necessary parameters and calls the appropriate helper
#' functions to obtain the selected gene lists according to the specified procedure.
#'
#' @param procedure The feature selection procedure to use.
#' @param expression_data The expression data frame.
#' @param target_var The target variable name.
#' @param threshold_value Threshold value for selection (if applicable).
#' @param disease_code Disease code for disease-related selections.
#' @param dea_genes List of differentially expressed genes (if applicable).
#' @param mrmr_genes List of mRMR selected genes (if applicable
#' @param alternative_genes List of alternative selected genes (if applicable).
#' @param directory_to_load Directory to load data from (if applicable).
#' @param directory_to_save Directory to save results to (if applicable).
#' @param mrmr_path Path to the mRMR executable.
#'
#' @return A character vector of selected gene names based on the specified procedure.
#'
#' @export
run_feature_selection <- function(procedure, expression_data, target_var, threshold_value = NULL, disease_code = NULL, dea_genes = NULL, mrmr_genes = NULL, alternative_genes = NULL, directory_to_load = NULL, directory_to_save = NULL, mrmr_path = "./mrmr") {
  
  switch(procedure,
         "mrmr" = obtain_mrmr(expression_data, target_var, threshold_value, directory_to_save, mrmr_path),
         "disease_related" = obtain_disease_related_curated_genes(disease_code, directory_to_load, directory_to_save),
         "disease_related_entire_list" = obtain_disease_related_entire_genes(disease_code, directory_to_load, directory_to_save),
         "data_driven" = obtain_data_driven_genes(dea_genes, mrmr_genes, expression_data, target_var, threshold_value, directory_to_load, directory_to_save, mrmr_path),
         "omnipath_disease_related" = {
           seeds <- obtain_disease_related_curated_genes(disease_code, directory_to_load, directory_to_save)
           obtain_omnipath_expansion(seeds, directory_to_save)
         },
         "omnipath_data_driven" = {
           seeds <- obtain_data_driven_genes(dea_genes, mrmr_genes, expression_data, target_var, threshold_value, directory_to_load, directory_to_save, mrmr_path)
           obtain_omnipath_expansion(seeds, directory_to_save)
         },
         "omnipath_intersection" = obtain_omnipath_combined(expression_data, target_var, disease_code, threshold_value, dea_genes, mrmr_genes, directory_to_load, directory_to_save, operation = "intersection", mrmr_path = mrmr_path),
         "omnipath_union" = obtain_omnipath_combined(expression_data, target_var, disease_code, threshold_value, dea_genes, mrmr_genes, directory_to_load, directory_to_save, operation = "union", mrmr_path = mrmr_path),
         "alternative" = obtain_alternative_genes(alternative_genes, directory_to_load, directory_to_save),
         stop("Invalid procedure specified.")
  )

}

#' Obtain alternative selection of genes
#'
#' @description
#' This function loads an alternative list of genes from a specified text file. The file is expected to
#' contain gene names separated by commas. The function also saves the loaded gene list to a specified
#' directory if the `directory_to_save` parameter is provided.
#'
#' @param alternative_genes Character string of the filename containing the alternative gene list.
#' @param directory_to_load Directory to load the gene list from.
#' @param directory_to_save Directory to save the loaded gene list to (optional).
#'
#' @return A character vector of gene names loaded from the specified file.
#'
#' @export
obtain_alternative_genes <- function(alternative_genes, directory_to_load, directory_to_save = NULL) {

  cat("Obtaining alternative list of genes:", alternative_genes, "\n")
  file_path <- file.path(directory_to_load, paste0(alternative_genes, ".txt"))
  genes_list <- scan(file_path, what = "character", sep = ",")
  if (!is.null(directory_to_save)) {
    dir_path <- file.path(directory_to_save, "feature_selection")
    if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
    write(genes_list, file = file.path(dir_path, paste0(alternative_genes, ".txt")), ncolumns = 1)
  }

  return(genes_list)
}

#' Obtain mRMR genes
#'
#' @description
#' This function runs the mRMR feature selection algorithm on the provided expression data and target variable.
#' It can also perform bootstrapping to determine a score threshold for selecting features. The selected mRMR 
#' genes are saved to a specified directory if the `directory_to_save` parameter is provided.
#'
#' @param expression_data Data frame containing the expression data.
#' @param target_var Character string of the target variable column name.
#' @param threshold_value Numeric value for the mRMR score threshold. If NULL, bootstrapping will be performed to determine the threshold.
#' @param directory_to_save Directory to save the selected mRMR genes to (optional).
#' @param mrmr_path Path to the mRMR executable.
#'
#' @return A character vector of selected mRMR gene names.
#'
#' @export
obtain_mrmr <- function(expression_data, target_var, threshold_value = NULL, directory_to_save = NULL, mrmr_path = "./mrmr") {

  cat("Obtaining mrmr genes\n")
  if(is.null(threshold_value)){
    cat("Start of bootstrapping process\n")
    all_mrmr_scores_df <- run_mrmr_bootstrapping(expression_data, target_var, n_iterations = 100, directory_to_save = directory_to_save, mrmr_path = mrmr_path)
    threshold_value <- extract_score_threshold(all_mrmr_scores_df, "score")
  }
  mrmr_raw <- run_mrmr(expression_data, target_var, variables = ncol(expression_data) - 1, samples = nrow(expression_data), mrmr_path = mrmr_path)
  mrmr_results <- extract_mrmr_features(mrmr_raw)
  mrmr_list <- trimws(mrmr_results[mrmr_results$Score >= threshold_value, ]$Name)
  if (!is.null(directory_to_save)) {
    save_helper_txt(mrmr_list, "mrmr.txt", directory_to_save, "feature_selection")
  }

  return(mrmr_list)
}

#' Obtain Data Driven Genes (DEA + mRMR)
#'
#' @description
#' This function combines the results from Differential Expression Analysis (DEA) and mRMR feature
#' selection to obtain a list of data-driven genes. It takes in the lists of DEA and mRMR genes, 
#' as well as the expression data and target variable, to determine the final set of genes. 
#' The resulting list can be saved to a specified directory if the `directory_to_save` parameter is provided.
#'
#' @param dea_genes List of differentially expressed genes (optional).
#' @param mrmr_genes List of mRMR selected genes (optional).
#' @param expression_data Data frame containing the expression data.
#' @param target_var Character string of the target variable column name.
#' @param threshold_value Numeric value for the mRMR score threshold. If NULL, bootstrapping will be performed to determine the threshold.
#' @param directory_to_load Directory to load the DEA results from.
#' @param directory_to_save Directory to save the selected data-driven genes to (optional).
#' @param mrmr_path Path to the mRMR executable.
#'
#' @return A character vector of selected data-driven gene names.
#'
#' @export
obtain_data_driven_genes <- function(dea_genes = NULL, mrmr_genes = NULL, expression_data, target_var, threshold_value = NULL, directory_to_load, directory_to_save = NULL, mrmr_path = "./mrmr"){

    cat("Obtaining data driven genes\n")
    if(is.null(mrmr_genes)) mrmr_genes <- obtain_mrmr(expression_data, target_var, threshold_value, directory_to_save, mrmr_path)
    if(is.null(dea_genes)){
      dea_genes <- scan(file.path(directory_to_load, "dea.txt"), what = "character", sep = ",")
    }
    data_driven <- union(dea_genes, mrmr_genes)
    if(!is.null(directory_to_save)) save_helper_txt(data_driven, "data_driven.txt", directory_to_save, "feature_selection")

    return(data_driven)
}

#' Combined Omnipath Expansion (Intersection/Union)
#' 
#' @description
#' This function combines the Omnipath expansions of both disease-related and data-driven gene lists 
#' using either intersection or union operations. It takes in the necessary parameters to obtain the 
#' disease-related and data-driven gene lists, performs the Omnipath expansion for each, and then 
#' combines the results based on the specified operation. The final list of genes can be saved to 
#' a specified directory if the `directory_to_save` parameter is provided.
#'
#' @param expression_data Data frame containing the expression data.
#' @param target_var Character string of the target variable column name.
#' @param disease_code Disease code for disease-related selections.
#' @param threshold_value Numeric value for the mRMR score threshold. If NULL, bootstrapping will be performed to determine the threshold.
#' @param dea_genes List of differentially expressed genes (optional).
#' @param mrmr_genes List of mRMR selected genes (optional).
#' @param directory_to_load Directory to load the necessary data from.
#' @param directory_to_save Directory to save the combined gene list to (optional).
#' @param operation Character string specifying the combination operation: "intersection" or "union".
#' @param mrmr_path Path to the mRMR executable.
#'
#' @return A character vector of combined gene names based on the specified operation.
#'
#' @export
obtain_omnipath_combined <- function(expression_data, target_var, disease_code, threshold_value, dea_genes, mrmr_genes, directory_to_load, directory_to_save = NULL, operation = "intersection", mrmr_path = "./mrmr") {
  disease_related <- obtain_disease_related_curated_genes(disease_code, directory_to_load, directory_to_save)
  data_driven <- obtain_data_driven_genes(dea_genes, mrmr_genes, expression_data, target_var, threshold_value, directory_to_load, directory_to_save, mrmr_path)
  
  exp_dis <- obtain_omnipath_expansion(disease_related)
  exp_dat <- obtain_omnipath_expansion(data_driven)
  
  genes_list <- if(operation == "intersection") intersect(exp_dat, exp_dis) else union(exp_dat, exp_dis)
  if (!is.null(directory_to_save)) save_helper_txt(genes_list, paste0("omnipath_", operation, ".txt"), directory_to_save, "feature_selection")
  return(genes_list)
}

#' Obtain disease-related entire genes
#'
#' @description
#' This function loads the entire list of disease-related genes from a specified file. The file is 
#' expected to be a TSV file containing a summary of gene-disease associations for the given disease code. 
#' The function extracts the gene names from the file and saves the list to a specified directory if the 
#' `directory_to_save` parameter is provided.
#'
#' @param disease_code Disease code for which to obtain the gene list.
#' @param directory_to_load Directory to load the gene list from.
#' @param directory_to_save Directory to save the loaded gene list to (optional).
#'
#' @return A character vector of gene names related to the specified disease.
#'
#' @export
obtain_disease_related_entire_genes <- function(disease_code, directory_to_load, directory_to_save = NULL){

    cat("Obtaining disease related entire genes\n")
    file_path <- file.path(directory_to_load, "disgenet_tables", paste0(disease_code, "_disease_gda_summary.tsv"))
    disgenet_entire <- utils::read.csv(file_path, sep ="\t")
    genes_list <- disgenet_entire$Gene
    if (!is.null(directory_to_save)) save_helper_txt(genes_list, "disease_related_entire_list.txt", directory_to_save, "feature_selection")
    
    return(genes_list)
}

#' Obtain curated disease genes
#'
#' @description 
#' This function loads a curated list of disease-related genes from a specified file. The file is expected
#' to be a TSV file containing gene-disease associations with a "curated" curation effort for the given disease code.
#' The function extracts the gene symbols from the file and saves the list to a specified directory if the 
#' `directory_to_save` parameter is provided.
#'
#' @param disease_code Disease code for which to obtain the curated gene list.
#' @param directory_to_load Directory to load the gene list from.
#' @param directory_to_save Directory to save the loaded gene list to (optional).
#'
#' @return A character vector of curated gene symbols related to the specified disease.
#'
#' @export
obtain_disease_related_curated_genes <- function(disease_code, directory_to_load, directory_to_save = NULL) {
cat("Obtaining disease related curated genes\n")
  disgenet <- utils::read.csv(file.path(directory_to_load, "disgenet_curated_gene_disease_associations.tsv"), sep = "\t")
  gene_list <- disgenet[disgenet$diseaseId == disease_code, ]$geneSymbol
  if (!is.null(directory_to_save)) save_helper_txt(gene_list, "disease_related.txt", directory_to_save, "feature_selection")
  return(gene_list)
}

#' Obtain Omnipath expansion
#'
#' @description
#' This function takes a list of seed genes and retrieves their interaction partners from the Omnipath database. 
#' It filters the interactions to include only those with a curation effort greater than 1 and returns a unique list 
#' of genes that includes both the original seed genes and their interaction partners.
#' The resulting list of genes can be saved to a specified directory if the `directory_to_save` parameter is provided.
#'
#' @param genes A character vector of seed gene symbols to expand.
#' @param directory_to_save Directory to save the expanded gene list to (optional).
#'
#' @return A character vector of gene symbols that includes the original seed genes and their Omnipath interaction partners.
#'
#' @export
obtain_omnipath_expansion <- function(genes, directory_to_save = NULL){

    cat("Obtaining omnipath expansion genes\n")
    omnipath <- OmnipathR::import_all_interactions(organism = 9606, directed  = 'no') %>%
        dplyr::filter(.data$source_genesymbol %in% genes | .data$target_genesymbol %in% genes) %>%
        dplyr::filter(.data$curation_effort > 1)
    expansion_omnipath <- unique(c(omnipath$target_genesymbol, omnipath$source_genesymbol, genes))
    if (!is.null(directory_to_save)) save_helper_txt(expansion_omnipath, "omnipath_expansion.txt", directory_to_save, "feature_selection")

    return(expansion_omnipath)
}

#' @keywords internal
run_mrmr <- function(data, target_var, variables, samples, mrmr_path) {

  data[[target_var]] <- as.integer(as.factor(data[[target_var]]))
  temp_file <- tempfile(fileext = ".csv")
  write.csv(data, temp_file, row.names = FALSE)
  command <- sprintf("%s -i %s -n 500 -t 1 -v %s -s %s", mrmr_path, temp_file, variables, samples)

  return(system(command, intern = TRUE))
}

#' @keywords internal
extract_mrmr_features <- function(results) {

  mrmr_start <- which(grepl("\\*\\*\\* mRMR features \\*\\*\\*", results)) + 2
  feature_lines <- results[mrmr_start:length(results)]
  feature_lines <- feature_lines[feature_lines != ""]
  end_idx <- which(grepl("^\\s*\\*{3}", feature_lines))[1]
  feature_lines <- feature_lines[1:(end_idx - 1)]
  df <- utils::read.table(text = feature_lines, header = FALSE, sep = "\t")
  colnames(df) <- c("Order", "Feature", "Name", "Score")

  return(df)
}

#' @keywords internal
run_mrmr_bootstrapping <- function(expression_data, target_var, n_iterations = 10, directory_to_save = NULL, mrmr_path = "./mrmr") {

  all_results <- data.frame()
  for (i in 1:n_iterations) {
    shuffled <- shuffling_sample_disease_categories(expression_data, target_var)
    raw <- run_mrmr(shuffled, target_var, ncol(shuffled)-1, nrow(shuffled), mrmr_path)
    df <- extract_mrmr_features(raw)
    df$iteration <- i
    df$score <- df$Score # for consistency with your threshold extractor
    all_results <- rbind(all_results, df)
  }
  
  return(all_results)
}

#' @keywords internal
shuffling_sample_disease_categories <- function(data, target_var){

  data[[target_var]] <- sample(data[[target_var]])

  return(data)
}

#' @keywords internal
extract_score_threshold <- function(data, score_col, threshold_fraction = 0.01) {

  vals <- data[[score_col]]
  ks <- stats::ks.test(vals, "pnorm", mean = mean(vals), sd = stats::sd(vals))

  if (ks$p.value < 0.05) {
    return(sort(vals, decreasing = TRUE)[ceiling(threshold_fraction * length(vals))])
  } else {
    return(mean(vals) + stats::qnorm(1 - threshold_fraction) * stats::sd(vals))
  }
}