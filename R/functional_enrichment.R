#' @importFrom enrichR enrichr
NULL

#' Run Functional Enrichment using EnrichR
#'
#' @description
#' This function performs functional enrichment analysis using the EnrichR package. 
#' It takes a list of gene symbols and a list of EnrichR databases, connects to the EnrichR API,
#' retrieves the enrichment results, filters them based on the adjusted P-value threshold, and returns the filtered results.
#' The function also prints the number of significant terms found for each database to the console.
#' If no databases are specified, it defaults to using "GO_Biological_Process_2023", "KEGG_2021_Human", and "Reactome_Pathways_2024".
#' 
#' @param genes Character vector of gene symbols.
#' @param dbs Character vector of Enrichr databases.
#' @param alpha Numeric. Adjusted P-value threshold (default 0.05).
#' 
#' @return A list of data frames containing the filtered enrichment results for each specified database.
#' Each data frame includes only the terms with an adjusted P-value less than or equal to the specified alpha threshold.
#' 
#' @export
run_functional_enrichment <- function(genes, dbs = NULL, alpha = 0.05) {

  # 1. Set default databases if none provided
  if (is.null(dbs)) {
    dbs <- c("GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_Pathways_2024")
  }
  # 2. Connect to EnrichR and retrieve results
  cat("Connecting to Enrichr for databases:", paste(dbs, collapse = ", "), "...\n")
  results <- enrichR::enrichr(genes, dbs)
  # 3. Filter results based on adjusted P-value
  filtered_results <- lapply(results, function(df) {
    if (nrow(df) == 0) return(df)
    df[df$Adjusted.P.value <= alpha, ]
  })
  # 4. Print number of significant terms found for each database
  for (db in names(filtered_results)) {
    cat(sprintf("- %s: %d significant terms found.\n", db, nrow(filtered_results[[db]])))
  }

  return(filtered_results)
}
