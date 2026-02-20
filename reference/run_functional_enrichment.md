# Run Functional Enrichment using EnrichR

This function performs functional enrichment analysis using the EnrichR
package. It takes a list of gene symbols and a list of EnrichR
databases, connects to the EnrichR API, retrieves the enrichment
results, filters them based on the adjusted P-value threshold, and
returns the filtered results. The function also prints the number of
significant terms found for each database to the console. If no
databases are specified, it defaults to using
"GO_Biological_Process_2023", "KEGG_2021_Human", and
"Reactome_Pathways_2024".

## Usage

``` r
run_functional_enrichment(genes, dbs = NULL, alpha = 0.05)
```

## Arguments

- genes:

  Character vector of gene symbols.

- dbs:

  Character vector of Enrichr databases.

- alpha:

  Numeric. Adjusted P-value threshold (default 0.05).

## Value

A list of data frames containing the filtered enrichment results for
each specified database. Each data frame includes only the terms with an
adjusted P-value less than or equal to the specified alpha threshold.
