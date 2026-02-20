# Obtain curated disease genes

This function loads a curated list of disease-related genes from a
specified file. The file is expected to be a TSV file containing
gene-disease associations with a "curated" curation effort for the given
disease code. The function extracts the gene symbols from the file and
saves the list to a specified directory if the `directory_to_save`
parameter is provided.

## Usage

``` r
obtain_disease_related_curated_genes(
  disease_code,
  directory_to_load,
  directory_to_save = NULL
)
```

## Arguments

- disease_code:

  Disease code for which to obtain the curated gene list.

- directory_to_load:

  Directory to load the gene list from.

- directory_to_save:

  Directory to save the loaded gene list to (optional).

## Value

A character vector of curated gene symbols related to the specified
disease.
