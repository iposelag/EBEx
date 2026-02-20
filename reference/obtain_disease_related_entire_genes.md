# Obtain disease-related entire genes

This function loads the entire list of disease-related genes from a
specified file. The file is expected to be a TSV file containing a
summary of gene-disease associations for the given disease code. The
function extracts the gene names from the file and saves the list to a
specified directory if the `directory_to_save` parameter is provided.

## Usage

``` r
obtain_disease_related_entire_genes(
  disease_code,
  directory_to_load,
  directory_to_save = NULL
)
```

## Arguments

- disease_code:

  Disease code for which to obtain the gene list.

- directory_to_load:

  Directory to load the gene list from.

- directory_to_save:

  Directory to save the loaded gene list to (optional).

## Value

A character vector of gene names related to the specified disease.
