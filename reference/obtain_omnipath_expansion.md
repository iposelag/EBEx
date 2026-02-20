# Obtain Omnipath expansion

This function takes a list of seed genes and retrieves their interaction
partners from the Omnipath database. It filters the interactions to
include only those with a curation effort greater than 1 and returns a
unique list of genes that includes both the original seed genes and
their interaction partners. The resulting list of genes can be saved to
a specified directory if the `directory_to_save` parameter is provided.

## Usage

``` r
obtain_omnipath_expansion(genes, directory_to_save = NULL)
```

## Arguments

- genes:

  A character vector of seed gene symbols to expand.

- directory_to_save:

  Directory to save the expanded gene list to (optional).

## Value

A character vector of gene symbols that includes the original seed genes
and their Omnipath interaction partners.
