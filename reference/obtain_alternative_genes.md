# Obtain alternative selection of genes

This function loads an alternative list of genes from a specified text
file. The file is expected to contain gene names separated by commas.
The function also saves the loaded gene list to a specified directory if
the `directory_to_save` parameter is provided.

## Usage

``` r
obtain_alternative_genes(
  alternative_genes,
  directory_to_load,
  directory_to_save = NULL
)
```

## Arguments

- alternative_genes:

  Character string of the filename containing the alternative gene list.

- directory_to_load:

  Directory to load the gene list from.

- directory_to_save:

  Directory to save the loaded gene list to (optional).

## Value

A character vector of gene names loaded from the specified file.
