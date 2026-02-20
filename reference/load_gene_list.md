# Load a gene list from a file

This function loads a gene list from a specified text file. The file is
expected to contain gene names separated by commas.

## Usage

``` r
load_gene_list(genes_list, base_path = "../COPD/feature_selection/")
```

## Arguments

- genes_list:

  Character string of the filename.

- base_path:

  Base directory. Defaults to "../COPD/raw_data/".

## Value

A character vector of gene names.
