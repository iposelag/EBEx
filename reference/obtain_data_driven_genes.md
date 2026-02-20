# Obtain Data Driven Genes (DEA + mRMR)

This function combines the results from Differential Expression Analysis
(DEA) and mRMR feature selection to obtain a list of data-driven genes.
It takes in the lists of DEA and mRMR genes, as well as the expression
data and target variable, to determine the final set of genes. The
resulting list can be saved to a specified directory if the
`directory_to_save` parameter is provided.

## Usage

``` r
obtain_data_driven_genes(
  dea_genes = NULL,
  mrmr_genes = NULL,
  expression_data,
  target_var,
  threshold_value = NULL,
  directory_to_load,
  directory_to_save = NULL,
  mrmr_path = NULL
)
```

## Arguments

- dea_genes:

  List of differentially expressed genes (optional).

- mrmr_genes:

  List of mRMR selected genes (optional).

- expression_data:

  Data frame containing the expression data.

- target_var:

  Character string of the target variable column name.

- threshold_value:

  Numeric value for the mRMR score threshold. If NULL, bootstrapping
  will be performed to determine the threshold.

- directory_to_load:

  Directory to load the DEA results from.

- directory_to_save:

  Directory to save the selected data-driven genes to (optional).

- mrmr_path:

  Path to the mRMR executable.

## Value

A character vector of selected data-driven gene names.
