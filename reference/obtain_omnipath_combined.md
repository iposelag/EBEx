# Combined Omnipath Expansion (Intersection/Union)

This function combines the Omnipath expansions of both disease-related
and data-driven gene lists using either intersection or union
operations. It takes in the necessary parameters to obtain the
disease-related and data-driven gene lists, performs the Omnipath
expansion for each, and then combines the results based on the specified
operation. The final list of genes can be saved to a specified directory
if the `directory_to_save` parameter is provided.

## Usage

``` r
obtain_omnipath_combined(
  expression_data,
  target_var,
  disease_code,
  threshold_value,
  dea_genes,
  mrmr_genes,
  directory_to_load,
  directory_to_save = NULL,
  operation = "intersection",
  mrmr_path = NULL
)
```

## Arguments

- expression_data:

  Data frame containing the expression data.

- target_var:

  Character string of the target variable column name.

- disease_code:

  Disease code for disease-related selections.

- threshold_value:

  Numeric value for the mRMR score threshold. If NULL, bootstrapping
  will be performed to determine the threshold.

- dea_genes:

  List of differentially expressed genes (optional).

- mrmr_genes:

  List of mRMR selected genes (optional).

- directory_to_load:

  Directory to load the necessary data from.

- directory_to_save:

  Directory to save the combined gene list to (optional).

- operation:

  Character string specifying the combination operation: "intersection"
  or "union".

- mrmr_path:

  Path to the mRMR executable.

## Value

A character vector of combined gene names based on the specified
operation.
