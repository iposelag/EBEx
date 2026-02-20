# Obtain mRMR genes

This function runs the mRMR feature selection algorithm on the provided
expression data and target variable. It can also perform bootstrapping
to determine a score threshold for selecting features. The selected mRMR
genes are saved to a specified directory if the `directory_to_save`
parameter is provided.

## Usage

``` r
obtain_mrmr(
  expression_data,
  target_var,
  threshold_value = NULL,
  directory_to_save = NULL,
  mrmr_path = NULL
)
```

## Arguments

- expression_data:

  Data frame containing the expression data.

- target_var:

  Character string of the target variable column name.

- threshold_value:

  Numeric value for the mRMR score threshold. If NULL, bootstrapping
  will be performed to determine the threshold.

- directory_to_save:

  Directory to save the selected mRMR genes to (optional).

- mrmr_path:

  Path to the mRMR executable.

## Value

A character vector of selected mRMR gene names.
