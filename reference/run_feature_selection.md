# Main Function for Feature Selection

This function serves as the main entry point for performing feature
selection based on various procedures. It takes in the necessary
parameters and calls the appropriate helper functions to obtain the
selected gene lists according to the specified procedure.

## Usage

``` r
run_feature_selection(
  procedure,
  expression_data,
  target_var,
  threshold_value = NULL,
  disease_code = NULL,
  dea_genes = NULL,
  mrmr_genes = NULL,
  alternative_genes = NULL,
  directory_to_load = NULL,
  directory_to_save = NULL,
  mrmr_path = NULL
)
```

## Arguments

- procedure:

  The feature selection procedure to use.

- expression_data:

  The expression data frame.

- target_var:

  The target variable name.

- threshold_value:

  Threshold value for selection (if applicable).

- disease_code:

  Disease code for disease-related selections.

- dea_genes:

  List of differentially expressed genes (if applicable).

- mrmr_genes:

  List of mRMR selected genes (if applicable

- alternative_genes:

  List of alternative selected genes (if applicable).

- directory_to_load:

  Directory to load data from (if applicable).

- directory_to_save:

  Directory to save results to (if applicable).

- mrmr_path:

  Path to the mRMR executable.

## Value

A character vector of selected gene names based on the specified
procedure.
