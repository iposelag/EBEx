# Generate SHAP Contributions for a Single Sample

This function generates SHAP contributions for a single sample using a
DALEX explainer object. It constructs a sample label based on the first
column of the data and the row name, computes the SHAP values, and
optionally saves the results as an RDS file in a structured directory.

## Usage

``` r
generate_shap_contributions(
  explainer,
  data_to_shap,
  target_var,
  sample_idx,
  directory_to_save = NULL
)
```

## Arguments

- explainer:

  DALEX explainer object.

- data_to_shap:

  Dataframe containing the samples.

- target_var:

  Character. Name of the target variable column.

- sample_idx:

  Integer. Index of the sample to explain.

- directory_to_save:

  Character. Optional directory to save the RDS file.

## Value

A list containing the sample label and the SHAP values dataframe.
