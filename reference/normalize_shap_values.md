# Normalize SHAP Values

This function takes a dataframe of SHAP values and normalizes them
within each sample. The normalization is done by dividing the mean SHAP
value of each variable by the sum of absolute mean SHAP values for that
sample, resulting in a 'norm_shap' value that represents the relative
importance of each gene within that sample.

## Usage

``` r
normalize_shap_values(shap_df)
```

## Arguments

- shap_df:

  A dataframe with columns variable_name, contribution, and sample.

## Value

A dataframe with an added 'norm_shap' column.
