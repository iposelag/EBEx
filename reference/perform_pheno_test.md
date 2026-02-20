# Perform Statistical Test for Gene-Phenotype Association

This function takes a dataframe containing mean SHAP values for a gene
and a phenotypic variable, checks for normality, and automatically
selects the appropriate statistical test (correlation for numeric
variables, t-test or ANOVA for normally distributed groups, and
non-parametric tests otherwise).

## Usage

``` r
perform_pheno_test(data, variable)
```

## Arguments

- data:

  Dataframe containing 'mean_shap' and the target variable.

- variable:

  Character. Name of the phenotypic variable column.

## Value

A tidy dataframe with test results, including p-values, test type, and
normality information.
