# Perform Statistical Test for Gene-Phenotype Association

This function takes a dataframe containing mean SHAP values for a gene
and a phenotypic variable, checks for normality, and automatically
selects the appropriate statistical test (correlation for numeric
variables, t-test or ANOVA for normally distributed groups, and
non-parametric tests otherwise). It also reports an effect size and,
where applicable, the direction of the association.

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

A tidy dataframe with test results, including p-values, test type,
normality information, `effect_size` (Pearson r / Spearman rho;
difference in means for t-test; Wilcoxon r for Mann-Whitney U;
eta-squared for ANOVA and Kruskal-Wallis) and `direction` (NA for ANOVA
and Kruskal-Wallis).
