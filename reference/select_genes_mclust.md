# Select Candidate Genes based on Mclust Clusters

This function identifies candidate genes based on the clusters defined
by the Mclust model. It determines a threshold based on the minimum
score of the selected clusters and filters the genes accordingly. The
resulting candidate genes are sorted by their SHAP scores and returned
as a dataframe. The function also prints the defined threshold and the
number of selected genes to the console.

## Usage

``` r
select_genes_mclust(mclust_list, sel_groups)
```

## Arguments

- mclust_list:

  List. Output from `run_mclust_analysis`.

- sel_groups:

  Numeric vector. The cluster IDs chosen as "high SHAP".

## Value

A sorted dataframe of selected genes.
