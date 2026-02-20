# Select Candidate Genes by Threshold

This function filters the original data based on a user-defined
threshold for the specified column. It returns a sorted dataframe of
selected genes and prints the number and percentage of genes selected.

## Usage

``` r
get_candidate_genes(analysis_list, threshold)
```

## Arguments

- analysis_list:

  List. The output from `compute_density_analysis`.

- threshold:

  Numeric. The threshold decided by the user.

## Value

A sorted dataframe of selected genes.
