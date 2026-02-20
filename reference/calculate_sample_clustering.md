# Perform Sample Clustering

This function computes a distance matrix based on the provided method
(either "euclidean" or "cosine") and performs hierarchical clustering on
the samples. It then cuts the resulting dendrogram into a specified
number of clusters (k) and returns both the hierarchical clustering
object and the cluster assignments for each sample.

## Usage

``` r
calculate_sample_clustering(mat, method = "cosine", k = 5)
```

## Arguments

- mat:

  Matrix. Expression or SHAP matrix (genes as rows, samples as columns).

- method:

  Character. "euclidean" or "cosine".

- k:

  Integer. Number of clusters to cut.

## Value

A list containing the hierarchical clustering object and the cluster
assignments.
