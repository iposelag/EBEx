# Create Preprocessing Recipe

This function creates a preprocessing recipe for machine learning models
using the `recipes` package. It includes steps for optional Box-Cox
transformation, correlation filtering to remove highly correlated
predictors, downsampling to address class imbalance, and optional
normalization of predictors.

## Usage

``` r
create_recipe(data, target_var, transformation = NULL, normalize = FALSE)
```

## Arguments

- data:

  Data frame containing the training data.

- target_var:

  Character string of the target variable column name.

- transformation:

  Character string specifying the type of transformation to apply (e.g.,
  "boxcox"). Default is NULL (no transformation).

- normalize:

  Logical. Whether to apply normalization to the predictors. Default is
  FALSE.

## Value

A `recipes` recipe object that can be used in a machine learning
workflow.
