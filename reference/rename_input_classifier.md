# Rename Composite Input-Classifier Names

Takes a name like "data_driven_rf" and converts it to "data-driven RF"
using the existing mapping functions in the package.

## Usage

``` r
rename_input_classifier(x, sep = "_")
```

## Arguments

- x:

  Character vector of composite names.

- sep:

  Character. The separator used in the original names (default "\_").

## Value

A character vector of pretty composite names.
