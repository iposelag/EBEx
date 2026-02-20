# Split Data into Training and Test Sets

This function loads an expression dataset, splits it into training and
test sets while maintaining class balance, and optionally saves the
resulting datasets as .Rda files. It also prints the class imbalance
percentages for both the training and test sets to the console.

## Usage

``` r
obtain_split_data(
  directory_to_load,
  file_name,
  target_var,
  directory_to_save = NULL,
  verbose = TRUE
)
```

## Arguments

- directory_to_load:

  Directory to load the expression data from.

- file_name:

  Name of the expression data .Rda file (without extension).

- target_var:

  Name of the target variable in the dataset.

- directory_to_save:

  Directory to save the split datasets (optional).

## Value

A list containing the training and test datasets as dataframes.
