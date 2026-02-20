# Run Differential Expression Analysis using Limma

This function performs differential expression analysis using the Limma
package. It constructs a design matrix based on the specified grouping
variable, fits a linear model to the expression data, and computes
contrasts to identify differentially expressed genes. The results are
saved as both RDS and CSV files in the specified output directory. The
function returns the fitted model object containing the DEA results.

## Usage

``` r
run_dea_limma(
  expr_matrix,
  metadata,
  gene_names,
  output_dir,
  group_col,
  contrast,
  design_formula = NULL,
  colnames_design_formula = NULL
)
```

## Arguments

- expr_matrix:

  Matrix or dataframe. Expression data (genes as rows, samples as
  columns).

- metadata:

  Dataframe. Sample metadata.

- gene_names:

  Dataframe or vector. Gene symbols or annotations.

- output_dir:

  Character. Directory to save results.

- group_col:

  Character. Name of the column in metadata for grouping.

- contrast:

  Character. Contrast strings (e.g., "COPD - CTRL").

- design_formula:

  Formula. Optional custom design formula.

- colnames_design_formula:

  Character vector. Optional names for design matrix columns.

## Value

MArrayLM object. The fitted model containing DEA results.
