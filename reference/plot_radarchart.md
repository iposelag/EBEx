# Plot Radar Chart

Generates a PDF radar chart with specific color logic for known and
unknown input lists.

## Usage

``` r
plot_radarchart(data, metric, plot_name, output_dir = NULL, min, max, seq)
```

## Arguments

- data:

  List. Output from `prepare_data_for_radarchart`.

- metric:

  Character. Metric name for the title.

- plot_name:

  Character. Name for the output file.

- output_dir:

  Character. Directory to save the PDF.

- min:

  Numeric. Minimum axis value.

- max:

  Numeric. Maximum axis value.

- seq:

  Numeric. Step for axis labels.

## Value

Invisibly returns the plot object.
