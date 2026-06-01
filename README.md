
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!--# <img src="man/figures/logo4.png" align="top" height="100" /> -->

# <img src="man/figures/logo.png" align="top" height="130" style="margin-rght: 20px;"/>

<!-- badges: start -->
<!-- badges: end -->

EBEx (**E**nsemble-**B**ased **Ex**plainable framework) is an
open-source R package implementing a multi-step machine learning
pipeline for disease-relevant gene prioritisation from transcriptomic
data. EBEx is designed for the high-dimensional, low-sample-size regime
that characterises most clinical transcriptomic cohorts, prioritising
interpretability and robustness over model complexity.

## Installation

You can install the development version of EBEx from
[GitHub](https://github.com/iposelag/EBEx) directly in R (\>= 4.1) with:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("iposelag/EBEx")
```

## Dependencies

EBEx integrates `CRAN` and `Bioconductor` packages. `Bioconductor`
dependencies include `Biobase`, `limma`, `sva`, `ComplexHeatmap`, and
`OmnipathR.` Core CRAN dependencies comprise the `tidymodels` ecosystem
(`parsnip`, `recipes`, `workflows`, `tune`, `rsample`, `yardstick`,
`dials`, `themis`), machine learning backends (`ranger`, `kernlab`,
`kknn`, `xgboost`, `glmnet`), explainability tools (`DALEX`,
`DALEXtra`), and additional utilities (`dplyr`, `ggplot2`, `mclust`,
`enrichR`, `RColorBrewer`).

## Getting Started

-   [Pipeline Overview](articles/pipeline_overview.html): conceptual
    description of the EBEx framework and its design rationale
-   [Step 01: Feature Selection](articles/feature_selection_guide.html):
    example of gene list generation from data-driven and knowledge-based
    strategies
-   [Step 02: Classification and Candidate Gene
    Selection](articles/candidate_genes_performance.html): example of
    classifier training, explainability-score aggregation, and candidate
    gene selection

## Citation

Comming soon!

## License

MIT + file LICENSE
