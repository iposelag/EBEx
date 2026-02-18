#!/usr/bin/env Rscript
######################################
## Differential Expression Analysis using limma (DEA)
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(MyPipelinePkg)
library(Biobase)
# Global parameters
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
# Relative paths
raw_data_dir <- "test/raw_data/"
output_dir <- "test/analysis/DEA/"

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. DATA PREPARATION (project-specific)
# Load ExpressionSet
cat("Loading ExpressionSet...\n")
eset_path <- file.path(raw_data_dir, "eset.Rda")
load(eset_path)  # Load the ExpressionSet object
# Relevel condition factor and extract expression, metadata, and gene names
eset$dis_condition <- relevel(factor(eset$dis_condition), ref = "CTRL")
expression <- Biobase::exprs(eset)
# Subset training samples
cat("Subsetting training samples...\n")
samples_train <- read.csv(file.path("test/analysis/preprocessing", "train_samples.csv"))
expression <- expression[,samples_train$sample_id]
metadata <- Biobase::pData(eset)
metadata   <- metadata[samples_train$sample_id, ]
fdata <- Biobase::fData(eset)
gene_names <- data.frame(gene_symbol = fdata$GENE_SYMBOL)
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. MODEL DEFINITION
cat("Defining design formula for DEA...\n")
design_formula <- formula(~ 0 + dis_condition + platform_id + sex)
colnames_design <- c("CTRL", "COPD", "GPL6480", "Female")
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. RUN DEA
cat("Running DEA...\n")
dea_fit <- run_dea_limma(
  expr_matrix = expression, 
  metadata = metadata, 
  gene_names = gene_names, 
  output_dir = output_dir, 
  group_col = "dis_condition",
  design_formula = design_formula, 
  colnames_design_formula = colnames_design,
  contrast = "COPD - CTRL"
)
# Generate Venn diagrams, P-value histograms, and MD plots for DEA results
cat("Generating DEA diagnostic plots...\n")
dea_results <- dea_analysis(
  dea_fit, 
  pv = 0.01, 
  lfc = log2(1.5), 
  coef = 1, 
  output_dir = output_dir
)
cat("________________________________\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. FINISH
cat("DEA analysis completed. Results and plots saved in:", output_dir, "\n")