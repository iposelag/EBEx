
#!/usr/bin/env Rscript
######################################
## EBEx Pipeline: ML training and SHAP explainability
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(devtools)
devtools::load_all()
devtools::document()
library(EBEx)
library(dplyr)
library(optparse)

# Define CLI options
option_list <- list(
  make_option(c("-p", "--procedure"), type="character", default="alternative", 
              help="Feature selection procedure [default %default]"),
  make_option(c("-g", "--genes"), type="character", default="candidate_genes", 
              help="Name of the gene list file [default %default]"),
  make_option(c("-t", "--target"), type="character", default="dis_condition", 
              help="Target variable [default %default]"),

  # Optional parameters for specific procedures (e.g., threshold for DEA, disease code for disease-related genes)
  make_option(c("--threshold"), type="numeric", default=NULL, help="Numeric threshold (for mRMR)"),
  make_option(c("--disease"), type="character", default=NULL, help="Disease code (DisGeNET)"),
  make_option(c("--dea_list"), type="character", default=NULL, help="DEA gene list name"),
  make_option(c("--mrmr_list"), type="character", default=NULL, help="mRMR gene list name"),
)

opt <- parse_args(OptionParser(option_list=option_list))

# Global parameters
set.seed(1234)
ini <- Sys.time()

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
# Parameters
procedure         <- opt$procedure
alternative_genes <- opt$genes
target_var        <- opt$target
threshold_value   <- opt$threshold # NULL
disease_code      <- opt$disease # NULL
dea_genes         <- opt$dea_list
mrmr_genes        <- opt$mrmr_list
# Relative paths
base_dir <- "test"
preprocessing_dir <- file.path(base_dir, "analysis/preprocessing")
features_dir <- file.path(base_dir, "feature_selection")
output_dir <- file.path(base_dir, paste0("results_", features_list))
# Model settings
ml_models_to_run_vector <- c("rf", "knn", "svm_r", "svm_p", "glm", "xgb")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

print_message("================================================================")
print_message("STARTING ANALYSIS FOR:", alternative_genes)
print_message("PROCEDURE:", procedure)
if(!is.null(threshold_value)) print_message("THRESHOLD:", threshold_value)
if(!is.null(disease_code))    print_message("DISEASE CODE:", disease_code)
print_message("OUTPUT DIRECTORY:", output_dir)
print_message("================================================================")

## ----------------------------------------------------------------------------------------------------------------------------------------
# # 2. PARALELIZACIÓN PENDIENTE!!!
# # Usaremos los cores asignados por Slurm
# num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
# cl <- makePSOCKcluster(num_cores)
# registerDoParallel(cl)
# cat("Corriendo con", num_cores, "núcleos...\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. PREPROCESSING: Split
# 2.1. Load train and test sets
print_message("Loading train and test sets...")
load(file.path(preprocessing_dir, "expression_train.Rda"))
load(file.path(preprocessing_dir, "expression_test.Rda"))
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. FEATURE SELECTION
genes_list <- run_feature_selection(
  procedure = procedure,
  expression_data = expression_train,
  target_var = target_var,
  threshold_value = threshold_value,
  disease_code = disease_code,
  dea_genes = dea_genes,
  mrmr_genes  = mrmr_genes,
  alternative_genes = alternative_genes,
  directory_to_load = features_dir,
  directory_to_save = features_dir
)
# Filter training data
col_selection <- c(target_var, genes_list)
expression_train <- expression_train[, colnames(expression_train) %in% col_selection]
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. CLASSIFICATION PERFORMANCE
# 4.1 Run ML models
print_message("Running ML models on candidate gene list...")
results <- ML_models(
  data = expression_train,
  target_var = target_var,
  models_to_run = ml_models_to_run_vector,
  directory_to_save = test_output
)
# 4.2 Extract results
print_message("Extracting model results...")
models_results <- extract_models_results(
  models_to_run = ml_models_to_run_vector,
  results_models = results,
  expression_train = expression_train,
  expression_test = expression_test,
  target_var = target_var
)
# Save list of results
saveRDS(models_results, file = file.path(output_dir,"models_results.rds"))
# 4.3 Process data for plotting
print_message("Processing results for plotting...")
models_results_cv_df <- process_cross_validation_metrics(
  results_models = models_results,
  classifiers = ml_models_to_run_vector
)
models_results_test_df <- process_test_metrics(
  results_models = models_results,
  classifiers = ml_models_to_run_vector
)
models_results_cv_df$classifier <- rename_classifier(models_results_cv_df$classifier)
models_results_cv_df$input_list <- rep("Candidate List", nrow(models_results_cv_df))
models_results_test_df$classifier <- rename_classifier(models_results_test_df$classifier)
models_results_test_df$input_list <- rep("Candidate List", nrow(models_results_test_df))
# 4.4. Plot classification performance
# Define minmax and normMCC for radar charts
minmax <- c(0.5,0.9,0.1)
metric <- "normMCC"
# Plot Test
print_message("Plotting classification performance radar charts...")
data_for_plot <- prepare_data_for_radarchart(models_results_test_df, metric, "estimate", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = paste0("radarchart_test_", procedure), ml_models_dir, minmax[1], minmax[2], minmax[3])
# Cross validation
data_for_plot <- prepare_data_for_radarchart(models_results_cv_df, metric, "mean", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = paste0("radarchart_cross_validation_", procedure), ml_models_dir, minmax[1], minmax[2], minmax[3])
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. EXPLAINABILITY
print_message("Calculating SHAP explainability...")
shap_results <- calculate_all_shap(ml_models_to_run = ml_models_to_run_vector,
                                   results_models = results,
                                   expression_data = expression_train,
                                   target_var = target_var,
                                   directory_to_save = test_output)
saveRDS(shap_results, file = file.path(test_output,"shap_results.rds"))
print_message("________________________________")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 6. FINISH
stopCluster(cl)
fin <- Sys.time()
print_message("Pipeline completed!")
print_message("Total time:", fin - ini)