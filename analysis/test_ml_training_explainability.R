
#!/usr/bin/env Rscript
######################################
## Ranking of genes based on aggregated SHAP values
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# Required Libraries
library(devtools)
devtools::load_all()
devtools::document()
library(MyPipelinePkg)
library(dplyr)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
set.seed(1234)
ini <- Sys.time()

test_dir <- "test"
target_var <- "dis_condition"
test_output <- file.path(test_dir, "results_data_driven")
ml_models_to_run_vector <- c("rf")
# ml_models_to_run_vector <- c("rf", "knn", "svm_r", "svm_p", "glm", "xgb")

if (!dir.exists(test_output)) dir.create(test_output, recursive = TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------------
# # 2. PARALELIZACIÓN (Tidymodels lo detecta automáticamente)
# # Usaremos los cores asignados por Slurm
# num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", 1))
# cl <- makePSOCKcluster(num_cores)
# registerDoParallel(cl)
# cat("Corriendo con", num_cores, "núcleos...\n")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. PREPROCESSING: Split

# Create train and test sets
cat("Splitting data into train and test sets...\n")
split_data <- obtain_split_data(directory_to_load = file.path(test_dir, "input_data"),
                                file_name = "expression",
                                target_var = target_var,
                                directory_to_save = test_output)

expression_train <- split_data$train
expression_test <- split_data$test

# Save training samples
train_samples <- rownames(expression_train)
write.csv(train_samples, file = file.path(test_output, "train_samples.csv"), row.names = FALSE, col.names = "sample_id")

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. FEATURE SELECTION
genes_list <- run_feature_selection(
  procedure = "alternative",
  alternative_genes = "data_driven",
  directory_to_load = file.path(test_dir, "input_data"),
  directory_to_save = test_dir
)

# Filter training data
col_selection <- c(target_var, genes_list)
expression_train <- expression_train[, colnames(expression_train) %in% col_selection]

## ----------------------------------------------------------------------------------------------------------------------------------------
# 4. CLASSIFICATION PERFORMANCE
# 4.1 Run ML models
cat("Running ML models on candidate gene list...\n")
results <- ML_models(
  data = expression_train,
  target_var = target_var,
  models_to_run = ml_models_to_run_vector,
  directory_to_save = test_output
)

# 4.2 Extract results
cat("Extracting model results...\n")
models_results <- extract_models_results(
  models_to_run = ml_models_to_run_vector,
  results_models = results,
  expression_train = expression_train,
  expression_test = expression_test,
  target_var = target_var
)

# Save list of results
saveRDS(models_results, file = file.path(test_output,"models_results.rds"))

# 4.3 Process data for plotting
cat("Processing results for plotting...\n")
models_results_cv_df <- process_cross_validation_metrics(
  models_results = models_results,
  classifiers = ml_models_to_run_vector
)
models_results_test_df <- process_test_metrics(
  models_results = models_results,
  classifiers = ml_models_to_run_vector
)
models_results_cv_df$classifier <- rename_classifier(models_results_cv_df$classifier)
models_results_cv_df$input_list <- rep("Candidate List", nrow(models_results_cv_df))
models_results_test_df$classifier <- rename_classifier(models_results_test_df$classifier)
models_results_test_df$input_list <- rep("Candidate List", nrow(models_results_test_df))

# Plot classification performance
# Define minmax and normMCC for radar charts
minmax <- c(0.5,0.9,0.1)
metric <- "normMCC"
# Plot Test
cat("Plotting classification performance radar charts...\n")
data_for_plot <- prepare_data_for_radarchart(models_results_test_df, metric, "estimate", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = "radarchart_test_candidate", test_output, minmax[1], minmax[2], minmax[3])
# Cross validation
data_for_plot <- prepare_data_for_radarchart(models_results_cv_df, metric, "mean", minmax[1], minmax[2])
plot_radarchart(data_for_plot, metric, plot_name = "radarchart_cross_validation_candidate", test_output, minmax[1], minmax[2], minmax[3])

## ----------------------------------------------------------------------------------------------------------------------------------------
# 5. EXPLAINABILITY
cat("Calculating SHAP explainability...\n")
shap_results <- calculate_all_shap(ml_models_to_run = ml_models_to_run_vector,
                                   results_models = results,
                                   expression_data = expression_train,
                                   target_var = target_var,
                                   directory_to_save = test_output)

saveRDS(shap_results, file = file.path(test_output,"shap_results.rds"))

## ----------------------------------------------------------------------------------------------------------------------------------------
# END
stopCluster(cl)
fin <- Sys.time()
cat("Pipeline completed!\n")
cat("Total time:", fin - ini, "\n")
