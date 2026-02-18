#!/usr/bin/env Rscript
######################################
## Preprocessing: Split data into train and test sets
## Iria Pose
## 2025
######################################

## ----------------------------------------------------------------------------------------------------------------------------------------
# 0. Setup & libraries
library(EBEx)
# Global parameters
set.seed(1234)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 1. CONFIGURATION
# Relative paths
test_dir <- "test"
target_var <- "dis_condition"
output_dir <- file.path(test_dir, "analysis/preprocessing/")
# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 2. PREPROCESSING: Split data into train and test sets
# Create train and test sets
cat("Splitting data into train and test sets...\n")
split_data <- obtain_split_data(directory_to_load = file.path(test_dir, "raw_data"),
                                file_name = "expression",
                                target_var = target_var,
                                directory_to_save = output_dir)
# Extract train and test sets
expression_train <- split_data$train
expression_test <- split_data$test
# Save training samples
train_samples <- data.frame(sample_id = rownames(expression_train))
write.csv(train_samples, file = file.path(output_dir, "train_samples.csv"), row.names = FALSE)
# Save test samples
test_samples <- data.frame(sample_id = rownames(expression_test))
write.csv(test_samples, file = file.path(output_dir, "test_samples.csv"), row.names = FALSE)

## ----------------------------------------------------------------------------------------------------------------------------------------
# 3. FINISH
cat(sprintf("Data splitting completed. Train and test samples saved. %s\n", output_dir))