#' @importFrom rsample initial_split training testing
NULL

#' Split Data into Training and Test Sets
#'
#' @description
#' This function loads an expression dataset, splits it into training and test sets while maintaining
#' class balance, and optionally saves the resulting datasets as .Rda files. It also prints the class imbalance
#' percentages for both the training and test sets to the console.
#'
#' @param directory_to_load Directory to load the expression data from.
#' @param file_name Name of the expression data .Rda file (without extension).
#' @param target_var Name of the target variable in the dataset.
#' @param directory_to_save Directory to save the split datasets (optional).
#'
#' @return A list containing the training and test datasets as dataframes.
#'
#' @export
obtain_split_data <- function(directory_to_load, file_name, target_var, directory_to_save = NULL, verbose = TRUE) {

  set.seed(1234)
  # Load the expression data
  file_path <- file.path(directory_to_load, paste0(file_name, ".Rda"))
  if (!file.exists(file_path)) stop("Data file not found: ", file_path)
  obj_name <- load(file_path)
  expression <- get(obj_name)
  # Ensure target variable is a factor
  expression[[target_var]] <- as.factor(expression[[target_var]])
  # Split the data using rsample
  expression_split <- rsample::initial_split(expression, strata = target_var)
  # Extract training and test sets
  expression_train <- rsample::training(expression_split)
  expression_test  <- rsample::testing(expression_split)
  # Use the shared helper for saving .Rda files
  if (!is.null(directory_to_save)) {
    save_helper_rda(expression_train, "expression_train", "expression_train.Rda", directory_to_save, "")
    save_helper_rda(expression_test, "expression_test", "expression_test.Rda", directory_to_save, "")
  }
  print_message("Class imbalance: train set", verbose = verbose)
  print_message(round(prop.table(table(expression_train[[target_var]])) * 100, 3), verbose = verbose)
  print_message("\nClass imbalance: test set", verbose = verbose)
  print_message(round(prop.table(table(expression_test[[target_var]])) * 100, 3), verbose = verbose)
  
  return(list(train = expression_train, test = expression_test))
}