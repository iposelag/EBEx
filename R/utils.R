#' @importFrom dplyr select mutate bind_rows case_when
#' @importFrom purrr map_dfr imap_dfr
#' @importFrom rlang .data
NULL

## --------------------------------------------------------------------------------------------------------------------
# FUNCTIONS TO SAVE DATA

#' Internal helper to handle directory creation and saving text files
#' @keywords internal
save_helper_txt <- function(data, filename, directory, subfolder = "") {

  out_dir <- file.path(directory, subfolder)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write(data, file = file.path(out_dir, filename), ncolumns = 1)
  
  return(invisible(NULL))
}

#' Internal helper to handle directory creation and saving .Rda files
#' @keywords internal
save_helper_rda <- function(data_obj, obj_name, filename, directory, subfolder = "") {
  out_dir <- file.path(directory, subfolder)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # To use save() inside a function, we assign the object to a temporary name
  # in a temporary environment to ensure the .Rda file contains the right object name.
  tmp_env <- new.env()
  assign(obj_name, data_obj, envir = tmp_env)
  save(list = obj_name, file = file.path(out_dir, filename), envir = tmp_env)

  return(invisible(NULL))
}

#' Internal helper to save multiple objects into one .Rda file
#' @keywords internal
save_helper_multi_rda <- function(obj_list, name_list, filename, directory, subfolder = "") {

  cat("Saving multiple .Rda objects to", file.path(directory, subfolder, filename), "...")
  out_dir <- file.path(directory, subfolder)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  # Create a temporary environment
  tmp_env <- new.env()
  # Assign each object its specific name within that environment
  for (i in seq_along(obj_list)) {
    assign(name_list[[i]], obj_list[[i]], envir = tmp_env)
  }
  # Save the list of names from that environment into one file
  save(list = unlist(name_list), file = file.path(out_dir, filename), envir = tmp_env)

  return(invisible(NULL))
}

#' Internal helper to handle directory creation and saving CSV files
#' @keywords internal
save_helper_csv <- function(data, filename, directory, subfolder = "", verbose = TRUE) {

  print_message("Saving CSV to", file.path(directory, subfolder, filename), "...", verbose = verbose)
  out_dir <- file.path(directory, subfolder)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(data, file = file.path(out_dir, filename), row.names = FALSE)

}

#' Internal helper to save .rds files
#' @keywords internal
save_helper_rds <- function(data, filename, directory, subfolder = "", verbose = TRUE) {

  print_message("Saving RDS to", file.path(directory, subfolder, filename), "...", verbose = verbose)
  out_dir <- file.path(directory, subfolder)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(data, file = file.path(out_dir, filename))

}

## --------------------------------------------------------------------------------------------------------------------
# RENAME FUNCTIONS

#' Rename classifier names for display
#' @keywords internal
rename_classifier <- function(x) {
  classifier_mapping <- c(
    rf = "RF",
    svm_r = "SVM-rad",
    svm_p = "SVM-poly",
    knn = "kNN",
    xgb = "XGB",
    xgb_bundle = "XGB",
    glm = "GLM"
  )
  res <- ifelse(x %in% names(classifier_mapping), 
                classifier_mapping[x],
                x)

  return(unname(res))
}

#' Rename input list names for display
#' @keywords internal
rename_input_list <- function(x) {

  # Add your specific mappings here
  input_mapping <- c(
    dea = "DEA",
    mrmr = "mRMR",
    mrmr_30 = "mRMR 30",
    mrmr_76 = "mRMR 76",
    data_driven = "data-driven",
    guildify_data_driven = "data-driven GUILDify",
    omnipath_data_driven = "data-driven OmniPath",
    disease_related = "COPD-related curated",
    guildify_disease_related = "COPD-related curated GUILDify",
    guildify_functional_based_disease_related = "COPD-related curated GUILDify functional based",
    omnipath_disease_related = "COPD-related curated OmniPath",
    disease_related_entire_list = "COPD-related entire list",
    omnipath_intersection = "OmniPath intersection",
    omnipath_union = "OmniPath union",
    candidate = "Candidate List"
  )
  res <- ifelse(x %in% names(input_mapping), 
                input_mapping[x], 
                x)

  return(unname(res))
}

#' Rename Composite Input-Classifier Names
#'
#' @description
#' Takes a name like "data_driven_rf" and converts it to "data-driven RF" 
#' using the existing mapping functions in the package.
#'
#' @param x Character vector of composite names.
#' @param sep Character. The separator used in the original names (default "_").
#' @return A character vector of pretty composite names.
#' @export
rename_input_classifier <- function(x, sep = "_") {
  
  # Lista de clasificadores conocidos para poder separar correctamente
  classifs <- c("rf", "svm_r", "svm_p", "knn", "xgb", "glm", "xgb_bundle")
  
  sapply(x, function(name) {
    # 1. Intentar encontrar qué clasificador está al final del string
    # Usamos una expresión regular para encontrar el sufijo
    found_clf <- NULL
    for (clf in classifs) {
      pattern <- paste0("_", clf, "$")
      if (grepl(pattern, name)) {
        found_clf <- clf
        break
      }
    }
    
    if (!is.null(found_clf)) {
      # 2. Separar la parte del input de la del clasificador
      # Quitamos el "_clasificador" del final
      input_part <- gsub(paste0("_", found_clf, "$"), "", name)
      
      # 3. Traducir ambas partes usando tus funciones
      pretty_input <- rename_input_list(input_part)
      pretty_clf   <- rename_classifier(found_clf)
      
      # 4. Unir con un espacio (ideal para UpSet plots)
      return(paste(pretty_input, pretty_clf))
    } else {
      # Si no detecta un clasificador al final, traduce el string completo
      return(rename_input_list(name))
    }
  }) %>% unname()
}

## --------------------------------------------------------------------------------------------------------------------
# METRICS PROCESSING

#' Process Cross-Validation Metrics from Results List
#' converts the nested list of cross-validation metrics into a long-format data frame
#'
#' @keywords internal
process_cross_validation_metrics <- function(results_models, classifiers) {

  all_metrics <- purrr::map_dfr(classifiers, function(cl) {
    if (is.null(results_models$cross_validation[[cl]])) return(NULL)
    df <- as.data.frame(results_models$cross_validation[[cl]]$model_metrics)
    # Detectar cuál de las dos columnas de error existe
    err_col <- if ("std_err" %in% colnames(df)) "std_err" else "sd"
    # Standardize Tidymodels names
    df %>%
      dplyr::select(metric = .data$.metric, mean = .data$mean, sd = .data[[err_col]]) %>%
      dplyr::mutate(classifier = cl)
  })
  # Ensure normMCC exists
  if (!"normMCC" %in% all_metrics$metric && "mcc" %in% all_metrics$metric) {
    mcc_rows <- all_metrics %>% dplyr::filter(.data$metric == "mcc")
    norm_mcc <- mcc_rows %>%
      dplyr::mutate(metric = "normMCC", mean = (.data$mean + 1) / 2, sd = NA)
    all_metrics <- dplyr::bind_rows(all_metrics, norm_mcc)
  }
  
  return(all_metrics)
}

#' Process Test Metrics from Results List: converts the nested list of test metrics into a long-format data frame
#' @keywords internal
process_test_metrics <- function(results_models, classifiers) {

  all_metrics <- purrr::map_dfr(classifiers, function(cl) {
    if (is.null(results_models$test[[cl]])) return(NULL)
    df <- as.data.frame(results_models$test[[cl]]$model_metrics)
    df %>%
      dplyr::select(metric = .data$.metric, estimate = .data$.estimate) %>%
      dplyr::mutate(classifier = cl)
  })
  # Ensure normMCC exists
  if (!"normMCC" %in% all_metrics$metric && "mcc" %in% all_metrics$metric) {
    mcc_rows <- all_metrics %>% dplyr::filter(.data$metric == "mcc")
    norm_mcc <- mcc_rows %>%
      dplyr::mutate(metric = "normMCC", estimate = (.data$estimate + 1) / 2)
    all_metrics <- dplyr::bind_rows(all_metrics, norm_mcc)
  }
  
  return(all_metrics)
}

#' Process Computation Times
#' @keywords internal
process_times <- function(results_models, unit = "mins") {

  valid_units <- c("secs", "mins", "hours", "days")
  if (!(unit %in% valid_units)) stop("Invalid unit. Choose from 'secs', 'mins', 'hours', or 'days'.")
  purrr::imap_dfr(results_models$times, function(val, cl) {
    data.frame(classifier = cl, time = as.numeric(val, units = unit))
  })

}

## --------------------------------------------------------------------------------------------------------------------
# VERBOSE MESSAGES
#' Print Timestamped Messages
#'
#' @param ... Arguments to be passed to the message.
#' @param verbose Logical. If FALSE, the message is suppressed.
#' Defaults to the global option "EBEx.verbose" or TRUE.
#' @export
print_message <- function(..., verbose = getOption("EBEx.verbose", default = TRUE)) {
  
  # If verbose is FALSE, exit the function immediately without printing
  if (!verbose) return(invisible(NULL))
  
  args <- list(...)
  
  # Your formatting logic
  args <- lapply(args, function(x) {
    if (is.null(x)) return("NULL")
    if (is.table(x)) return(paste(paste(names(x), x, sep = ": "), collapse = ", "))
    if (is.vector(x) && length(x) > 1) return(paste(x, collapse = ", "))
    return(as.character(x))
  })
  
  message_text <- paste(unlist(args), collapse = " ")
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ":    ", message_text, "\n")
  
  return(invisible(NULL))
}

## ----------------------------------------------------------------------------------------------------------------------------------------
# Vectors to define consistent colors
genes_list_colors <- c("DEA" = "#bc6c25",
                  "mRMR" = "#f6bd60",
                  "mRMR & DEA" = "#19787F",
                  "mRMR 30" = "#e9d8a6",
                  "mRMR 76" = "#ee9b00",
                  "data-driven" = "#19787F",
                  "data-driven OmniPath" = "#8EBCBF",
                  "data-driven GUILDify" = "#d2e7d6",
                  "data-driven GUILDify functional based" = "#9b2226",
                  "COPD-related curated" = "#325486",
                  "COPD-related curated GUILDify" = "#cfe2f3",
                  "COPD-related entire list & COPD-related curated GUILDify" = "#cfe2f3",
                  "COPD-related curated & COPD-related entire list & COPD-related curated GUILDify" = "#325486",
                  "COPD-related curated GUILDify functional based" = "#FF6600",
                  "COPD-related curated OmniPath" = "#9AAAC3",
                  "COPD-related entire list" = "#bea9de",
                  "OmniPath intersection" = "#FF9999",
                  "OmniPath union" = "#ffbaba")

ml_models_colors <- c("RF" = "#6b3e26",
               "SVM-rad" = "#ffc5d9",
               "SVM-poly" = "#c2f2d0",
               "kNN" = "#ffcb85",
               "GLM" = "#fdf5c9",
               "XGB" = "#ff6f69")

dis_condition <- c("CTRL" = "#ffdcdb", "COPD" = "#91a8d0")

GOLD_stage <- c("0-At Risk" = "#f6e0b5", "1-Mild COPD" = "#eea990",
                "2-Moderate COPD" = "#aa6f73", "3-Severe COPD" = "#a39193",
                "4-Very Severe COPD" = "#66545e")

sex <- c("1-Male" = "#e1f7d5", "2-Female" = "#c9c9ff")

smoker <- c("1-Current" = "#83adb5", 
            "2-Ever (>100)" = "#c7bbc9",
            "3-Never"="#5e3c58")

age <-  c("(27,35]" = "#ece6ff", "(35,45]" = "#efbbff",
          "(45,55]" = "#d896ff", "(55,65]" = "#be29ec",
         "(65,75]" = "#800080", "(75,91]" = "#660066")

pneumocystis_colonization <- c("Negative" = "#b1cc74",
                              "Positive" = "#eec643")

platform_id <- c("GPL14550" = "#918450",
                 "GPL6480" = "#a41623")

gtex_group_colors <- c(
  Brain = "#8da0cb", Adipose = "#66c2a5", Cardiovascular = "#fc8d62",
  Digestive = "#a6d854", Endocrine = "#ffd92f", Immune = "#e78ac3",
  Musculoskeletal = "#b3b3b3", Skin = "#e5c494", Kidney = "#80b1d3",
  "Reproductive (F)" = "#fb9a99", "Reproductive (M)" = "#cab2d6",
  Bladder = "#fdbf6f", Lung = "#ffff99")

gene_status_colors <- c(
  "DEG" = "#bc6c25",
  "Extended" = "#c1121f",
  "Canonical" = "#0072B2",
  "Disease Related curated" = "#325486",
  "Other" = "grey80"
)

## ----------------------------------------------------------------------------------------------------------------------------------------
# MAPPING VECTORS

#' Default Gene Name Mapping for GTEx
#' @keywords internal
gtex_gene_map <- c(
  "CEP170B" = "KIAA0284", "EVA1A" = "FAM176A", "NIBAN2" = "FAM129B",
  "ASIC4" = "ACCN4", "TENM2" = "ODZ2", "CEMIP" = "KIAA1199",
  "FAM30A" = "KIAA0125", "CCAR2" = "DBC1", "JCHAIN" = "IGJ",
  "HIF1A" = "MOP-1", "MROH5" = "FLJ43860"
)
tissue_fine_map <- c(
  # ---- Adipose ----
  "Adipose_Subcutaneous" = "Adipose - Subcutaneous",
  "Adipose_Visceral_Omentum" = "Adipose - Visceral",
  # ---- Adrenal ----
  "Adrenal_Gland" = "Adrenal Gland",
  # ---- Bladder ----
  "Bladder" = "Bladder",
  # ---- Blood vessels ----
  "Artery_Aorta" = "Artery - Aorta",
  "Artery_Coronary" = "Artery - Coronary",
  "Artery_Tibial" = "Artery - Tibial",
  "Cells_EBV-transformed_lymphocytes" = "EBV-transformed lymphocytes",
  # ---- Brain ----
  "Brain_Amygdala" = "Brain - Amygdala",
  "Brain_Anterior_cingulate_cortex_BA24" = "Brain - Anterior cingulate cortex (BA24)",
  "Brain_Caudate_basal_ganglia" = "Brain - Caudate (basal ganglia)",
  "Brain_Cerebellar_Hemisphere" = "Brain - Cerebellar Hemisphere",
  "Brain_Cerebellum" = "Brain - Cerebellum",
  "Brain_Cortex" = "Brain - Cortex",
  "Brain_Frontal_Cortex_BA9" = "Brain - Frontal Cortex (BA9)",
  "Brain_Hippocampus" = "Brain - Hippocampus",
  "Brain_Hypothalamus" = "Brain - Hypothalamus",
  "Brain_Nucleus_accumbens_basal_ganglia" = "Brain - Nucleus accumbens (basal ganglia)",
  "Brain_Putamen_basal_ganglia" = "Brain - Putamen (basal ganglia)",
  "Brain_Spinal_cord_cervical_c-1" = "Brain - Spinal cord (cervical c-1)",
  "Brain_Substantia_nigra" = "Brain - Substantia nigra",
  # ---- Breast ----
  "Breast_Mammary_Tissue" = "Breast - Mammary Tissue",
  # ---- Cervix ----
  "Cervix_Ectocervix" = "Cervix - Ectocervix",
  "Cervix_Endocervix" = "Cervix - Endocervix",
  # ---- Colon ----
  "Colon_Sigmoid" = "Colon - Sigmoid",
  "Colon_Transverse" = "Colon - Transverse",
  # "Colon_Transverse_Mixed_Cell" = "Colon - Transverse",
  # "Colon_Transverse_Mucosa" = "Colon - Transverse",
  # "Colon_Transverse_Muscularis" = "Colon - Transverse",
  # ---- Esophagus ----
  "Esophagus_Gastroesophageal_Junction" = "Esophagus - Gastroesophageal Junction",
  "Esophagus_Mucosa" = "Esophagus - Mucosa",
  "Esophagus_Muscularis" = "Esophagus - Muscularis",
  # ---- Heart ----
  "Heart_Atrial_Appendage" = "Heart - Atrial Appendage",
  "Heart_Left_Ventricle" = "Heart - Left Ventricle",
  # ---- Kidney ----
  "Kidney_Cortex" = "Kidney - Cortex",
  "Kidney_Medulla" = "Kidney - Medulla",
  # ---- Liver ----
  "Liver" = "Liver",
  # "Liver_Hepatocyte" = "Liver Hepatocyte",
  # "Liver_Mixed_Cell" = "Liver Mixed Cell",
  # "Liver_Portal_Tract" = "Liver Portal Tract",
  # ---- Pancreas ----
  "Pancreas" = "Pancreas",
  # "Pancreas_Acini" = "Pancreas Acini",
  # "Pancreas_Islets" = "Pancreas Islets",
  # "Pancreas_Mixed_Cell" = "Pancreas Mixed Cell",
  # ---- Skin ----
  "Cells_Cultured_fibroblasts" = "Cells - Cultured fibroblasts",
  "Skin_Not_Sun_Exposed_Suprapubic" = "Skin - Not Sun Exposed",
  "Skin_Sun_Exposed_Lower_leg" = "Skin - Sun Exposed",
  # ---- Intestine / stomach ----
  "Small_Intestine_Terminal_Ileum" = "Small Intestine - Terminal Ileum",
  # "Small_Intestine_Terminal_Ileum_Lymphode_Aggregate" = "Small Intestine - Terminal Ileum Lymphode Aggregate",
  # "Small_Intestine_Terminal_Ileum_Mixed_Cell" = "Small Intestine - Terminal Ileum Mixed Cell",
  # ---- Stomach ----
  "Stomach" = "Stomach",
  # "Stomach_Mixed_Cell" = "Stomach",
  # "Stomach_Mucosa" = "Stomach",
  # "Stomach_Muscularis" = "Stomach",
  # ---- Other single tissues ----
  "Fallopian_Tube" = "Fallopian Tube",
  "Lung" = "Lung",
  "Minor_Salivary_Gland" = "Minor Salivary Gland",
  "Muscle_Skeletal" = "Muscle - Skeletal",
  "Nerve_Tibial" = "Nerve - Tibial",
  "Ovary" = "Ovary",
  "Pituitary" = "Pituitary",
  "Prostate" = "Prostate",
  "Spleen" = "Spleen",
  "Testis" = "Testis",
  "Thyroid" = "Thyroid",
  "Uterus" = "Uterus",
  "Vagina" = "Vagina",
  "Whole_Blood" = "Whole Blood"
)

tissue_macro_map <- c(
  # ---- Brain ---- todo el SNC
  "Brain - Amygdala" = "Brain",
  "Brain - Anterior cingulate cortex (BA24)" = "Brain",
  "Brain - Caudate (basal ganglia)" = "Brain",
  "Brain - Cerebellar Hemisphere" = "Brain",
  "Brain - Cerebellum" = "Brain",
  "Brain - Cortex" = "Brain",
  "Brain - Frontal Cortex (BA9)" = "Brain",
  "Brain - Hippocampus" = "Brain",
  "Brain - Hypothalamus" = "Brain",
  "Brain - Nucleus accumbens (basal ganglia)" = "Brain",
  "Brain - Putamen (basal ganglia)" = "Brain",
  "Brain - Spinal cord (cervical c-1)" = "Brain",
  "Brain - Substantia nigra" = "Brain",
  # ---- Adipose ---- tejido adiposo
  "Adipose - Subcutaneous" = "Adipose",
  "Adipose - Visceral" = "Adipose",
  # ---- Cardiovascular ---- corazón + vasos
  "Artery - Aorta" = "Cardiovascular",
  "Artery - Coronary" = "Cardiovascular",
  "Artery - Tibial" = "Cardiovascular",
  "Heart - Atrial Appendage" = "Cardiovascular",
  "Heart - Left Ventricle" = "Cardiovascular",
  # ---- Digestive ---- GI + hígado + páncreas + glándulas salivares
  "Colon - Sigmoid" = "Digestive",
  "Colon - Transverse" = "Digestive",
  "Esophagus - Gastroesophageal Junction" = "Digestive",
  "Esophagus - Mucosa" = "Digestive",
  "Esophagus - Muscularis" = "Digestive",
  "Stomach" = "Digestive",
  "Small Intestine - Terminal Ileum" = "Digestive",
  "Liver" = "Digestive",
  # "Liver Hepatocyte" = "Digestive",
  # "Liver Mixed Cell" = "Digestive",
  # "Liver Portal Tract" = "Digestive",
  "Pancreas" = "Digestive",
  # "Pancreas Acini" = "Digestive",
  # "Pancreas Islets" = "Digestive",
  # "Pancreas Mixed Cell" = "Digestive",
  "Minor Salivary Gland" = "Digestive",
  # ---- Endocrine ---- glandulas hormonales
  "Adrenal Gland" = "Endocrine",
  "Pituitary" = "Endocrine",
  "Thyroid" = "Endocrine",
  # ---- Immune / Blood ---- sangre + linfocitos
  "Whole Blood" = "Immune",
  "EBV-transformed lymphocytes" = "Immune",
  "Spleen" = "Immune",
  # ---- Musculoskeletal ---- músculo + nervios
  "Muscle - Skeletal" = "Musculoskeletal",
  "Nerve - Tibial" = "Musculoskeletal",
  # ---- Reproductive (F) ---- mujer
  "Breast - Mammary Tissue" = "Reproductive (F)",
  "Cervix - Ectocervix" = "Reproductive (F)",
  "Cervix - Endocervix" = "Reproductive (F)",
  "Fallopian Tube" = "Reproductive (F)",
  "Ovary" = "Reproductive (F)",
  "Uterus" = "Reproductive (F)",
  "Vagina" = "Reproductive (F)",
  # ---- Reproductive (M) ---- hombre
  "Prostate" = "Reproductive (M)",
  "Testis" = "Reproductive (M)",
  # ---- Bladder
  "Bladder" = "Bladder",
  # ---- Skin ---- piel
  "Cells - Cultured fibroblasts" = "Skin",
  "Skin - Not Sun Exposed" = "Skin",
  "Skin - Sun Exposed" = "Skin",
  # Kidney
  "Kidney - Cortex" = "Kidney",
  "Kidney - Medulla" = "Kidney",
  # ----  Lung ----
  "Lung" = "Lung"
)