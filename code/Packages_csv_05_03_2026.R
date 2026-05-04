library(here)

packages_used <- tibble::tribble(
  ~Package,           ~Version,    ~Purpose,
  
  # Core data manipulation
  "tidyverse",        "2.0.0",     "Meta-package loading core data manipulation libraries (dplyr, tidyr, purrr, stringr, ggplot2, readr, forcats, tibble)",
  "dplyr",            "1.1.4",     "Data manipulation: filter, mutate, group_by, summarise, joins",
  "tidyr",            "1.3.1",     "Data reshaping: pivot_longer, pivot_wider, separate, unite",
  "purrr",            "1.0.2",     "Functional iteration: map, map_dfr, walk for applying functions across lists and outcomes",
  "stringr",          "1.5.1",     "String manipulation: str_detect, str_remove, str_squish, str_to_title for cleaning text fields",
  "readr",            "2.1.5",     "Fast CSV reading and writing",
  "forcats",          "1.0.0",     "Factor manipulation for ordered taxonomic levels in figures",
  "tibble",           "3.2.1",     "Modern data frames; tribble() used for inline table construction",
  "rlang",            "1.1.4",     "Tidy evaluation helpers (ensym, as_name, .data) for dynamic column references",
  "here",             "1.0.1",     "Project-relative path management for reproducibility across working directories",
  
  # Soil-specific
  "aqp",              "2.0.3",     "Algorithms for Quantitative Pedology; Munsell-to-sRGB color conversion via parseMunsell()",
  
  # Modeling framework
  "tidymodels",       "1.2.0",     "Meta-package loading modeling workflow libraries (recipes, parsnip, workflows, tune, rsample, yardstick)",
  "recipes",          "1.1.0",     "Feature preprocessing: median imputation, type coercion, factor handling",
  "parsnip",          "1.2.1",     "Unified model specification interface across RF, XGBoost, and other engines",
  "workflows",        "1.1.4",     "Bundles recipes and model specifications for tuning and fitting",
  "tune",             "1.2.1",     "Hyperparameter tuning grid search with cross-validation",
  "rsample",          "1.2.1",     "Resampling: nested cross-validation (5x5) and train/test splits",
  "yardstick",        "1.3.1",     "Classification metrics: accuracy, kap, bal_accuracy, f_meas, top_k_accuracy",
  "dials",            "1.3.0",     "Hyperparameter range specifications for tuning grids",
  
  # Model engines
  "ranger",           "0.16.0",    "Random forest engine used via parsnip::rand_forest(engine = 'ranger')",
  "xgboost",          "1.7.7.1",   "Gradient boosting engine used via parsnip::boost_tree(engine = 'xgboost')",
  
  # Dimensionality reduction and clustering
  "uwot",             "0.2.2",     "UMAP implementation; trains and saves dimensionality reduction models for the umap_hdb model variants",
  "dbscan",           "1.2-0",     "HDBSCAN clustering applied to UMAP-projected feature space",
  
  # Visualization
  "ggplot2",          "3.5.1",     "All publication figures: bar charts, density plots, confusion matrices, hierarchical accuracy panels",
  "patchwork",        "1.2.0",     "Multi-panel figure assembly (e.g., 3x2 confusion matrix grid, master EDA figures)",
  "scales",           "1.3.0",     "Axis scale formatting and color palette utilities",
  
  # Utility
  "beepr",            "2.0.0",     "Audible completion notification at end of long-running model fits"
)

cat("Total packages documented:", nrow(packages_used), "\n")

write.csv(packages_used,
          here("outputs", "reports_05_02", "packages_used.csv"),
          row.names = FALSE)

cat("Saved to:", here("outputs", "reports_05_02", "packages_used.csv"), "\n")