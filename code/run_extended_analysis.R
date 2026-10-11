## run_extended_analysis.R
## Runs summarise_extended() on the completed extended_main_results.rds.
## Called AFTER run_extended.R has completed and saved extended_main_results.rds.
## Usage: cd code && Rscript run_extended_analysis.R [results_dir]

args <- commandArgs(trailingOnly = TRUE)
source("load_functions.R")
r_dir <- normalizePath("R")
load_all_definitions(r_dir)

results_dir <- if (length(args)) args[1] else file.path("..", "results")
rds_path <- file.path(results_dir, "extended_main_results.rds")

if (!file.exists(rds_path)) {
  stop("extended_main_results.rds not found in: ", results_dir,
       "\nRun run_extended.R first.")
}

cat("Loading extended_main_results.rds ...\n")
ext_main <- readRDS(rds_path)
cat(sprintf("Loaded: %d rows\n", nrow(ext_main)))

setwd(results_dir)
ext_specs <- make_ext_specs()
summarise_extended(ext_main, ext_specs)
