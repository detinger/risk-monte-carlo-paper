## Entry point for the 50-instance extension experiment.
## Run from the code/ directory:
##   Rscript run_extended.R ../results
##
## Results land in results/checkpoints_extended/ and results/extended_*.csv

args <- commandArgs(trailingOnly = TRUE)
source("load_functions.R")
r_dir <- normalizePath("R")
load_all_definitions(r_dir)

results_dir <- if (length(args)) args[1] else file.path("..", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
setwd(results_dir)
cat("Writing results to", getwd(), "\n")

run_toplevel(file.path(r_dir, "12_extended_instances.R"))
