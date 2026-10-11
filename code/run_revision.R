# Runs only the revision experiments (R/11, DEVELOPMENT_PLAN.md Phase 2),
# without the rest of the pipeline.
#
# Usage, from the code/ directory:
#   Rscript run_revision.R [results_dir]
# results_dir defaults to ../results. Put the original Design A results
# (full_production_results.rds, or the checkpoints_production/ folder) there
# first to run on the paper's instances; see design_a_specs() in R/06.
# About 7-8 hours on 6 cores. Checkpointed: rerun the same command to resume.

args <- commandArgs(trailingOnly = TRUE)
source("load_functions.R")
r_dir <- normalizePath("R")
load_all_definitions(r_dir)

results_dir <- if (length(args)) args[1] else file.path("..", "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
setwd(results_dir)
cat("Writing results to", getwd(), "\n")

## Make full_production available for summarise_scale_balanced in R/11.
## If full_production_results.rds exists here (put the original Design A results
## in results_dir before running), the raw-vs-normalized comparison is included
## in scale_balanced_summary.rds.  Otherwise the summary reports normalized
## results only and the raw comparison can be added later.
if (file.exists("full_production_results.rds")) {
  full_production <- readRDS("full_production_results.rds")
  cat("Loaded full_production from full_production_results.rds\n")
} else {
  full_production <- NULL
  cat("Note: full_production_results.rds not found; raw comparison in scale-balanced summary will be skipped.\n")
}

run_toplevel(file.path(r_dir, "11_revision_experiments.R"))
