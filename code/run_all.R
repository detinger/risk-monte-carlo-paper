# Entry point: sources all scripts in dependency order. The full pipeline
# takes days of compute (see README.md); every long run is checkpointed and
# resumes when rerun.
#
# Run from the code/ directory:  Rscript run_all.R
# (in RStudio, source this file; it finds its own directory).
#
# Output (RDS, CSV, PDF, checkpoint folders) goes to ../results. External
# data for R/10 is read from ../data. Original result files placed in
# ../results are picked up instead of recomputed.
#
# Alternatives: smoke_test.R checks that everything runs (a few minutes);
# run_revision.R runs only the revision experiments (R/11).

if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))  # RStudio

CODE_DIR    <- normalizePath(".")
RESULTS_DIR <- file.path(dirname(CODE_DIR), "results")
DATA_DIR    <- file.path(dirname(CODE_DIR), "data")
stopifnot(file.exists(file.path(CODE_DIR, "R", "01_core_functions.R")))
dir.create(RESULTS_DIR, showWarnings = FALSE)
setwd(RESULTS_DIR)
run_script <- function(f) source(file.path(CODE_DIR, "R", f))

# --- Core definitions (always needed) ---
run_script("01_core_functions.R")
run_script("03_additional_strategies.R")
run_script("04_additional_algorithms.R")

# --- Illustration and worked example on one small instance (§8–9.5) ---
# 05 must follow 02 directly: its rescoring draws continue the random number
# stream that 02 leaves, so any other order changes the §9.4 results.
run_script("02_illustration.R")
run_script("05_worked_example.R")

# --- Full-scale Design A run (§10) ---
# Checkpointed parallel run; produces full_production_results.rds
run_script("06_full_scale_run.R")

# --- OFAT Design B sweep (§11.1–11.3) ---
# Produces ofat_results.rds. Runs before 07, which needs config_tag() from here.
run_script("08_ofat_sweep.R")

# --- Design A analysis and validation (§10.7–10.8) ---
# Requires full_production (06) and ofat_results (08)
run_script("07_full_scale_analysis.R")

# --- Targeted checks (§11.4–11.7) ---
# τ sweep, config check, scalarization check, instance-level inference
run_script("09_targeted_checks.R")

# --- Real-world data comparison (§12) ---
# Requires PSPLIB j30/j60, MineLib KD and ReleasePlanner data files in
# ../data (see R/10_real_world_data.R), and the readxl package
run_script("10_real_world_data.R")

# --- Revision experiments (DEVELOPMENT_PLAN.md, Phase 2) ---
# Placebo, critical-path proxy, equal wall-clock, probe-size sweep, zero-range audit
run_script("11_revision_experiments.R")
