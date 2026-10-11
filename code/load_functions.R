# Helpers to use the code in R/ without running its experiments.
#
# Every script in R/ defines functions and then runs its experiment at full
# size when sourced. load_all_definitions() loads only the function
# definitions, constants (ALL_CAPS names), and library() calls of every
# script, so a single experiment can be run on its own (run_revision.R) or
# at a small size (smoke_test.R). Loading draws no random numbers.

R_FILES <- c("01_core_functions.R", "02_illustration.R", "03_additional_strategies.R",
             "04_additional_algorithms.R", "05_worked_example.R", "06_full_scale_run.R",
             "07_full_scale_analysis.R", "08_ofat_sweep.R", "09_targeted_checks.R",
             "10_real_world_data.R", "11_revision_experiments.R", "12_extended_instances.R")

is_definition <- function(e, libraries = TRUE) {
  if (!is.call(e)) return(FALSE)
  if (identical(e[[1]], as.name("library"))) return(libraries)
  if (!(identical(e[[1]], as.name("<-")) || identical(e[[1]], as.name("=")))) return(FALSE)
  rhs_is_fn <- is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))
  lhs_is_const <- is.name(e[[2]]) && grepl("^[A-Z][A-Z0-9_]+$", as.character(e[[2]]))
  rhs_is_fn || lhs_is_const
}

## Function definitions, constants, and (optionally) library() calls of one file.
load_definitions <- function(path, libraries = TRUE) {
  for (e in parse(path, keep.source = FALSE)) if (is_definition(e, libraries)) eval(e, globalenv())
}

## All files, in order. R/10's library(readxl) is skipped if readxl is not installed.
load_all_definitions <- function(r_dir = "R") {
  for (f in R_FILES)
    load_definitions(file.path(r_dir, f),
                     libraries = f != "10_real_world_data.R" || requireNamespace("readxl", quietly = TRUE))
}

## Everything except the definitions, in file order: the script's runs and
## summaries. Optionally restricted to a line range, skipping expressions
## whose source contains one of `skip`.
run_toplevel <- function(path, from = 1, to = Inf, skip = character()) {
  ex <- parse(path, keep.source = TRUE)
  refs <- attr(ex, "srcref")
  for (k in seq_along(ex)) {
    line <- refs[[k]][1]
    if (line < from || line > to || is_definition(ex[[k]])) next
    txt <- paste(as.character(refs[[k]]), collapse = "\n")
    if (any(vapply(skip, grepl, logical(1), x = txt, fixed = TRUE))) next
    tryCatch(eval(ex[[k]], globalenv()), error = function(e)
      stop(sprintf("%s, line %d: %s", basename(path), line, conditionMessage(e)), call. = FALSE))
  }
}
