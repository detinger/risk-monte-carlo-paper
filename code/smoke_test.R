# Smoke test: runs every stage of the pipeline at a tiny size, in the order
# of run_all.R, to check that the code parses, each file's dependencies are
# defined before it runs, and the analysis scripts run end to end on output
# of the right shape. Sections 8 and 9.4 run at full size and are checked
# against the outputs stored in the ODD HTML (regression check). Other
# numbers it prints mean nothing. R/10 is only partly run, since it needs
# external data.
#
# Usage, from the code/ directory (about 6 minutes):
#   Rscript smoke_test.R [output_dir]
# Output files (CSV, RDS, checkpoints, plots.pdf) go to output_dir, by
# default a fresh temporary directory.

args <- commandArgs(trailingOnly = TRUE)
code_dir <- normalizePath(".")
out_dir <- if (length(args)) args[1] else file.path(tempdir(), "smoke")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- normalizePath(out_dir)
src <- function(f) file.path(code_dir, "R", f)
stopifnot(file.exists(src("01_core_functions.R")))
source(file.path(code_dir, "load_functions.R"))  # load_definitions(), run_toplevel()
setwd(out_dir)
pdf("plots.pdf")
n_cores <- max(1, parallel::detectCores() - 1)

stage_times <- c()
stage <- function(name, expr) {
  cat("\n\n########", name, "########\n")
  t0 <- proc.time()[3]
  force(expr)
  stage_times[name] <<- round(proc.time()[3] - t0, 1)
  cat(sprintf("######## %s: OK (%.1f s)\n", name, stage_times[name]))
}
expect_equal <- function(got, want, what, tol = 1e-6) {
  if (!isTRUE(all.equal(unname(got), unname(want), tolerance = tol)))
    stop(sprintf("%s: got %s, expected %s", what, paste(signif(got, 7), collapse = " "),
                 paste(signif(want, 7), collapse = " ")), call. = FALSE)
  cat("check passed:", what, "\n")
}

## Tiny Design A settings, shared by every step that has to reproduce Design A.
TINY <- list(n = 20, must_k = 3, budget_fraction = 0.5, deadline_quantile = 0.4,
             eval_budget = 200, reps = 2, risk_alpha = 0.25, nsim_probe = 10, nsim_eval = 20,
             nsim_posthoc = 100, ga_pop = 20, tabu_neigh = 5, moead_K = 20)
tiny_specs <- lapply(1:4, function(i) list(seed = 100 + i, density = c(0.06, 0.12, 0.18, 0.24)[i]))


stage("parse all files", {
  for (f in list.files(file.path(code_dir, "R"), pattern = "\\.R$")) invisible(parse(src(f)))
})

stage("01, 03, 04 definitions", {
  for (f in c("01_core_functions.R", "03_additional_strategies.R", "04_additional_algorithms.R"))
    load_definitions(src(f))
})

stage("repair trace draws no random numbers", {
  inst <- list(acts = generate_activities(30, 2, seed = 3), dag = generate_dag(30, 0.1, seed = 3))
  topo <- topo_order(inst$dag$preds); must <- which(inst$acts$must)
  budget <- 0.5 * sum(inst$acts$cost)
  rep_fn <- function(strategy, trace) { set.seed(9)
    apply_repair_strategy(rep(1L, 30), inst$acts, inst$dag$preds, inst$dag$succ, topo, budget, must,
                          deadline = 120, strategy = strategy, nsim_probe = 10, trace = trace) }
  for (s in c("risk_aware", "risk_shuffled", "critical_path")) {
    tr <- new_cut_trace()
    stopifnot(identical(rep_fn(s, NULL), rep_fn(s, tr)), tr$cuts > 0,
              cost_roi(rep_fn(s, NULL), inst$acts)["cost"] <= budget)
    cat(sprintf("%-14s cuts %d, zero range %d, differs from EFF %d\n", s, tr$cuts, tr$zero_range, tr$differs_from_eff))
  }
})

stage("02 illustration (as is, checked against the ODD output)", {
  run_toplevel(src("02_illustration.R"))
  hv <- aggregate(hypervolume ~ repair + objective, data = design, FUN = mean)$hypervolume
  expect_equal(hv, c(9.790873, 10.413667, 0.011115, 0.011115), "§8.4 mean hypervolumes")
  expect_equal(hypervolume_2d(ref$pareto_front), 11.2851, "§8.3 reference hypervolume", tol = 1e-5)
})

stage("05 worked example (as is, checked against the ODD output)", {
  load_definitions(src("05_worked_example.R"))
  run_toplevel(src("05_worked_example.R"))
  hv <- with(aggregate(hypervolume ~ algorithm + repair, data = full_res[full_res$placement == "multi", ], FUN = mean),
             setNames(hypervolume, paste(algorithm, repair)))
  expect_equal(hv[c("NSGA2 efficiency_only", "SA risk_aware", "MOEAD risk_aware_random", "Random penalty")],
               c(9.960644, 2.458733, 6.565733, 2.1783), "§9.4 mean hypervolumes")
})

stage("06 definitions and top level (full-size runs replaced by tiny Design A)", {
  load_definitions(src("06_full_scale_run.R"))
  run_toplevel(src("06_full_scale_run.R"), skip = c("for (algo in names(ALGO_SCOPE))",
                                                "full_production <- run_production_batch(",
                                                "saveRDS(full_production"))
  full_production <<- do.call(run_production_batch, c(
    list(instance_specs = tiny_specs, n_cores = n_cores, checkpoint_dir = "checkpoints_production"), TINY))
  combos <- sum(lengths(ALGO_SCOPE)) * length(REPAIRS)
  stopifnot(nrow(full_production) == 4 * TINY$reps * combos)
  repaired <- full_production[full_production$repair != "penalty", ]
  stopifnot(all(repaired$feasible_fraction == 1), all(repaired$hypervolume[repaired$placement == "multi"] > 0))
  cat("rows:", nrow(full_production), "| all repair-based runs feasible\n")
})

stage("08 Design B OFAT, tiny", {
  load_definitions(src("08_ofat_sweep.R"))
  OFAT_FACTORS <<- list(
    n                 = list(baseline = 20,   levels = c(14, 30)),
    density           = list(baseline = 0.10, levels = c(0.05, 0.10)),
    budget_fraction   = list(baseline = 0.50, levels = c(0.35, 0.50)),
    deadline_quantile = list(baseline = 0.40, levels = c(0.40, 0.55)),
    risk_alpha        = list(baseline = 0.25, levels = c(0, 0.25)))
  ofat_results <<- run_full_ofat(runs_per_setting = 3, eval_budget = 150, n_cores = n_cores)
  stopifnot(nrow(ofat_results) == 10 * 3 * length(OFAT_ALGOS) * length(OFAT_REPAIRS))
  ## The two baseline cells that are levels of their sweeps are the same runs.
  cell <- function(f, lv) { d <- ofat_results[ofat_results$factor == f & ofat_results$level == lv, ]
                            d$hypervolume[order(d$run, d$algorithm, d$repair)] }
  stopifnot(identical(cell("density", 0.10), cell("budget_fraction", 0.50)))
  cat("check passed: the density 0.10 and budget_fraction 0.50 cells are identical runs\n")
  run_toplevel(src("08_ofat_sweep.R"), skip = c("run_full_ofat(runs_per_setting", "saveRDS(ofat_results"))
})

stage("07 Design A analysis (validation at tiny size)", {
  load_definitions(src("07_full_scale_analysis.R"))
  run_toplevel(src("07_full_scale_analysis.R"),
               skip = c("validate_against_reference(seeds", "saveRDS(validation", "summarise_recovery(validation"))
  validation <- validate_against_reference(seeds = 8001:8002, n = 12, reps = 1, eval_budget = 200,
                                           ga_pop = 20, tabu_neigh = 5, moead_K = 20, nsim_probe = 10,
                                           nsim_eval = 20, nsim_ref = 100, n_cores = 2)
  summarise_recovery(validation)
})

stage("09 targeted checks, tiny (incl. reproduction of Design A rows)", {
  load_definitions(src("09_targeted_checks.R"))
  inst <- make_instance(TINY$n, seed = 7001, density = 0.10)
  p <- derive_algo_params(150, tabu_neigh = 5)
  for (tau in c(0.1, 8)) {
    res <- run_tabu_with_temperature("efficiency_random", 2L, inst, p, risk_alpha = 0.25,
                                     nsim_probe = 10, nsim_eval = 20, seed = 1000, temperature = tau)
    cat("tau =", tau, " front size:", length(res$pareto_idx), "\n")
  }
  config_check <<- run_config_check(full_production, reps = 1, n_cores = n_cores, eval_budget = 150,
                                    ga_pop = 20, tabu_neigh = 5, moead_K = 20, nsim_probe = 10,
                                    nsim_eval = 20, nsim_posthoc = 100, n = TINY$n)
  scalarization_check <<- do.call(run_scalarization_check, c(
    list(full_production = full_production, reps = TINY$reps, verify_reps = TINY$reps, n_cores = n_cores),
    TINY[setdiff(names(TINY), "reps")]))
  ## Runs the summaries, including the check that raw tabu runs reproduce the
  ## stored Design A rows exactly; the full-size runs are skipped.
  run_toplevel(src("09_targeted_checks.R"), skip = c("for (tau in tau_values)",
                                                 "config_check <- run_config_check(full_production, reps = 10",
                                                 "saveRDS(config_check",
                                                 "scalarization_check <- run_scalarization_check(full_production, reps = 10",
                                                 "saveRDS(scalarization_check"))
})

stage("10 real-world data (functions only; data files not available)", {
  DATA_DIR <- "."   # default needed before load_definitions evaluates RP_DIR
  load_definitions(src("10_real_world_data.R"), libraries = FALSE)
  st <- t(sapply(tiny_specs, function(s) network_stats(generate_dag(50, s$density, s$seed)$preds)))
  stopifnot(all(st[, "arcs"] >= 49))  # every non-terminal node has a successor
  print(round(st, 3))
})

stage("11 revision experiments, tiny", {
  load_definitions(src("11_revision_experiments.R"))
  arms <- REVISION_MAIN_ARMS; arms$nsim_probe <- TINY$nsim_probe
  common <- list(eval_budget = TINY$eval_budget, nsim_eval = TINY$nsim_eval, nsim_posthoc = TINY$nsim_posthoc,
                 n = TINY$n, tabu_neigh = TINY$tabu_neigh, n_cores = n_cores)
  revision_main <<- do.call(run_revision_experiments, c(list(tiny_specs, arms, reps = TINY$reps), common))
  revision_probe <<- do.call(run_revision_experiments, c(list(tiny_specs, REVISION_PROBE_ARMS, reps = 1), common))
  stopifnot(nrow(revision_main) == 4 * TINY$reps * nrow(arms))
  ## Arms that repeat a Design A strategy follow the Design A search: the best
  ## value found does not depend on the rescoring draws, so it must match.
  for (a in c("EFF", "EFF-SM", "RISK")) for (pl in c("multi", "single")) {
    rv <- revision_main[revision_main$arm == a & revision_main$placement == pl, ]
    if (!nrow(rv)) next
    da <- full_production[full_production$algorithm == "Tabu" & full_production$placement == pl &
                          full_production$repair == rv$repair[1], ]
    mm <- merge(rv[, c("instance", "rep", "best_roi")], da[, c("instance", "rep", "best_roi")],
                by = c("instance", "rep"))
    stopifnot(nrow(mm) == nrow(rv))
    expect_equal(mm$best_roi.x, mm$best_roi.y, sprintf("%s %s follows the Design A search", a, pl), tol = 1e-12)
  }
  summarise_revision(revision_main, revision_probe, n_boot = 200)
})

stage("11 scale-balanced experiments, tiny (§13.5)", {
  ## run_sa_normalized, run_moead_normalized, run_one_budgeted_normalized,
  ## run_scale_balanced_block, run_scale_balanced_experiments, summarise_scale_balanced.
  sb_common <- list(eval_budget = TINY$eval_budget, nsim_eval = TINY$nsim_eval,
                    nsim_posthoc = TINY$nsim_posthoc, n = TINY$n, must_k = TINY$must_k,
                    tabu_neigh = TINY$tabu_neigh, ga_pop = TINY$ga_pop, moead_K = TINY$moead_K,
                    nsim_probe = TINY$nsim_probe, n_cores = n_cores,
                    checkpoint_dir = "checkpoints_scale_balanced_smoke")
  sb_results <<- do.call(run_scale_balanced_experiments, c(list(tiny_specs, reps = TINY$reps), sb_common))
  sb_algos <- setdiff(names(ALGO_SCOPE), "GA_single")
  sb_algos <- sb_algos[vapply(sb_algos, function(a) "multi" %in% ALGO_SCOPE[[a]], logical(1))]
  expected_rows <- length(tiny_specs) * TINY$reps * length(sb_algos) * length(REPAIRS)
  stopifnot(nrow(sb_results) == expected_rows)
  stopifnot(all(sb_results$feasible_fraction[sb_results$repair != "penalty"] == 1))
  cat("rows:", nrow(sb_results), "| algorithms:", paste(sort(unique(sb_results$algorithm)), collapse = " "), "\n")
  summarise_scale_balanced(sb_results, full_production, n_boot = 50)
})

stage("load_all_definitions() (used by run_revision.R)", load_all_definitions(file.path(code_dir, "R")))

invisible(dev.off())
cat("\n\nAll stages passed. Output in", out_dir, "\n")
print(stage_times)
