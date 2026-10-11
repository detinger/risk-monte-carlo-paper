# Full-scale Design A run (ODD §10.1–10.6)
# derive_algo_params(), make_instance(), parallel checkpointed execution.
# Requires: source("04_additional_algorithms.R")


## 10.1 Equal evaluation budgets across algorithms

derive_algo_params <- function(eval_budget, ga_pop = 100, tabu_neigh = 10, moead_K = 100) {
  list(
    eval_budget = eval_budget,
    pop_size    = ga_pop,
    n_gen_ga    = max(1, floor((eval_budget - ga_pop) / ga_pop)),
    random_eval = eval_budget,
    sa_iter     = eval_budget,                              # 1 evaluation / iteration
    tabu_iter   = max(1, floor(eval_budget / tabu_neigh)),   # tabu_neigh evaluations / iteration
    tabu_neigh  = tabu_neigh,
    moead_K     = moead_K,
    moead_gen   = max(1, floor(eval_budget / moead_K))       # K evaluations / generation
  )
}


## 10.2 A portable deadline

## make_instance() is defined in R/01 (section 7.1), since section 9.5 uses it too.

run_one_budgeted <- function(algo, strategy, m, inst, p, risk_alpha, nsim_probe, nsim_eval, seed) {
  with(inst, {
    repfn <- function(x) apply_repair_strategy(x, acts, preds, succ, topo, budget, must_idx,
                                                deadline, strategy = strategy,
                                                risk_alpha = risk_alpha, nsim_probe = nsim_probe)
    evalfn <- make_evaluate(m, acts, preds, topo, budget, deadline, nsim_eval = nsim_eval,
                             penalize = (strategy == "penalty"))
    n <- nrow(acts)
    switch(algo,
      Random    = run_random_search(n, p$random_eval, repfn, evalfn, seed = seed),
      GA_single = run_evolution(n, p$pop_size, p$n_gen_ga, 0.7, 0.15, repfn, evalfn, seed = seed),
      NSGA2     = run_evolution(n, p$pop_size, p$n_gen_ga, 0.7, 0.15, repfn, evalfn, seed = seed),
      SA        = run_sa(n, p$sa_iter, repfn, evalfn, m = m, seed = seed),
      Tabu      = run_tabu(n, p$tabu_iter, repfn, evalfn, m = m,
                            neighbourhood_size = p$tabu_neigh, seed = seed),
      MOEAD     = run_moead(n, p$moead_K, p$moead_gen, repfn, evalfn, seed = seed)
    )
  })
}


## 10.3 Checkpointed, parallel execution

library(parallel)

## lapply() over X on n_cores, used by every batch runner in R/06-R/11.
## A worker that fails returns a "try-error" from mclapply() instead of
## stopping the batch; this stops with the failed tasks' messages, so a
## failure cannot be silently bound into the results.
parallel_map <- function(X, FUN, n_cores) {
  if (.Platform$OS.type == "windows" && n_cores > 1) {
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))
    clusterEvalQ(cl, library(parallel))
    ## Exports every object currently defined (all functions from sections
    ## 7, 9 and 10 that the worker depends on, transitively) rather than
    ## hand-maintaining a dependency list that silently goes stale the next
    ## time a function in that chain changes.
    clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
    results <- parLapply(cl, X, FUN)
  } else {
    results <- mclapply(X, FUN, mc.cores = n_cores, mc.preschedule = FALSE)
  }
  failed <- vapply(results, function(r) inherits(r, "try-error") || is.null(r), logical(1))
  if (any(failed))
    stop(sprintf("%d of %d parallel tasks failed:\n%s", sum(failed), length(results),
                 paste(unique(vapply(results[failed], function(r)
                   if (is.null(r)) "no result (worker killed?)" else as.character(r), character(1))),
                   collapse = "\n")), call. = FALSE)
  results
}

## Writes one checkpoint file per instance. If that file already exists,
## the instance is skipped entirely and its saved result is returned, so a
## multi-hour or multi-day batch can be killed (or crash) and resumed
## without losing finished work or re-spending computation on it.
run_instance_production <- function(inst_id, seed, density, eval_budget = 15000,
                                     reps = 15, risk_alpha = 0.25, nsim_probe = 30,
                                     nsim_eval = 60, nsim_posthoc = 500,
                                     n = 50, must_k = 3, budget_fraction = 0.5,
                                     deadline_quantile = 0.4,
                                     ga_pop = 100, tabu_neigh = 10, moead_K = 100,
                                     checkpoint_dir = "checkpoints_production") {
  dir.create(checkpoint_dir, showWarnings = FALSE)
  ckpt <- file.path(checkpoint_dir, sprintf("inst_%03d.rds", inst_id))
  if (file.exists(ckpt)) {
    message(sprintf("[instance %d] checkpoint found, skipping", inst_id))
    return(readRDS(ckpt))
  }

  inst <- make_instance(n, seed, density, must_k, budget_fraction, deadline_quantile)
  p <- derive_algo_params(eval_budget, ga_pop = ga_pop, tabu_neigh = tabu_neigh, moead_K = moead_K)
  rows <- list()

  for (r in seq_len(reps)) {
    scen_big <- matrix(rtri(nsim_posthoc * n, inst$acts$optimistic, inst$acts$pessimistic,
                             inst$acts$realistic), nrow = nsim_posthoc, byrow = TRUE)
    for (algo in names(ALGO_SCOPE)) for (placement in ALGO_SCOPE[[algo]]) {
      m <- if (placement == "single") 1L else 2L
      for (strategy in REPAIRS) {
        seed_run <- inst_id * 1e6 + r * 1e4 +
          which(names(ALGO_SCOPE) == algo) * 1e2 + which(REPAIRS == strategy) * 10 + m
        t0 <- proc.time()[3]
        res <- tryCatch(
          run_one_budgeted(algo, strategy, m, inst, p, risk_alpha, nsim_probe, nsim_eval, seed_run),
          error = function(e) { message(sprintf("FAILED inst=%d algo=%s strategy=%s: %s",
                                                  inst_id, algo, strategy, conditionMessage(e))); NULL }
        )
        wall <- proc.time()[3] - t0
        if (is.null(res)) next
        sc <- score_result(res, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen_big)
        rows[[length(rows) + 1]] <- data.frame(
          instance = inst_id, seed = seed, density = density, rep = r,
          algorithm = algo, placement = placement, repair = strategy,
          hypervolume = sc$hypervolume, front_size = sc$front_size, best_roi = sc$best_roi,
          risk_at_best_roi = sc$risk_at_best_roi, feasible_fraction = sc$feasible_fraction,
          eval_budget = eval_budget, wall_sec = wall)
      }
    }
  }
  out <- do.call(rbind, rows)
  saveRDS(out, ckpt)
  out
}

## Parallel batch over many instances. 
## we use the commented PSOCK-cluster block below.
## Every parameter run_instance_production() accepts is exposed here too
## and passed straight through.
run_production_batch <- function(instance_specs, eval_budget = 15000, reps = 15,
                                  risk_alpha = 0.25, nsim_probe = 30, nsim_eval = 60, nsim_posthoc = 500,
                                  n = 50, must_k = 3, budget_fraction = 0.5, deadline_quantile = 0.4,
                                  ga_pop = 100, tabu_neigh = 10, moead_K = 100,
                                  n_cores = max(1, parallel::detectCores() - 1),
                                  checkpoint_dir = "checkpoints_production") {
  worker <- function(i) {
    spec <- instance_specs[[i]]
    run_instance_production(inst_id = i, seed = spec$seed, density = spec$density,
                             eval_budget = eval_budget, reps = reps, risk_alpha = risk_alpha,
                             nsim_probe = nsim_probe, nsim_eval = nsim_eval, nsim_posthoc = nsim_posthoc,
                             n = n, must_k = must_k, budget_fraction = budget_fraction,
                             deadline_quantile = deadline_quantile, ga_pop = ga_pop,
                             tabu_neigh = tabu_neigh, moead_K = moead_K, checkpoint_dir = checkpoint_dir)
  }

  do.call(rbind, parallel_map(seq_along(instance_specs), worker, n_cores))
}


## 10.4 Statistical procedure

cliffs_delta <- function(x, y) (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
delta_label <- function(d) {
  ad <- abs(d)
  if (ad < 0.147) "negligible" else if (ad < 0.33) "small" else if (ad < 0.474) "medium" else "large"
}

analyse_repair_effect <- function(res, placement, outcome = "hypervolume") {
  d <- res[res$placement == placement, ]
  out <- list()
  for (a in unique(d$algorithm)) {
    da <- d[d$algorithm == a, ]
    da$block <- interaction(da$instance, da$rep, drop = TRUE)
    wide <- reshape(da[, c("block", "repair", outcome)], idvar = "block", timevar = "repair", direction = "wide")
    wide <- wide[stats::complete.cases(wide), , drop = FALSE]
    mat <- as.matrix(wide[, -1, drop = FALSE])
    colnames(mat) <- sub(paste0("^", outcome, "\\."), "", colnames(mat))
    if (nrow(mat) < 3 || ncol(mat) < 2) { out[[a]] <- list(note = "insufficient data"); next }
    fr <- tryCatch(friedman.test(mat), error = function(e) NULL)
    strategies <- colnames(mat); deltas <- list()
    for (i in seq_along(strategies)) for (j in seq_along(strategies)) if (i < j) {
      dd <- cliffs_delta(mat[, i], mat[, j])
      deltas[[paste(strategies[i], "vs", strategies[j])]] <- c(delta = round(dd, 3), magnitude = delta_label(dd))
    }
    out[[a]] <- list(n_blocks = nrow(mat), friedman = fr, cliffs_delta = deltas)
  }
  out
}

## Shared by sections 10.5 and 10.6, so pilot data and a full production
## run are always reported the same way. Prints all pairwise Cliff's
## deltas since which pair matters most depends on the algorithm.
print_repair_analysis <- function(analysis_list, label = "") {
  if (nzchar(label)) cat("===", label, "===\n")
  for (a in names(analysis_list)) {
    r <- analysis_list[[a]]
    if (!is.null(r$note)) { cat(sprintf("%-10s %s\n", a, r$note)); next }
    cat(sprintf("%-10s n_blocks=%d  Friedman chi-sq=%.1f (df=%d) p=%.3g\n",
                a, r$n_blocks, r$friedman$statistic, r$friedman$parameter, r$friedman$p.value))
    for (nm in names(r$cliffs_delta)) {
      dd <- r$cliffs_delta[[nm]]
      cat(sprintf("           %-40s delta=%-7s (%s)\n", nm, dd["delta"], dd["magnitude"]))
    }
  }
  invisible(analysis_list)
}


## 10.5 Cost accounting

inst <- make_instance(50, seed = 11, density = 0.15)
p <- derive_algo_params(500, ga_pop = 15, tabu_neigh = 6, moead_K = 12)
for (algo in names(ALGO_SCOPE)) for (placement in ALGO_SCOPE[[algo]]) {
  m <- if (placement == "single") 1L else 2L
  for (strategy in c(REPAIRS)) {
    t0 <- proc.time()[3]
    run_one_budgeted(algo, strategy, m, inst, p, risk_alpha = 0.25, nsim_probe = 30, nsim_eval = 60, seed = 1)
    cat(algo, placement, strategy, round(proc.time()[3] - t0, 2), "s\n")
  }
}


## 10.6 Running full-scale comparison

## Instance i has generator seed 100 + i and a density drawn from U(0.05, 0.25).
## In the original run the densities came from runif() without a seed of
## their own, so they depended on whatever ran before in the R session, and
## the code alone cannot regenerate them (replaying sections 8-10 in order
## gives other values than the ones the paper reports). They are stored in
## the original results, columns `seed` and `density`, one row per instance.
## design_a_specs() reuses them from full_production_results.rds or from the
## per-instance checkpoints when either is available, and otherwise draws a
## new, seeded set of 20 instances (reproducible, but not the paper's).
design_a_specs <- function(results_file = "full_production_results.rds",
                           checkpoint_dir = "checkpoints_production", n_inst = 20, seed = 2026) {
  stored <- NULL
  if (file.exists(results_file)) {
    stored <- readRDS(results_file)
  } else {
    ckpts <- file.path(checkpoint_dir, sprintf("inst_%03d.rds", seq_len(n_inst)))
    if (all(file.exists(ckpts))) stored <- do.call(rbind, lapply(ckpts, readRDS))
  }
  if (!is.null(stored)) {
    sp <- unique(stored[, c("instance", "seed", "density")])
    sp <- sp[order(sp$instance), ]
    stopifnot(nrow(sp) == n_inst, identical(as.integer(sp$instance), seq_len(n_inst)))
    message("Design A instances: reusing the stored seeds and densities")
    return(lapply(seq_len(n_inst), function(i) list(seed = sp$seed[i], density = sp$density[i])))
  }
  warning("Original Design A results not found: drawing a new, seeded set of instances. ",
          "Results will not match the paper's Design A.", call. = FALSE)
  set.seed(seed)
  dens <- runif(n_inst, 0.05, 0.25)
  lapply(seq_len(n_inst), function(i) list(seed = 100 + i, density = dens[i]))
}

instance_specs <- design_a_specs()
## Production setting: nsim_probe = 30 Monte Carlo draws per cut (Design B,
## Section 11, uses 15).
full_production <- run_production_batch(instance_specs, eval_budget = 3000, reps = 30,
                                         nsim_probe = 30, n_cores = 6)
saveRDS(full_production, "full_production_results.rds")
