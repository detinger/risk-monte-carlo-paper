# Revision experiments (DEVELOPMENT_PLAN.md, Phase 2)
# 2.1 probe-size sweep, 2.2 zero-range audit, 2.3 placebo (RAND-RISK),
# 2.4 deterministic proxy (CRIT), 2.6 equal wall-clock (EFF with 3x budget).
# Tabu search on Design A's instances and configuration (Section 10).
# Requires: source("06_full_scale_run.R"), and compare_pair() from
# 09_targeted_checks.R for the summary.


## 13.1 Arms

## An arm is one repair strategy under one placement, evaluation budget
## (as a multiple of Design A's 3000), and probe count. `seed_as` names the
## strategy whose run seed the arm uses. Arms that repeat a Design A
## strategy use Design A's seed formula, so on the same instances their
## search follows the Design A runs exactly (the rescoring draws differ, see
## run_revision_block()). EFF-x3 shares EFF's seeds, and the probe arms share
## RISK's: each is the same search given more evaluations, or fewer or more
## probe draws. RAND-RISK and CRIT take strategy indices 8 and 9, after the
## seven of REPAIRS, which keeps their seeds apart from every Design A run.
REVISION_REPAIRS <- c(REPAIRS, "risk_shuffled", "critical_path")

revision_arm <- function(label, repair, placement, budget_mult = 1, nsim_probe = 30, seed_as = repair)
  data.frame(label = label, repair = repair, placement = placement, budget_mult = budget_mult,
             nsim_probe = nsim_probe, seed_as = seed_as, stringsAsFactors = FALSE)

REVISION_MAIN_ARMS <- rbind(
  revision_arm("EFF",       "efficiency_only",   "multi"),
  revision_arm("EFF-SM",    "efficiency_random", "multi"),
  revision_arm("RISK",      "risk_aware",        "multi"),
  revision_arm("RAND-RISK", "risk_shuffled",     "multi"),
  revision_arm("CRIT",      "critical_path",     "multi"),
  revision_arm("EFF-x3",    "efficiency_only",   "multi", budget_mult = 3),
  revision_arm("EFF",       "efficiency_only",   "single"),
  revision_arm("RISK",      "risk_aware",        "single"),
  revision_arm("RAND-RISK", "risk_shuffled",     "single"),
  revision_arm("CRIT",      "critical_path",     "single"))

## Probe-size sweep; RISK at nsim_probe = 30 is in the main arms.
REVISION_PROBE_ARMS <- do.call(rbind, lapply(c(5, 10, 60), function(k)
  revision_arm(paste0("RISK-p", k), "risk_aware", "multi", nsim_probe = k)))

arms_tag <- function(arms) {
  s <- paste(apply(arms, 1, paste, collapse = ":"), collapse = "|")
  sprintf("%05d", sum(utf8ToInt(s) * seq_along(utf8ToInt(s))) %% 100000)
}


## 13.2 Running one block

## Same repair and evaluation set-up as run_one_budgeted(), tabu branch,
## plus the trace counters of the repair operator.
run_tabu_arm <- function(inst, arm, p, risk_alpha, nsim_eval, seed, trace = NULL) {
  m <- if (arm$placement == "single") 1L else 2L
  with(inst, {
    repfn <- function(x) apply_repair_strategy(x, acts, preds, succ, topo, budget, must_idx,
                                                deadline, strategy = arm$repair,
                                                risk_alpha = risk_alpha, nsim_probe = arm$nsim_probe,
                                                trace = trace)
    evalfn <- make_evaluate(m, acts, preds, topo, budget, deadline, nsim_eval = nsim_eval)
    run_tabu(nrow(acts), p$tabu_iter, repfn, evalfn, m = m,
             neighbourhood_size = p$tabu_neigh, seed = seed)
  })
}

## One block: one instance and one replication, every arm. All arms of a
## block are rescored on the same scenario matrix, seeded per block.
run_revision_block <- function(inst_id, seed, density, r, arms, eval_budget = 3000,
                               risk_alpha = 0.25, nsim_eval = 60, nsim_posthoc = 500,
                               n = 50, must_k = 3, budget_fraction = 0.5, deadline_quantile = 0.4,
                               tabu_neigh = 10, checkpoint_dir = "checkpoints_revision") {
  ## The file name carries the instance's seed and density, so checkpoints
  ## from another instance set (design_a_specs() falls back to a new one when
  ## the original results are missing) are never picked up by mistake.
  dir.create(checkpoint_dir, showWarnings = FALSE)
  ckpt <- file.path(checkpoint_dir, sprintf("rev_inst_%03d_s%d_d%.6f_rep_%02d_E%d_%s.rds",
                                            inst_id, seed, density, r, eval_budget, arms_tag(arms)))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  inst <- make_instance(n, seed, density, must_k, budget_fraction, deadline_quantile)
  set.seed(seed * 100 + r)
  scen_big <- matrix(rtri(nsim_posthoc * n, inst$acts$optimistic, inst$acts$pessimistic,
                          inst$acts$realistic), nrow = nsim_posthoc, byrow = TRUE)
  rows <- list()
  for (a in seq_len(nrow(arms))) {
    arm <- arms[a, ]
    m <- if (arm$placement == "single") 1L else 2L
    p <- derive_algo_params(eval_budget * arm$budget_mult, tabu_neigh = tabu_neigh)
    seed_run <- inst_id * 1e6 + r * 1e4 + which(names(ALGO_SCOPE) == "Tabu") * 1e2 +
      which(REVISION_REPAIRS == arm$seed_as) * 10 + m
    trace <- new_cut_trace()
    t0 <- proc.time()[3]
    res <- run_tabu_arm(inst, arm, p, risk_alpha, nsim_eval, seed_run, trace)
    wall <- proc.time()[3] - t0
    sc <- score_result(res, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen_big)
    rows[[a]] <- data.frame(
      instance = inst_id, seed = seed, density = density, rep = r, algorithm = "Tabu",
      placement = arm$placement, arm = arm$label, repair = arm$repair,
      eval_budget = eval_budget * arm$budget_mult, nsim_probe = arm$nsim_probe,
      hypervolume = sc$hypervolume, front_size = sc$front_size, best_roi = sc$best_roi,
      risk_at_best_roi = sc$risk_at_best_roi, feasible_fraction = sc$feasible_fraction,
      wall_sec = wall, risk_cuts = trace$cuts,
      mean_leaves = if (trace$cuts) trace$leaves / trace$cuts else NA_real_,
      zero_range_share = if (trace$cuts) trace$zero_range / trace$cuts else NA_real_,
      differs_from_eff_share = if (trace$cuts) trace$differs_from_eff / trace$cuts else NA_real_)
  }
  out <- do.call(rbind, rows)
  saveRDS(out, ckpt)
  out
}

## Every (instance, replication) block in parallel; one checkpoint per block,
## so an interrupted run resumes where it stopped. Arguments in ... go to
## run_revision_block().
run_revision_experiments <- function(instance_specs, arms, reps = 30,
                                     n_cores = max(1, parallel::detectCores() - 1), ...) {
  blocks <- expand.grid(rep = seq_len(reps), inst_id = seq_along(instance_specs))
  worker <- function(k) {
    i <- blocks$inst_id[k]
    run_revision_block(i, instance_specs[[i]]$seed, instance_specs[[i]]$density,
                       blocks$rep[k], arms, ...)
  }
  do.call(rbind, parallel_map(seq_len(nrow(blocks)), worker, n_cores))
}


## 13.3 Summary

## Instance-level comparisons go through compare_pair() (Section 11.7), which
## matches strategies by its `repair` column; here that is the arm label.
summarise_revision <- function(main, probe = NULL, n_boot = 2000) {
  by_arm <- function(d) { d$repair <- d$arm; d }
  means <- do.call(rbind, lapply(split(main, list(main$placement, main$arm), drop = TRUE), function(d)
    data.frame(placement = d$placement[1], arm = d$arm[1], eval_budget = d$eval_budget[1],
               hypervolume = mean(d$hypervolume), front_size = mean(d$front_size),
               best_roi = mean(d$best_roi, na.rm = TRUE), risk_at_best_roi = mean(d$risk_at_best_roi, na.rm = TRUE),
               wall_sec = mean(d$wall_sec))))
  eff_wall <- setNames(means$wall_sec[means$arm == "EFF"], means$placement[means$arm == "EFF"])
  means$wall_vs_eff <- round(means$wall_sec / eff_wall[means$placement], 2)
  means <- means[order(means$placement, means$arm), ]
  cat("=== Means by placement and arm (tabu search, Design A instances) ===\n")
  print(means, row.names = FALSE, digits = 3)

  dm <- by_arm(main[main$placement == "multi", ])
  multi_pairs <- list(c("RISK", "EFF"), c("EFF-SM", "EFF"), c("RAND-RISK", "EFF"), c("RAND-RISK", "RISK"),
                      c("CRIT", "EFF"), c("CRIT", "RISK"), c("EFF-x3", "EFF"), c("EFF-x3", "RISK"))
  multi <- do.call(rbind, lapply(multi_pairs, function(pr) compare_pair(dm, pr[1], pr[2], n_boot = n_boot)))
  cat("\n=== Multi-objective placement, hypervolume: instance-level comparisons ===\n")
  cat("RAND-RISK ~ RISK: the gain comes from perturbing the cut, not from risk information (2.3)\n")
  cat("CRIT ~ RISK: a deterministic critical-path rule does as well as Monte Carlo probes (2.4)\n")
  cat("EFF-x3 ~ RISK: given RISK's wall-clock time, efficiency-only repair closes the gap (2.6)\n\n")
  print(multi, row.names = FALSE)

  ds <- by_arm(main[main$placement == "single", ])
  single <- do.call(rbind, lapply(list(c("RISK", "EFF"), c("RAND-RISK", "EFF"), c("CRIT", "EFF"), c("CRIT", "RISK")),
    function(pr) rbind(cbind(outcome = "risk_at_best_roi", compare_pair(ds, pr[1], pr[2], "risk_at_best_roi", n_boot)),
                       cbind(outcome = "best_roi", compare_pair(ds, pr[1], pr[2], "best_roi", n_boot)))))
  cat("\n=== Single-objective placement: risk and value at the best-value solution ===\n")
  print(single, row.names = FALSE)

  ## Per-instance effect sizes (DEVELOPMENT_PLAN 3.3): Cliff's delta over the
  ## replications of each instance, against EFF.
  per_inst <- do.call(rbind, lapply(c("RISK", "EFF-SM", "RAND-RISK", "CRIT", "EFF-x3"), function(a)
    do.call(rbind, lapply(sort(unique(dm$instance)), function(i) {
      di <- dm[dm$instance == i, ]
      data.frame(arm = a, instance = i,
                 delta = cliff_rank(di$hypervolume[di$arm == a], di$hypervolume[di$arm == "EFF"]))
    }))))
  dist <- do.call(rbind, lapply(split(per_inst, per_inst$arm), function(d)
    data.frame(arm = d$arm[1], at_1 = sum(d$delta == 1), large = sum(d$delta >= 0.474),
               negligible = sum(abs(d$delta) < 0.147), negative = sum(d$delta <= -0.147),
               median = round(median(d$delta), 3), min = round(min(d$delta), 3), max = round(max(d$delta), 3))))
  cat("\n=== Per-instance Cliff's delta against EFF, hypervolume (counts of instances) ===\n")
  print(dist, row.names = FALSE)

  ## Zero-range audit (2.2), risk-scored arms only.
  scored <- main[main$risk_cuts > 0, ]
  audit <- do.call(rbind, lapply(split(scored, list(scored$placement, scored$arm), drop = TRUE), function(d)
    data.frame(placement = d$placement[1], arm = d$arm[1], cuts_per_run = round(mean(d$risk_cuts)),
               mean_leaves = round(mean(d$mean_leaves), 2),
               zero_range_share = round(mean(d$zero_range_share), 3),
               differs_from_eff_share = round(mean(d$differs_from_eff_share), 3))))
  cat("\n=== Risk-scored cuts: share with a zero-range risk signal, share differing from EFF's cut ===\n")
  print(audit, row.names = FALSE)

  out <- list(means = means, multi = multi, single = single, per_instance = per_inst,
              per_instance_summary = dist, audit = audit)

  if (!is.null(probe)) {
    reps_p <- sort(unique(probe$rep))
    base <- main[main$placement == "multi" & main$arm == "RISK" & main$rep %in% reps_p, ]
    dp <- by_arm(rbind(base[, names(probe)], probe))
    probe_means <- aggregate(cbind(hypervolume, risk_at_best_roi, wall_sec, zero_range_share) ~ nsim_probe,
                             dp, mean)
    probe_cmp <- do.call(rbind, lapply(setdiff(unique(probe$arm), "RISK"), function(a)
      compare_pair(dp, a, "RISK", n_boot = n_boot)))
    cat("\n=== Probe-size sweep (2.1), RISK, multi-objective placement, replications",
        min(reps_p), "to", max(reps_p), "===\n")
    print(probe_means, row.names = FALSE, digits = 3)
    cat("\n")
    print(probe_cmp, row.names = FALSE)
    out$probe_means <- probe_means; out$probe <- probe_cmp
  }

  for (nm in names(out)) write.csv(out[[nm]], sprintf("revision_%s.csv", nm), row.names = FALSE)
  invisible(out)
}


## 13.4 Running the revision experiments

## Main arms: 20 instances x 30 replications x 10 arms. About 40 core-hours
## (tabu search, from the Design A cost table: EFF ~13 s, RISK ~41 s per run),
## roughly 6 to 7 hours on 6 cores. Probe arms, 10 replications: about 6 core-hours.
instance_specs <- design_a_specs()
revision_main <- run_revision_experiments(instance_specs, REVISION_MAIN_ARMS, reps = 30, n_cores = 6)
saveRDS(revision_main, "revision_main_results.rds")
revision_probe <- run_revision_experiments(instance_specs, REVISION_PROBE_ARMS, reps = 10, n_cores = 6)
saveRDS(revision_probe, "revision_probe_results.rds")

revision_summary <- summarise_revision(revision_main, revision_probe)
saveRDS(revision_summary, "revision_summary.rds")


## 13.5 Scale-balanced runs for all algorithms (Phase 2.5)
##
## Extends the tabu-only scalarization check (R/09 §11.6) to all six algorithms
## and all seven repair strategies, 30 replications each, on Design A's 20
## instances.  The objective scale is set per-instance from the ROI of the
## efficiency-only reference solution, matching §11.6.
##
## Algorithms that do not use a scalarized objective (NSGA-II, random search)
## are scale-invariant by construction; they are still included so that the
## table of results covers the full algorithm set and makes the non-effect of
## normalization on those algorithms explicit.  GA_single is single-objective
## only and is also omitted from the multi-objective comparison here, as in
## Design A.
##
## Estimated cost: all algorithms × 7 strategies × 30 reps × 20 instances.
## Wall-clock is dominated by RISK repair.  At Design A rates (EFF ~13 s,
## RISK ~41 s per tabu run), this is roughly 3–4× the revision main-arm cost,
## or about 20–25 hours on 6 cores.  Checkpointed: rerun to resume.


## 13.5.1 Scale-balanced SA
##
## SA selects moves by a randomly re-weighted sum of the two objectives.
## With raw objectives (value ~tens, risk in [0,1]) the risk term contributes
## little to that sum.  Dividing both objectives by obj_scale before forming
## the delta puts them on comparable units.  Everything else (archive, final
## evaluation) is unchanged.
run_sa_normalized <- function(n_genes, n_iter, repair_fn, evaluate_fn, obj_scale,
                               t0 = 1, cool = 0.995, k_flip = 2, seed = 1) {
  set.seed(seed)
  x <- repair_fn(sample(c(0L, 1L), n_genes, replace = TRUE)); obj <- evaluate_fn(x)
  archive <- list(list(x = x, obj = obj)); T <- t0
  for (it in seq_len(n_iter)) {
    xn <- repair_fn(flip_k_bits(x, k_flip)); objn <- evaluate_fn(xn)
    w <- runif(1)
    delta <- (w * objn[1] / obj_scale[1] + (1 - w) * objn[2] / obj_scale[2]) -
             (w * obj[1]  / obj_scale[1] + (1 - w) * obj[2]  / obj_scale[2])
    if (delta < 0 || runif(1) < exp(-delta / max(T, 1e-6))) { x <- xn; obj <- objn }
    archive <- archive_insert(archive, xn, objn)
    T <- T * cool
  }
  objs <- matrix(vapply(archive, function(p) p$obj, numeric(2)), ncol = 2, byrow = TRUE)
  list(population = archive, objectives = objs, pareto_idx = fast_nondominated_sort(objs)[[1]])
}


## 13.5.2 Scale-balanced MOEA/D
##
## MOEA/D's Tchebycheff aggregation computes max(wk * |obj - z|) where z is
## the running ideal point.  With unscaled objectives the large value dimension
## dominates.  Normalizing by obj_scale before computing Tchebycheff distances
## and updating z gives the two objectives equal unit weight.
run_moead_normalized <- function(n_genes, K, n_gen, repair_fn, evaluate_fn, obj_scale,
                                  neighbour_T = 4, pcx = 0.7, pmut = 0.15, seed = 1) {
  set.seed(seed)
  w <- cbind(seq(0, 1, length.out = K), 1 - seq(0, 1, length.out = K))
  nb <- lapply(seq_len(K), function(k) {
    d <- rowSums((w - matrix(w[k, ], K, 2, byrow = TRUE))^2)
    order(d)[seq_len(min(neighbour_T, K))]
  })
  pop <- vector("list", K)
  for (k in seq_len(K)) {
    x <- repair_fn(sample(c(0L, 1L), n_genes, replace = TRUE))
    pop[[k]] <- list(x = x, obj = evaluate_fn(x))
  }
  ## Ideal point tracked in normalized space.
  z <- apply(do.call(rbind, lapply(pop, function(p) p$obj / obj_scale)), 2, min)
  tcheby <- function(obj, wk) max(wk * abs(obj / obj_scale - z))
  crossover <- function(a, b) {
    if (runif(1) > pcx) return(a)
    pt <- sample(2:(n_genes - 1), 1); c(a[1:pt], b[(pt + 1):n_genes])
  }
  mutate <- function(a) { flip <- runif(n_genes) < pmut; a[flip] <- 1L - a[flip]; a }
  for (g in seq_len(n_gen)) for (k in seq_len(K)) {
    parents <- sample(nb[[k]], 2)
    child_x <- repair_fn(mutate(crossover(pop[[parents[1]]]$x, pop[[parents[2]]]$x)))
    child_obj <- evaluate_fn(child_x)
    z <- pmin(z, child_obj / obj_scale)
    for (j in nb[[k]]) if (tcheby(child_obj, w[j, ]) <= tcheby(pop[[j]]$obj, w[j, ]))
      pop[[j]] <- list(x = child_x, obj = child_obj)
  }
  objs <- matrix(vapply(pop, function(p) p$obj, numeric(2)), ncol = 2, byrow = TRUE)
  list(population = pop, objectives = objs, pareto_idx = fast_nondominated_sort(objs)[[1]])
}


## 13.5.3 Dispatch for all algorithms with scale-balanced objectives
##
## For NSGA-II and random search, which do not use a scalarized objective,
## this falls through to run_one_budgeted() (obj_scale has no effect on their
## move acceptance; it is computed but ignored).
run_one_budgeted_normalized <- function(algo, strategy, m, inst, p, risk_alpha, nsim_probe,
                                         nsim_eval, seed, obj_scale) {
  with(inst, {
    repfn  <- function(x) apply_repair_strategy(x, acts, preds, succ, topo, budget, must_idx,
                                                 deadline, strategy = strategy,
                                                 risk_alpha = risk_alpha, nsim_probe = nsim_probe)
    evalfn <- make_evaluate(m, acts, preds, topo, budget, deadline, nsim_eval = nsim_eval,
                             penalize = (strategy == "penalty"))
    n <- nrow(acts)
    switch(algo,
      Random    = run_random_search(n, p$random_eval, repfn, evalfn, seed = seed),
      GA_single = run_evolution(n, p$pop_size, p$n_gen_ga, 0.7, 0.15, repfn, evalfn, seed = seed),
      NSGA2     = run_evolution(n, p$pop_size, p$n_gen_ga, 0.7, 0.15, repfn, evalfn, seed = seed),
      SA        = run_sa_normalized(n, p$sa_iter, repfn, evalfn, obj_scale = obj_scale, seed = seed),
      Tabu      = run_tabu_normalized(n, p$tabu_iter, repfn, evalfn, m = m, obj_scale = obj_scale,
                                      neighbourhood_size = p$tabu_neigh, seed = seed),
      MOEAD     = run_moead_normalized(n, p$moead_K, p$moead_gen, repfn, evalfn,
                                       obj_scale = obj_scale, seed = seed)
    )
  })
}


## 13.5.4 One block: one instance × one replication, all multi-objective
##         algorithms × all repair strategies, scale-balanced.
##         obj_scale uses the same reference as §11.6: ROI of EFF repair on the
##         full activity set, which is a budget-feasible upper bound for this
##         instance.
run_scale_balanced_block <- function(inst_id, seed, density, r, eval_budget = 3000,
                                      risk_alpha = 0.25, nsim_eval = 60, nsim_posthoc = 500,
                                      n = 50, must_k = 3, budget_fraction = 0.5,
                                      deadline_quantile = 0.4, ga_pop = 100, tabu_neigh = 10,
                                      moead_K = 100, nsim_probe = 30,
                                      checkpoint_dir = "checkpoints_scale_balanced") {
  dir.create(checkpoint_dir, showWarnings = FALSE)
  ckpt <- file.path(checkpoint_dir,
                    sprintf("sb_inst_%03d_s%d_d%.6f_rep_%02d.rds", inst_id, seed, density, r))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  inst <- make_instance(n, seed, density, must_k, budget_fraction, deadline_quantile)

  ## Scale from the EFF reference solution (identical to §11.6).
  ref_sel   <- repair_activity_network(rep(1L, nrow(inst$acts)),
                                       inst$acts, inst$preds, inst$succ, inst$topo,
                                       inst$budget, inst$must_idx, inst$deadline, alpha = 1)
  obj_scale <- c(unname(cost_roi(ref_sel, inst$acts)["roi"]), 1)

  p <- derive_algo_params(eval_budget, ga_pop = ga_pop, tabu_neigh = tabu_neigh, moead_K = moead_K)

  set.seed(seed * 100 + r)
  scen_big <- matrix(rtri(nsim_posthoc * n, inst$acts$optimistic, inst$acts$pessimistic,
                          inst$acts$realistic), nrow = nsim_posthoc, byrow = TRUE)

  ## Multi-objective placement only; GA_single is excluded.
  sb_algos <- setdiff(names(ALGO_SCOPE), "GA_single")
  sb_algos <- sb_algos[vapply(sb_algos, function(a) "multi" %in% ALGO_SCOPE[[a]], logical(1))]

  rows <- list()
  for (algo in sb_algos) {
    m <- 2L
    for (strategy in REPAIRS) {
      seed_run <- inst_id * 1e6 + r * 1e4 + which(names(ALGO_SCOPE) == algo) * 1e2 +
        which(REPAIRS == strategy) * 10 + m
      t0 <- proc.time()[3]
      res <- tryCatch(
        run_one_budgeted_normalized(algo, strategy, m, inst, p, risk_alpha, nsim_probe,
                                    nsim_eval, seed_run, obj_scale),
        error = function(e) {
          message(sprintf("FAILED sb inst=%d algo=%s strategy=%s: %s",
                          inst_id, algo, strategy, conditionMessage(e))); NULL
        }
      )
      wall <- proc.time()[3] - t0
      if (is.null(res)) next
      sc <- score_result(res, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen_big)
      rows[[length(rows) + 1]] <- data.frame(
        instance = inst_id, seed = seed, density = density, rep = r,
        algorithm = algo, placement = "multi", repair = strategy,
        hypervolume = sc$hypervolume, front_size = sc$front_size, best_roi = sc$best_roi,
        risk_at_best_roi = sc$risk_at_best_roi, feasible_fraction = sc$feasible_fraction,
        eval_budget = eval_budget, wall_sec = wall, obj_scale1 = obj_scale[1])
    }
  }
  out <- do.call(rbind, rows)
  saveRDS(out, ckpt)
  out
}


## 13.5.5 Parallel batch over all Design A instances.
run_scale_balanced_experiments <- function(instance_specs, reps = 30,
                                            n_cores = max(1, parallel::detectCores() - 1), ...) {
  blocks <- expand.grid(rep = seq_len(reps), inst_id = seq_along(instance_specs))
  worker <- function(k) {
    i <- blocks$inst_id[k]
    run_scale_balanced_block(i, instance_specs[[i]]$seed, instance_specs[[i]]$density,
                              blocks$rep[k], ...)
  }
  do.call(rbind, parallel_map(seq_len(nrow(blocks)), worker, n_cores))
}


## 13.5.6 Summary: effect of normalization on each algorithm × strategy.
##         Also compares raw (Design A) vs. scale-balanced for each algorithm.
summarise_scale_balanced <- function(sb, full_production = NULL, n_boot = 2000) {
  ## Restrict to the instances and replications present in both datasets.
  reps_sb <- sort(unique(sb$rep))

  ## Mean hypervolume by algorithm × strategy for each scalarization.
  means_sb  <- aggregate(hypervolume ~ algorithm + repair, sb,  mean)
  means_sb$scalarization  <- "normalized"
  if (!is.null(full_production)) {
    raw <- full_production[full_production$placement == "multi" &
                           full_production$rep %in% reps_sb &
                           full_production$algorithm %in% unique(sb$algorithm), ]
    means_raw <- aggregate(hypervolume ~ algorithm + repair, raw, mean)
    means_raw$scalarization <- "raw"
    means <- rbind(means_sb, means_raw)
  } else {
    cat("Note: full_production is NULL; raw-aggregation comparison skipped.\n")
    means <- means_sb
  }
  means <- means[order(means$algorithm, means$repair, means$scalarization), ]
  cat("=== Mean hypervolume by algorithm, strategy, and scalarization ===\n")
  print(means, row.names = FALSE, digits = 3)

  ## Cliff's delta: RISK vs EFF under normalized objectives, per algorithm.
  pairs <- list(c("risk_aware", "efficiency_only"),
                c("efficiency_random", "efficiency_only"),
                c("risk_aware", "efficiency_random"))
  by_algo <- do.call(rbind, lapply(sort(unique(sb$algorithm)), function(a) {
    da <- sb[sb$algorithm == a, ]; da$repair <- da$repair
    do.call(rbind, lapply(pairs, function(pr) {
      res <- tryCatch(compare_pair(da, pr[1], pr[2], n_boot = n_boot),
                      error = function(e) NULL)
      if (is.null(res)) return(NULL)
      cbind(algorithm = a, res)
    }))
  }))
  cat("\n=== Cliff's delta, normalized objectives, per algorithm ===\n")
  print(by_algo, row.names = FALSE)

  ## Per-instance delta distribution (DEVELOPMENT_PLAN §3.3) under normalized objectives.
  per_inst <- do.call(rbind, lapply(sort(unique(sb$algorithm)), function(a) {
    da <- sb[sb$algorithm == a, ]
    do.call(rbind, lapply(c("risk_aware", "efficiency_random"), function(arm) {
      do.call(rbind, lapply(sort(unique(da$instance)), function(i) {
        di <- da[da$instance == i, ]
        data.frame(algorithm = a, arm = arm, instance = i,
                   delta = cliff_rank(di$hypervolume[di$repair == arm],
                                      di$hypervolume[di$repair == "efficiency_only"]))
      }))
    }))
  }))
  dist <- do.call(rbind, lapply(split(per_inst, list(per_inst$algorithm, per_inst$arm), drop = TRUE), function(d)
    data.frame(algorithm = d$algorithm[1], arm = d$arm[1],
               at_1 = sum(d$delta == 1), large = sum(d$delta >= 0.474),
               negligible = sum(abs(d$delta) < 0.147), negative = sum(d$delta <= -0.147),
               median_delta = round(median(d$delta), 3))))
  cat("\n=== Per-instance Cliff's delta vs EFF, normalized objectives ===\n")
  print(dist, row.names = FALSE)

  out <- list(means = means, by_algo = by_algo, per_instance = per_inst, per_instance_summary = dist)
  for (nm in names(out)) write.csv(out[[nm]], sprintf("scale_balanced_%s.csv", nm), row.names = FALSE)
  invisible(out)
}


## 13.5.7 Running the scale-balanced experiments
##
## 20 instances × 30 reps × 5 algorithms × 7 strategies = 21 000 runs.
## Cost estimate: dominated by RISK repair for Tabu (~41 s) and SA (~200 s).
## Allow 20–25 hours on 6 cores; checkpointed so an interrupted run resumes.
sb_results <- run_scale_balanced_experiments(instance_specs, reps = 30, n_cores = 16)
saveRDS(sb_results, "scale_balanced_results.rds")

sb_summary <- summarise_scale_balanced(sb_results, full_production)
saveRDS(sb_summary, "scale_balanced_summary.rds")
