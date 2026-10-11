# Targeted checks (ODD §11.4–11.7)
# τ sweep, config check, scalarization check, instance-level inference.
# Requires: full_production, ofat_results, and source("04_additional_algorithms.R")


## 11.4 Temperature effect

## Minimal wrapper for the tau sensitivity sweep only. Mirrors
## run_one_budgeted's Tabu branch exactly, 
## forwarding temperature through to apply_repair_strategy.
## Does not modify run_one_budgeted itself, so nothing else in the paper's
## pipeline is touched by this.
run_tabu_with_temperature <- function(strategy, m, inst, p, risk_alpha,
                                       nsim_probe, nsim_eval, seed, temperature) {
  with(inst, {
    repfn <- function(x) apply_repair_strategy(x, acts, preds, succ, topo, budget, must_idx,
                                                deadline, strategy = strategy,
                                                risk_alpha = risk_alpha, nsim_probe = nsim_probe,
                                                temperature = temperature)
    evalfn <- make_evaluate(m, acts, preds, topo, budget, deadline, nsim_eval = nsim_eval,
                             penalize = (strategy == "penalty"))
    n <- nrow(acts)
    run_tabu(n, p$tabu_iter, repfn, evalfn, m = m,
             neighbourhood_size = p$tabu_neigh, seed = seed)
  })
}


tau_values <- c(0.1, 0.5, 1, 2, 4, 8, 16, 32, 64)
n_reps <- 20
results <- data.frame()

for (tau in tau_values) {
  hv_eff_only   <- numeric(n_reps)
  hv_eff_random <- numeric(n_reps)
  hv_risk_aware <- numeric(n_reps)

  for (r in seq_len(n_reps)) {
    inst <- make_instance(50, seed = 7000 + r, density = 0.10, must_k = 3,
                           budget_fraction = 0.5, deadline_quantile = 0.4)
    p <- derive_algo_params(3000, tabu_neigh = 8)
    scen <- matrix(rtri(300 * 50, inst$acts$optimistic, inst$acts$pessimistic, inst$acts$realistic),
                    nrow = 300, byrow = TRUE)

    res_eo <- run_one_budgeted("Tabu", "efficiency_only", 2L, inst, p, risk_alpha = 0.25,
                                nsim_probe = 30, nsim_eval = 60, seed = r * 1000)
    res_er <- run_tabu_with_temperature("efficiency_random", 2L, inst, p, risk_alpha = 0.25,
                                nsim_probe = 30, nsim_eval = 60, seed = r * 1000, temperature = tau)
    res_ra <- run_one_budgeted("Tabu", "risk_aware", 2L, inst, p, risk_alpha = 0.25,
                                nsim_probe = 30, nsim_eval = 60, seed = r * 1000)

    hv_eff_only[r]   <- score_result(res_eo, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen)$hypervolume
    hv_eff_random[r] <- score_result(res_er, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen)$hypervolume
    hv_risk_aware[r] <- score_result(res_ra, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen)$hypervolume
  }

  results <- rbind(results, data.frame(
    tau = tau,
    delta_vs_risk_aware = cliffs_delta(hv_eff_random, hv_risk_aware),   # >0: randomization beats risk-aware
    delta_vs_eff_only   = cliffs_delta(hv_eff_random, hv_eff_only)      # randomization's own benefit at this tau
  ))
}

write.csv(results, "tau_sensitivity.csv", row.names = FALSE)

print(results, row.names = FALSE)


## 11.5 Configuration check: Design B’s search configuration on Design A’s instances

## Design A and Design B differ in instance sampling and in search
## configuration. This check applies Design B's configuration (Section 11.2)
## to Design A's 20 instances (Section 10.6), for tabu search and the three
## strategies on which the designs disagree. Comparing its result with
## Design A (same instances, other configuration) and with the Design B
## baseline cell (same configuration, other instances) separates the two
## sources. Instances are rebuilt from the seeds and densities stored in
## full_production, with Design A's make_instance() arguments. Run seeds
## follow Design A's formula; rescoring draws are seeded per block.

CONFIG_CHECK_REPAIRS <- c("efficiency_only", "efficiency_random", "risk_aware")

config_check_instance <- function(inst_id, seed, density, reps = 10, algo = "Tabu",
                                  eval_budget = 6000, ga_pop = 40, tabu_neigh = 8, moead_K = 30,
                                  risk_alpha = 0.25, nsim_probe = 15, nsim_eval = 40,
                                  nsim_posthoc = 300, n = 50, must_k = 3,
                                  budget_fraction = 0.5, deadline_quantile = 0.4,
                                  checkpoint_dir = "checkpoints_config_check") {
  dir.create(checkpoint_dir, showWarnings = FALSE)
  ckpt <- file.path(checkpoint_dir, sprintf("config_inst_%03d_E%d_reps%d.rds", inst_id, eval_budget, reps))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  inst <- make_instance(n, seed, density, must_k, budget_fraction, deadline_quantile)
  p <- derive_algo_params(eval_budget, ga_pop = ga_pop, tabu_neigh = tabu_neigh, moead_K = moead_K)
  rows <- list()
  for (r in seq_len(reps)) {
    set.seed(seed * 100 + r)
    scen_big <- matrix(rtri(nsim_posthoc * n, inst$acts$optimistic, inst$acts$pessimistic,
                            inst$acts$realistic), nrow = nsim_posthoc, byrow = TRUE)
    for (strategy in CONFIG_CHECK_REPAIRS) {
      seed_run <- inst_id * 1e6 + r * 1e4 + which(names(ALGO_SCOPE) == algo) * 1e2 +
        which(REPAIRS == strategy) * 10 + 2L
      res <- run_one_budgeted(algo, strategy, 2L, inst, p, risk_alpha, nsim_probe, nsim_eval, seed_run)
      sc <- score_result(res, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen_big)
      rows[[length(rows) + 1]] <- data.frame(
        instance = inst_id, seed = seed, density = density, rep = r, algorithm = algo,
        repair = strategy, hypervolume = sc$hypervolume, front_size = sc$front_size,
        risk_at_best_roi = sc$risk_at_best_roi, eval_budget = eval_budget)
    }
  }
  out <- do.call(rbind, rows)
  saveRDS(out, ckpt)
  out
}

run_config_check <- function(full_production, reps = 10,
                             n_cores = max(1, parallel::detectCores() - 1), ...) {
  specs <- unique(full_production[, c("instance", "seed", "density")])
  specs <- specs[order(specs$instance), ]
  worker <- function(k) config_check_instance(specs$instance[k], specs$seed[k], specs$density[k],
                                              reps = reps, ...)
  do.call(rbind, parallel_map(seq_len(nrow(specs)), worker, n_cores))
}

summarise_config_check <- function(config_check, full_production, ofat_results, algo = "Tabu") {
  pooled <- function(d) {
    hv <- function(r) d$hypervolume[d$repair == r]
    c(risk_aware_vs_eff_only       = cliffs_delta(hv("risk_aware"), hv("efficiency_only")),
      eff_random_vs_eff_only       = cliffs_delta(hv("efficiency_random"), hv("efficiency_only")),
      eff_random_vs_risk_aware     = cliffs_delta(hv("efficiency_random"), hv("risk_aware")))
  }
  per_instance <- function(d) {
    dl <- sapply(sort(unique(d$instance)), function(i) {
      di <- d[d$instance == i, ]
      cliffs_delta(di$hypervolume[di$repair == "efficiency_random"], di$hypervolume[di$repair == "risk_aware"])
    })
    sprintf("%d / %d / %d", sum(dl > 0), sum(dl < 0), sum(dl == 0))
  }
  dA <- full_production[full_production$placement == "multi" & full_production$algorithm == algo &
                        full_production$repair %in% CONFIG_CHECK_REPAIRS, ]
  dB <- ofat_results[ofat_results$factor == "density" & ofat_results$level == 0.10 &
                     ofat_results$algorithm == algo & ofat_results$repair %in% CONFIG_CHECK_REPAIRS, ]
  dB$instance <- dB$run
  tab <- data.frame(
    source = c("Design A: configuration A, instances A",
               "This section: configuration B, instances A",
               "Design B baseline: configuration B, fresh instances"),
    rbind(pooled(dA), pooled(config_check), pooled(dB)),
    random_vs_risk_per_instance = c(per_instance(dA), per_instance(config_check), per_instance(dB)),
    check.names = FALSE)
  tab[, 2:4] <- round(tab[, 2:4], 3)
  cat("Cliff's delta, tabu search, multi-objective placement.\n")
  cat("Last column: instances favouring randomization / favouring risk-aware / no difference.\n\n")
  print(tab, row.names = FALSE)
  write.csv(tab, "config_check_summary.csv", row.names = FALSE)
  invisible(tab)
}

## Tabu search, 20 instances x 10 replications x 3 strategies = 600 runs at
## eval_budget = 6000. One block (three strategies) takes about two minutes on
## one core, so the run takes roughly 1.5 to 2.5 hours on six cores.
config_check <- run_config_check(full_production, reps = 10, n_cores = 6)
saveRDS(config_check, "config_check_results.rds")

summarise_config_check(config_check, full_production, ofat_results)


## 11.6. Scalarization effect check

SCALARIZATION_REPAIRS <- c("efficiency_only", "efficiency_random", "risk_aware")

run_tabu_normalized <- function(n_genes, n_iter, repair_fn, evaluate_fn, m, obj_scale,
                                neighbourhood_size = 10, tenure = 5, seed = 1) {
  set.seed(seed)
  x <- repair_fn(sample(c(0L, 1L), n_genes, replace = TRUE)); obj <- evaluate_fn(x)
  best_obj <- obj; archive <- list(list(x = x, obj = obj)); tabu_until <- integer(n_genes)
  for (it in seq_len(n_iter)) {
    cand_idx <- sample(n_genes, min(neighbourhood_size, n_genes))
    best_delta <- Inf; best_move <- NULL; best_move_x <- NULL; best_move_obj <- NULL
    for (bit in cand_idx) {
      xn <- x; xn[bit] <- 1L - xn[bit]; xn <- repair_fn(xn); objn <- evaluate_fn(xn)
      is_tabu <- tabu_until[bit] >= it
      aspires <- if (m == 1L) objn[1] < best_obj[1] else dominates(objn, best_obj)
      if (is_tabu && !aspires) next
      scal_n <- if (m == 1L) objn[1] else sum(objn / obj_scale)     # the only change
      if (scal_n < best_delta) {
        best_delta <- scal_n; best_move <- bit; best_move_x <- xn; best_move_obj <- objn
      }
    }
    if (!is.null(best_move)) {
      x <- best_move_x; obj <- best_move_obj
      tabu_until[best_move] <- it + tenure
      if (m == 1L) { if (obj[1] < best_obj[1]) best_obj <- obj } else if (dominates(obj, best_obj)) best_obj <- obj
      archive <- archive_insert(archive, x, obj)
    }
  }
  objs <- matrix(vapply(archive, function(p) p$obj, numeric(m)), ncol = m, byrow = TRUE)
  list(population = archive, objectives = objs, pareto_idx = fast_nondominated_sort(objs)[[1]])
}

## Same repair and evaluation set-up as run_one_budgeted(), tabu branch only.
run_tabu_normalized_budgeted <- function(strategy, inst, p, risk_alpha, nsim_probe, nsim_eval, seed) {
  with(inst, {
    ref_sel <- repair_activity_network(rep(1L, nrow(acts)), acts, preds, succ, topo, budget, must_idx,
                                       deadline = deadline, alpha = 1)
    obj_scale <- c(unname(cost_roi(ref_sel, acts)["roi"]), 1)
    repfn <- function(x) apply_repair_strategy(x, acts, preds, succ, topo, budget, must_idx,
                                                deadline, strategy = strategy,
                                                risk_alpha = risk_alpha, nsim_probe = nsim_probe)
    evalfn <- make_evaluate(2L, acts, preds, topo, budget, deadline, nsim_eval = nsim_eval,
                            penalize = (strategy == "penalty"))
    run_tabu_normalized(nrow(acts), p$tabu_iter, repfn, evalfn, m = 2L, obj_scale = obj_scale,
                        neighbourhood_size = p$tabu_neigh, seed = seed)
  })
}

## One Design A instance. Rescoring draws for replication r > 1 are reproduced
## by replaying the last run of replication r - 1 in Design A (MOEA/D,
## multi-objective, penalty). In the first verify_reps replications, raw
## efficiency-only tabu search is also run, as a check that instances, seeds,
## and draws match the stored Design A rows.
scalarization_check_instance <- function(inst_id, seed, density, reps = 10, verify_reps = 2,
                                         eval_budget = 3000, risk_alpha = 0.25, nsim_probe = 30,
                                         nsim_eval = 60, nsim_posthoc = 500, n = 50, must_k = 3,
                                         budget_fraction = 0.5, deadline_quantile = 0.4,
                                         ga_pop = 100, tabu_neigh = 10, moead_K = 100,
                                         checkpoint_dir = "checkpoints_scalarization") {
  dir.create(checkpoint_dir, showWarnings = FALSE)
  ckpt <- file.path(checkpoint_dir, sprintf("scal_inst_%03d_reps%d.rds", inst_id, reps))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  inst <- make_instance(n, seed, density, must_k, budget_fraction, deadline_quantile)
  p <- derive_algo_params(eval_budget, ga_pop = ga_pop, tabu_neigh = tabu_neigh, moead_K = moead_K)
  seed_of <- function(r, algo, strategy, m) inst_id * 1e6 + r * 1e4 +
    which(names(ALGO_SCOPE) == algo) * 1e2 + which(REPAIRS == strategy) * 10 + m
  last_algo <- tail(names(ALGO_SCOPE), 1)
  last_m <- if (tail(ALGO_SCOPE[[last_algo]], 1) == "single") 1L else 2L
  last_strategy <- tail(REPAIRS, 1)

  score_row <- function(res, scen_big, r, strategy, scalarization) {
    sc <- score_result(res, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen_big)
    data.frame(instance = inst_id, rep = r, algorithm = "Tabu", placement = "multi",
               repair = strategy, scalarization = scalarization,
               hypervolume = sc$hypervolume, front_size = sc$front_size,
               best_roi = sc$best_roi, risk_at_best_roi = sc$risk_at_best_roi)
  }
  rows <- list()
  for (r in seq_len(reps)) {
    if (r > 1) invisible(tryCatch(run_one_budgeted(last_algo, last_strategy, last_m, inst, p, risk_alpha,
                                                   nsim_probe, nsim_eval, seed_of(r - 1, last_algo, last_strategy, last_m)),
                                  error = function(e) NULL))
    scen_big <- matrix(rtri(nsim_posthoc * n, inst$acts$optimistic, inst$acts$pessimistic,
                            inst$acts$realistic), nrow = nsim_posthoc, byrow = TRUE)
    for (strategy in SCALARIZATION_REPAIRS) {
      res <- run_tabu_normalized_budgeted(strategy, inst, p, risk_alpha, nsim_probe, nsim_eval,
                                          seed_of(r, "Tabu", strategy, 2L))
      rows[[length(rows) + 1]] <- score_row(res, scen_big, r, strategy, "normalized")
    }
    if (r <= verify_reps) {
      res <- run_one_budgeted("Tabu", "efficiency_only", 2L, inst, p, risk_alpha, nsim_probe, nsim_eval,
                              seed_of(r, "Tabu", "efficiency_only", 2L))
      rows[[length(rows) + 1]] <- score_row(res, scen_big, r, "efficiency_only", "raw_check")
    }
  }
  out <- do.call(rbind, rows)
  saveRDS(out, ckpt)
  out
}

run_scalarization_check <- function(full_production, reps = 10,
                                    n_cores = max(1, parallel::detectCores() - 1), ...) {
  specs <- unique(full_production[, c("instance", "seed", "density")])
  specs <- specs[order(specs$instance), ]
  worker <- function(k) scalarization_check_instance(specs$instance[k], specs$seed[k], specs$density[k],
                                                     reps = reps, ...)
  do.call(rbind, parallel_map(seq_len(nrow(specs)), worker, n_cores))
}

summarise_scalarization_check <- function(scal, full_production) {
  ## Check: raw efficiency-only runs reproduce the stored Design A rows.
  chk <- merge(scal[scal$scalarization == "raw_check", c("instance", "rep", "hypervolume")],
               full_production[full_production$algorithm == "Tabu" & full_production$placement == "multi" &
                               full_production$repair == "efficiency_only", c("instance", "rep", "hypervolume")],
               by = c("instance", "rep"), suffixes = c("_check", "_stored"))
  if (!nrow(chk) || !isTRUE(all.equal(chk$hypervolume_check, chk$hypervolume_stored, tolerance = 1e-10)))
    stop("Raw efficiency-only runs do not reproduce Design A; the pairing is not valid.")
  cat("Check passed:", nrow(chk), "raw runs reproduce the stored Design A rows exactly.\n\n")

  reps_used <- sort(unique(scal$rep))
  raw <- full_production[full_production$algorithm == "Tabu" & full_production$placement == "multi" &
                         full_production$repair %in% SCALARIZATION_REPAIRS & full_production$rep %in% reps_used &
                         full_production$instance %in% unique(scal$instance),
                         c("instance", "rep", "repair", "hypervolume", "front_size", "risk_at_best_roi")]
  raw$scalarization <- "raw"
  norm <- scal[scal$scalarization == "normalized", names(raw)]
  both <- rbind(raw, norm)

  effects <- function(d) {
    hv <- function(r) d$hypervolume[d$repair == r]
    im <- aggregate(hypervolume ~ instance + repair, d, mean)
    g <- function(r) im$hypervolume[im$repair == r][order(im$instance[im$repair == r])]
    dd <- g("risk_aware") - g("efficiency_only")
    data.frame(risk_aware_vs_eff_only   = cliffs_delta(hv("risk_aware"), hv("efficiency_only")),
               eff_random_vs_eff_only   = cliffs_delta(hv("efficiency_random"), hv("efficiency_only")),
               eff_random_vs_risk_aware = cliffs_delta(hv("efficiency_random"), hv("risk_aware")),
               risk_vs_eff_instances    = sprintf("%d / %d / %d", sum(dd > 1e-9), sum(abs(dd) <= 1e-9), sum(dd < -1e-9)))
  }
  tab <- rbind(cbind(scalarization = "raw", effects(both[both$scalarization == "raw", ])),
               cbind(scalarization = "normalized", effects(both[both$scalarization == "normalized", ])))
  tab[, 2:4] <- round(tab[, 2:4], 3)
  cat("Cliff's delta, tabu search, Design A instances, replications", min(reps_used), "to", max(reps_used), "\n")
  cat("Last column: instances where risk-aware repair wins / ties / loses against efficiency-only.\n")
  print(tab, row.names = FALSE)

  desc <- aggregate(cbind(hypervolume, front_size, risk_at_best_roi, single_point = front_size == 1) ~
                      scalarization + repair, both, mean)
  cat("\nMeans by scalarization and strategy:\n")
  print(desc[order(desc$repair, desc$scalarization), ], row.names = FALSE, digits = 3)

  own <- sapply(SCALARIZATION_REPAIRS, function(r)
    cliffs_delta(norm$hypervolume[norm$repair == r], raw$hypervolume[raw$repair == r]))
  cat("\nEffect of normalization itself, normalized vs raw, same strategy (Cliff's delta):\n")
  print(round(own, 3))
  write.csv(tab, "scalarization_check_summary.csv", row.names = FALSE)
  invisible(list(effects = tab, means = desc, normalization_effect = own))
}

## Tabu search, 20 instances x 10 replications x 3 strategies = 600 runs at
## eval_budget = 3000, plus the MOEA/D replays and 40 raw check runs.
scalarization_check <- run_scalarization_check(full_production, reps = 10, n_cores = 6)
saveRDS(scalarization_check, "scalarization_check_results.rds")

summarise_scalarization_check(scalarization_check, full_production)

n6 <- scalarization_check[scalarization_check$scalarization == "normalized", ]
im <- aggregate(hypervolume ~ instance + repair, n6, mean)
g <- function(r) im$hypervolume[im$repair == r][order(im$instance[im$repair == r])]
{cat(sprintf("Normalized: per-instance difference in mean hypervolume, risk_aware minus efficiency_only: %.2f to %.2f\n",
            min(g("risk_aware") - g("efficiency_only")), max(g("risk_aware") - g("efficiency_only"))))
cat(sprintf("Normalized: risk at best value, risk_aware vs efficiency_only, delta = %.3f\n",
            cliffs_delta(n6$risk_at_best_roi[n6$repair == "risk_aware"], n6$risk_at_best_roi[n6$repair == "efficiency_only"])))
cat(sprintf("Normalized: best value, risk_aware vs efficiency_only, delta = %.3f\n",
            cliffs_delta(n6$best_roi[n6$repair == "risk_aware"], n6$best_roi[n6$repair == "efficiency_only"])))}


## 11.7 Instance-level inference

## Cliff's delta through ranks: identical to P(X > Y) - P(X < Y), with ties
## counted as zero, and fast enough for the bootstrap.
cliff_rank <- function(x, y) {
  nx <- length(x); ny <- length(y)
  r <- rank(c(x, y))
  2 * (sum(r[seq_len(nx)]) - nx * (nx + 1) / 2) / (nx * ny) - 1
}

## Matched-pairs rank-biserial correlation of paired differences d.
## Zero differences are dropped, as in the Wilcoxon signed-rank test;
## their number is reported separately as ties.
rank_biserial_paired <- function(d, tol = 1e-9) {
  d <- d[abs(d) > tol]
  if (!length(d)) return(0)
  r <- rank(abs(d))
  (sum(r[d > 0]) - sum(r[d < 0])) / sum(r)
}

wilcox_p <- function(d, tol = 1e-9) {
  if (sum(abs(d) > tol) < 1) return(NA_real_)
  suppressWarnings(wilcox.test(d[abs(d) > tol], mu = 0, exact = FALSE)$p.value)
}

## Paired comparison of strategy a against strategy b for one algorithm.
## d must hold one row per instance, run, and strategy.
compare_pair <- function(d, a, b, outcome = "hypervolume", n_boot = 2000, tol = 1e-9) {
  da <- d[d$repair == a, ]; da <- da[order(da$instance, da$rep), ]
  db <- d[d$repair == b, ]; db <- db[order(db$instance, db$rep), ]
  stopifnot(identical(da$instance, db$instance), identical(da$rep, db$rep))
  xa <- da[[outcome]]; xb <- db[[outcome]]
  ok <- !is.na(xa) & !is.na(xb)
  xa <- xa[ok]; xb <- xb[ok]; inst <- da$instance[ok]

  ## Instance level
  ma <- tapply(xa, inst, mean); mb <- tapply(xb, inst, mean)
  dif <- ma - mb

  ## Cluster bootstrap of the pooled Cliff's delta
  ids <- sort(unique(inst)); rows <- split(seq_along(inst), inst)
  boot <- replicate(n_boot, {
    pick <- unlist(rows[as.character(sample(ids, length(ids), replace = TRUE))], use.names = FALSE)
    cliff_rank(xa[pick], xb[pick])
  })

  data.frame(comparison = paste(a, "vs", b),
             pooled_delta = round(cliff_rank(xa, xb), 3),
             boot_lo = round(unname(quantile(boot, 0.025)), 3),
             boot_hi = round(unname(quantile(boot, 0.975)), 3),
             inst_won = sum(dif > tol), inst_tied = sum(abs(dif) <= tol), inst_lost = sum(dif < -tol),
             mean_diff = round(mean(dif), 3),
             wilcoxon_p = signif(wilcox_p(dif, tol), 3),
             rank_biserial = round(rank_biserial_paired(dif, tol), 3))
}

table1_pairs <- list(c("risk_aware", "efficiency_only"),
                     c("efficiency_random", "efficiency_only"),
                     c("risk_aware_random", "risk_aware"),
                     c("marginal_random", "repair_marginal"),
                     c("penalty", "efficiency_only"),
                     c("efficiency_random", "risk_aware"))

set.seed(2026)

## Friedman's test: 600 run-level blocks against 20 instance-level blocks
dA <- full_production[full_production$placement == "multi", ]
friedman_both <- do.call(rbind, lapply(c("Tabu", "SA", "NSGA2", "MOEAD", "Random"), function(a) {
  da <- dA[dA$algorithm == a, ]
  blk <- aggregate(hypervolume ~ instance + rep + repair, da, mean)
  blk$block <- paste(blk$instance, blk$rep)
  f_run <- friedman.test(hypervolume ~ repair | block, data = blk)
  im <- aggregate(hypervolume ~ instance + repair, da, mean)
  f_inst <- friedman.test(hypervolume ~ repair | instance, data = im)
  f_nopen <- friedman.test(hypervolume ~ repair | instance, data = im[im$repair != "penalty", ])
  data.frame(algorithm = a,
             chisq_600_blocks = round(unname(f_run$statistic), 1),
             chisq_20_instances = round(unname(f_inst$statistic), 1),
             p_20_instances = signif(f_inst$p.value, 3),
             chisq_20_no_penalty = round(unname(f_nopen$statistic), 1),
             p_20_no_penalty = signif(f_nopen$p.value, 3))
}))
{cat("=== Design A, multi-objective placement: Friedman's test across strategies ===\n")
cat("600 run-level blocks and 20 instance means, all seven strategies; last two columns without penalty\n")
print(friedman_both, row.names = FALSE)}

##  Paired comparisons at the instance level (Table 1 pairs, plus
##  randomized efficiency against risk-aware repair)
instance_level_A <- do.call(rbind, lapply(c("Tabu", "SA", "NSGA2", "MOEAD", "Random"), function(a) {
  da <- dA[dA$algorithm == a, ]
  cbind(algorithm = a, do.call(rbind, lapply(table1_pairs, function(p) compare_pair(da, p[1], p[2]))))
}))
{cat("\n=== Design A, multi-objective placement: hypervolume, instance-level comparisons ===\n")
cat("pooled_delta: Cliff's delta over the 600 runs (unpaired, descriptive)\n")
cat("boot_lo, boot_hi: 95% cluster-bootstrap interval, instances resampled\n")
cat("inst_won/tied/lost: instances whose mean favours the first strategy / tie exactly / favour the second\n")
cat("wilcoxon_p, rank_biserial: paired over the 20 instance means (zero differences dropped)\n\n")
print(instance_level_A, row.names = FALSE)}

## Single-objective placement, penalty excluded: risk and value of the
## best solution, risk-aware against efficiency-only repair
dS <- full_production[full_production$placement == "single" & full_production$repair != "penalty", ]
instance_level_single <- do.call(rbind, lapply(c("Tabu", "SA", "GA_single", "Random"), function(a) {
  da <- dS[dS$algorithm == a, ]
  rbind(cbind(algorithm = a, outcome = "risk_at_best_roi",
              compare_pair(da, "risk_aware", "efficiency_only", "risk_at_best_roi")),
        cbind(algorithm = a, outcome = "best_roi",
              compare_pair(da, "risk_aware", "efficiency_only", "best_roi")))
}))
{cat("\n=== Design A, single-objective placement, penalty excluded: risk_aware vs efficiency_only ===\n")
print(instance_level_single, row.names = FALSE)}

## Design B. Each instance is run once, so the paired Wilcoxon tests of
## Section 11.3 are already at the instance level. Added here: the
## matched-pairs rank-biserial correlation next to Cliff's delta.
design_b_paired <- do.call(rbind, lapply(split(ofat_results, list(ofat_results$factor, ofat_results$level, ofat_results$algorithm), drop = TRUE), function(d) {
  one <- function(a, b) {
    xa <- d[d$repair == a, ]; xa <- xa[order(xa$run), ]
    xb <- d[d$repair == b, ]; xb <- xb[order(xb$run), ]
    dif <- xa$hypervolume - xb$hypervolume
    c(delta = cliff_rank(xa$hypervolume, xb$hypervolume),
      rank_biserial = rank_biserial_paired(dif),
      won = sum(dif > 1e-9), tied = sum(abs(dif) <= 1e-9), lost = sum(dif < -1e-9))
  }
  ra <- one("risk_aware", "efficiency_only"); er <- one("efficiency_random", "risk_aware")
  data.frame(factor = d$factor[1], level = d$level[1], algorithm = d$algorithm[1],
             ra_delta = round(ra["delta"], 3), ra_rank_biserial = round(ra["rank_biserial"], 3),
             ra_won = ra["won"], ra_tied = ra["tied"], ra_lost = ra["lost"],
             er_delta = round(er["delta"], 3), er_rank_biserial = round(er["rank_biserial"], 3),
             er_won = er["won"], er_tied = er["tied"], er_lost = er["lost"])
}))
design_b_paired <- design_b_paired[order(design_b_paired$algorithm, design_b_paired$factor, design_b_paired$level), ]
rownames(design_b_paired) <- NULL
{cat("\n=== Design B, Tabu: risk_aware vs efficiency_only (ra_) and efficiency_random vs risk_aware (er_) ===\n")
print(design_b_paired[design_b_paired$algorithm == "Tabu", -3], row.names = FALSE)
}

## Randomized efficiency against risk-aware repair, Tabu, the 26 conditions
## of the four instance factors (the alpha sweep excluded).
rq2 <- design_b_paired[design_b_paired$algorithm == "Tabu" & design_b_paired$factor != "risk_alpha", ]
{cat(sprintf("\nDesign B, Tabu, %d conditions, efficiency_random vs risk_aware:\n", nrow(rq2)))
cat(sprintf("  Cliff's delta: risk-aware ahead with at least a small effect in %d, negligible in %d, randomized ahead with at least a small effect in %d\n",
            sum(rq2$er_delta <= -0.147), sum(abs(rq2$er_delta) < 0.147), sum(rq2$er_delta >= 0.147)))
cat("  negligible at:", paste(rq2$factor, rq2$level)[abs(rq2$er_delta) < 0.147], sep = "\n    ")
cat(sprintf("  instances: risk-aware ahead on more instances in %d, randomized ahead on more in %d, equal in %d\n",
            sum(rq2$er_lost > rq2$er_won), sum(rq2$er_won > rq2$er_lost), sum(rq2$er_won == rq2$er_lost)))}

## Scalarization check (Section 11.6): instance-level comparisons,
## replications 1 to 10

if (exists("scalarization_check")) {
  sc_norm <- scalarization_check[scalarization_check$scalarization == "normalized", ]
  sc_raw  <- full_production[full_production$algorithm == "Tabu" & full_production$placement == "multi" &
                             full_production$rep <= 10 &
                             full_production$repair %in% c("efficiency_only", "efficiency_random", "risk_aware"), ]
  scal_level <- rbind(
    cbind(aggregation = "raw",        compare_pair(sc_raw,  "risk_aware", "efficiency_only")),
    cbind(aggregation = "raw",        compare_pair(sc_raw,  "efficiency_random", "efficiency_only")),
    cbind(aggregation = "normalized", compare_pair(sc_norm, "risk_aware", "efficiency_only")),
    cbind(aggregation = "normalized", compare_pair(sc_norm, "efficiency_random", "efficiency_only")),
    cbind(aggregation = "normalized", compare_pair(sc_norm, "risk_aware", "efficiency_only", "risk_at_best_roi")),
    cbind(aggregation = "normalized", compare_pair(sc_norm, "risk_aware", "efficiency_only", "best_roi")))
  scal_level$outcome <- c("hypervolume", "hypervolume", "hypervolume", "hypervolume", "risk_at_best_roi", "best_roi")
  cat("\n=== Scalarization check, Tabu, replications 1 to 10: instance-level comparisons ===\n")
  print(scal_level, row.names = FALSE)
}

instance_level_inference <- list(friedman = friedman_both, design_a = instance_level_A,
                                 single = instance_level_single, design_b = design_b_paired,
                                 scalarization = if (exists("scal_level")) scal_level else NULL)
saveRDS(instance_level_inference, "instance_level_inference.rds")
cat("\nSaved instance_level_inference.rds\n")
