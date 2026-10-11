# Full-scale analysis: Design A results and validation (ODD §10.7–10.8)
# Descriptive summaries, Table 1 construction, cost, instance-level heterogeneity,
# small-instance validation against exhaustive enumeration.
# Requires: full_production and ofat_results loaded (run 06 and 08 first)


## 10.7 Summary of the full-scale comparison

## Descriptive summary: mean outcome per algorithm x repair x placement.

## Each metric is averaged SEPARATELY per algorithm x repair cell, not via a
## single aggregate() call across multiple response columns. 

## aggregate() with cbind(a, b, c) performs listwise deletion: a row is dropped 
## from EVERY column's mean if ANY of the cbind'd columns is NA for that row,
## even columns with no NAs of their own. 

## best_roi and risk_at_best_roi are NA whenever no feasible solution was found
## (all of penalty's failed runs, for instance);
## aggregating them together with hypervolume would
## silently drop those rows from the hypervolume mean too, producing a
## survivorship-biased number for penalty specifically. Meaning each
## column independently, with na.rm = TRUE, avoids this entirely.
mean_by_cell <- function(data, cols) {
  out <- data.frame()
  for (a in unique(data$algorithm)) for (r in unique(data$repair)) {
    dd <- data[data$algorithm == a & data$repair == r, ]
    row <- data.frame(algorithm = a, repair = r)
    for (col in cols) row[[col]] <- mean(dd[[col]], na.rm = TRUE)
    out <- rbind(out, row)
  }
  out
}

summary_multi  <- mean_by_cell(full_production[full_production$placement == "multi", ],
                                c("hypervolume", "front_size", "best_roi", "risk_at_best_roi", "feasible_fraction"))
summary_single <- mean_by_cell(full_production[full_production$placement == "single", ],
                                c("best_roi", "risk_at_best_roi", "feasible_fraction"))
write.csv(summary_multi,  "summary_multi_placement.csv",  row.names = FALSE)
write.csv(summary_single, "summary_single_placement.csv", row.names = FALSE)
cat("--- multi-objective placement: mean hypervolume, front size, best value, risk at best value, feasible fraction ---\n"); print(summary_multi)

cat("\n--- single-objective placement: mean best value, risk at best value, feasible fraction ---\n"); print(summary_single)

## Density-binned breakdown: Design A is not a factorial design over
## density the way the OFAT sweep (Section 11) is - each instance draws
## its own density independently, 0.05 to 0.25. Tertile binning turns this
## into a natural three-level comparison without a separate experiment.
full_production$density_bin <- cut(full_production$density,
  breaks = quantile(full_production$density, c(0, 1/3, 2/3, 1)),
  labels = c("low", "mid", "high"), include.lowest = TRUE)
cliffs_delta <- function(x, y) (sum(outer(x, y, ">")) - sum(outer(x, y, "<"))) / (length(x) * length(y))
d_multi <- full_production[full_production$placement == "multi", ]
delta_by_density <- data.frame()
for (a in unique(d_multi$algorithm)) for (b in c("low","mid","high")) {
  dd <- d_multi[d_multi$algorithm == a & d_multi$density_bin == b, ]
  eff  <- dd$hypervolume[dd$repair == "efficiency_only"]
  risk <- dd$hypervolume[dd$repair == "risk_aware"]
  delta_by_density <- rbind(delta_by_density, data.frame(algorithm = a, density_bin = b,
                                                          delta = round(cliffs_delta(risk, eff), 3)))
}
write.csv(delta_by_density, "delta_by_density_multi.csv", row.names = FALSE)
cat("\n--- risk_aware vs efficiency_only, Cliff's delta by density tertile ---\n")

print(reshape(delta_by_density, idvar = "algorithm", timevar = "density_bin", direction = "wide"), row.names = FALSE)

## Inferential statistics: printed via the shared helper from section 10.4,
## and the analysis objects themselves saved 
res_multi_prod  <- analyse_repair_effect(full_production, "multi",  "hypervolume")
res_single_prod <- analyse_repair_effect(full_production, "single", "best_roi")
saveRDS(list(multi = res_multi_prod, single = res_single_prod), "full_production_analysis.rds")
print_repair_analysis(res_multi_prod,  "PRODUCTION: multi-objective placement, hypervolume")

print_repair_analysis(res_single_prod, "PRODUCTION: single-objective placement, best_roi")

## Diagnostic plots, saved to PDF -.
#pdf("production_boxplots.pdf", width = 10, height = 6)
for (pl in c("multi", "single")) {
  outcome <- if (pl == "multi") "hypervolume" else "best_roi"
  d <- full_production[full_production$placement == pl, ]
  for (a in unique(d$algorithm)) {
    da <- d[d$algorithm == a, ]
    boxplot(da[[outcome]] ~ da$repair, main = paste(a, "-", pl, "placement"),
            ylab = outcome, xlab = "repair strategy")
  }
}

#dev.off()

cat("\nSaved: full_production_results.rds, full_production_analysis.rds,\n",
    "       summary_multi_placement.csv, summary_single_placement.csv,\n",
    "       delta_by_density_multi.csv, production_boxplots.pdf\n")

res <- full_production
d <- res[res$placement == "multi", ]

mag <- function(x) ifelse(abs(x) < 0.147, "negligible",
                    ifelse(abs(x) < 0.33, "small",
                    ifelse(abs(x) < 0.474, "medium", "large")))

pairs <- list(
  c("efficiency_only", "risk_aware"),
  c("efficiency_only", "efficiency_random"),
  c("risk_aware", "risk_aware_random"),
  c("repair_marginal", "marginal_random"),
  c("efficiency_only", "penalty")
)
pair_names <- c("risk_aware vs eff_only", "eff_random vs eff_only",
                 "risk_aware_random vs risk_aware", "marginal_random vs marginal", "penalty vs eff_only")

algos <- c("Tabu", "Random", "NSGA2", "SA", "MOEAD")
out <- data.frame(algorithm = algos)
for (i in seq_along(pairs)) {
  p <- pairs[[i]]
  vals <- sapply(algos, function(a) {
    da <- d[d$algorithm == a, ]
    x <- da$hypervolume[da$repair == p[2]]
    y <- da$hypervolume[da$repair == p[1]]
    dd <- cliffs_delta(x, y)
    sprintf("%.3f (%s)", dd, mag(dd))
  })
  out[[pair_names[i]]] <- vals
}
print(out, row.names = FALSE)

write.csv(out, "table1_complete.csv", row.names = FALSE)

## front_size vs hypervolume correlation, per algorithm, across its own
## seven strategies. Confirms whether a collapsed front is a precondition
## for low hypervolume specifically for Tabu, or a general property.

summary_for_corr <- mean_by_cell(d_multi, c("hypervolume", "front_size"))
front_size_hv_correlation <- data.frame()
for (a in unique(summary_for_corr$algorithm)) {
  da <- summary_for_corr[summary_for_corr$algorithm == a, ]
  front_size_hv_correlation <- rbind(front_size_hv_correlation, data.frame(
    algorithm = a, correlation = round(cor(da$front_size, da$hypervolume), 3),
    front_size_min = round(min(da$front_size), 2), front_size_max = round(max(da$front_size), 2)))
}
write.csv(front_size_hv_correlation, "front_size_hypervolume_correlation.csv", row.names = FALSE)
cat("--- front_size vs hypervolume correlation, within each algorithm's 7 strategies ---\n")

print(front_size_hv_correlation, row.names = FALSE)

## Computational cost: risk-aware repair's overhead by algorithm.

cost_table <- data.frame()
for (a in unique(d_multi$algorithm)) {
  da <- d_multi[d_multi$algorithm == a, ]
  eff_sec  <- mean(da$wall_sec[da$repair == "efficiency_only"], na.rm = TRUE)
  risk_sec <- mean(da$wall_sec[da$repair == "risk_aware"], na.rm = TRUE)
  cost_table <- rbind(cost_table, data.frame(algorithm = a, eff_sec = round(eff_sec, 2),
                                              risk_sec = round(risk_sec, 2), ratio = round(risk_sec / eff_sec, 1)))
}
write.csv(cost_table, "computational_cost_by_algorithm.csv", row.names = FALSE)
cat("\n--- wall_sec: efficiency_only vs risk_aware, multi placement ---\n")

print(cost_table, row.names = FALSE)

## Risk-reduction side effect in single-objective placement: does
## risk_aware repair produce a lower-risk solution even when nothing in
## the search objective asks for it? risk_at_best_roi has NAs only for
## penalty (confirmed separately); efficiency_only and risk_aware have
## none, so this comparison is not exposed to the listwise-deletion issue.

d_single <- full_production[full_production$placement == "single", ]
single_risk_table <- data.frame()
for (a in unique(d_single$algorithm)) {
  da <- d_single[d_single$algorithm == a, ]
  eff_r  <- mean(da$risk_at_best_roi[da$repair == "efficiency_only"], na.rm = TRUE)
  risk_r <- mean(da$risk_at_best_roi[da$repair == "risk_aware"], na.rm = TRUE)
  single_risk_table <- rbind(single_risk_table, data.frame(algorithm = a, eff_risk = round(eff_r, 3),
                                                            risk_aware_risk = round(risk_r, 3), diff = round(eff_r - risk_r, 3)))
}
write.csv(single_risk_table, "single_placement_risk_reduction.csv", row.names = FALSE)
cat("\n--- risk_at_best_roi, single placement: efficiency_only vs risk_aware ---\n")

print(single_risk_table, row.names = FALSE)

## Full feasible_fraction table, all algorithms x repairs -- confirms
## penalty is the only strategy with any infeasibility issue.

feas_table <- mean_by_cell(d_multi, "feasible_fraction")
write.csv(feas_table, "feasible_fraction_full_table.csv", row.names = FALSE)
cat("\n--- feasible_fraction, multi placement, all combinations ---\n")

print(reshape(feas_table, idvar = "algorithm", timevar = "repair", direction = "wide"), row.names = FALSE)

## Instance-level heterogeneity in Tabu's multi-objective effect: is
## the aggregate delta (0.576) a moderate, consistent improvement, or does
## it hide a more extreme underlying pattern?

dt <- d_multi[d_multi$algorithm == "Tabu", ]
inst_delta_table <- data.frame()
for (i in sort(unique(dt$instance))) {
  di <- dt[dt$instance == i, ]
  eff  <- di$hypervolume[di$repair == "efficiency_only"]
  risk <- di$hypervolume[di$repair == "risk_aware"]
  inst_delta_table <- rbind(inst_delta_table, data.frame(instance = i, delta = round(cliffs_delta(risk, eff), 3)))
}
write.csv(inst_delta_table, "tabu_instance_level_deltas.csv", row.names = FALSE)
cat("\n--- per-instance Cliff's delta, Tabu, risk_aware vs efficiency_only ---\n")

print(inst_delta_table, row.names = FALSE)

n_inst <- nrow(inst_delta_table)
cat(sprintf("\n%d/%d instances at delta=1.000 (complete dominance), %d/%d at delta=0.000 (no effect), %d/%d elsewhere\n",
            sum(inst_delta_table$delta == 1), n_inst, sum(inst_delta_table$delta == 0), n_inst,
            sum(inst_delta_table$delta != 1 & inst_delta_table$delta != 0), n_inst))

cat("\nSaved: front_size_hypervolume_correlation.csv, computational_cost_by_algorithm.csv,\n",
    "       single_placement_risk_reduction.csv, feasible_fraction_full_table.csv,\n",
    "       tabu_instance_level_deltas.csv\n")


## 10.8 Small-instance validation against exhaustive enumeration

validate_one_instance <- function(seed, n = 18, density = 0.10, must_k = 2,
                                  budget_fraction = 0.5, deadline_quantile = 0.4,
                                  algos = c("Tabu", "SA", "NSGA2", "MOEAD", "Random"),
                                  strategies = c("efficiency_only", "efficiency_random", "risk_aware"),
                                  reps = 3, eval_budget = 3000, ga_pop = 100, tabu_neigh = 10,
                                  moead_K = 100, risk_alpha = 0.25, nsim_probe = 30,
                                  nsim_eval = 60, nsim_ref = 300,
                                  checkpoint_dir = "checkpoints_validation") {
  dir.create(checkpoint_dir, showWarnings = FALSE)
  ckpt <- file.path(checkpoint_dir, sprintf("validation_%d_%s.rds", seed, config_tag()))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  inst <- make_instance(n, seed = seed, density = density, must_k = must_k,
                        budget_fraction = budget_fraction, deadline_quantile = deadline_quantile)
  set.seed(seed + 1e5)
  scen <- matrix(rtri(nsim_ref * n, inst$acts$optimistic, inst$acts$pessimistic, inst$acts$realistic),
                 nrow = nsim_ref, byrow = TRUE)
  ref <- enumerate_reference(inst$acts, inst$preds, inst$succ, inst$budget, inst$must_idx,
                             inst$topo, inst$deadline, scen, max_n = 20)
  ref_hv <- hypervolume_2d(ref$pareto_front)
  p <- derive_algo_params(eval_budget, ga_pop = ga_pop, tabu_neigh = tabu_neigh, moead_K = moead_K)

  rows <- list()
  for (r in seq_len(reps)) for (a in algos) for (st in strategies) {
    seed_run <- seed * 1000 + r * 100 + which(algos == a) * 10 + which(strategies == st)
    res <- run_one_budgeted(a, st, 2L, inst, p, risk_alpha, nsim_probe, nsim_eval, seed_run)
    sc <- score_result(res, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen)
    rows[[length(rows) + 1]] <- data.frame(
      seed = seed, rep = r, algorithm = a, repair = st,
      feasible_sets = ref$n_feasible_closed_sets, ref_front = nrow(ref$pareto_front),
      ref_hv = ref_hv, hv = sc$hypervolume, front_size = sc$front_size,
      recovery = sc$hypervolume / ref_hv)
  }
  out <- do.call(rbind, rows)
  saveRDS(out, ckpt)
  out
}

validate_against_reference <- function(seeds = 8001:8005, ...,
                                       n_cores = max(1, parallel::detectCores() - 1)) {
  worker <- function(s) validate_one_instance(s, ...)
  do.call(rbind, parallel_map(seeds, worker, n_cores))
}

summarise_recovery <- function(v) {
  inst <- unique(v[, c("seed", "feasible_sets", "ref_front", "ref_hv")])
  cat("Reference instances:\n"); print(inst, row.names = FALSE)
  s <- aggregate(recovery ~ algorithm + repair, v, function(x)
    c(min = min(x), median = median(x), exact = mean(x > 0.999)))
  s <- do.call(data.frame, s)
  names(s) <- c("algorithm", "repair", "min_recovery", "median_recovery", "share_exact")
  cat("\nRecovery of the reference hypervolume, by algorithm and repair strategy:\n")
  print(s[order(s$algorithm, s$repair), ], row.names = FALSE, digits = 3)
  cat(sprintf("\nOverall: %d runs, minimum recovery %.3f, %.0f%% of runs recover the reference exactly.\n",
              nrow(v), min(v$recovery), 100 * mean(v$recovery > 0.999)))
  write.csv(s, "validation_recovery_summary.csv", row.names = FALSE)
  invisible(s)
}

validation <- validate_against_reference(seeds = 8001:8005, n_cores = 5)
saveRDS(validation, "validation_results.rds")

summarise_recovery(validation)

## Single-objective placement with penalty excluded. analyse_repair_effect()
## keeps complete blocks only. A penalty run with no feasible solution has
## best_roi = NA, so the previous single-objective analysis is restricted to
## the 120 blocks (4 instances) where penalty succeeds, and all strategies
## tie on those blocks. Excluding penalty keeps all 600 blocks.
res_single_nopen <- analyse_repair_effect(full_production[full_production$repair != "penalty", ],
                                          "single", "best_roi")
print_repair_analysis(res_single_nopen, "PRODUCTION: single-objective placement, best_roi, penalty excluded")

## Share of blocks in which best value differs across the six repair strategies.
d_single_rep <- full_production[full_production$placement == "single" &
                                full_production$repair != "penalty", ]
varies <- aggregate(best_roi ~ algorithm + instance + rep, d_single_rep,
                    function(x) diff(range(x)) > 1e-9)
cat("\n--- share of single-objective blocks where best value differs across repair strategies ---\n")

print(aggregate(best_roi ~ algorithm, varies, mean))

## Effect size for risk at the best-value solution, risk_aware vs efficiency_only.
single_risk_delta <- do.call(rbind, lapply(unique(d_single_rep$algorithm), function(a) {
  da <- d_single_rep[d_single_rep$algorithm == a, ]
  data.frame(algorithm = a,
             delta = round(cliffs_delta(da$risk_at_best_roi[da$repair == "risk_aware"],
                                        da$risk_at_best_roi[da$repair == "efficiency_only"]), 3))
}))
cat("\n--- risk at best value, risk_aware vs efficiency_only, Cliff's delta ---\n")

print(single_risk_delta, row.names = FALSE)

write.csv(single_risk_delta, "single_placement_risk_delta.csv", row.names = FALSE)
