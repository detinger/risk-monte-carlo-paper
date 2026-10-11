# Illustration: repair divergence and near-exact validation (ODD §8)
# Requires: source("01_core_functions.R")


## 8. Illustration: repair divergence and near-exact validation on one instance

n <- 15
acts <- generate_activities(n, must_k = 2, seed = 11)
dag  <- generate_dag(n, density = 0.15, seed = 11)
topo <- topo_order(dag$preds)
must_idx <- which(acts$must)
floor_cost <- must_closure_cost(acts, dag$preds, must_idx)
budget <- floor_cost + 0.5 * (sum(acts$cost) - floor_cost)
deadline <- 45
cat("must-closure floor cost:", floor_cost, " | chosen budget:", budget, "\n")


## 8.1 The precondition from section 5, made concrete

## An under-provisioned budget: lower than the must-closure floor. No
## repair strategy can succeed here, and the existing repair_to_budget_graph
## function in the earlier .Rmd files does not warn when this happens; it
## just returns an over-budget "repaired" solution.
tight_budget <- floor_cost - 5
stuck <- repair_activity_network(rep(1L, n), acts, dag$preds, dag$succ, topo,
                                  tight_budget, must_idx, deadline = deadline, alpha = 1)
cat("cost achieved:", cost_roi(stuck, acts)["cost"], " vs. infeasible budget:", tight_budget, "\n")


## 8.2 Efficiency-only vs. risk-aware repair, same instance

eff_only <- repair_activity_network(rep(1L, n), acts, dag$preds, dag$succ, topo,
                                     budget, must_idx, deadline = deadline, alpha = 1)
set.seed(77)
risk_aware <- repair_activity_network(rep(1L, n), acts, dag$preds, dag$succ, topo,
                                       budget, must_idx, deadline = deadline,
                                       alpha = 0.25, nsim_probe = 60)

scen_big <- matrix(rtri(600 * n, acts$optimistic, acts$pessimistic, acts$realistic),
                    nrow = 600, byrow = TRUE)
mc_eff  <- project_duration_mc(eff_only,   acts, dag$preds, topo, deadline = deadline, scenario_matrix = scen_big)
mc_risk <- project_duration_mc(risk_aware, acts, dag$preds, topo, deadline = deadline, scenario_matrix = scen_big)

cat("efficiency-only : kept =", sum(eff_only),   " cost =", cost_roi(eff_only,   acts)["cost"],
    " value =", cost_roi(eff_only,   acts)["roi"], " P(overrun) =", mc_eff$p_over,  "\n")

cat("risk-aware      : kept =", sum(risk_aware), " cost =", cost_roi(risk_aware, acts)["cost"],
    " value =", cost_roi(risk_aware, acts)["roi"], " P(overrun) =", mc_risk$p_over, "\n")

cat("selections differ at activities:", which(eff_only != risk_aware), "\n")


## 8.3 Check against the near-exact reference

set.seed(42)
scen_ref <- matrix(rtri(300 * n, acts$optimistic, acts$pessimistic, acts$realistic), nrow = 300, byrow = TRUE)
t0 <- proc.time()[3]
ref <- enumerate_reference(acts, dag$preds, dag$succ, budget, must_idx, topo, deadline, scen_ref, max_n = 20)
cat("enumeration wall time (s):", round(proc.time()[3] - t0, 1),
    " | distinct feasible closed sets:", ref$n_feasible_closed_sets,
    " | Pareto-optimal points:", nrow(ref$pareto_front), "\n")

cat("reference hypervolume:", hypervolume_2d(ref$pareto_front), "\n")


## 8.4 The 2x2 design on the same instance

design <- run_factorial_design(acts, dag$preds, dag$succ, topo, budget, must_idx, deadline,
                                risk_alpha = 0.25, pop_size = 50, n_gen = 10, reps = 10, base_seed = 500)
print(aggregate(cbind(hypervolume, best_roi, risk_at_best_roi) ~ repair + objective,
                 data = design, FUN = mean))

d <- design[design$objective == "multi_ROI_risk", ]
wt <- wilcox.test(d$hypervolume[d$repair == "risk_aware"],
                   d$hypervolume[d$repair == "efficiency_only"],
                   paired = TRUE, alternative = "greater")
cat("paired Wilcoxon, risk-aware > efficiency-only hypervolume, p =", round(wt$p.value, 4),
    " (n =", nrow(d)/2, "pairs -- far too few to claim significance; see section 9)\n")

best_multi_hv <- max(aggregate(hypervolume ~ repair, data = d, FUN = mean)$hypervolume)
cat("best factorial cell recovers", round(100 * best_multi_hv / hypervolume_2d(ref$pareto_front), 1),
    "% of the reference hypervolume\n")
