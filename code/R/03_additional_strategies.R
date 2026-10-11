# Additional repair strategies (ODD §9.1)
# MARG and MARG-SM cutting criteria.
# Requires: source("01_core_functions.R")


## 9.1 Additional repair strategies

## repair_marginal cuts by absolute value, not value/cost. Because value is
## additive over the selected set, the marginal value lost by removing leaf i
## is exactly roi_i regardless of what else is selected, so this is a
## genuinely different, still-simple rule from the ratio-based one.
repair_marginal <- function(sel, acts, preds, succ, budget, must_idx) {
  sel <- closure_predecessors(sel, preds, must_idx)
  repeat {
    if (cost_roi(sel, acts)["cost"] <= budget) break
    L <- setdiff(leaves_in_set(sel, succ), must_idx)
    if (!length(L)) break
    cut <- L[which.min(acts$roi[L])]
    sel[cut] <- 0L
  }
  sel
}

## Stochastic version of repair_marginal, same pattern as
## repair_activity_network_randomized above.
repair_marginal_randomized <- function(sel, acts, preds, succ, budget, must_idx = integer(0), temperature = 8) {
  sel <- closure_predecessors(sel, preds, must_idx)
  repeat {
    if (cost_roi(sel, acts)["cost"] <= budget) break
    L <- setdiff(leaves_in_set(sel, succ), must_idx)
    if (!length(L)) break
    roi_rng <- range(acts$roi)
    roi_norm <- if (diff(roi_rng) > 0) (acts$roi - roi_rng[1]) / diff(roi_rng) else rep(0.5, length(acts$roi))
    score <- 1 - roi_norm[L]
    w <- exp(score * temperature); w <- w / sum(w)
    cut <- L[sample.int(length(L), 1, prob = w)]
    sel[cut] <- 0L
  }
  sel
}

## Seven strategies: three deterministic criteria (efficiency, risk-aware,
## marginal), each with a randomized twin, plus penalty (which has no
## selection criterion to randomize). Two more controls for the revision
## (R/11): risk_shuffled (RAND-RISK, placebo) and critical_path (CRIT,
## deterministic proxy); see leaf_risk_signal() in R/01. They are kept out
## of REPAIRS so that the seeds of the original runs do not change.
apply_repair_strategy <- function(sel, acts, preds, succ, topo, budget, must_idx,
                                   deadline, strategy, risk_alpha = 0.25, nsim_probe = 20, temperature = 8,
                                   trace = NULL) {
  switch(strategy,
    efficiency_only    = repair_activity_network(sel, acts, preds, succ, topo, budget, must_idx,
                                                  deadline = deadline, alpha = 1),
    efficiency_random  = repair_activity_network_randomized(sel, acts, preds, succ, topo, budget, must_idx,
                                                  deadline = deadline, alpha = 1, temperature = temperature),
    risk_aware         = repair_activity_network(sel, acts, preds, succ, topo, budget, must_idx,
                                                  deadline = deadline, alpha = risk_alpha, nsim_probe = nsim_probe,
                                                  trace = trace),
    risk_shuffled      = repair_activity_network(sel, acts, preds, succ, topo, budget, must_idx,
                                                  deadline = deadline, alpha = risk_alpha, nsim_probe = nsim_probe,
                                                  risk_signal = "shuffled", trace = trace),
    critical_path      = repair_activity_network(sel, acts, preds, succ, topo, budget, must_idx,
                                                  deadline = deadline, alpha = risk_alpha,
                                                  risk_signal = "critical_path", trace = trace),
    risk_aware_random  = repair_activity_network_randomized(sel, acts, preds, succ, topo, budget, must_idx,
                                                  deadline = deadline, alpha = risk_alpha, nsim_probe = nsim_probe,
                                                  temperature = temperature),
    repair_marginal    = repair_marginal(sel, acts, preds, succ, budget, must_idx),
    marginal_random    = repair_marginal_randomized(sel, acts, preds, succ, budget, must_idx, temperature = temperature),
    penalty            = closure_predecessors(sel, preds, must_idx),
    stop("unknown repair strategy: ", strategy)
  )
}

## For "penalty", precedence stays a hard constraint (closure), but budget
## is priced into fitness rather than repaired. Every other strategy passes
## penalize = FALSE: their repair step already guarantees feasibility, so
## there is nothing left to price in.
make_evaluate <- function(m, acts, preds, topo, budget, deadline,
                           nsim_eval = 40, penalize = FALSE, penalty_coef = 0.08) {
  pen <- function(cost) if (penalize) penalty_coef * max(0, cost - budget) else 0
  if (m == 1L) {
    function(x) { cr <- cost_roi(x, acts); c(-cr["roi"] + pen(cr["cost"])) }
  } else {
    function(x) {
      cr <- cost_roi(x, acts)
      mc <- project_duration_mc(x, acts, preds, topo, nsim = nsim_eval, deadline = deadline)
      c(-cr["roi"] + pen(cr["cost"]), mc$p_over)
    }
  }
}
