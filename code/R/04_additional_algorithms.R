# Additional algorithms and full comparison harness (ODD §9.2–9.3)
# Random search, GA/NSGA-II engine, SA, tabu search, MOEA/D, run_full_comparison().
# Requires: source("01_core_functions.R"), source("03_additional_strategies.R")


## 9.2 Additional algorithms

run_random_search <- function(n_genes, n_eval, repair_fn, evaluate_fn, seed = 1) {
  set.seed(seed)
  archive <- list()
  seen <- new.env(parent = emptyenv())
  for (k in seq_len(n_eval)) {
    x <- repair_fn(sample(c(0L, 1L), n_genes, replace = TRUE))
    key <- paste(x, collapse = "")
    if (!is.null(seen[[key]])) next
    seen[[key]] <- TRUE
    obj <- evaluate_fn(x)
    archive <- archive_insert(archive, x, obj)
  }
  m <- length(archive[[1]]$obj)
  objs <- matrix(vapply(archive, function(p) p$obj, numeric(m)), ncol = m, byrow = TRUE)
  list(population = archive, objectives = objs, pareto_idx = seq_along(archive))
}

## Builds up from nothing: at each step add the not-yet-selected activity
## with the best score (predecessors pulled in via closure) that still fits
## the budget, until nothing more fits. weight_risk > 0 blends a Monte
## Carlo risk term into the construction score; sweeping it traces an
## approximate front for the multi-objective placement.
run_greedy <- function(acts, preds, succ, topo, budget, must_idx, deadline, weight_risk = 0) {
  n <- nrow(acts)
  sel <- closure_predecessors(rep(0L, n), preds, must_idx)
  repeat {
    cr <- cost_roi(sel, acts)
    remaining <- setdiff(seq_len(n), which(sel == 1L))
    if (!length(remaining)) break
    best_gain <- -Inf; best_j <- NA_integer_; best_cand <- NULL
    for (j in remaining) {
      cand <- closure_predecessors(replace(sel, j, 1L), preds, must_idx)
      cand_cost <- cost_roi(cand, acts)["cost"]
      if (cand_cost > budget) next
      added_roi  <- cost_roi(cand, acts)["roi"] - cr["roi"]
      added_cost <- cand_cost - cr["cost"]
      score <- added_roi / max(added_cost, 1e-9)
      if (weight_risk > 0 && !is.null(deadline)) {
        mc <- project_duration_mc(cand, acts, preds, topo, nsim = 25, deadline = deadline)
        score <- (1 - weight_risk) * score - weight_risk * mc$p_over
      }
      if (score > best_gain) { best_gain <- score; best_j <- j; best_cand <- cand }
    }
    if (is.na(best_j)) break
    sel <- best_cand
  }
  sel
}

flip_k_bits <- function(x, k = 1) { idx <- sample(length(x), k); x[idx] <- 1L - x[idx]; x }
dominates <- function(a, b) all(a <= b) && any(a < b)

archive_insert <- function(archive, x, obj, max_size = 200) {
  keep <- vector("logical", length(archive)); dominated_by_new <- FALSE
  for (i in seq_along(archive)) {
    if (dominates(archive[[i]]$obj, obj)) dominated_by_new <- TRUE
    keep[i] <- !dominates(obj, archive[[i]]$obj)
  }
  archive <- archive[keep]
  if (!dominated_by_new) archive[[length(archive) + 1]] <- list(x = x, obj = obj)
  if (length(archive) > max_size) {
    m <- length(archive[[1]]$obj)
    objs <- matrix(vapply(archive, function(p) p$obj, numeric(m)), ncol = m, byrow = TRUE)
    cd <- crowding_distance(objs, seq_len(nrow(objs)))
    keep_idx <- order(-cd[as.character(seq_len(nrow(objs)))])[seq_len(max_size)]
    archive <- archive[keep_idx]
  }
  archive
}

## m == 1: classic accept-better/accept-worse-with-probability SA.
## m == 2: the walk accepts by a randomly re-weighted scalarisation each
## step (so it wanders the front rather than collapsing to one point),
## while every point seen is offered to a separate Pareto archive, which is
## what gets reported.
run_sa <- function(n_genes, n_iter, repair_fn, evaluate_fn, m,
                    t0 = 1, cool = 0.995, k_flip = 2, seed = 1) {
  set.seed(seed)
  x <- repair_fn(sample(c(0L, 1L), n_genes, replace = TRUE)); obj <- evaluate_fn(x)
  archive <- list(list(x = x, obj = obj)); T <- t0
  for (it in seq_len(n_iter)) {
    xn <- repair_fn(flip_k_bits(x, k_flip)); objn <- evaluate_fn(xn)
    if (m == 1L) {
      delta <- objn[1] - obj[1]
    } else {
      w <- runif(1)
      delta <- (w * objn[1] + (1 - w) * objn[2]) - (w * obj[1] + (1 - w) * obj[2])
    }
    if (delta < 0 || runif(1) < exp(-delta / max(T, 1e-6))) { x <- xn; obj <- objn }
    archive <- archive_insert(archive, xn, objn)
    T <- T * cool
  }
  objs <- matrix(vapply(archive, function(p) p$obj, numeric(m)), ncol = m, byrow = TRUE)
  list(population = archive, objectives = objs, pareto_idx = fast_nondominated_sort(objs)[[1]])
}

## Neighbourhood = a random sample of single-bit-flip moves each iteration.
## A flipped bit is tabu for `tenure` iterations unless the move improves on
## the best (or archive) seen so far (aspiration).
run_tabu <- function(n_genes, n_iter, repair_fn, evaluate_fn, m,
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
      scal_n <- if (m == 1L) objn[1] else sum(objn)
      if (scal_n < best_delta) {
        best_delta <- scal_n; best_move <- bit; best_move_x <- xn; best_move_obj <- objn
      }
    }
    if (!is.null(best_move)) {
      ## Keep the repaired candidate that was evaluated. Repairing the raw move
      ## a second time can return a different solution whenever repair is
      ## stochastic (the softmax strategies, and risk-aware repair through its
      ## Monte Carlo probe), which would pair one solution with the objective
      ## of another.
      x <- best_move_x; obj <- best_move_obj
      tabu_until[best_move] <- it + tenure
      if (m == 1L) { if (obj[1] < best_obj[1]) best_obj <- obj } else if (dominates(obj, best_obj)) best_obj <- obj
      archive <- archive_insert(archive, x, obj)
    }
  }
  objs <- matrix(vapply(archive, function(p) p$obj, numeric(m)), ncol = m, byrow = TRUE)
  list(population = archive, objectives = objs, pareto_idx = fast_nondominated_sort(objs)[[1]])
}

## K evenly spaced weight vectors decompose the two objectives into K
## scalar subproblems solved via Tchebycheff scalarisation against a
## running ideal point z. Reproduction draws parents from each
## subproblem's neighbourhood (its closest weight vectors); a child can
## replace any neighbour it scalarises better for.
run_moead <- function(n_genes, K, n_gen, repair_fn, evaluate_fn,
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
  z <- apply(do.call(rbind, lapply(pop, function(p) p$obj)), 2, min)
  tcheby <- function(obj, wk) max(wk * abs(obj - z))
  crossover <- function(a, b) {
    if (runif(1) > pcx) return(a)
    pt <- sample(2:(n_genes - 1), 1); c(a[1:pt], b[(pt + 1):n_genes])
  }
  mutate <- function(a) { flip <- runif(n_genes) < pmut; a[flip] <- 1L - a[flip]; a }
  for (g in seq_len(n_gen)) for (k in seq_len(K)) {
    parents <- sample(nb[[k]], 2)
    child_x <- repair_fn(mutate(crossover(pop[[parents[1]]]$x, pop[[parents[2]]]$x)))
    child_obj <- evaluate_fn(child_x)
    z <- pmin(z, child_obj)
    for (j in nb[[k]]) if (tcheby(child_obj, w[j, ]) <= tcheby(pop[[j]]$obj, w[j, ]))
      pop[[j]] <- list(x = child_x, obj = child_obj)
  }
  objs <- matrix(vapply(pop, function(p) p$obj, numeric(2)), ncol = 2, byrow = TRUE)
  list(population = pop, objectives = objs, pareto_idx = fast_nondominated_sort(objs)[[1]])
}


## 9.3 The comparison harness

ALGO_SCOPE <- list(Random = c("single", "multi"), GA_single = c("single"), NSGA2 = c("multi"),
                    SA = c("single", "multi"), Tabu = c("single", "multi"), MOEAD = c("multi"))
REPAIRS <- c("efficiency_only", "efficiency_random", "risk_aware", "risk_aware_random",
             "repair_marginal", "marginal_random", "penalty")
run_one <- function(algo, strategy, m, acts, preds, succ, topo, budget, must_idx, deadline,
                     risk_alpha, nsim_probe, nsim_eval, seed,
                     pop_size, n_gen, sa_iter, tabu_iter, random_eval, moead_K) {
  n <- nrow(acts)
  repfn <- function(x) apply_repair_strategy(x, acts, preds, succ, topo, budget, must_idx,
                                              deadline, strategy = strategy,
                                              risk_alpha = risk_alpha, nsim_probe = nsim_probe)
  evalfn <- make_evaluate(m, acts, preds, topo, budget, deadline, nsim_eval = nsim_eval,
                           penalize = (strategy == "penalty"))
  switch(algo,
    Random    = run_random_search(n, random_eval, repfn, evalfn, seed = seed),
    GA_single = run_evolution(n, pop_size, n_gen, 0.7, 0.15, repfn, evalfn, seed = seed),
    NSGA2     = run_evolution(n, pop_size, n_gen, 0.7, 0.15, repfn, evalfn, seed = seed),
    SA        = run_sa(n, sa_iter, repfn, evalfn, m = m, seed = seed),
    Tabu      = run_tabu(n, tabu_iter, repfn, evalfn, m = m, seed = seed),
    MOEAD     = run_moead(n, moead_K, n_gen, repfn, evalfn, seed = seed)
  )
}

## Reduces a run to a fair, comparable outcome: the final population/archive
## is re-scored on TRUE roi/risk against a large shared scenario matrix,
## restricted to the FEASIBLE subset -- this matters specifically for
## "penalty", whose raw individuals need not respect budget; every other
## strategy is feasible by construction already.
## value/risk tradeoff endpoints alongside
## the existing hypervolume/front_size/feasible_fraction. 
## min_risk and roi_at_min_risk are added specifically to capture the OTHER end of
## the tradeoff, not just the value-maximizing corner.
score_result <- function(res, acts, preds, topo, budget, deadline, scen_big) {
  xs <- unique(lapply(res$population, function(p) p$x))
  rows <- lapply(xs, function(x) {
    mc <- project_duration_mc(x, acts, preds, topo, deadline = deadline, scenario_matrix = scen_big)
    c(roi = unname(cost_roi(x, acts)["roi"]), risk = mc$p_over,
      feasible = as.numeric(cost_roi(x, acts)["cost"] <= budget))
  })
  df <- as.data.frame(do.call(rbind, rows))
  feas_frac <- mean(df$feasible)
  df <- df[df$feasible == 1, c("roi", "risk"), drop = FALSE]
  if (!nrow(df)) return(list(hypervolume = 0, front_size = 0, feasible_fraction = feas_frac,
                              best_roi = NA, risk_at_best_roi = NA, min_risk = NA, roi_at_min_risk = NA))
  dom <- vapply(seq_len(nrow(df)), function(i)
    any(df$roi >= df$roi[i] & df$risk <= df$risk[i] & (df$roi[i] < df$roi | df$risk[i] > df$risk)), logical(1))
  front <- df[!dom, , drop = FALSE]
  best_idx <- which.max(front$roi)
  safest_idx <- which.min(front$risk)
  list(hypervolume = hypervolume_2d(as.matrix(front)), front_size = nrow(front), feasible_fraction = feas_frac,
       best_roi = front$roi[best_idx], risk_at_best_roi = front$risk[best_idx],
       min_risk = front$risk[safest_idx], roi_at_min_risk = front$roi[safest_idx])
}

run_full_comparison <- function(acts, preds, succ, topo, budget, must_idx, deadline,
                                 risk_alpha = 0.25, nsim_probe = 20, nsim_eval = 40, nsim_posthoc = 300,
                                 reps = 3, base_seed = 1000, pop_size = 20, n_gen = 10,
                                 sa_iter = 120, tabu_iter = 60, random_eval = 150, moead_K = 10) {
  n <- nrow(acts); rows <- list()
  for (r in seq_len(reps)) {
    scen_big <- matrix(rtri(nsim_posthoc * n, acts$optimistic, acts$pessimistic, acts$realistic),
                        nrow = nsim_posthoc, byrow = TRUE)
    for (algo in names(ALGO_SCOPE)) for (placement in ALGO_SCOPE[[algo]]) {
      m <- if (placement == "single") 1L else 2L
      for (strategy in REPAIRS) {
        seed <- base_seed + r * 10000 + which(names(ALGO_SCOPE) == algo) * 100 +
          which(REPAIRS == strategy) * 10 + m
        t0 <- proc.time()[3]
        res <- run_one(algo, strategy, m, acts, preds, succ, topo, budget, must_idx, deadline,
                        risk_alpha, nsim_probe, nsim_eval, seed,
                        pop_size, n_gen, sa_iter, tabu_iter, random_eval, moead_K)
        wall <- proc.time()[3] - t0
        sc <- score_result(res, acts, preds, topo, budget, deadline, scen_big)
        rows[[length(rows) + 1]] <- data.frame(
          rep = r, algorithm = algo, placement = placement, repair = strategy,
          hypervolume = sc$hypervolume, front_size = sc$front_size, best_roi = sc$best_roi,
          risk_at_best_roi = sc$risk_at_best_roi, feasible_fraction = sc$feasible_fraction, wall_sec = wall)
      }
    }
  }
  do.call(rbind, rows)
}
