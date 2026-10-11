# Core functions (ODD §7)
# Triangular sampling, critical-path, closure, repair operator,
# evolutionary engine, near-exact enumeration, pilot 2x2 harness.


## 7.1 Triangular sampling, activity and network generation

rtri <- function(n, a, b, c) {
  stopifnot(all(a <= c), all(c <= b))
  u  <- runif(n)
  Fc <- (c - a) / (b - a)
  ifelse(u < Fc,
         a + sqrt(u * (b - a) * (c - a)),
         b - sqrt((1 - u) * (b - a) * (b - c)))
}

generate_activities <- function(n, must_k = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  must <- rep(FALSE, n)
  must[sample(seq_len(n), must_k)] <- TRUE
  data.frame(
    id          = seq_len(n),
    cost        = sample(5:40, n, replace = TRUE),
    optimistic  = sample(1:10,  n, replace = TRUE),
    realistic   = sample(10:20, n, replace = TRUE),
    pessimistic = sample(20:35, n, replace = TRUE),
    roi         = round(runif(n, 1.0, 3.0), 2), # value
    must        = must
  )
}

## Edges only run i -> j for i < j, so the graph is acyclic by construction.
generate_dag <- function(n, density, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  succ <- vector("list", n); preds <- vector("list", n)
  for (i in seq_len(n)) preds[[i]] <- integer(0)
  for (i in 1:(n - 1)) {
    cand <- (i + 1):n
    choose <- cand[runif(length(cand)) < density]
    if (!length(choose)) choose <- i + 1L   # guarantee every non-terminal node has a successor
    succ[[i]] <- sort(unique(choose))
  }
  succ[[n]] <- integer(0)
  for (i in 1:n) for (j in succ[[i]]) preds[[j]] <- sort(unique(c(preds[[j]], i)))
  list(preds = preds, succ = succ)
}

topo_order <- function(preds) {
  n <- length(preds)
  indeg <- vapply(preds, length, integer(1))
  S <- which(indeg == 0L); order <- integer(0)
  succ <- vector("list", n); for (i in seq_len(n)) succ[[i]] <- integer(0)
  for (j in seq_len(n)) if (length(preds[[j]])) for (p in preds[[j]]) succ[[p]] <- c(succ[[p]], j)
  while (length(S)) {
    v <- S[1]; S <- S[-1]; order <- c(order, v)
    for (w in succ[[v]]) { indeg[w] <- indeg[w] - 1L; if (indeg[w] == 0L) S <- c(S, w) }
  }
  if (length(order) != n) stop("Graph is not a DAG.")
  order
}


## A complete instance (ODD section 10.2): activities, network, budget, and a
## deadline at a quantile of the efficiency-only reference solution's duration.
## Uses repair_activity_network() and project_duration_mc(), defined below.
make_instance <- function(n, seed, density, must_k = 3, budget_fraction = 0.5,
                           deadline_quantile = 0.4, nsim_calib = 500) {
  acts <- generate_activities(n, must_k = must_k, seed = seed)
  dag  <- generate_dag(n, density = density, seed = seed)
  topo <- topo_order(dag$preds)
  must_idx <- which(acts$must)
  floor_cost <- must_closure_cost(acts, dag$preds, must_idx)
  budget <- floor_cost + budget_fraction * (sum(acts$cost) - floor_cost)

  ref_sel <- repair_activity_network(rep(1L, n), acts, dag$preds, dag$succ, topo,
                                      budget, must_idx, deadline = NULL, alpha = 1)
  ref_mc <- project_duration_mc(ref_sel, acts, dag$preds, topo, nsim = nsim_calib)
  deadline <- as.numeric(quantile(ref_mc$draws, deadline_quantile))

  list(acts = acts, preds = dag$preds, succ = dag$succ, topo = topo, must_idx = must_idx,
       budget = budget, deadline = deadline, floor_cost = floor_cost,
       seed = seed, density = density, budget_fraction = budget_fraction,
       deadline_quantile = deadline_quantile)
}


## 7.2 Critical-path forward pass and Monte Carlo risk

project_duration_once <- function(dur, preds, keep, topo) {
  n <- length(dur); keep <- as.integer(keep > 0)
  idx_keep <- which(keep == 1L)
  if (!length(idx_keep)) return(0)
  ES <- numeric(n); EF <- rep(NA_real_, n)
  for (i in topo) {
    if (keep[i] == 0L) next
    p <- preds[[i]]
    if (length(p)) {
      pk <- p[keep[p] == 1L]
      ES[i] <- if (length(pk)) suppressWarnings(max(EF[pk], na.rm = TRUE)) else 0
      if (!is.finite(ES[i])) ES[i] <- 0
    } else ES[i] <- 0
    EF[i] <- ES[i] + dur[i]
  }
  mx <- suppressWarnings(max(EF[idx_keep], na.rm = TRUE))
  if (!is.finite(mx)) 0 else mx
}

## `scenario_matrix`, if supplied, fixes the random duration draws so that
## different repair calls or different methods being compared see the same
## randomness (common random numbers) instead of independently resampled
## noise.
project_duration_mc <- function(sel, acts, preds, topo, nsim = 100,
                                 deadline = NULL, scenario_matrix = NULL) {
  n <- length(sel)
  if (is.null(scenario_matrix)) {
    scenario_matrix <- matrix(rtri(nsim * n, acts$optimistic, acts$pessimistic, acts$realistic),
                               nrow = nsim, byrow = TRUE)
  }
  nsim <- nrow(scenario_matrix)
  draws <- vapply(seq_len(nsim), function(s) {
    project_duration_once(scenario_matrix[s, ], preds, sel, topo)
  }, numeric(1))
  p_over <- if (is.null(deadline)) NA_real_ else mean(draws > deadline)
  list(mean = mean(draws), sd = stats::sd(draws), p_over = p_over, draws = draws)
}


## 7.3 Precedence closure and feasibility bookkeeping

closure_predecessors <- function(sel, preds, must_idx = integer(0)) {
  n <- length(sel); sel <- as.integer(sel)
  if (length(must_idx)) sel[must_idx] <- 1L
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    for (i in seq_len(n)) if (sel[i] == 1L && length(preds[[i]])) {
      need <- preds[[i]][sel[preds[[i]]] == 0L]
      if (length(need)) { sel[need] <- 1L; changed <- TRUE }
    }
  }
  sel
}

leaves_in_set <- function(sel, succ) {
  which(sel == 1L & vapply(seq_along(sel), function(i) {
    if (sel[i] == 0L) return(FALSE)
    s <- succ[[i]]
    if (!length(s)) return(TRUE)
    !any(sel[s] == 1L)
  }, logical(1)))
}

cost_roi <- function(sel, acts) {
  idx <- which(sel == 1L)
  if (!length(idx)) return(c(cost = 0, roi = 0))
  c(cost = sum(acts$cost[idx]), roi = sum(acts$roi[idx]))
}

## Feasibility precondition described in section 5: the bare must-closure has
## to fit the budget, or no repair strategy can ever succeed.
must_closure_cost <- function(acts, preds, must_idx) {
  cost_roi(closure_predecessors(rep(0L, nrow(acts)), preds, must_idx), acts)["cost"]
}


## 7.4 The repair operator

## Risk signal of each removable leaf, before normalization over L.
##   "mc"            the Monte Carlo risk reduction delta_p (the paper's rule)
##   "shuffled"      delta_p randomly permuted across the leaves: same values
##                   and cost, no link to which leaf carries the risk (placebo)
##   "critical_path" deterministic proxy: the leaf's finish time in a forward
##                   pass with PERT mean durations (o + 4m + p) / 6; no sampling
leaf_risk_signal <- function(L, sel, acts, preds, topo, deadline, nsim_probe, risk_signal) {
  if (risk_signal == "critical_path") {
    pert <- (acts$optimistic + 4 * acts$realistic + acts$pessimistic) / 6
    ES <- numeric(length(sel)); EF <- numeric(length(sel))
    for (i in topo) if (sel[i] == 1L) {
      pk <- preds[[i]][sel[preds[[i]]] == 1L]
      ES[i] <- if (length(pk)) max(EF[pk]) else 0
      EF[i] <- ES[i] + pert[i]
    }
    return(EF[L])
  }
  scen <- matrix(rtri(nsim_probe * length(sel), acts$optimistic, acts$pessimistic, acts$realistic),
                  nrow = nsim_probe, byrow = TRUE)
  p_before <- project_duration_mc(sel, acts, preds, topo, deadline = deadline,
                                   scenario_matrix = scen)$p_over
  delta_p <- vapply(L, function(i) {
    cand <- sel; cand[i] <- 0L
    p_after <- project_duration_mc(cand, acts, preds, topo, deadline = deadline,
                                    scenario_matrix = scen)$p_over
    p_before - p_after
  }, numeric(1))
  if (risk_signal == "shuffled" && length(L) > 1) delta_p <- delta_p[sample.int(length(L))]
  delta_p
}

## `trace`, if supplied, is an environment from new_cut_trace() that counts,
## over all risk-scored cuts: the cuts, the leaves they chose from, the cuts
## where the risk signal had zero range over L (so the score fell back to
## efficiency alone), and the cuts that removed a different leaf than
## efficiency-only repair would have. It draws no random numbers, so results
## are the same with or without it.
new_cut_trace <- function() {
  tr <- new.env()
  tr$cuts <- 0; tr$leaves <- 0; tr$zero_range <- 0; tr$differs_from_eff <- 0
  tr
}
trace_cut <- function(trace, L, signal, cut, eff_norm) {
  if (is.null(trace)) return(invisible())
  trace$cuts <- trace$cuts + 1
  trace$leaves <- trace$leaves + length(L)
  trace$zero_range <- trace$zero_range + (diff(range(signal)) == 0)
  trace$differs_from_eff <- trace$differs_from_eff + (cut != L[which.min(eff_norm[L])])
}

repair_activity_network <- function(sel, acts, preds, succ, topo,
                                     budget, must_idx = integer(0),
                                     deadline = NULL, alpha = 1,
                                     nsim_probe = 40, risk_signal = "mc", trace = NULL) {
  sel <- closure_predecessors(sel, preds, must_idx)
  cr  <- cost_roi(sel, acts) # value/cost
  if (cr["cost"] <= budget) return(sel)

  eff <- acts$roi / pmax(acts$cost, 1e-9)
  eff_rng <- range(eff)
  eff_norm <- if (diff(eff_rng) > 0) (eff - eff_rng[1]) / diff(eff_rng) else rep(0.5, length(eff))

  repeat {
    cr <- cost_roi(sel, acts)
    if (cr["cost"] <= budget) break
    L <- setdiff(leaves_in_set(sel, succ), must_idx)
    if (!length(L)) break  # at the must-closure floor; see section 5

    if (alpha >= 1 || is.null(deadline)) {
      cut <- L[which.min(eff_norm[L])]
    } else {
      delta_p <- leaf_risk_signal(L, sel, acts, preds, topo, deadline, nsim_probe, risk_signal)
      dp_rng <- range(delta_p)
      dp_norm <- if (diff(dp_rng) > 0) (delta_p - dp_rng[1]) / diff(dp_rng) else rep(0, length(L))
      ## cut_score is high for leaves that are both low-value and whose
      ## removal buys a meaningful drop in deadline-overrun risk.
      cut_score <- alpha * (1 - eff_norm[L]) + (1 - alpha) * dp_norm
      cut <- L[which.max(cut_score)]
      trace_cut(trace, L, delta_p, cut, eff_norm)
    }
    sel[cut] <- 0L
  }
  sel
}

## Stochastic twin of the function above, added in section 12 once the
## collapse investigation pointed at deterministic repair's diversity loss
## (Martins & Ribas, 2021) as a contributing factor. Same cut_score logic;
## the leaf actually cut is drawn with probability proportional to
## exp(score * temperature) instead of chosen by strict argmax. temperature
## controls how close to deterministic (high) vs uniform-random (low) this
## is. Every repair strategy below has a _random counterpart built on this.
repair_activity_network_randomized <- function(sel, acts, preds, succ, topo,
                                                 budget, must_idx = integer(0),
                                                 deadline = NULL, alpha = 1,
                                                 nsim_probe = 40, temperature = 8) {
  sel <- closure_predecessors(sel, preds, must_idx)
  cr  <- cost_roi(sel, acts)
  if (cr["cost"] <= budget) return(sel)

  eff <- acts$roi / pmax(acts$cost, 1e-9)
  eff_rng <- range(eff)
  eff_norm <- if (diff(eff_rng) > 0) (eff - eff_rng[1]) / diff(eff_rng) else rep(0.5, length(eff))

  repeat {
    cr <- cost_roi(sel, acts)
    if (cr["cost"] <= budget) break
    L <- setdiff(leaves_in_set(sel, succ), must_idx)
    if (!length(L)) break

    if (alpha >= 1 || is.null(deadline)) {
      score <- 1 - eff_norm[L]
    } else {
      scen <- matrix(rtri(nsim_probe * length(sel), acts$optimistic, acts$pessimistic, acts$realistic),
                      nrow = nsim_probe, byrow = TRUE)
      p_before <- project_duration_mc(sel, acts, preds, topo, deadline = deadline, scenario_matrix = scen)$p_over
      delta_p <- vapply(L, function(i) {
        cand <- sel; cand[i] <- 0L
        p_before - project_duration_mc(cand, acts, preds, topo, deadline = deadline, scenario_matrix = scen)$p_over
      }, numeric(1))
      dp_rng <- range(delta_p)
      dp_norm <- if (diff(dp_rng) > 0) (delta_p - dp_rng[1]) / diff(dp_rng) else rep(0, length(L))
      score <- alpha * (1 - eff_norm[L]) + (1 - alpha) * dp_norm
    }
    w <- exp(score * temperature); w <- w / sum(w)
    cut <- L[sample.int(length(L), 1, prob = w)]
    sel[cut] <- 0L
  }
  sel
}


## 7.5 A single evolutionary engine for both objective placements

fast_nondominated_sort <- function(objs) {
  n <- nrow(objs)
  dom_count <- integer(n); dom_by <- vector("list", n)
  fronts <- list(integer(0))
  for (p in seq_len(n)) {
    dom_by[[p]] <- integer(0)
    for (q in seq_len(n)) if (p != q) {
      if (all(objs[p,] <= objs[q,]) && any(objs[p,] < objs[q,])) {
        dom_by[[p]] <- c(dom_by[[p]], q)
      } else if (all(objs[q,] <= objs[p,]) && any(objs[q,] < objs[p,])) {
        dom_count[p] <- dom_count[p] + 1L
      }
    }
    if (dom_count[p] == 0L) fronts[[1]] <- c(fronts[[1]], p)
  }
  k <- 1L
  while (length(fronts[[k]])) {
    nxt <- integer(0)
    for (p in fronts[[k]]) for (q in dom_by[[p]]) {
      dom_count[q] <- dom_count[q] - 1L
      if (dom_count[q] == 0L) nxt <- c(nxt, q)
    }
    k <- k + 1L
    fronts[[k]] <- nxt
  }
  fronts[seq_len(k - 1L)]
}

crowding_distance <- function(objs, idx) {
  m <- ncol(objs); k <- length(idx)
  if (k <= 2) { d <- rep(Inf, k); names(d) <- idx; return(d) }
  d <- setNames(numeric(k), idx)
  for (j in seq_len(m)) {
    ord <- idx[order(objs[idx, j])]
    d[as.character(ord[1])] <- Inf
    d[as.character(ord[k])] <- Inf
    rng <- objs[ord[k], j] - objs[ord[1], j]
    if (rng <= 0) next
    for (t in 2:(k - 1)) {
      d[as.character(ord[t])] <- d[as.character(ord[t])] +
        (objs[ord[t + 1], j] - objs[ord[t - 1], j]) / rng
    }
  }
  d
}

crowded_tournament <- function(rank, crowd, i, j) {
  if (rank[i] != rank[j]) return(if (rank[i] < rank[j]) i else j)
  if (crowd[i] != crowd[j]) return(if (crowd[i] > crowd[j]) i else j)
  sample(c(i, j), 1)
}

run_evolution <- function(n_genes, pop_size, n_gen, pcx, pmut,
                           repair_fn, evaluate_fn, seed = 1) {
  set.seed(seed)
  init <- function() sample(c(0L, 1L), n_genes, replace = TRUE, prob = c(0.45, 0.55))
  repair_and_eval <- function(x) { x <- repair_fn(x); list(x = x, obj = evaluate_fn(x)) }

  pop <- lapply(seq_len(pop_size), function(i) repair_and_eval(init()))
  m <- length(pop[[1]]$obj)
  ## matrix(..., byrow = TRUE), not t(vapply(...)), because the latter
  ## returns a 1 x pop_size matrix rather than pop_size x 1 when m == 1.
  objs <- matrix(vapply(pop, function(p) p$obj, numeric(m)), ncol = m, byrow = TRUE)

  crossover <- function(a, b) {
    if (runif(1) > pcx) return(a)
    pt <- sample(2:(n_genes - 1), 1)
    c(a[1:pt], b[(pt + 1):n_genes])
  }
  mutate <- function(a) { flip <- runif(n_genes) < pmut; a[flip] <- 1L - a[flip]; a }

  for (g in seq_len(n_gen)) {
    fronts <- fast_nondominated_sort(objs)
    rank <- integer(nrow(objs)); crowd <- numeric(nrow(objs))
    for (f in seq_along(fronts)) {
      rank[fronts[[f]]] <- f
      cd <- crowding_distance(objs, fronts[[f]])
      crowd[fronts[[f]]] <- cd[as.character(fronts[[f]])]
    }
    offspring <- vector("list", pop_size)
    for (k in seq_len(pop_size)) {
      i <- sample(pop_size, 1); j <- sample(pop_size, 1); p1 <- crowded_tournament(rank, crowd, i, j)
      i <- sample(pop_size, 1); j <- sample(pop_size, 1); p2 <- crowded_tournament(rank, crowd, i, j)
      offspring[[k]] <- repair_and_eval(mutate(crossover(pop[[p1]]$x, pop[[p2]]$x)))
    }
    off_objs <- matrix(vapply(offspring, function(p) p$obj, numeric(m)), ncol = m, byrow = TRUE)
    all_ind <- c(pop, offspring); all_objs <- rbind(objs, off_objs)
    fronts2 <- fast_nondominated_sort(all_objs)
    next_idx <- integer(0)
    for (f in fronts2) {
      if (length(next_idx) + length(f) <= pop_size) {
        next_idx <- c(next_idx, f)
      } else {
        cd <- crowding_distance(all_objs, f)
        ord <- f[order(-cd[as.character(f)])]
        next_idx <- c(next_idx, ord[seq_len(pop_size - length(next_idx))])
        break
      }
    }
    pop <- all_ind[next_idx]; objs <- all_objs[next_idx, , drop = FALSE]
  }
  fronts <- fast_nondominated_sort(objs)
  list(population = pop, objectives = objs, pareto_idx = fronts[[1]])
}


## 7.6 Small-instance near-exact reference

enumerate_reference <- function(acts, preds, succ, budget, must_idx, topo,
                                 deadline, scenario_matrix, max_n = 20) {
  n <- nrow(acts)
  if (n > max_n) stop(sprintf("n = %d exceeds the exhaustive-enumeration limit (%d).", n, max_n))
  seen <- new.env(parent = emptyenv())
  rows <- list()
  for (code in 0:(2^n - 1)) {
    raw <- as.integer(intToBits(code))[1:n]
    closed <- closure_predecessors(raw, preds, must_idx)
    key <- paste(closed, collapse = "")
    if (!is.null(seen[[key]])) next
    seen[[key]] <- TRUE
    cr <- cost_roi(closed, acts)
    if (cr["cost"] > budget) next
    mc <- project_duration_mc(closed, acts, preds, topo, deadline = deadline,
                               scenario_matrix = scenario_matrix)
    rows[[length(rows) + 1]] <- data.frame(roi = unname(cr["roi"]), cost = unname(cr["cost"]),
                                            risk = mc$p_over, n_kept = sum(closed),
                                            stringsAsFactors = FALSE)
  }
  df <- do.call(rbind, rows)
  dom <- vapply(seq_len(nrow(df)), function(i) {
    any(df$roi >= df$roi[i] & df$risk <= df$risk[i] & (df$roi[i] < df$roi | df$risk[i] > df$risk))
  }, logical(1))
  list(all_feasible = df, pareto_front = df[!dom, ], n_feasible_closed_sets = nrow(df))
}

hypervolume_2d <- function(pf, ref = c(0, 1)) {
  ord <- order(-pf[, "roi"]); pf <- pf[ord, , drop = FALSE]
  hv <- 0; prev_risk <- ref[2]
  for (i in seq_len(nrow(pf))) {
    width <- pf[i, "roi"] - ref[1]; height <- prev_risk - pf[i, "risk"]
    if (width > 0 && height > 0) hv <- hv + width * height
    prev_risk <- min(prev_risk, pf[i, "risk"])
  }
  hv
}


## 7.7 The 2x2 factorial harness (pilot design)

run_factorial_design <- function(acts, preds, succ, topo, budget, must_idx, deadline,
                                  risk_alpha = 0.3, pop_size = 24, n_gen = 12,
                                  pcx = 0.7, pmut = 0.15, nsim_probe = 20,
                                  nsim_eval = 40, nsim_posthoc = 400,
                                  reps = 5, base_seed = 100) {
  n <- nrow(acts)
  cells <- expand.grid(repair = c("efficiency_only", "risk_aware"),
                        objective = c("single_ROI", "multi_ROI_risk"),
                        stringsAsFactors = FALSE)
  make_repair <- function(kind) {
    a <- if (kind == "efficiency_only") 1 else risk_alpha
    function(x) repair_activity_network(x, acts, preds, succ, topo, budget, must_idx,
                                         deadline = deadline, alpha = a, nsim_probe = nsim_probe)
  }
  eval_single <- function(x) c(-cost_roi(x, acts)["roi"])
  eval_multi  <- function(x) {
    mc <- project_duration_mc(x, acts, preds, topo, nsim = nsim_eval, deadline = deadline)
    c(-cost_roi(x, acts)["roi"], mc$p_over)
  }
  rows <- list()
  for (r in seq_len(reps)) {
    posthoc_scen <- matrix(rtri(nsim_posthoc * n, acts$optimistic, acts$pessimistic, acts$realistic),
                            nrow = nsim_posthoc, byrow = TRUE)
    for (c_i in seq_len(nrow(cells))) {
      kind <- cells$repair[c_i]; obj_kind <- cells$objective[c_i]
      m <- if (obj_kind == "single_ROI") 1L else 2L
      eval_fn <- if (m == 1L) eval_single else eval_multi
      seed <- base_seed + r * 1000 + c_i
      t0 <- proc.time()[3]
      res <- run_evolution(n, pop_size, n_gen, pcx, pmut,
                            repair_fn = make_repair(kind), evaluate_fn = eval_fn, seed = seed)
      wall <- proc.time()[3] - t0
      final_x <- unique(lapply(res$population, function(p) p$x))
      pts <- do.call(rbind, lapply(final_x, function(x) {
        mc <- project_duration_mc(x, acts, preds, topo, deadline = deadline, scenario_matrix = posthoc_scen)
        c(roi = unname(cost_roi(x, acts)["roi"]), risk = mc$p_over)
      }))
      dom <- vapply(seq_len(nrow(pts)), function(i) {
        any(pts[,"roi"] >= pts[i,"roi"] & pts[,"risk"] <= pts[i,"risk"] &
            (pts[i,"roi"] < pts[,"roi"] | pts[i,"risk"] > pts[,"risk"]))
      }, logical(1))
      front <- pts[!dom, , drop = FALSE]
      hv <- hypervolume_2d(front)
      best_i <- which.max(front[, "roi"])
      rows[[length(rows) + 1]] <- data.frame(
        rep = r, repair = kind, objective = obj_kind, hypervolume = hv,
        front_size = nrow(front), best_roi = front[best_i, "roi"],
        risk_at_best_roi = front[best_i, "risk"], wall_sec = wall, stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
