# OFAT sensitivity sweep (ODD §11.1–11.3)
# Design B: n, density, budget_fraction, deadline_q, risk_alpha.
# Requires: source("04_additional_algorithms.R")


## 11.2 Running sweep

OFAT_FACTORS <- list(
  n                 = list(baseline = 50,   levels = c(14, 30, 60, 75, 150, 300)),
  density           = list(baseline = 0.1, levels = c(0.02, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30)),
  budget_fraction   = list(baseline = 0.50, levels = seq(0.05, 0.95, by = 0.15)),
  deadline_quantile = list(baseline = 0.40, levels = seq(0.20, 0.95, by = 0.15)),
  risk_alpha        = list(baseline = 0.25, levels = seq(0, 1, by = 0.2))
)
## The baseline values of density (0.10) and budget_fraction (0.50) are also
## levels of their own sweeps. Those two settings are the same configuration,
## run on the same instance seeds (seed_base + 1..20) with the same run seeds,
## so they return identical results: they are one measurement shown twice,
## not two independent confirmations. The baselines of n, deadline_quantile
## and risk_alpha are not levels of their sweeps; ofat_baseline_rows() below
## takes the baseline cell from the density sweep.
OFAT_ALGOS <- c("NSGA2", "SA", "Tabu", "MOEAD")
OFAT_REPAIRS <- REPAIRS

config_tag <- function() {
  s <- paste(c(OFAT_ALGOS, OFAT_REPAIRS), collapse = "|")
  sprintf("%05d", sum(utf8ToInt(s)) %% 100000)
}

## Rows of the baseline configuration, taken from the first factor whose
## levels include its baseline value.
ofat_baseline_rows <- function(ofat_results) {
  for (f in names(OFAT_FACTORS)) {
    b <- OFAT_FACTORS[[f]]$baseline
    if (any(abs(OFAT_FACTORS[[f]]$levels - b) < 1e-9))
      return(ofat_results[ofat_results$factor == f & abs(ofat_results$level - b) < 1e-9, ])
  }
  stop("no factor has its baseline among its levels")
}

## One (factor, level) setting: RUNS_PER_SETTING fresh random instances,
## each run once (reps = 1 per instance here; the instances themselves are
## the replication). Checkpointed the same way as
## section 10, one file per (factor, level).

run_ofat_setting <- function(factor_name, level_value, runs_per_setting = 20,
                              eval_budget = 6000, nsim_probe = 15, nsim_eval = 40,
                              nsim_posthoc = 300, seed_base = 5000,
                              checkpoint_dir = "checkpoints_ofat") {
  dir.create(checkpoint_dir, showWarnings = FALSE)
  ckpt <- file.path(checkpoint_dir, sprintf("ofat_%s_%s_%s.rds", factor_name, gsub("\\.", "p", level_value), config_tag()))
  if (file.exists(ckpt)) return(readRDS(ckpt))

  base <- lapply(OFAT_FACTORS, `[[`, "baseline")
  base[[factor_name]] <- level_value
  n <- base$n
  must_k <- max(2, round(n / 15))
  p <- derive_algo_params(eval_budget, ga_pop = 40, tabu_neigh = 8, moead_K = 30)

  rows <- list()
  for (u in seq_len(runs_per_setting)) {
    seed <- seed_base + u
    inst <- make_instance(n, seed = seed, density = base$density, must_k = must_k,
                           budget_fraction = base$budget_fraction,
                           deadline_quantile = base$deadline_quantile)
    scen_big <- matrix(rtri(nsim_posthoc * n, inst$acts$optimistic, inst$acts$pessimistic,
                             inst$acts$realistic), nrow = nsim_posthoc, byrow = TRUE)
    for (algo in OFAT_ALGOS) for (strategy in OFAT_REPAIRS) {
      seed_run <- seed * 1000 + which(OFAT_REPAIRS == strategy)
      res <- tryCatch(
        run_one_budgeted(algo, strategy, 2L, inst, p, base$risk_alpha, nsim_probe, nsim_eval, seed_run),
        error = function(e) { message(sprintf("FAILED factor=%s level=%s algo=%s strategy=%s: %s",
                                               factor_name, level_value, algo, strategy, conditionMessage(e))); NULL })
      if (is.null(res)) next
      sc <- score_result(res, inst$acts, inst$preds, inst$topo, inst$budget, inst$deadline, scen_big)
      rows[[length(rows) + 1]] <- data.frame(
        factor = factor_name, level = level_value, run = u, algorithm = algo, repair = strategy,
        hypervolume = sc$hypervolume, feasible_fraction = sc$feasible_fraction,
        front_size = sc$front_size, best_roi = sc$best_roi, risk_at_best_roi = sc$risk_at_best_roi,
        min_risk = sc$min_risk, roi_at_min_risk = sc$roi_at_min_risk)
    }
  }
  out <- do.call(rbind, rows)
  saveRDS(out, ckpt)
  out
}

## Full sweep: every factor x every level. Runs OFAT_ALGOS x OFAT_REPAIRS
## (4 x 7 = 28 combinations) x runs_per_setting instances, once per level.
## Parallelised across (factor, level) settings, the natural independent
## unit here.

run_full_ofat <- function(runs_per_setting = 20, eval_budget = 6000, n_cores = max(1, detectCores() - 1)) {
  settings <- do.call(rbind, lapply(names(OFAT_FACTORS), function(f) {
    data.frame(factor = f, level = OFAT_FACTORS[[f]]$levels)
  }))
  worker <- function(i) {
    run_ofat_setting(settings$factor[i], settings$level[i],
                      runs_per_setting = runs_per_setting, eval_budget = eval_budget)
  }
  do.call(rbind, parallel_map(seq_len(nrow(settings)), worker, n_cores))
}

ofat_results <- run_full_ofat(runs_per_setting = 20, eval_budget = 6000, n_cores = 6)
saveRDS(ofat_results, "ofat_results.rds")

expected_settings <- sum(vapply(OFAT_FACTORS, function(f) length(f$levels), integer(1)))
got_settings <- length(unique(interaction(ofat_results$factor, ofat_results$level, drop = TRUE)))
cat("OFAT sweep collected", nrow(ofat_results), "rows across", got_settings,
    "of", expected_settings, "(factor, level) settings\n")

if (got_settings < expected_settings) {
  warning(sprintf("%d settings produced no data -- check checkpoints_ofat/ for missing or empty files",
                   expected_settings - got_settings))
}


## 11.3 Summary of the sweep

summary_tab <- aggregate(hypervolume ~ factor + level + algorithm + repair, data = ofat_results, FUN = mean)
write.csv(summary_tab, "ofat_summary.csv", row.names = FALSE)
cat("--- OFAT summary (also saved to ofat_summary.csv) ---\n")

print(summary_tab)

## One panel per factor x algorithm: x = level, y = mean hypervolume,
## coloured by repair strategy. 

repair_style <- data.frame(
  repair = c("efficiency_only", "efficiency_random", "risk_aware", "risk_aware_random",
             "repair_marginal", "marginal_random", "penalty"),
  col    = c("steelblue", "steelblue", "firebrick", "firebrick",
             "forestgreen", "forestgreen", "grey40"),
  lty    = c(1, 2, 1, 2, 1, 2, 1),
  stringsAsFactors = FALSE
)

#pdf("ofat_plots.pdf", width = 8, height = 6)
for (f in unique(summary_tab$factor)) {
  d <- summary_tab[summary_tab$factor == f, ]
  ## n spans 14 to 300, unevenly spaced -- on a linear axis the 14-75 range
  ## (where Tabu's effect jumps from small to large) gets squeezed into
  ## about a fifth of the plot width; log-scale gives it roughly half.
  log_x <- (f == "n")
  for (a in unique(d$algorithm)) {
    da <- d[d$algorithm == a, ]
    plot(NULL, xlim = range(da$level), ylim = range(da$hypervolume), main = paste(f, "-", a),
         xlab = f, ylab = "mean hypervolume", log = if (log_x) "x" else "")
    for (rp in unique(da$repair)) {
      dd <- da[da$repair == rp, ]; dd <- dd[order(dd$level), ]
      sty <- repair_style[repair_style$repair == rp, ]
      lines(dd$level, dd$hypervolume, type = "b", pch = 16, col = sty$col, lty = sty$lty)
    }
    legend("topleft", legend = repair_style$repair, col = repair_style$col, lty = repair_style$lty,
           bty = "n", cex = 0.7)
  }
}

#dev.off()
#cat("Saved: ofat_summary.csv, ofat_plots.pdf\n")

## Figure 1 in the paper: risk-aware vs. efficiency-only repair, Cliff's delta by n,
## one panel per algorithm. 

## One row per (algorithm, n level): Cliff's delta for risk_aware vs
## efficiency_only, computed from the raw paired runs, not from
## pre-aggregated means.

d_n <- ofat_results[ofat_results$factor == "n", ]
algos <- c("NSGA2", "MOEAD", "SA", "Tabu")   # panel order: flattest to steepest
levels_n <- sort(unique(d_n$level))

delta_tab <- do.call(rbind, lapply(algos, function(a) {
  do.call(rbind, lapply(levels_n, function(lv) {
    dd <- d_n[d_n$algorithm == a & d_n$level == lv, ]
    eff  <- dd$hypervolume[dd$repair == "efficiency_only"]
    risk <- dd$hypervolume[dd$repair == "risk_aware"]
    data.frame(algorithm = a, n = lv, delta = cliffs_delta(risk, eff))
  }))
}))

## Baseline n (not a level of the n sweep), from the baseline cell, drawn
## as an open marker so it is not read as one of the sweep's own levels.
d_base <- ofat_baseline_rows(ofat_results)
n_base <- OFAT_FACTORS$n$baseline
delta_base <- do.call(rbind, lapply(algos, function(a) {
  dd <- d_base[d_base$algorithm == a, ]
  data.frame(algorithm = a, n = n_base,
             delta = cliffs_delta(dd$hypervolume[dd$repair == "risk_aware"],
                                  dd$hypervolume[dd$repair == "efficiency_only"]))
}))

## Publication figure: 2x3 grid (one panel unused), log-x, effect-size
## reference bands from Section 2.3's thresholds (0.147 / 0.33 / 0.474).
## y-axis range is taken from the data itself, not hardcoded, since a few
## points (NSGA-II, MOEA/D) dip slightly negative and a fixed positive-only
## range would silently clip them off the plot.
y_range <- range(c(delta_tab$delta, delta_base$delta, 0, 0.474), na.rm = TRUE)
y_range <- y_range + c(-1, 1) * 0.05 * diff(y_range)

## Resets the interactive plot pane; skipped in batch runs, where it would
## close a pdf() device opened by the caller.
if (interactive()) grDevices::graphics.off()
#png("figure1_n_sweep.png", width = 900, height = 300)

par(mfrow = c(1,4), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 2, 0))

for (a in algos) {
  dd <- delta_tab[delta_tab$algorithm == a, ]
  dd <- dd[order(dd$n), ]

  plot(NULL, xlim = range(c(levels_n, n_base)), ylim = y_range, log = "x",
       xaxt = "n", xlab = "n", ylab = "Cliff's delta",
       main = a, cex.main = 1.1, font.main = 1)
  axis(1, at = sort(c(levels_n, n_base)), labels = sort(c(levels_n, n_base)))

  ## effect-size reference bands (Romano et al., 2006), drawn first so the
  ## data line sits on top
  ## With log = "x", par("usr")[1:2] are log10 of the x limits.
  ux <- 10^par("usr")[1:2]
  rect(ux[1], 0.474, ux[2], par("usr")[4], col = "grey92", border = NA)
  rect(ux[1], 0.33,  ux[2], 0.474,          col = "grey96", border = NA)
  abline(h = c(0.147, 0.33, 0.474), col = "grey75", lty = 3, lwd = 0.8)
  abline(h = 0, col = "grey50", lty = 1, lwd = 0.8)
  box()

  lines(dd$n, dd$delta, type = "b", pch = 16, col = "firebrick", lwd = 1.8)
  points(n_base, delta_base$delta[delta_base$algorithm == a], pch = 5, col = "firebrick", cex = 1.2)
  ## Labels for the threshold lines
  text(ux[1], c(0.147, 0.33, 0.474), c("small", "medium", "large"),
       adj = c(-0.1, -0.3), cex = 0.6, col = "grey40", xpd = NA)
}

#mtext("Grey bands: medium (>= 0.33) and large (>= 0.474) effect size", side = 3, line = 0, outer = TRUE, cex = 0.8)
#dev.off()

d_n <- ofat_results[ofat_results$factor == "n", ]
algos <- c("NSGA2", "MOEAD", "SA", "Tabu")
levels_n <- sort(unique(d_n$level))

results <- do.call(rbind, lapply(algos, function(a) {
  do.call(rbind, lapply(levels_n, function(lv) {
    dd <- d_n[d_n$algorithm == a & d_n$level == lv, ]
    eff  <- dd$hypervolume[dd$repair == "efficiency_only"]
    risk <- dd$hypervolume[dd$repair == "risk_aware"]
    wt <- tryCatch(wilcox.test(risk, eff, paired = TRUE)$p.value, error = function(e) NA)
    data.frame(algorithm = a, n = lv, n_pairs = length(eff), wilcox_p = wt)
  }))
}))
print(results, row.names = FALSE)

cat("\nSummary for the text: of", nrow(results), "(algorithm, n) cells,",
    sum(results$wilcox_p < 0.05, na.rm = TRUE), "show p < 0.05 (uncorrected).\n")

## Statistical significance across the full OFAT sweep, all five factors

algos <- sort(unique(ofat_results$algorithm))  # read from data

results <- data.frame()
for (fac in unique(ofat_results$factor)) {
  d_fac <- ofat_results[ofat_results$factor == fac, ]
  for (a in algos) {
    for (lv in sort(unique(d_fac$level))) {
      dd <- d_fac[d_fac$algorithm == a & d_fac$level == lv, ]
      eff  <- dd$hypervolume[dd$repair == "efficiency_only"]
      risk <- dd$hypervolume[dd$repair == "risk_aware"]
      wt <- tryCatch(wilcox.test(risk, eff, paired = TRUE)$p.value, error = function(e) NA)
      results <- rbind(results, data.frame(factor = fac, algorithm = a, level = lv,
                                            n_pairs = length(eff), wilcox_p = wt))
    }
  }
}

write.csv(results, "ofat_wilcoxon_all_factors.csv", row.names = FALSE)
{cat("=== per-cell results, by factor ===\n")
for (fac in unique(results$factor)) {
  cat(sprintf("\n--- %s ---\n", fac))
  print(results[results$factor == fac, c("algorithm","level","n_pairs","wilcox_p")], row.names = FALSE)
}}

{cat("\n=== summary, all factors combined ===\n")
cat(sprintf("of %d (factor, algorithm, level) cells, %d show p < 0.05 (uncorrected)\n",
            nrow(results), sum(results$wilcox_p < 0.05, na.rm = TRUE)))
}

{cat("\n=== summary, per factor ===\n")
for (fac in unique(results$factor)) {
  r <- results[results$factor == fac, ]
  cat(sprintf("%-18s %d/%d cells p<0.05\n", fac, sum(r$wilcox_p < 0.05, na.rm=TRUE), nrow(r)))
}}

comparisons <- list(
  c("risk_aware", "efficiency_only"),
  c("efficiency_random", "efficiency_only"),
  c("risk_aware_random", "risk_aware"),
  c("marginal_random", "repair_marginal")
)
 
## One row per (factor, level, algorithm, comparison, metric). Long
## format, deliberately -- easy to filter/pivot for any specific figure
## or table without re-running anything.
results <- data.frame()
for (fac in unique(ofat_results$factor)) {
  d_fac <- ofat_results[ofat_results$factor == fac, ]
  for (lv in sort(unique(d_fac$level))) {
    d_lv <- d_fac[d_fac$level == lv, ]
    for (algo in unique(d_lv$algorithm)) {
      d_algo <- d_lv[d_lv$algorithm == algo, ]
      for (cmp in comparisons) {
        x_repair <- cmp[1]; y_repair <- cmp[2]
        for (metric in c("hypervolume", "risk_at_best_roi")) {
          x <- d_algo[[metric]][d_algo$repair == x_repair]
          y <- d_algo[[metric]][d_algo$repair == y_repair]
          x <- x[!is.na(x)]; y <- y[!is.na(y)]  # per-column NA handling, not listwise
          if (length(x) == 0 || length(y) == 0) next
          delta <- cliffs_delta(x, y)
          results <- rbind(results, data.frame(
            factor = fac, level = lv, algorithm = algo,
            comparison = paste(x_repair, "vs", y_repair),
            metric = metric, n_x = length(x), n_y = length(y),
            mean_x = round(mean(x), 3), mean_y = round(mean(y), 3),
            delta = round(delta, 3)
          ))
        }
      }
    }
  }
}
 
write.csv(results, "ofat_full_analysis.csv", row.names = FALSE)
cat("Saved ofat_full_analysis.csv --", nrow(results), "rows.\n\n")

## Quick orientation: the risk_aware vs efficiency_only comparison on
## hypervolume, one row per factor x level x algorithm, the single most
## load-bearing comparison in the whole design.
{cat("=== risk_aware vs efficiency_only, hypervolume, all factors/levels/algorithms ===\n")
headline <- results[results$comparison == "risk_aware vs efficiency_only" & results$metric == "hypervolume", ]
print(headline[order(headline$factor, headline$level, headline$algorithm),
                c("factor","level","algorithm","delta")], row.names = FALSE)}

## risk_alpha needs separate handling, not the same treatment as the
## other four factors: efficiency_random and efficiency_only both
## hardcode alpha=1 internally, so eff_random_vs_eff_only cannot move
## across risk_alpha levels by construction, and at risk_alpha=1,
## risk_aware becomes mathematically identical to efficiency_only (the
## validity-check identity already established elsewhere), so any
## apparent "substitutability" near alpha=1 reflects risk_aware
## degrading toward efficiency_only, not a genuine condition-dependent
## finding comparable to the other four factors.

results <- data.frame()
for (fac in unique(ofat_results$factor)) {
  d_fac <- ofat_results[ofat_results$factor == fac, ]
  for (algo in unique(d_fac$algorithm)) {
    for (lv in sort(unique(d_fac$level))) {
      dd <- d_fac[d_fac$algorithm == algo & d_fac$level == lv, ]
      eff_only   <- dd$hypervolume[dd$repair == "efficiency_only"]
      eff_random <- dd$hypervolume[dd$repair == "efficiency_random"]
      risk_aware <- dd$hypervolume[dd$repair == "risk_aware"]
      if (!length(eff_only) || !length(eff_random) || !length(risk_aware)) next
      results <- rbind(results, data.frame(
        factor = fac, algorithm = algo, level = lv,
        eff_random_vs_risk_aware = round(cliffs_delta(eff_random, risk_aware), 3),
        eff_random_vs_eff_only   = round(cliffs_delta(eff_random, eff_only), 3)
      ))
    }
  }
}
results <- results[order(results$factor, results$algorithm, results$level), ]

## Flag risk_alpha rows explicitly rather than let them read as
## equivalent to the other four factors.
results$note <- ifelse(results$factor == "risk_alpha",
                        "mechanical: eff_random/eff_only never use risk_alpha; risk_aware -> eff_only at alpha=1",
                        "")

write.csv(results, "substitutability_by_condition.csv", row.names = FALSE)

{cat("=== Tabu, all factors except risk_alpha ===\n")
tabu_main <- results[results$algorithm == "Tabu" & results$factor != "risk_alpha",
                      c("factor", "level", "eff_random_vs_risk_aware", "eff_random_vs_eff_only")]
print(tabu_main, row.names = FALSE)

cat(sprintf("\nlevels where eff_random > risk_aware, excluding risk_alpha: %d / %d\n",
            sum(tabu_main$eff_random_vs_risk_aware > 0), nrow(tabu_main)))

cat("\n=== Tabu, risk_alpha (reported separately -- see note column) ===\n")
print(results[results$algorithm == "Tabu" & results$factor == "risk_alpha",
              c("level", "eff_random_vs_risk_aware", "eff_random_vs_eff_only")], row.names = FALSE)

cat("\nSaved substitutability_by_condition.csv -- all algorithms, all factors.\n")}

d75 <- ofat_results[ofat_results$factor == "n" & ofat_results$level == 75 & ofat_results$algorithm == "Tabu", ]
ra <- d75[d75$repair == "risk_aware", ];      ra <- ra[order(ra$run), ]
eo <- d75[d75$repair == "efficiency_only", ]; eo <- eo[order(eo$run), ]
diffs <- ra$hypervolume - eo$hypervolume
{cat(sprintf("n = 75, Tabu: %d of %d instances show no difference (|difference| < 1), %d favour risk-aware, %d favour efficiency-only\n",
            sum(abs(diffs) < 1), length(diffs), sum(diffs >= 1), sum(diffs <= -1)))}
