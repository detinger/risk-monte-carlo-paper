# Worked example on one instance (ODD §9.4–9.5)
# Runs run_full_comparison() on the §8 instance and plots the network.
# Requires: source("04_additional_algorithms.R")


## 9.4 Running the comparison on the worked-example instance

t0 <- proc.time()[3]
full_res <- run_full_comparison(acts, dag$preds, dag$succ, topo, budget, must_idx, deadline, reps = 3)
cat("wall time (s):", round(proc.time()[3] - t0, 1), " | rows:", nrow(full_res), "\n")

{cat("--- multi-objective placement: mean hypervolume, algorithm x repair ---\n")
mres <- full_res[full_res$placement == "multi", ]
tab_m <- aggregate(hypervolume ~ algorithm + repair, data = mres, FUN = mean)
print(reshape(tab_m, idvar = "algorithm", timevar = "repair", direction = "wide"))}

{cat("--- single-objective placement: mean best value and risk-at-best-value ---\n")
sres <- full_res[full_res$placement == "single", ]
tab_s <- aggregate(cbind(best_roi, risk_at_best_roi) ~ algorithm + repair, data = sres, FUN = mean)
print(tab_s[order(tab_s$algorithm, tab_s$repair), ])}

{cat("--- feasible fraction under 'penalty' specifically (0 = never found a feasible point) ---\n")
pen <- full_res[full_res$repair == "penalty", ]
print(aggregate(feasible_fraction ~ algorithm + placement, data = pen, FUN = mean))}

{cat("--- Greedy, standalone (not crossed with repair) ---\n")
gf <- lapply(c(0, 0.4, 0.8), function(wgt) {
  x <- run_greedy(acts, dag$preds, dag$succ, topo, budget, must_idx, deadline, weight_risk = wgt)
  mc <- project_duration_mc(x, acts, dag$preds, topo, nsim = 300, deadline = deadline)
  data.frame(weight_risk = wgt, roi = cost_roi(x, acts)["roi"], risk = mc$p_over)
})
print(do.call(rbind, gf))}


## 9.5. Visualization

plot_activity_network <- function(inst, title = NULL) {
  n <- nrow(inst$acts)
  preds <- inst$preds
 
  ## Longest-path level of each node: 0 for sources, else 1 + max(level of
  ## predecessors). Computed in topological order so every predecessor's
  ## level is already known by the time a node is processed.
  level <- integer(n)
  for (i in inst$topo) {
    if (length(preds[[i]]) == 0) {
      level[i] <- 0
    } else {
      level[i] <- 1 + max(level[preds[[i]]])
    }
  }
 
  ## x-position: spread nodes within each level evenly, sorted by node
  ## index within the level purely for a stable, reproducible layout.
  x <- numeric(n)
  for (lv in unique(level)) {
    idx <- which(level == lv)
    x[idx] <- seq(-1, 1, length.out = length(idx) + 2)[-c(1, length(idx) + 2)]
    if (length(idx) == 1) x[idx] <- 0
  }
  y <- -level  # top-to-bottom: sources at the top
 
  plot(NULL, xlim = range(x) + c(-0.3, 0.3), ylim = range(y) + c(-0.5, 0.5),
       xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
       main = title, cex.main = 1.0, font.main = 1)
 
  ## Edges first, so nodes draw on top and stay legible.
  for (i in seq_len(n)) {
    for (j in inst$succ[[i]]) {
      arrows(x[i], y[i], x[j], y[j], length = 0.05, angle = 20,
             col = "grey65", lwd = 0.7)
    }
  }
 
  ## Nodes: mandatory activities distinguished by fill colour, 
  ## so the must-closure story from Section 2.1 is visible directly
  ## in the picture, not only in the accompanying text.
  is_must <- seq_len(n) %in% inst$must_idx
  points(x, y, pch = 21, cex = 1.6,
         bg = ifelse(is_must, "firebrick", "steelblue2"),
         col = "grey30", lwd = 0.8)
 
  legend("topright", legend = c("mandatory", "optional"), pch = 21,
         pt.bg = c("firebrick", "steelblue2"), col = "grey30",
         bty = "n", cex = 0.8)
 
  invisible(list(x = x, y = y, level = level))
}


inst <- make_instance(50, seed = 11, density = 0.10, must_k = 3, budget_fraction = 0.5)

#png("network_examples.png", width = 1500, height = 550)
par(mfrow = c(1, 3), mar = c(1, 1, 3, 1))
for (d in c(0.02, 0.10, 0.25)) {
  inst <- make_instance(50, seed = 11, density = d, must_k = 3, budget_fraction = 0.5)
  n_edges <- sum(sapply(inst$preds, length))
  plot_activity_network(inst, title = sprintf("n=50, density=%.2f (%d edges)", d, n_edges))
}
#dev.off()
