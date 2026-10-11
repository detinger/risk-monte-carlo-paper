# Real-world data comparison (ODD §12)
# PSPLIB (j30/j60), MineLib KD, ReleasePlanner datasets.
# Requires: external data files (see comments in code for download URLs), in
# DATA_DIR (run_all.R sets it to ../data; default: the working directory).
# Also requires full_production (R/06) for the Design A comparisons.

DATA_DIR <- if (exists("DATA_DIR")) DATA_DIR else "."


## 12.2.1 PSPLib

## read from already-extracted folders directly 
read_instance_set <- function(dir) {
  files <- list.files(dir, pattern = "\\.sm$", full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("No .sm files found in ", dir,
                            " -- check the folder actually contains PSPLIB .sm files.")
  files
}

j30_files <- read_instance_set(file.path(DATA_DIR, "unzipped_j30"))
j60_files <- read_instance_set(file.path(DATA_DIR, "unzipped_j60"))
cat("j30 files found:", length(j30_files), "\n")

cat("j60 files found:", length(j60_files), "\n")

## parsing
parse_sm <- function(path) {
  lines <- readLines(path, warn = FALSE)
  n_jobs <- as.integer(sub(".*:\\s*", "", grep("^jobs", lines, value = TRUE)))

  prec_start <- grep("^PRECEDENCE RELATIONS:", lines) + 2
  prec_end <- grep("^\\*{10,}", lines)
  prec_end <- prec_end[prec_end > prec_start][1] - 1
  prec_lines <- lines[prec_start:prec_end]

  edges_from <- integer(0); edges_to <- integer(0)
  for (ln in prec_lines) {
    parts <- as.integer(strsplit(trimws(ln), "\\s+")[[1]])
    if (length(parts) < 3) next
    job <- parts[1]; n_succ <- parts[3]
    if (n_succ > 0) {
      succs <- parts[4:(3 + n_succ)]
      edges_from <- c(edges_from, rep(job, n_succ))
      edges_to <- c(edges_to, succs)
    }
  }

  dur_start <- grep("^REQUESTS/DURATIONS:", lines)
  durations <- NA_real_
  if (length(dur_start)) {
    dur_start <- dur_start + 3
    dur_end <- grep("^\\*{10,}", lines); dur_end <- dur_end[dur_end > dur_start][1] - 1
    dur_lines <- lines[dur_start:dur_end]
    durations <- vapply(dur_lines, function(ln) {
      parts <- as.numeric(strsplit(trimws(ln), "\\s+")[[1]])
      if (length(parts) >= 3) parts[3] else NA_real_
    }, numeric(1), USE.NAMES = FALSE)
    durations <- durations[!is.na(durations)]
  }

  n_possible_edges <- n_jobs * (n_jobs - 1) / 2
  list(file = basename(path), n_jobs = n_jobs, n_edges = length(edges_from),
       density = length(edges_from) / n_possible_edges,
       mean_duration = mean(durations, na.rm = TRUE),
       max_duration = suppressWarnings(max(durations, na.rm = TRUE)))
}

summarize_instance_set <- function(files, n_sample = NULL, seed = 1) {
  if (!is.null(n_sample) && n_sample < length(files)) {
    set.seed(seed)
    files <- sample(files, n_sample)
  }
  parsed <- lapply(files, function(f) tryCatch(parse_sm(f), error = function(e) NULL))
  failed <- sum(vapply(parsed, is.null, logical(1)))
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  df <- do.call(rbind, lapply(parsed, function(x) data.frame(
    file = x$file, n_jobs = x$n_jobs, density = x$density,
    mean_duration = x$mean_duration, max_duration = x$max_duration)))
  list(n_requested = length(files), n_parsed = nrow(df), n_failed = failed, per_instance = df)
}
cat("Retrieval and parsing pipeline defined.\n")

j30_result <- summarize_instance_set(j30_files)
cat("j30: requested", j30_result$n_requested, " parsed", j30_result$n_parsed,
    " failed", j30_result$n_failed, "\n")

j60_result <- summarize_instance_set(j60_files)
cat("j60: requested", j60_result$n_requested, " parsed", j60_result$n_parsed,
    " failed", j60_result$n_failed, "\n")

cat("\nj30 n_jobs distribution (expect all 32):\n"); print(table(j30_result$per_instance$n_jobs))

cat("j60 n_jobs distribution (expect all 62):\n"); print(table(j60_result$per_instance$n_jobs))

saveRDS(j30_result, "psplib_j30_summary.rds")
saveRDS(j60_result, "psplib_j60_summary.rds")
cat("\nSaved psplib_j30_summary.rds and psplib_j60_summary.rds.\n")

## Compares PSPLIB's real density distribution against Design A's actual
## generated instances. The `density` column of full_production is the arc
## probability passed to generate_dag(), which also adds an arc i -> i + 1
## wherever a node drew no successor, so the realized density of a network
## is higher. Both sides are compared on realized density (arcs / (n(n-1)/2)).
## PSPLIB counts include the two dummy jobs and their arcs.

designA_specs <- unique(full_production[, c("instance", "seed", "density")])
designA_specs <- designA_specs[order(designA_specs$instance), ]
designA_preds <- lapply(seq_len(nrow(designA_specs)), function(k)
  generate_dag(50, designA_specs$density[k], designA_specs$seed[k])$preds)
designA_density <- vapply(designA_preds, function(p) sum(lengths(p)) / (50 * 49 / 2), numeric(1))

j30_density <- readRDS("psplib_j30_summary.rds")$per_instance$density
j60_density <- readRDS("psplib_j60_summary.rds")$per_instance$density

cat("=== Design A (n=50, 20 real generated instances) ===\n")

cat(sprintf("n=%d  mean=%.4f  sd=%.4f  range=%.4f-%.4f\n",
            length(designA_density), mean(designA_density), sd(designA_density),
            min(designA_density), max(designA_density)))

cat("\n=== PSPLIB j30 (n=32, 480 real instances) ===\n")

cat(sprintf("n=%d  mean=%.4f  sd=%.4f  range=%.4f-%.4f\n",
            length(j30_density), mean(j30_density), sd(j30_density),
            min(j30_density), max(j30_density)))

cat("\n=== PSPLIB j60 (n=62, 480 real instances) ===\n")

cat(sprintf("n=%d  mean=%.4f  sd=%.4f  range=%.4f-%.4f\n",
            length(j60_density), mean(j60_density), sd(j60_density),
            min(j60_density), max(j60_density)))

cat("\n=== Design A vs. each PSPLIB set, Wilcoxon two-sample ===\n")

wt_j30 <- wilcox.test(designA_density, j30_density)
wt_j60 <- wilcox.test(designA_density, j60_density)
cat(sprintf("vs j30 (n=32): p = %.4f\n", wt_j30$p.value))

cat(sprintf("vs j60 (n=62): p = %.4f\n", wt_j60$p.value))

interp_density_n50 <- mean(j30_density) + (mean(j60_density) - mean(j30_density)) *
  (50 - 32) / (62 - 32)
cat(sprintf("\nLinear interpolation, PSPLIB-consistent density at n=50: %.4f\n", interp_density_n50))

cat(sprintf("Design A's actual mean: %.4f (%s the interpolated value, and %s even the n=32 anchor)\n",
            mean(designA_density),
            ifelse(mean(designA_density) > interp_density_n50, "above", "below"),
            ifelse(mean(designA_density) > mean(j30_density), "above", "below")))

designA_psplib_comparison <- list(
  designA = c(n = length(designA_density), mean = mean(designA_density), sd = sd(designA_density)),
  j30 = c(n = length(j30_density), mean = mean(j30_density), sd = sd(j30_density)),
  j60 = c(n = length(j60_density), mean = mean(j60_density), sd = sd(j60_density)),
  wilcox_p_vs_j30 = wt_j30$p.value,
  wilcox_p_vs_j60 = wt_j60$p.value,
  interpolated_n50 = interp_density_n50
)
saveRDS(designA_psplib_comparison, "designA_psplib_density_comparison.rds")
cat("\nSaved designA_psplib_density_comparison.rds\n")


## 12.2.2 MineLib’s Arizona’s Copper Deposit (KD)

## Structural comparison: MineLib KD (real Arizona copper deposit, CPIT
## restricted to period 0) vs. this paper's synthetic generator.
blocks <- read.table(file.path(DATA_DIR, "kd.blocks", "kd.blocks"), header = FALSE,
                      col.names = c("id","x","y","z","tonn","value","dest","cu_pct","proc_profit"))
prec_lines <- readLines(file.path(DATA_DIR, "kd.prec", "kd.prec"))

n_blocks <- nrow(blocks)

## Precedence: total edges, mean predecessors per block
n_edges <- 0L
for (ln in prec_lines) {
  parts <- as.integer(strsplit(trimws(ln), "\\s+")[[1]])
  n_pred <- parts[2]
  n_edges <- n_edges + n_pred
}

kd_mean_preds <- n_edges / n_blocks

{cat("=== MineLib KD (n=14,153, real block/precedence data) ===\n")
cat("total precedence edges:", n_edges, "\n")
cat("mean predecessors per block:", round(kd_mean_preds, 3), "\n")
 
cat("\n=== Value structure ===\n")
cat("fraction negative value (waste):", round(mean(blocks$value < 0), 3), "\n")
cat("fraction positive value (ore):", round(mean(blocks$value > 0), 3), "\n")}

## Design A comparison: mean predecessors of Design A's own networks,
## regenerated from their seeds and arc probabilities (section 12.2.1).
designA_mean_preds <- vapply(designA_preds, function(p) mean(lengths(p)), numeric(1))
 
{cat("\n=== Design A (n=50, 20 real generated instances) ===\n")
cat(sprintf("mean predecessors per activity: mean=%.3f  sd=%.3f  range=%.3f-%.3f\n",
            mean(designA_mean_preds), sd(designA_mean_preds),
            min(designA_mean_preds), max(designA_mean_preds)))
 
cat(sprintf("\nKD's mean predecessors is %.1fx Design A's mean (%.3f vs %.3f)\n",
            kd_mean_preds / mean(designA_mean_preds), kd_mean_preds, mean(designA_mean_preds)))
cat("KD's precedence is geometric (spatial adjacency); Design A's is\n")
cat("arbitrary-pair, so this is a portable statistic across mechanisms,\n")
cat("not a claim that the two generating processes are otherwise alike.\n")}


## 12.2.3 ReleasePlanner Dataset 2

# The dataset files have been downloaded from https://sites.google.com/site/mrkarim/data-sets?authuser=0

library(readxl)
RP_DIR   <- file.path(DATA_DIR, "ReleasePlanner Dataset 2")
RP_FILES <- c(dataset_1 = "dataset-1.xls", dataset_2 = "dataset-2.xls")

## One sheet as text, with column names taken from the given header row.
rp_sheet <- function(file, sheet, header_row = 2) {
  raw <- as.data.frame(suppressMessages(read_excel(file, sheet = sheet, col_names = FALSE, col_types = "text")))
  hdr <- trimws(as.character(unlist(raw[header_row, ])))
  hdr[is.na(hdr) | !nzchar(hdr)] <- paste0("unnamed_", which(is.na(hdr) | !nzchar(hdr)))
  body <- raw[-seq_len(header_row), , drop = FALSE]
  names(body) <- hdr
  rownames(body) <- NULL
  body
}
rp_empty <- function(x) is.na(x) | !nzchar(trimws(x)) | trimws(x) == "FALSE"
rp_split <- function(x) {
  if (rp_empty(x)) return(character(0))
  out <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  out[nzchar(out)]
}

## Precedence network of one dataset, as a predecessor list over its features.
rp_network <- function(file) {
  req <- rp_sheet(file, "Requirements")
  req <- req[!rp_empty(req$ID), ]
  ids <- trimws(req$ID)
  arcs <- unique(do.call(rbind, lapply(seq_along(ids), function(k) {
    to <- rp_split(req$Precedes[k])
    if (length(to)) data.frame(from = ids[k], to = to) else NULL
  })))
  coupled <- unique(do.call(rbind, lapply(seq_along(ids), function(k) {
    to <- rp_split(req$Couplings[k])
    if (length(to)) data.frame(a = pmin(ids[k], to), b = pmax(ids[k], to)) else NULL
  })))
  unknown <- setdiff(c(arcs$from, arcs$to), ids)
  if (length(unknown)) stop("Unknown feature IDs in 'Precedes': ", paste(unknown, collapse = ", "))
  preds <- lapply(ids, function(i) match(arcs$from[arcs$to == i], ids))
  list(ids = ids, preds = preds, n_coupled_pairs = if (is.null(coupled)) 0L else nrow(coupled), req = req)
}

## Structural statistics of a network given as a predecessor list. Works for
## any node order: a topological order is computed first, and a cycle stops.
network_stats <- function(preds) {
  n <- length(preds)
  indeg <- lengths(preds); succ <- vector("list", n)
  for (j in seq_len(n)) for (i in preds[[j]]) succ[[i]] <- c(succ[[i]], j)
  order <- integer(0); queue <- which(indeg == 0)
  while (length(queue)) {
    u <- queue[1]; queue <- queue[-1]; order <- c(order, u)
    for (v in succ[[u]]) { indeg[v] <- indeg[v] - 1; if (indeg[v] == 0) queue <- c(queue, v) }
  }
  if (length(order) < n) stop("The precedence network contains a cycle.")
  level <- integer(n); anc <- vector("list", n)
  for (u in order) {
    p <- preds[[u]]
    level[u] <- if (length(p)) 1 + max(level[p]) else 0
    anc[[u]] <- unique(c(p, unlist(anc[p])))
  }
  involved <- lengths(preds) > 0 | lengths(succ) > 0
  n_arcs <- sum(lengths(preds))
  c(n = n, arcs = n_arcs, density = n_arcs / (n * (n - 1) / 2), mean_predecessors = n_arcs / n,
    share_independent = mean(!involved), longest_chain = max(level),
    mean_ancestors = mean(lengths(anc)), max_ancestors = max(lengths(anc)))
}

## Evidence that the datasets are candidate sets before selection: total
## effort demanded by all features against each resource's capacity.
rp_capacity <- function(file) {
  res <- rp_sheet(file, "Resources")
  req <- rp_sheet(file, "Requirements"); req <- req[!rp_empty(req$ID), ]
  releases <- names(res)[4:ncol(res)]
  do.call(rbind, lapply(seq_len(nrow(res)), function(k) {
    r <- trimws(res[k, 1])
    if (rp_empty(r) || !r %in% names(req)) return(NULL)
    cap <- as.numeric(unlist(res[k, releases]))
    demand <- sum(as.numeric(req[[r]]), na.rm = TRUE)
    data.frame(resource = r, demand = demand, capacity_one_release = cap[1],
               capacity_all_releases = sum(cap, na.rm = TRUE),
               demand_over_one_release = round(demand / cap[1], 2),
               demand_over_all_releases = round(demand / sum(cap, na.rm = TRUE), 2))
  }))
}

## Cost against value per feature. Value is the stakeholder-weighted mean
## score on the "Value" criterion; cost is "Cost of development" where
## present, otherwise total effort across resources.
rp_cost_value <- function(file) {
  req <- rp_sheet(file, "Requirements"); req <- req[!rp_empty(req$ID), ]
  effort_cols <- names(req)[(which(names(req) == "Group") + 1):ncol(req)]
  cost_cols <- if ("Cost of development" %in% names(req)) "Cost of development" else effort_cols
  cost <- rowSums(sapply(cost_cols, function(cc) as.numeric(req[[cc]])), na.rm = TRUE)
  sh <- rp_sheet(file, "Stakeholders"); sh <- sh[!rp_empty(sh[[1]]), ]
  w <- setNames(as.numeric(sh$Weight), trimws(sh[[1]]))
  val_sheet <- rp_sheet(file, "Value")
  value <- sapply(trimws(req$ID), function(fid) {
    row <- val_sheet[trimws(val_sheet[[1]]) == fid, , drop = FALSE]
    if (!nrow(row)) return(NA)
    scores <- unlist(row[1, -1]); keep <- !rp_empty(scores) & names(scores) %in% names(w)
    if (!any(keep)) return(NA)
    weighted.mean(as.numeric(scores[keep]), w[names(scores)[keep]])
  })
  ok <- is.finite(cost) & cost > 0 & is.finite(value)
  c(features = sum(ok), spearman = cor(cost[ok], value[ok], method = "spearman"),
    cost_max_over_min = max(cost[ok]) / min(cost[ok]))
}

## ---- Run ------------------------------------------------------------------
rp <- lapply(RP_FILES, function(f) rp_network(file.path(RP_DIR, f)))
rp_stats <- t(sapply(rp, function(x) network_stats(x$preds)))

specs <- unique(full_production[, c("instance", "seed", "density")])
specs <- specs[order(specs$instance), ]
designA_stats <- t(sapply(seq_len(nrow(specs)), function(k)
  network_stats(generate_dag(50, specs$density[k], specs$seed[k])$preds)))

## Design B at n = 50, lowest density levels and baseline, for a comparison
## at matched density (same seeds as the OFAT sweep, Section 11.2).
designB_stats <- t(sapply(c(0.02, 0.05, 0.10), function(rho)
  colMeans(t(sapply(5000 + 1:20, function(sd) network_stats(generate_dag(50, rho, sd)$preds))))))
rownames(designB_stats) <- c("Design B, density 0.02 (mean)", "Design B, density 0.05 (mean)",
                             "Design B, density 0.10 (mean)")

rp_comparison <- rbind(
  "ReleasePlanner Dataset 1"  = rp_stats["dataset_1", ],
  "ReleasePlanner Dataset 2"  = rp_stats["dataset_2", ],
  "Design A (mean)"           = colMeans(designA_stats),
  "Design A (min)"            = apply(designA_stats, 2, min),
  "Design A (max)"            = apply(designA_stats, 2, max),
  designB_stats)
cat("=== Precedence structure: ReleasePlanner datasets and the networks used in this paper ===\n")

print(round(rp_comparison, 3))

## Dataset 1 has n = 50, the size of Design A's instances: where does it fall
## relative to the range of Design A's 20 networks?
d1 <- rp_stats["dataset_1", ]
cols <- c("density", "mean_predecessors", "share_independent", "longest_chain", "mean_ancestors", "max_ancestors")
position <- sapply(cols, function(cc) {
  if (d1[cc] < min(designA_stats[, cc])) "below Design A's range"
  else if (d1[cc] > max(designA_stats[, cc])) "above Design A's range" else "within Design A's range"
})
cat("\n=== Dataset 1 (n = 50) relative to Design A's 20 networks ===\n")

print(data.frame(statistic = cols, dataset_1 = round(d1[cols], 3), design_A_min = round(apply(designA_stats[, cols], 2, min), 3),
                 design_A_max = round(apply(designA_stats[, cols], 2, max), 3), position = position), row.names = FALSE)

cat("\n=== Coupled feature pairs (not represented in this paper's model) ===\n")

print(sapply(rp, `[[`, "n_coupled_pairs"))

cat("\n=== Candidate sets before selection: demand of all features against capacity ===\n")

rp_capacity_tab <- lapply(RP_FILES, function(f) rp_capacity(file.path(RP_DIR, f)))
for (k in names(rp_capacity_tab)) { cat("\n", k, "\n"); print(rp_capacity_tab[[k]], row.names = FALSE) }

cat("\n=== Cost against stakeholder value per feature ===\n")

rp_cost_value_tab <- t(sapply(RP_FILES, function(f) rp_cost_value(file.path(RP_DIR, f))))
print(round(rp_cost_value_tab, 3))

saveRDS(list(comparison = rp_comparison, designA = designA_stats, capacity = rp_capacity_tab,
             cost_value = rp_cost_value_tab), "releaseplanner_comparison.rds")
