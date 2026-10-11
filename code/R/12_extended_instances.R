## 12_extended_instances.R
##
## Extends the bimodal-response analysis from 20 to 70 instances by running
## 50 additional Design A configurations (n=50, same parameters) with tabu
## search only, comparing RISK vs EFF.
##
## Purpose: strengthen the density-predictor claim (Spearman rs = -0.705 on
## n=20) with a 70-instance combined sample.
##
## Run via:   cd code && Rscript run_extended.R ../results

## ---- Constants ----

EXT_N_NEW   <- 50L
EXT_SEED    <- 2027L    ## separate from Design A's fallback seed (2026)
EXT_SEEDS   <- 200L + seq_len(EXT_N_NEW)   ## 201-250

## Arms: RISK and EFF only, multi-objective placement
EXT_ARMS <- rbind(
  revision_arm("EFF",  "efficiency_only", "multi"),
  revision_arm("RISK", "risk_aware",      "multi")
)


## ---- Functions ----

make_ext_specs <- function(n = EXT_N_NEW, rng_seed = EXT_SEED) {
  set.seed(rng_seed)
  dens <- runif(n, 0.05, 0.25)
  lapply(seq_len(n), function(i) list(seed = EXT_SEEDS[i], density = dens[i]))
}

## Compute Cliff's delta (standalone, no codebase dependency for reporting)
cliff_rank_ext <- function(x, y) {
  sum(outer(x, y, ">")) / (length(x) * length(y)) -
  sum(outer(x, y, "<")) / (length(x) * length(y))
}

summarise_extended <- function(ext_main, ext_specs) {
  n <- length(ext_specs)
  per_inst <- do.call(rbind, lapply(seq_len(n), function(i) {
    di <- ext_main[ext_main$instance == i, ]
    x <- di$hypervolume[di$arm == "RISK"]
    y <- di$hypervolume[di$arm == "EFF"]
    data.frame(
      instance_ext = i,
      seed         = EXT_SEEDS[i],
      density      = ext_specs[[i]]$density,
      risk_delta   = round(cliff_rank_ext(x, y), 6)
    )
  }))

  ## Combine with the original 20-instance bimodal data
  bp_file <- "bimodal_predictor.csv"
  ri_file <- "revision_per_instance.csv"
  if (file.exists(bp_file) && file.exists(ri_file)) {
    orig_risk <- read.csv(ri_file, stringsAsFactors = FALSE)
    orig_risk <- orig_risk[orig_risk$arm == "RISK", c("instance", "delta")]
    names(orig_risk) <- c("instance_ext", "risk_delta")
    bimodal <- read.csv(bp_file, stringsAsFactors = FALSE)
    orig_comb <- merge(orig_risk, bimodal[, c("instance", "density")],
                       by.x = "instance_ext", by.y = "instance", all.x = TRUE)
    orig_comb$source <- "original"
    new_comb <- per_inst[, c("instance_ext", "density", "risk_delta")]
    new_comb$source <- "extended"
    all_inst <- rbind(orig_comb[, c("density", "risk_delta", "source")],
                      new_comb[, c("density", "risk_delta", "source")])
    ct <- cor.test(all_inst$density, all_inst$risk_delta, method = "spearman", exact = FALSE)
    cat(sprintf("\n=== Combined density vs RISK delta (%d instances) ===\n", nrow(all_inst)))
    cat(sprintf("Spearman rho = %.3f, p = %.6f\n", ct$estimate, ct$p.value))
    write.csv(all_inst, "extended_density_correlation.csv", row.names = FALSE)
    cat("Saved to extended_density_correlation.csv\n")
  } else {
    cat("Note: bimodal_predictor.csv or revision_per_instance.csv not found; combined analysis skipped.\n")
  }

  write.csv(per_inst, "extended_per_instance.csv", row.names = FALSE)
  cat("\nPer-instance results (sorted by density):\n")
  print(per_inst[order(per_inst$density), ], row.names = FALSE)
  invisible(per_inst)
}


## ---- Top-level run (executed by run_toplevel, not by load_definitions) ----

ext_specs <- make_ext_specs()

cat(sprintf("Running %d new instances x 30 reps x 2 arms (RISK, EFF) on %d cores\n",
            EXT_N_NEW, 16L))

ext_main <- run_revision_experiments(ext_specs, EXT_ARMS, reps = 30L, n_cores = 16L,
                                     checkpoint_dir = "checkpoints_extended")
saveRDS(ext_main, "extended_main_results.rds")

summarise_extended(ext_main, ext_specs)
