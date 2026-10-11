# Bimodal response analysis: what predicts which instances show delta=1 vs delta=0?
# Run from the results/ directory (or adjust paths).

## ---- 1. Instance metadata from checkpoint filenames ----
fnames <- list.files("checkpoints_revision", pattern = "\\.rds$")
parsed <- regmatches(fnames, regexpr("rev_inst_(\\d+)_s(\\d+)_d([0-9.]+)", fnames))
meta <- do.call(rbind, lapply(parsed, function(x) {
  m <- regmatches(x, regexec("rev_inst_(\\d+)_s(\\d+)_d([0-9.]+)", x))[[1]]
  data.frame(instance = as.integer(m[2]), seed = as.integer(m[3]),
             density = as.numeric(m[4]), stringsAsFactors = FALSE)
}))
meta <- unique(meta[order(meta$instance), ])

## ---- 2. Per-instance delta for RISK vs EFF ----
per_inst <- read.csv("revision_per_instance.csv", stringsAsFactors = FALSE)
risk_delta <- per_inst[per_inst$arm == "RISK", c("instance", "delta")]
names(risk_delta)[2] <- "risk_delta"

## ---- 3. Per-instance zero_range_share from checkpoint files ----
## Each checkpoint (one instance × one replication × one arm) stores a data frame
## with columns including zero_range_share if the audit code ran.
## Fall back to computing from raw checkpoints if available.
ckpt_files <- list.files("checkpoints_revision", pattern = "\\.rds$", full.names = TRUE)

## Extract per-instance mean zero_range_share for RISK arm only.
## Checkpoint filename encodes arm via the numeric suffix (E value).
## RISK arm: seed_as = "risk_aware" → strategy index 3 → suffix 09977 (inst*1e6 + rep*1e4 + 1*1e2 + 3*10 + 2)
## Actually, easiest: load a sample and check column names.
sample_ckpt <- readRDS(ckpt_files[1])
cat("Checkpoint columns:", names(sample_ckpt), "\n")

## If zero_range_share is available per row, aggregate.
if ("zero_range_share" %in% names(sample_ckpt)) {
  zr_data <- do.call(rbind, lapply(ckpt_files, function(f) {
    d <- readRDS(f)
    d[d$repair == "risk_aware", c("instance", "rep", "zero_range_share")]
  }))
  zr_by_inst <- aggregate(zero_range_share ~ instance, zr_data, mean)
} else {
  cat("zero_range_share not in checkpoint; skipping per-instance ZR computation.\n")
  zr_by_inst <- data.frame(instance = meta$instance,
                           zero_range_share = NA_real_)
}

## ---- 4. Merge and analyse ----
d <- merge(meta, risk_delta, by = "instance")
d <- merge(d, zr_by_inst, by = "instance", all.x = TRUE)
d$response <- ifelse(d$risk_delta >= 0.474, "large", "negligible")

cat("\n=== Per-instance data ===\n")
print(d[order(d$instance), ], row.names = FALSE)

cat("\n=== Spearman correlation: density vs RISK delta ===\n")
ct <- cor.test(d$density, d$risk_delta, method = "spearman", exact = FALSE)
cat(sprintf("rho = %.3f, p = %.4f\n", ct$estimate, ct$p.value))

if (!all(is.na(d$zero_range_share))) {
  cat("\n=== Spearman correlation: zero_range_share vs RISK delta ===\n")
  ct2 <- cor.test(d$zero_range_share, d$risk_delta, method = "spearman", exact = FALSE)
  cat(sprintf("rho = %.3f, p = %.4f\n", ct2$estimate, ct2$p.value))
}

cat("\n=== Mean density by response group ===\n")
print(aggregate(density ~ response, d, mean), row.names = FALSE)

cat("\n=== Wilcoxon test: density differs between large and negligible response? ===\n")
wt <- wilcox.test(density ~ response, data = d, exact = FALSE)
cat(sprintf("W = %.0f, p = %.4f\n", wt$statistic, wt$p.value))

## Save for reference
write.csv(d, "bimodal_predictor.csv", row.names = FALSE)
cat("\nSaved to bimodal_predictor.csv\n")
