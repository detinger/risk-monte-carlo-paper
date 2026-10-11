# When Does Risk-Aware Repair Pay Off?

Replication code for the paper submitted to *MDPI Applied Sciences*.

The central question: does a repair operator that accounts for Monte Carlo deadline-overrun risk outperform a pure efficiency (ROI/cost) rule when trimming over-budget project activity networks?

---

## Repository layout

```
code/
  R/                       # numbered scripts, sourced in order
    01_core_functions.R    # sampling, DAG, critical path, repair operators, NSGA-II engine
    02_illustration.R      # §8 repair-divergence demo and near-exact validation
    03_additional_strategies.R  # MARG and MARG-SM repair rules
    04_additional_algorithms.R  # Random, Greedy, SA, Tabu, MOEA/D
    05_worked_example.R    # §9.4–9.5 multi-algorithm comparison on one instance
    06_full_scale_run.R    # §10 Design A: parallel, checkpointed production run
    07_full_scale_analysis.R    # §10.7–10.8 analysis and validation
    08_ofat_sweep.R        # §11.1–11.3 Design B OFAT sensitivity sweep
    09_targeted_checks.R   # §11.4–11.7 τ sweep, config check, scalarization check
    10_real_world_data.R   # §12 PSPLIB, MineLib KD, ReleasePlanner comparison
    11_revision_experiments.R   # §13 placebo, CRIT proxy, equal-budget, probe sweep
    12_extended_instances.R     # 50-instance extension (70 instances total)
    bimodal_analysis.R     # density-predictor analysis
  load_functions.R         # load definitions without running experiments
  run_all.R                # full pipeline (days of compute)
  run_revision.R           # revision experiments only (~7–8 h on 6 cores)
  run_extended.R           # extended-instances experiment
  run_extended_analysis.R  # post-hoc summary of extended run
  smoke_test.R             # end-to-end sanity check (~6 minutes)
paper/
  paper-MDPI-AppliedSciences-merged.docx
results/                   # output CSVs, RDS files, checkpoint folders
data/                      # external real-world data (not in the repo; see §10)
media/
  image1.png
```

---

## Requirements

- **R ≥ 4.1** — only base R and the `parallel` package (built-in) are required for all experiments except the real-world comparison.
- `readxl` — needed only for `R/10_real_world_data.R` (ReleasePlanner `.xlsx` file). The script skips that file gracefully if the package is absent.
- `rstudioapi` — auto-detected when sourcing inside RStudio; not needed otherwise.

No CRAN packages need to be installed before running the smoke test or the main experiments.

---

## Quick start — smoke test

Verifies the entire pipeline end-to-end at tiny settings in about 6 minutes. Run from the `code/` directory:

```sh
cd code
Rscript smoke_test.R            # output to a temp directory
Rscript smoke_test.R /my/outdir # or choose a directory
```

The smoke test runs every stage (01–12) at reduced sizes, checks §8 and §9.4 outputs against the values in the ODD appendix, and prints a per-stage timing table.

---

## Running the experiments

All entry points are in `code/`. They write output to `../results` by default; pass a path as the first argument to override.

### Full pipeline (days)

```sh
cd code
Rscript run_all.R
```

Runs R/01 → R/11 in dependency order. Long-running stages (Design A, OFAT) are checkpointed in `results/checkpoints_production/` and resume automatically if the command is rerun.

### Revision experiments only (~7–8 h on 6 cores)

```sh
cd code
Rscript run_revision.R [results_dir]
```

Requires `full_production_results.rds` to already exist in `results_dir` (produced by Design A, R/06). If it is absent, the scale-balanced summary is produced without the raw-vs-normalised comparison.

### Extended instances (50 new + analysis)

```sh
cd code
Rscript run_extended.R [results_dir]        # produces extended_main_results.rds
Rscript run_extended_analysis.R [results_dir]
```

### Real-world data (§12)

Download the following datasets and place them in `data/` before running `run_all.R` (or source R/10 directly):

| Dataset | Expected path in `data/` |
|---------|--------------------------|
| PSPLIB j30 (unzipped) | `data/unzipped_j30/` (`.sm` files) |
| PSPLIB j60 (unzipped) | `data/unzipped_j60/` (`.sm` files) |
| MineLib KD | see comments in `R/10_real_world_data.R` |
| ReleasePlanner | see comments in `R/10_real_world_data.R` |

R/10 will error with a descriptive message if a dataset is missing.

---

## Key concepts

**Activity network.** A project is a DAG of activities, each with a cost, ROI (value), and a triangular duration distribution (optimistic / realistic / pessimistic). A selected subset must satisfy precedence closure and fit a budget.

**Repair operator.** When a candidate solution exceeds the budget, the operator greedily removes leaf activities until the budget is met. Three cutting criteria are compared:

| Label | Rule |
|-------|------|
| `EFF` | Lowest efficiency (ROI/cost) first — the baseline |
| `RISK` | Blend of efficiency and Monte Carlo risk-reduction signal (deadline-overrun probability drop per leaf removed) |
| `MARG` | Lowest marginal ROI first |
| `RAND-RISK` | Placebo: same risk signal, randomly shuffled across leaves |
| `CRIT` | Deterministic proxy: PERT critical-path finish time instead of Monte Carlo |

Stochastic (`_random`) variants of each rule draw the cut from a softmax over the score instead of argmax, reducing diversity collapse.

**Objective placement.** Each algorithm is run under two objective functions: single (maximise ROI only) and multi (maximise ROI, minimise P(overrun)), yielding 2 × N condition cells.

**Algorithms.** Random search, GA/NSGA-II, Simulated Annealing, Tabu search, MOEA/D — all share equal evaluation budgets (Design A: 3 000 evaluations).

**Output metrics.** Hypervolume of the Pareto front (multi-objective placement), best ROI, risk at best ROI, and feasible fraction.

---

## Output files (in `results/`)

| File | Contents |
|------|----------|
| `full_production_results.rds` | Design A raw rows (produced by R/06) |
| `revision_main_results.rds` | Revision experiment rows (R/11) |
| `revision_means.csv` | Means by arm and placement |
| `revision_per_instance.csv` | Per-instance Cliff's delta (RISK vs EFF) |
| `revision_summary.rds` | Bootstrap confidence intervals |
| `revision_audit.csv` | Zero-range and cut-trace diagnostics |
| `extended_main_results.rds` | 50-instance extension rows (R/12) |
| `extended_per_instance.csv` | Per-instance delta for extended set |
| `extended_density_correlation.csv` | Spearman rs vs. density for 70 instances |
| `bimodal_predictor.csv` | 20-instance bimodal analysis |
| `combined70_predictor.csv` | Combined 70-instance density predictor |
| `checkpoints_extended/` | Per-instance checkpoint RDS files |

---

## Reproducing specific paper numbers

| Paper location | Script | Notes |
|---------------|--------|-------|
| §8 hypervolumes | `R/02_illustration.R` | Checked in `smoke_test.R` stage "02" |
| §9.4 algorithm comparison | `R/05_worked_example.R` | Checked in smoke test stage "05" |
| Table 1 (Design A) | `R/07_full_scale_analysis.R` | Needs `full_production_results.rds` |
| §11 revision results | `R/11_revision_experiments.R` | Via `run_revision.R` |
| §13 density predictor | `R/12_extended_instances.R` | Via `run_extended.R` + `run_extended_analysis.R` |

---

## Developer notes

- Each R file in `R/` is self-contained: it defines functions at the top and runs its experiment below. `load_functions.R` provides `load_definitions()` and `run_toplevel()` to separate these two concerns — use them to load just the functions without triggering a full run.
- `parallel_map()` in R/06 falls back to `parLapply` on Windows and `mclapply` on Unix/macOS. Set `n_cores` to `1` to disable parallelism.
- All long runs are checkpointed. Rerunning the same command resumes from the last completed instance without recomputing.
- The random number stream in §8 (R/02) and §9.4 (R/05) is fixed: `run_all.R` sources them consecutively and the two stages share a stream. Sourcing R/05 without R/02 first will produce different numbers (but `smoke_test.R` handles this correctly).
