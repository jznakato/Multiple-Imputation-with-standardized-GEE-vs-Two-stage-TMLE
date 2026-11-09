## =============================================================================
## 04_summarize.R
## Stack the per-trial results (res_001.rds ... res_500.rds) and compute
## performance metrics and computation times for each estimator.
## =============================================================================
rm(list = ls())
library(dplyr)

## ---- paths (same as 03_run_sims.R) -------------------------------------------------
base_dir <- "/home/jznakato/MI_vs_TMLE_DGP2"
scen     <- "dgp2_main_J100_b0-1_bA1.5"
#run_label <- "m5_maxit10_AL-passive_slfull_exch"   # the results folder name from 03_run_sims.R
#run_label <- "m5_maxit10_AL-passive_slfull_exch_v2"
run_label<-"tmle_glm_only"
res_dir  <- file.path(base_dir, "sim_results", scen, run_label)
out_dir  <- file.path(res_dir, "summary")
dir.create(out_dir, showWarnings = FALSE)

## ---- 1. read and stack every trial ----------------------------------------------------
files <- list.files(res_dir, pattern = "^res_\\d+\\.rds$", full.names = TRUE)
message(length(files), " result files found")

all_out <- lapply(files, readRDS)
results <- bind_rows(lapply(all_out, `[[`, "res"))
timings <- bind_rows(lapply(all_out, `[[`, "timings"))

saveRDS(list(results = results, timings = timings),
        file.path(out_dir, "stacked_results.rds"))

## ---- 2. failures ----------------------------------------------------------------------
## trial-level failures (whole trial errored; no res file was written)
list.files(res_dir, pattern = "^ERROR_")

## estimator-level failures within trials
fail_tab <- timings %>%
  group_by(estimator) %>%
  summarise(n_trials = n(),
            n_failed = sum(status != "ok"),
            pct_failed = 100 * mean(status != "ok"),
            .groups = "drop")
fail_tab

## ---- 3. performance metrics -----------------------------------------------------------
## One row per estimator x estimand, over all successful trials.
##   bias      : mean(est - truth)
##   emp_se    : SD of the point estimates across trials        (sigma)
##   mean_se   : average estimated SE                             (sigma-hat)
##   se_ratio  : mean_se / emp_se
##   coverage  : % of 95% CIs containing the truth
##   power     : % of trials with p < 0.05 (only meaningful for RD/RR rows)
##   rmse      : sqrt(bias^2 + emp_se^2)
##   mcse_*    : Monte Carlo SEs of bias and coverage
## TMLE returns the RR SE on the log scale; put it on the ratio scale
## (delta method: SE(RR) = RR * SE(log RR)). CIs, coverage, p-values unchanged.
results <- results %>%
  mutate(se = ifelse(estimator == "tmle" & grepl("^RR_", estimand), se * est, se))
results <- results %>%
  mutate(reject_fix = ifelse(grepl("^RR_", estimand),
                             as.numeric(CI.lo > 1 | CI.hi < 1),
                             reject))
perf <- results %>%
  filter(!is.na(estimand), !is.na(est)) %>%
  group_by(estimator, estimand) %>%
  summarise(
    n_reps    = n(),
    mean_est  = mean(est),
    mean_CIlo = mean(CI.lo),
    mean_CIhi = mean(CI.hi),
    bias      = mean(bias),
    emp_se    = sd(est),
    mean_se   = mean(se),
    se_ratio  = mean(se) / sd(est),
    coverage  = 100 * mean(cover),
    power     = 100 * mean(reject_fix),
    rmse      = sqrt(bias^2 + emp_se^2),
    mcse_bias = sd(est) / sqrt(n()),
    mcse_cov  = 100 * sqrt(mean(cover) * (1 - mean(cover)) / n()),
    .groups   = "drop"
  ) %>%
  arrange(estimand, estimator)

perf
write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)

## cluster-level rows only, in a fixed estimator order (for the cATE table)
est_order <- c("mi_fcs_nopassive", "mi_to_nopassive",
               "mi_fcs_passive", "mi_to_passive", "tmle")
perf_cATE <- perf %>%
  filter(grepl("_cATE$", estimand)) %>%
  mutate(estimator = factor(estimator, levels = est_order)) %>%
  arrange(estimand, estimator)
perf_cATE

## ---- 4. computation time --------------------------------------------------------------
time_tab <- timings %>%
  filter(status == "ok") %>%
  group_by(estimator) %>%
  summarise(
    n_trials    = n(),
    mean_sec    = mean(elapsed_sec),
    sd_sec      = sd(elapsed_sec),
    median_sec  = median(elapsed_sec),
    mean_min    = mean(elapsed_sec) / 60,
    total_hours = sum(elapsed_sec) / 3600,
    .groups = "drop"
  ) %>%
  mutate(ratio_to_tmle = mean_sec / mean_sec[estimator == "tmle"]) %>%
  arrange(mean_sec)

time_tab
write.csv(time_tab, file.path(out_dir, "timings.csv"), row.names = FALSE)

## total compute time per trial (all 5 estimators)
timings %>% group_by(rep) %>% summarise(total_min = sum(elapsed_sec) / 60) %>%
  summarise(mean_min = mean(total_min), median_min = median(total_min))
