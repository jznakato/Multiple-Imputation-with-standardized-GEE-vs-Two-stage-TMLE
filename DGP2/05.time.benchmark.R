## =============================================================================
## 05_timing_benchmark.R
## Timing-only benchmark with the OFF-THE-SHELF MRStdCRT package GEE.
##
## Purpose: report computation times that reflect what a typical user would
## experience, i.e. calling MRStdCRT::MRStdCRT_fit() from the package once per
## scale (RD, RR) on each imputed dataset, instead of the faster shared-
## jackknife version used in the main simulation.
##
## Runs on a SUBSET of the saved trials. Does not touch the main results.
## Run it AFTER the main simulation has finished, so the machine is not shared
## with the 100 main workers.
## =============================================================================
rm(list = ls())
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)
library(future.apply)

## ---- paths and settings (same as 03_run_sims.R) -------------------------------
base_dir <- "/home/jznakato/MI_vs_TMLE_DGP2"
code_dir <- base_dir
scen     <- "dgp2_main_J100_b0-1_bA1.5"
data_dir <- file.path(base_dir, "sim_data", scen)

s <- list(SL.library             = c("SL.glm"),
          sl_name                = "slGLM",
          m                      = 5,
          maxit_fcs              = 10,
          nonpassive_interaction = "passive",
          corstr                 = "exchangeable")

n_bench   <- 20                     # number of trials to time
bench_reps <- 1:n_bench
n_workers <- 20                     # one trial per worker, a single wave

bench_dir <- file.path(base_dir, "sim_results", scen,
                       sprintf("timing_benchmark_package_m%d_maxit%d_AL-%s_%s",
                               s$m, s$maxit_fcs, s$nonpassive_interaction, s$sl_name))
dir.create(bench_dir, recursive = TRUE, showWarnings = FALSE)
message("Benchmark results go to: ", bench_dir)

truth <- readRDS(file.path(data_dir, "truth.rds"))
psi <- c(RD_indiv  = truth$RD.ind,    RR_indiv  = truth$RR.ind,
         RD_clust  = truth$RD.clust,  RR_clust  = truth$RR.clust,
         mu1_indiv = truth$int.ind,   mu0_indiv = truth$con.ind,
         mu1_clust = truth$int.clust, mu0_clust = truth$con.clust)

## ---- load code (same files as the main run) -------------------------------------
load_code <- function() {
  suppressPackageStartupMessages({
    library(dplyr); library(mice); library(miceadds); library(lme4)
    library(ltmle); library(MRStdCRT); library(SuperLearner)
  })
  files <- c("generate_data_v3.R", "do.stage1.R", "stage2/stage2.R",
             "stage2/tmle.R", "aps.update.R", "run.estimators.R",
             "do.mi.R", "MRStdCRT_fix.R", "02_estimators.R")
  for (f in files) source(file.path(code_dir, f), local = globalenv())
}
load_code()

## The PACKAGE's own function, not the modified copy in MRStdCRT_fix.R
pkg_fit <- MRStdCRT::MRStdCRT_fit
print(args(pkg_fit))     # check the argument names below match

## ---- off-the-shelf analysis of ONE imputed dataset -------------------------------
## As a user would run it: same data prep and trtprob as run_gee_std(), then the
## package function once for RD and once for RR. Each package call returns the
## cluster- and individual-level effects together.
gee_package_once <- function(d, corstr) {
  data_std <- d[d$Y1 == 1, ]
  prob <- data_std %>%
    group_by(id) %>% mutate(first_trt = first(A)) %>% ungroup() %>%
    mutate(prob_A_1 = mean(first_trt == 1, na.rm = TRUE),
           prob_A_0 = mean(first_trt == 0, na.rm = TRUE)) %>%
    mutate(assigned_value = ifelse(A == 1, prob_A_1, prob_A_0)) %>%
    pull(assigned_value)
  data_std$cluster <- data_std$id
  
  fit_args <- list(formula = Y2 ~ A + W1 + W2 + W3 + E1 + E2 + n,
                   data = data_std, cluster = "cluster", trt = "A",
                   trtprob = prob, family = binomial(link = "logit"),
                   method = "GEE", corstr = corstr)
  if ("jack" %in% names(formals(pkg_fit))) fit_args$jack <- 1
  
  fit_rd <- do.call(pkg_fit, c(fit_args, list(scale = "RD")))
  fit_rr <- do.call(pkg_fit, c(fit_args, list(scale = "RR")))
  invisible(list(fit_rd, fit_rr))
}

## ---- time one MI version: imputation and analysis separately ----------------------
time_mi <- function(Obs, time_ordered, passive, s) {
  fcs   <- !time_ordered
  maxit <- if (fcs) s$maxit_fcs else 1
  
  t_imp <- system.time({
    if (passive) {
      mi_spec <- make_mi_spec_passive(Obs, fcs = fcs)
      imps <- do.imputation.passive(data = mi_spec$dat, mi_spec = mi_spec, m = s$m,
                                    maxit = maxit,
                                    visitSequence = c("Y1", "L", "A_L", "Y2_star", "Y2"))
    } else {
      mi_spec <- make_mi_spec(Obs, fcs = fcs, interaction = s$nonpassive_interaction)
      imps <- do.imputation(data = mi_spec$dat, mi_spec = mi_spec, m = s$m,
                            maxit = maxit, visitSequence = c("Y1", "L", "A_L", "Y2"))
    }
  })[["elapsed"]]
  
  t_ana <- system.time(
    for (d in imps) gee_package_once(d, corstr = s$corstr)
  )[["elapsed"]]
  
  c(impute_sec = t_imp, analysis_sec = t_ana)
}

## ---- one benchmark trial -----------------------------------------------------------------
bench_one <- function(r) {
  out_file <- file.path(bench_dir, sprintf("time_%03d.rds", r))
  if (file.exists(out_file)) return("skipped")
  if (!exists("run_tmle", envir = globalenv(), mode = "function")) load_code()
  
  tryCatch({
    Obs <- observed_only(readRDS(file.path(data_dir, sprintf("rep_%03d.rds", r))))
    set.seed(200000 + r)
    
    versions <- list(
      mi_fcs_nopassive = c(time_ordered = FALSE, passive = FALSE),
      mi_to_nopassive  = c(time_ordered = TRUE,  passive = FALSE),
      mi_fcs_passive   = c(time_ordered = FALSE, passive = TRUE),
      mi_to_passive    = c(time_ordered = TRUE,  passive = TRUE)
    )
    
    t_tmle <- system.time(suppressMessages(suppressWarnings(
      run_tmle(Obs, psi, SL.library = s$SL.library))))[["elapsed"]]
    
    mi_rows <- lapply(names(versions), function(nm) {
      v  <- versions[[nm]]
      tt <- suppressMessages(suppressWarnings(
        time_mi(Obs, v[["time_ordered"]], v[["passive"]], s)))
      data.frame(estimator = nm, impute_sec = tt[["impute_sec"]],
                 analysis_sec = tt[["analysis_sec"]])
    })
    
    out <- rbind(
      data.frame(estimator = "tmle", impute_sec = NA, analysis_sec = NA),
      do.call(rbind, mi_rows)
    )
    out$total_sec <- ifelse(out$estimator == "tmle", t_tmle,
                            out$impute_sec + out$analysis_sec)
    out$rep <- r
    saveRDS(out, out_file)
    "done"
  }, error = function(e) paste("error:", conditionMessage(e)))
}

## =============================================================================
## BENCHMARK: all n_bench trials, one per worker (set up to run unattended)
## =============================================================================
todo <- bench_reps[!file.exists(file.path(bench_dir, sprintf("time_%03d.rds", bench_reps)))]
message(length(todo), " benchmark trials to run")
plan(multicore, workers = n_workers)
status <- future_lapply(todo, bench_one, future.seed = TRUE, future.chunk.size = 1)
plan(sequential)
print(table(unlist(status)))

## =============================================================================
## SUMMARY
## =============================================================================
bt <- dplyr::bind_rows(lapply(list.files(bench_dir, "^time_\\d+\\.rds$", full.names = TRUE), readRDS))
bench_summary <- bt %>%
  dplyr::group_by(estimator) %>%
  dplyr::summarise(n = dplyr::n(),
                   impute_min   = mean(impute_sec) / 60,
                   analysis_min = mean(analysis_sec) / 60,
                   total_min    = mean(total_sec) / 60,
                   sd_total_min = sd(total_sec) / 60,
                   .groups = "drop") %>%
  dplyr::mutate(ratio_to_tmle = total_min / total_min[estimator == "tmle"])
print(bench_summary)
write.csv(bench_summary, file.path(bench_dir, "benchmark_summary.csv"), row.names = FALSE)
