## =============================================================================
## 03_run_sims.R
## Read each saved trial, run all 5 estimators, save one result file per trial.
## Resumable: trials whose result file already exists are skipped.
## =============================================================================
rm(list = ls())


## One thread per process, so timings are not distorted by multithreaded BLAS
Sys.setenv(OMP_NUM_THREADS = 1, OPENBLAS_NUM_THREADS = 1, MKL_NUM_THREADS = 1)

library(future.apply)

## ---- paths -------------------------------------------------------------------
base_dir <- "/home/jznakato/MI_vs_TMLE_DGP2"
code_dir <- base_dir                         # folder holding your .R files
scen     <- "dgp2_main_J100_b0-1_bA1.5"
data_dir <- file.path(base_dir, "sim_data", scen)

## ---- estimator settings ----------------------------------------------------------
s <- list(SL.library             = c("SL.glm", "SL.earth", "SL.mean","SL.gam"),
          sl_name                = "slfull",
          m                      = 5,
          maxit_fcs              = 10,
          nonpassive_interaction = "passive")   # "passive" (clean 2x2) or "jav"

s$SL.library <- "SL.glm"
s$sl_name    <- "glm"
## new run_label, e.g. paste0(..., "_tmle_glm_only"), so it gets its own folder

# run_label <- sprintf("m%d_maxit%d_AL-%s_%s_exch_v2", s$m, s$maxit_fcs,
#                      s$nonpassive_interaction, paste(s$sl_name, collapse = "+"))

run_label <- "tmle_glm_only" 


res_dir <- file.path(base_dir, "sim_results", scen, run_label)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(s, file.path(res_dir, "estimator_settings.rds"))
message("Results go to: ", normalizePath(res_dir))

## ---- truth -------------------------------------------------------------------------
truth <- readRDS(file.path(data_dir, "truth.rds"))
psi <- c(RD_indiv  = truth$RD.ind,    RR_indiv  = truth$RR.ind,
         RD_clust  = truth$RD.clust,  RR_clust  = truth$RR.clust,
         mu1_indiv = truth$int.ind,   mu0_indiv = truth$con.ind,
         mu1_clust = truth$int.clust, mu0_clust = truth$con.clust)



## ---- code loader (runs on each worker) ---------------------------------------------
load_code <- function() {
  suppressPackageStartupMessages({
    library(dplyr); library(mice); library(miceadds); library(lme4)
    library(ltmle); library(MRStdCRT); library(SuperLearner)
  })
  files <- c("generate_data_v3.R", "do.stage1.R", "stage2/stage2.R",
             "stage2/tmle.R", "aps.update.R", "run.estimators.R",
             "do.mi.R", "MRStdCRT_fix.R",
             "02_estimators_v2.R")                 # always last
  for (f in files) source(file.path(code_dir, f), local = globalenv())
}
load_code()

## ---- one trial -----------------------------------------------------------------------
run_one_rep <- function(r) {
  out_file <- file.path(res_dir, sprintf("res_%03d.rds", r))
  if (file.exists(out_file)) return("skipped")
  
  ## load code once per worker
  if (!exists("run_all_estimators", envir = globalenv(), mode = "function")) load_code()
  
  status <- tryCatch({
    O   <- readRDS(file.path(data_dir, sprintf("rep_%03d.rds", r)))
    out <- do.call("run_all_estimators",
                   list(O = O, psi = psi, s = s, rep = r,estimators = "tmle", seed = 100000 + 10 * r))
    saveRDS(out, out_file)
    "done"
  }, error = function(e) {
    saveRDS(list(error = conditionMessage(e), rep = r),
            file.path(res_dir, sprintf("ERROR_%03d.rds", r)))
    paste("error:", conditionMessage(e))
  })
  status
}

## =============================================================================
## A. TEST FIRST: 2 trials, sequential (errors print directly)
## =============================================================================
# plan(sequential)
# test_status <- sapply(1:2, run_one_rep)
# test_status
# readRDS(file.path(res_dir, "res_001.rds"))$timings
# 

# ## timings from the two test trials
# tim <- dplyr::bind_rows(lapply(1:2, function(r)
#   readRDS(file.path(res_dir, sprintf("res_%03d.rds", r)))$timings))
# tim
# 
# ## per estimator: average over the 2 trials
# tim %>%
#   dplyr::group_by(estimator) %>%
#   dplyr::summarise(mean_sec = mean(elapsed_sec),
#                    mean_min = mean(elapsed_sec) / 60,
#                    status   = paste(unique(status), collapse = ", "),
#                    .groups  = "drop") %>%
#   dplyr::arrange(dplyr::desc(mean_sec))
# 
# ## total per trial, and projected full-run time
# per_trial_min <- tim %>% dplyr::group_by(rep) %>%
#   dplyr::summarise(total = sum(elapsed_sec) / 60) %>% dplyr::pull(total)
# per_trial_min
# 
# workers <- c(50, 100, 150)
# data.frame(workers,
#            est_hours = round(498 * mean(per_trial_min) / workers / 60, 1))

## =============================================================================
## B. FULL RUN: uncomment once the test looks right
## =============================================================================
R    <- 500
todo <- (1:R)[!file.exists(file.path(res_dir, sprintf("res_%03d.rds", 1:R)))]
message(length(todo), " trials to run")

plan(multicore, workers = 100)
t0 <- Sys.time()
status <- future_lapply(todo, run_one_rep, future.seed = TRUE,future.chunk.size = 1)
plan(sequential)
print(Sys.time() - t0)

table(unlist(status))
list.files(res_dir, pattern = "^ERROR_")     # trials that failed outright

benchmark_summary <- read_csv("~/MI_vs_TMLE_DGP2/sim_results/dgp2_main_J100_b0-1_bA1.5/timing_benchmark_package_m5_maxit10_AL-passive_slGLM/benchmark_summary.csv")

