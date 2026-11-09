## =============================================================================
## 02_estimators.R
## Estimator functions for the DGP 2 simulation:
##   - Two-stage TMLE
##   - 4 versions of MI + standardized GEE (2 x 2: time-ordering x structural zero)
## Every estimator returns rows with the SAME columns:
##   estimand, est, se, CI.lo, CI.hi, pval, bias, cover, reject
##
## This file REPLACES the old looper(). Do not source the old looper file after
## this one (duplicate function definitions overwrite each other).
##
## Requires (sourced beforehand): getYc(), get.weights2(), Stage2(),
##   make_mi_spec(), make_mi_spec_passive(), do.imputation(),
##   do.imputation.passive(), MRStdCRT_point() (from MRStdCRT_fix.R), pool_rubin()
## =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(mice)
  library(miceadds)   # 2l.binary
})

## -----------------------------------------------------------------------------
## 0. Helpers
## -----------------------------------------------------------------------------

## Keep only OBSERVED-data columns, so no estimator can see counterfactuals
## (Y1_star, L_1, Y2star_1, the true Y2_star, ...).
OBS_COLS <- c("id", "n", "E1", "E2", "W1", "W2", "W3",
              "A", "Delta", "Y1", "L", "Delta2", "Y2")

observed_only <- function(O) {
  missing_cols <- setdiff(OBS_COLS, names(O))
  if (length(missing_cols) > 0)
    stop("observed_only(): missing columns: ", paste(missing_cols, collapse = ", "))
  O[, OBS_COLS]
}

## Run one estimator: returns value, elapsed seconds, status, warning count.
## `expr` is evaluated lazily inside, so the clock starts before it runs.
timed_run <- function(expr) {
  t0  <- proc.time()
  val <- tryCatch(suppressWarnings(expr), error = function(e) e)
  elapsed <- (proc.time() - t0)[["elapsed"]]
  list(value   = if (inherits(val, "error")) NULL else val,
       elapsed = elapsed,
       status  = if (inherits(val, "error")) paste("error:", conditionMessage(val)) else "ok")
}

## -----------------------------------------------------------------------------
## 1. Two-stage TMLE
## -----------------------------------------------------------------------------

## Pull Txt (mu1) or Con (mu0) from a Stage2 output as one row.
extract_tmle_mean <- function(stage2_out, side = c("Txt", "Con"),
                              true_mean, df_for_t, use_cv = TRUE) {
  side   <- match.arg(side)
  prefix <- if (use_cv) paste0(side, ".CV") else side
  est   <- stage2_out[[paste0(prefix, ".est")]]
  CI.lo <- stage2_out[[paste0(prefix, ".CI.lo")]]
  CI.hi <- stage2_out[[paste0(prefix, ".CI.hi")]]
  se    <- stage2_out[[paste0(prefix, ".se")]]
  pval  <- 2 * pt(abs(est / se), df = df_for_t, lower.tail = FALSE)
  data.frame(est = est, se = se, CI.lo = CI.lo, CI.hi = CI.hi, pval = pval,
             bias   = est - true_mean,
             cover  = as.numeric(CI.lo <= true_mean & true_mean <= CI.hi),
             reject = as.numeric(pval < 0.05))
}

## Contrast row (RD or RR) from a Stage2 output, CV inference.
extract_tmle_contrast <- function(stage2_out) {
  x <- stage2_out[, c("CV.est", "CV.se", "CV.CI.lo", "CV.CI.hi",
                      "CV.pval", "CV.bias", "CV.cover", "CV.reject")]
  names(x) <- c("est", "se", "CI.lo", "CI.hi", "pval", "bias", "cover", "reject")
  x$cover  <- as.numeric(x$cover)
  x$reject <- as.numeric(x$reject)
  x
}

run_tmle <- function(Obs, psi, SL.library = "SL.glm") {
  
  indv_cov  <- c("W1", "W2", "W3")
  clust_cov <- c("E1", "E2")
  cand_cov  <- c("U", "E1", "E2")
  
  ## Stage 1
  Yc <- getYc(O = Obs, indv_cov = indv_cov, clust_cov = clust_cov, SL.library = SL.library)
  data.input <- cbind(Yc[, c("id", "n", "alpha", "A", "den_tmle", cand_cov)], Y = Yc$tmle)
  df_t <- length(unique(data.input$id)) - 2
  
  ## Weighted copies made from the SAME unweighted data (never re-weighting
  ## an already-weighted dataset), and passed explicitly to Stage2.
  data_indv  <- get.weights2(data = data.input, weighting = "indv")
  data_clust <- get.weights2(data = data.input, weighting = "clust")
  
  stage2 <- function(goal, target, dat, psi_val) {
    Stage2(goal = goal, target = target, data.input = dat,
           do.data.adapt = TRUE, do.cv.variance = TRUE,
           cand.QAdj = cand_cov, cand.Qform = "glm",
           cand.gAdj = cand_cov, cand.gform = "glm",
           verbose = FALSE, psi = psi_val)
  }
  
  ## Individual-level (iATE)
  rd_i <- stage2("RD",  "indv", data_indv, psi[["RD_indiv"]])
  rr_i <- stage2("aRR", "indv", data_indv, psi[["RR_indiv"]])
  
  ## Cluster-level (cATE)
  rd_c <- stage2("RD",  "clust", data_clust, psi[["RD_clust"]])
  rr_c <- stage2("aRR", "clust", data_clust, psi[["RR_clust"]])
  
  dplyr::bind_rows(
    cbind(estimand = "RD_iATE",  extract_tmle_contrast(rd_i)),
    cbind(estimand = "RR_iATE",  extract_tmle_contrast(rr_i)),
    cbind(estimand = "mu1_iATE", extract_tmle_mean(rd_i, "Txt", psi[["mu1_indiv"]], df_t)),
    cbind(estimand = "mu0_iATE", extract_tmle_mean(rd_i, "Con", psi[["mu0_indiv"]], df_t)),
    cbind(estimand = "RD_cATE",  extract_tmle_contrast(rd_c)),
    cbind(estimand = "RR_cATE",  extract_tmle_contrast(rr_c)),
    cbind(estimand = "mu1_cATE", extract_tmle_mean(rd_c, "Txt", psi[["mu1_clust"]], df_t)),
    cbind(estimand = "mu0_cATE", extract_tmle_mean(rd_c, "Con", psi[["mu0_clust"]], df_t))
  )
}

## -----------------------------------------------------------------------------
## 2. MI helpers
## Uses YOUR existing functions (sourced from do.mi.R before this file):
##   make_mi_spec()           non-passive spec (interaction = "jav" or "passive")
##   make_mi_spec_passive()   passive spec (Y2 through working variable Y2_star)
##   do.imputation()          mice for the non-passive spec
##   do.imputation.passive()  mice for the passive spec (passes where_mat)
##   pool_rubin()             Rubin's rules
## -----------------------------------------------------------------------------

## -----------------------------------------------------------------------------
## MRStdCRT_fit_all(): RD, RR and log RR at both levels from ONE jackknife
## -----------------------------------------------------------------------------
## A new helper; MRStdCRT_point() and MRStdCRT_fit() are NOT modified.
## It does what MRStdCRT_fit() does (one MRStdCRT_point() call on the full data
## plus one per left-out cluster, jack = 1 SE, t(m-1) CIs and p-values), but it
## keeps the four arm means (m1c, m0c, m1i, m0i) from every jackknife replicate
## and computes every contrast from them. RD and RR therefore come from the
## same 1 + J GEE fits instead of two separate MRStdCRT_fit() calls.
## `scale` in MRStdCRT_point() only affects elements 1:3 of its output; the
## means (elements 4:7) are the same for any scale.
MRStdCRT_fit_all <- function(formula, data, cluster, trt, trtprob,
                             method = "GEE", family = binomial(link = "logit"),
                             corstr = "exchangeable", alpha = 0.05) {
  
  get_means <- function(d, p) {
    MRStdCRT_point(formula, d, cluster, trt, p, family, corstr, method, scale = "RD")[[3]][4:7]
  }
  contrasts <- function(x) {
    x <- unname(x)                               # m1c, m0c, m1i, m0i
    c(mu1_clust = x[1], mu0_clust = x[2], mu1_indiv = x[3], mu0_indiv = x[4],
      RD_clust  = x[1] - x[2],      RD_indiv  = x[3] - x[4],
      RR_clust  = x[1] / x[2],      RR_indiv  = x[3] / x[4],
      lRR_clust = log(x[1] / x[2]), lRR_indiv = log(x[3] / x[4]))
  }
  
  ## full-data estimate
  full <- get_means(data, trtprob)
  
  ## leave-one-cluster-out jackknife (same loop as MRStdCRT_fit)
  cluster_names <- sort(unique(data[[cluster]]))
  m <- length(cluster_names)
  jack <- t(vapply(cluster_names, function(i) {
    keep <- data[[cluster]] != i
    get_means(droplevels(data[keep, , drop = FALSE]), trtprob[keep])
  }, numeric(4)))
  
  est <- contrasts(full)
  jk  <- t(apply(jack, 1, contrasts))
  se  <- sqrt(apply(jk, 2, var) * (m - 1)^2 / m)        # jack = 1 in MRStdCRT_fit
  thres <- qt(1 - alpha / 2, df = m - 1)
  
  data.frame(estimand = names(est),
             Estimate = unname(est),
             SE       = unname(se),
             CI.lo    = unname(est - thres * se),
             CI.hi    = unname(est + thres * se),
             pval     = unname(2 * pt(abs(est / se), df = m - 1, lower.tail = FALSE)),
             m        = m,
             row.names = NULL)
}

## -----------------------------------------------------------------------------
## run_gee_std_both(): your run_gee_std(), both levels from one call
## -----------------------------------------------------------------------------
## Same data prep, SAME trtprob calculation, same formula and GEE settings as
## run_gee_std(). Uses MRStdCRT_fit_all() instead of two MRStdCRT_fit() calls.
## Returns list(indiv = one wide row, clust = one wide row) with blocks
## RD_*, RR_*, mu1_*, mu0_*, lRR_* and the same columns as run_gee_std().
run_gee_std_both <- function(O, psi) {
  
  needed <- c("RD_clust","RD_indiv","RR_clust","RR_indiv",
              "mu1_clust","mu1_indiv","mu0_clust","mu0_indiv")
  missing <- setdiff(needed, names(psi))
  if (length(missing)) {
    stop("psi is missing required entries: ", paste(missing, collapse = ", "))
  }
  
  # --- data prep (unchanged from run_gee_std) ---
  data_std <- O[O$Y1 == 1, ]
  data_std_prob <- data_std %>%
    group_by(id) %>% mutate(first_trt = first(A)) %>% ungroup() %>%
    mutate(prob_A_1 = mean(first_trt == 1, na.rm = TRUE),
           prob_A_0 = mean(first_trt == 0, na.rm = TRUE)) %>%
    mutate(assigned_value = ifelse(A == 1, prob_A_1, prob_A_0))
  prob <- data_std_prob$assigned_value
  data_std$cluster <- data_std$id
  
  # --- one fit + one jackknife for everything ---
  fit <- MRStdCRT_fit_all(formula = Y2 ~ A + W1 + W2 + W3 + E1 + E2 + n,
                          data    = data_std,
                          cluster = "cluster",
                          trt     = "A",
                          trtprob = prob,
                          method  = "GEE",
                          family  = binomial(link = "logit"),
                          corstr  = "exchangeable")
  
  # --- same columns as pull_row() in run_gee_std ---
  pull_row <- function(name, target) {
    e <- fit[fit$estimand == name, ]
    data.frame(
      est    = e$Estimate,
      CI.lo  = e$CI.lo,
      CI.hi  = e$CI.hi,
      se     = e$SE,
      var    = e$SE^2,
      pval   = e$pval,
      bias   = e$Estimate - target,
      cover  = as.numeric(e$CI.lo <= target & target <= e$CI.hi),
      reject = as.numeric(e$pval < 0.05)
    )
  }
  prefix <- function(df, p) setNames(df, paste0(p, "_", names(df)))
  
  one_level <- function(suffix) {
    cbind(
      prefix(pull_row(paste0("RD_",  suffix), psi[[paste0("RD_",  suffix)]]), "RD"),
      prefix(pull_row(paste0("RR_",  suffix), psi[[paste0("RR_",  suffix)]]), "RR"),
      prefix(pull_row(paste0("mu1_", suffix), psi[[paste0("mu1_", suffix)]]), "mu1"),
      prefix(pull_row(paste0("mu0_", suffix), psi[[paste0("mu0_", suffix)]]), "mu0"),
      prefix(pull_row(paste0("lRR_", suffix), log(psi[[paste0("RR_", suffix)]])), "lRR")
    )
  }
  
  list(indiv = one_level("indiv"), clust = one_level("clust"))
}

## Pool one estimand block (e.g. all "RD_*" columns) with Rubin's rules.
pool_block <- function(res_imp, prefix, true_est, estimand_label) {
  cols  <- grep(paste0("^", prefix, "_"), names(res_imp), value = TRUE)
  block <- res_imp[, c("imp", cols)]
  names(block) <- sub(paste0("^", prefix, "_"), "", names(block))
  out <- pool_rubin(block, true_est = true_est)
  data.frame(estimand = estimand_label,
             est = out$est, se = out$se, CI.lo = out$CI.lo, CI.hi = out$CI.hi,
             pval = out$pval, bias = out$bias,
             cover = as.numeric(out$cover), reject = as.numeric(out$reject))
}

## -----------------------------------------------------------------------------
## 3. MI + standardized GEE (one function, 4 versions)
## -----------------------------------------------------------------------------
##   time_ordered = FALSE, passive = FALSE  -> MI 1: ignores both
##   time_ordered = TRUE,  passive = FALSE  -> MI 2: time-ordering only
##   time_ordered = FALSE, passive = TRUE   -> MI 3: structural zero only
##   time_ordered = TRUE,  passive = TRUE   -> MI 4: both
run_mi <- function(Obs, psi, time_ordered, passive, m = 5, maxit_fcs = 10,
                   nonpassive_interaction = "jav") {
  
  fcs   <- !time_ordered
  maxit <- if (fcs) maxit_fcs else 1          # time-ordered: single pass
  
  if (passive) {
    mi_spec <- make_mi_spec_passive(Obs, fcs = fcs)
    vs      <- c("Y1", "L", "A_L", "Y2_star", "Y2")
    imps    <- do.imputation.passive(data = mi_spec$dat, mi_spec = mi_spec,
                                     m = m, maxit = maxit, visitSequence = vs)
  } else {
    mi_spec <- make_mi_spec(Obs, fcs = fcs, interaction = nonpassive_interaction)
    vs      <- c("Y1", "L", "A_L", "Y2")
    imps    <- do.imputation(data = mi_spec$dat, mi_spec = mi_spec,
                             m = m, maxit = maxit, visitSequence = vs)
  }
  
  imps <- lapply(imps, function(d) d[which(d$Y1 == 1), ])   # focus population
  
  gee   <- lapply(imps, run_gee_std_both, psi = psi)          # one call per imputation, both levels
  res_i <- dplyr::bind_rows(lapply(gee, `[[`, "indiv"), .id = "imp")
  res_c <- dplyr::bind_rows(lapply(gee, `[[`, "clust"), .id = "imp")
  
  dplyr::bind_rows(
    pool_block(res_i, "RD",  psi[["RD_indiv"]],  "RD_iATE"),
    pool_block(res_i, "RR",  psi[["RR_indiv"]],  "RR_iATE"),
    pool_block(res_i, "mu1", psi[["mu1_indiv"]], "mu1_iATE"),
    pool_block(res_i, "mu0", psi[["mu0_indiv"]], "mu0_iATE"),
    pool_block(res_c, "RD",  psi[["RD_clust"]],  "RD_cATE"),
    pool_block(res_c, "RR",  psi[["RR_clust"]],  "RR_cATE"),
    pool_block(res_c, "mu1", psi[["mu1_clust"]], "mu1_cATE"),
    pool_block(res_c, "mu0", psi[["mu0_clust"]], "mu0_cATE")
  )
}

## -----------------------------------------------------------------------------
## 4. All five estimators on one dataset, each timed
## -----------------------------------------------------------------------------
## s = settings list with: SL.library, m, maxit_fcs, nonpassive_interaction
ESTIMATORS <- list(
  tmle             = function(Obs, psi, s) run_tmle(Obs, psi, SL.library = s$SL.library),
  mi_fcs_nopassive = function(Obs, psi, s) run_mi(Obs, psi, time_ordered = FALSE, passive = FALSE, m = s$m, maxit_fcs = s$maxit_fcs, nonpassive_interaction = s$nonpassive_interaction),
  mi_to_nopassive  = function(Obs, psi, s) run_mi(Obs, psi, time_ordered = TRUE,  passive = FALSE, m = s$m, maxit_fcs = s$maxit_fcs, nonpassive_interaction = s$nonpassive_interaction),
  mi_fcs_passive   = function(Obs, psi, s) run_mi(Obs, psi, time_ordered = FALSE, passive = TRUE,  m = s$m, maxit_fcs = s$maxit_fcs),
  mi_to_passive    = function(Obs, psi, s) run_mi(Obs, psi, time_ordered = TRUE,  passive = TRUE,  m = s$m, maxit_fcs = s$maxit_fcs)
)

## seed: if given, each estimator gets its own seed (seed + k), so results for
## one estimator do not change if another estimator is added or removed.
run_all_estimators <- function(O, psi, s, rep, estimators = names(ESTIMATORS),
                               seed = NULL) {
  Obs <- observed_only(O)
  out <- lapply(seq_along(estimators), function(k) {
    nm <- estimators[k]
    if (!is.null(seed)) set.seed(seed + k)
    tr  <- timed_run(ESTIMATORS[[nm]](Obs, psi, s))
    res <- if (is.null(tr$value)) data.frame(estimand = NA_character_)
    else tr$value
    res <- cbind(estimator = nm, res)
    tim <- data.frame(estimator = nm, elapsed_sec = tr$elapsed,
                      status = tr$status)
    list(res = res, tim = tim)
  })
  res <- dplyr::bind_rows(lapply(out, `[[`, "res")); res$rep <- rep
  tim <- dplyr::bind_rows(lapply(out, `[[`, "tim")); tim$rep <- rep
  list(res = res, timings = tim)
}

