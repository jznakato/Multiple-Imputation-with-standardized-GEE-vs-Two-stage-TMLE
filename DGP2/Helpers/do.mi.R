

##UPDATE 10/23/2026 A NEW STANDARD MI SPECIFICATION FOR NON PASSIVE FCS
make_mi_spec <- function(dat, cluster_var = "id", fcs = FALSE,
                         interaction = c("passive", "jav")) {
  interaction <- match.arg(interaction)
  if (!is.integer(dat[[cluster_var]]))
    dat[[cluster_var]] <- as.integer(as.factor(dat[[cluster_var]]))
  dat$A_L <- dat$A * dat$L
  if (interaction == "jav") dat$A_L[dat$A == 0] <- 0
  
  vars <- names(dat)
  base <- c("E1", "E2", "W1", "W2", "W3")
  pred <- matrix(0, length(vars), length(vars), dimnames = list(vars, vars))
  
  if (!fcs) {                                   # time-ordered (sequential)
    pred["Y1", c(base, "A")]                    <- 1
    pred["L",  c(base, "A", "Y1")]              <- 1
    pred["Y2", c(base, "A", "L", "A_L", "Y1")]  <- 1
  } else {                                      # naive FCS: ignores time order
    pred["Y1", c(base, "A", "L", "A_L", "Y2")]  <- 1
    pred["L",  c(base, "A", "Y1", "Y2")]        <- 1   # not A_L (built from L)
    pred["Y2", c(base, "A", "L", "A_L", "Y1")]  <- 1
  }
  pred[c("Y1", "L", "Y2"), cluster_var] <- -2
  
  meth <- setNames(rep("", length(vars)), vars)
  meth[c("Y1", "L", "Y2")] <- "2l.binary"
  
  if (interaction == "passive") {
    meth["A_L"] <- "~ I(A * L)"
  } else {                                      # JAV
    meth["A_L"] <- "2l.lmer"
    pred["A_L", c(base, "A", "Y1","L", if (fcs) "Y2")] <- 1   # not L (perfect prediction)
    pred["A_L", cluster_var] <- -2
  }
  list(dat = dat, impmeth = meth, predmat = pred)
}

# ============================================================================
# make_mi_spec_passive()
#
# Constructs a MICE specification that respects the structural relationship
# Y2 = Y1 * Y2_star, where Y2_star is the latent "potential PrEP uptake among
# the at-risk" working variable.
#
# Why passive imputation:
#   The standard MICE setup imputes Y1 and Y2 as two binary outcomes related
#   by a fitted logistic model. In the OPAL-style design Y2 is *deterministi-
#   cally* zero whenever Y1 = 0, which creates a near-absorbing (0,0) state
#   on the doubly-missing (Delta = 0) rows and prevents the chain from mixing.
#   Passive imputation removes this pathology by encoding the deterministic
#   structure directly:
#
#       Y2 := ifelse(Y1 == 0, 0, Y2_star)
#
#   where Y2_star is imputed only on rows that are at risk (Y1 = 1). The
#   chain on Y1 conditions on baseline covariates only -- never on the
#   downstream Y2 or Y2_star -- so Y1 cannot be dragged toward 0 by the
#   structural-zero attractor.
#
# What the function returns:
#   $dat        -- the augmented dataset with Y2_star added and Y2 replaced
#                  by a placeholder that mice() will overwrite passively
#   $impmeth    -- method vector for mice()
#   $predmat    -- predictor matrix for mice()
#   $where_mat  -- where matrix for mice() (restricts Y2_star imputation
#                  away from Y1 == 0 rows)
#
# Usage (drop-in replacement for make_mi_spec()):
#
#   spec <- make_mi_spec_passive(dat)
#   imp  <- mice(spec$dat,
#                m               = M,
#                maxit           = maxit,
#                method          = spec$impmeth,
#                predictorMatrix = spec$predmat,
#                where           = spec$where_mat,
#                printFlag       = FALSE)
#
#   completed <- mice::complete(imp, action = "long", include = FALSE)
#   # completed$Y2 is correctly reconstructed; ignore completed$Y2_star
#   # downstream and run your GEE / pooling on Y2 as usual.
# ============================================================================


make_mi_spec_passive <- function(dat,
                                 cluster_var      = "id",
                                 baseline_cluster = c("W1", "W2", "W3", "E1", "E2"),
                                 baseline_indv    = character(0),
                                 A_var            = "A",
                                 L_var            = "L",
                                 method_y1        = "2l.binary",
                                 method_y2_star   = "2l.binary",
                                 fcs              = FALSE) {
  
  # ----- sanity checks ------------------------------------------------------
  required <- c(cluster_var, "Delta", "Y1", "Delta2", "Y2")
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop("make_mi_spec_passive(): missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  baseline_all <- c(baseline_cluster, baseline_indv, L_var)
  missing_bl   <- setdiff(baseline_all, names(dat))
  if (length(missing_bl) > 0) {
    stop("make_mi_spec_passive(): missing baseline covariates: ",
         paste(missing_bl, collapse = ", "))
  }
  
  # cluster ID must be integer or factor for 2l methods
  if (!is.factor(dat[[cluster_var]]) && !is.integer(dat[[cluster_var]])) {
    dat[[cluster_var]] <- as.integer(as.factor(dat[[cluster_var]]))
  }
  dat$A_L <- dat$A * dat$L 
  


  
  # NEW
  dat$Y2_star <- with(dat,
                      ifelse(!is.na(Y1) & Y1 == 1, Y2, NA_real_)
  )
  # Initial placeholder for Y2 (mice will overwrite this passively at every
  # iteration). Set to 0 wherever Y1 is known to be 0 (structural), and to
  # Y2_star otherwise -- any remaining NAs are fine, mice's passive equation
  # will fill them once Y1 and Y2_star are imputed.
  dat$Y2 <- with(dat,
                 ifelse(!is.na(Y1) & Y1 == 0, 0, Y2_star)
  )
  
  # ----- 2. method vector ---------------------------------------------------
  imp_vars <- names(dat)
  meth <- setNames(rep("", length(imp_vars)), imp_vars)
  meth["Y1"]      <- method_y1
  meth["Y2_star"] <- method_y2_star
  meth["Y2"]      <- "~ I(ifelse(Y1 == 0, 0, Y2_star))"
  
  meth["L"]   <- "2l.binary"
  meth["A_L"] <- "~ I(A * L)" ##
  # ----- 3. predictor matrix ------------------------------------------------
  
  pred <- matrix(0, nrow = length(imp_vars), ncol = length(imp_vars),
                 dimnames = list(imp_vars, imp_vars))
  if (!fcs) {
    # Sequential: each model conditions only on variables preceding it in time
    pred["Y1",      c(baseline_cluster, A_var)]               <- 1
    pred["L",       c(baseline_cluster, A_var, "Y1")]         <- 1
    pred["Y2_star", c(baseline_all, A_var, "A_L")]            <- 1
  } else {
    # FCS: each model conditions on all other substantive variables
    pred["Y1",      c(baseline_cluster, A_var, "L", "A_L", "Y2")] <- 1
    pred["L",       c(baseline_cluster, A_var, "Y1", "Y2")]       <- 1  # no A_L (derived from L)
    pred["Y2_star", c(baseline_all, A_var, "A_L")]                <- 1  # no Y2 (derived from Y2_star)
  }
  # Cluster random intercept for all three 2l.binary models (both versions)
  pred[c("Y1", "L", "Y2_star"), cluster_var] <- -2
  
  wh[, "L"]   <- is.na(dat$L)
  wh[, "A_L"] <- is.na(dat$A_L)
  
  # For Delta = 0 rows we do not know Y1 a priori, so Y2_star starts NA and
  # the where matrix leaves it TRUE. At each iteration, after Y1 is drawn:
  #   * rows currently imputed as Y1 = 1 will get a Y2_star draw and the
  #     passive equation propagates it into Y2.
  #   * rows currently imputed as Y1 = 0 will *also* get a Y2_star draw
  #     (because mice's `where` is fixed at setup, not re-evaluated), but
  #     the passive equation discards it: Y2 := 0. The wasted draw does not
  #     contaminate Y2, so this is harmless.
  # We do NOT impute Delta or Delta2; their NAs (if any) are not used by the
  # Y1 or Y2_star imputation models because they are not in the predictor
  # matrix.
  
  list(
    dat       = dat,
    impmeth   = meth,
    predmat   = pred,
    where_mat = wh
  )
}

############################################################
## 3. HELPER: DO IMPUTATION
############################################################


do.imputation <- function(data, mi_spec, m = 5, maxit = 1,
                          visitSequence = c("Y1", "L", "A_L", "Y2"),
                          return_mids = FALSE) {
  
  imp <- mice::mice(
    data            = data,
    method          = mi_spec$impmeth,
    predictorMatrix = mi_spec$predmat,
    maxit           = maxit,
    m               = m,
    visitSequence   = visitSequence,
    printFlag       = FALSE
  )
  
  if (return_mids) return(imp)        # for plot(imp) / convergence checks
  mice::complete(imp, "all")          # list of m completed datasets
}




do.imputation.passive<-function(data,mi_spec,m=5,maxit=1,
                                visitSequence = c("Y1", "L", "A_L", "Y2_star", "Y2"),
                                return_mids = FALSE){
  
  # Run MI
  imp <- mice::mice(
    data = data,  # better to keep as data.frame
    method = mi_spec$impmeth,
    predictorMatrix = mi_spec$predmat,
    maxit = maxit,
    m = m,
    visitSequence =visitSequence,
    where = mi_spec$where_mat,
    printFlag = FALSE
  )
  if (return_mids) return(imp)
  ##returns list of 5 imputed datasets
  imp_datasets <- complete(imp, "all")

  return(imp_datasets)
}




############################################################
## 4. HELPER: RUBIN'S RULES FOR ONE MONTE CARLO REPLICATE
# Since we have m imputed completed datasets, we need to pool results via Rubin's rules to get one final
#estimate and standard errors for that
############################################################


pool_rubin <- function(res_imp, true_est, alpha = 0.05) {
  
  m <- nrow(res_imp)
  
  Q_bar <- mean(res_imp$est) ## mean point estimate
  U_bar <- mean(res_imp$var) ## intra cluster variance
  B <- var(res_imp$est) ## inter cluster variance
  
  T_var <- U_bar + (1 + 1 / m) * B ## total variance
  se_pool <- sqrt(T_var) ## pooled standard error
  
  
  ## get degrees of freedom for t-distribution
  if (is.na(B) || B < .Machine$double.eps) {
    df <- Inf
  } else {
    df <- (m - 1) * (1 + U_bar / ((1 + 1 / m) * B))^2
  }
  
  crit<-qt(0.025,df=df,lower.tail=F) 
  
  ci_lo <- Q_bar - crit * se_pool
  ci_hi <- Q_bar + crit * se_pool
  
  test_stat <- Q_bar / se_pool
  pval <- if (is.finite(df)) {
    2 * pt(abs(test_stat), df = df, lower.tail = FALSE)
  } else {
    2 * pnorm(abs(test_stat), lower.tail = FALSE)
  }
  
  
  data.frame(
    est    = Q_bar,
    se     = se_pool,
    CI.lo  = ci_lo,
    CI.hi  = ci_hi,
    pval   = pval,
    bias   = Q_bar - true_est,
    cover  = as.logical(ci_lo < true_est & true_est < ci_hi),
    reject = as.logical(pval < alpha),
    U_bar  = U_bar,
    B      = B,
    T_var  = T_var,
    df     = df
  )
}

### Helper function to run MI twice on one dataset ####

