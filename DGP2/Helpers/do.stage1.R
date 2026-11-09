# Code by Joy Nakato & Laura B. Balzer
# jznakato@berkeley.edu & laura.balzer@berkeley.edu

##=================================================================================
# getYc: Function to get cluster-level endpoint in Stage 1 
## Inputs:
# O: observed data
# ind_cov: Individual-level covariates
# clust_cov: Cluster-level covariates 
# SL.library: SuperLearner libary (default=NULL uses GLM)
## Output: 
# Returns a cluster-level data set after iterating over all clusters
##=================================================================================
getYc <- function(O, indv_cov, clust_cov, SL.library = NULL) {
  
  clusters <- unique(O$id)
  J <- length(clusters)
  stage1_out <- data.frame(matrix(NA, nrow = J, ncol = 18))
  
  for (j in clusters) {
    Oc <- O[O$id == j, ]
    
    
## UPDATED CNODE implementation of LTMLE
    # ---- Data construction ----
    Oc_num <- data.frame(
      # Baseline covariates (fully observed)
      W1 = Oc$W1,
      W2 = Oc$W2,
      W3 = Oc$W3,
      
      # C1: censored when Delta=0 (Y1 not observed)
      C1 = BinaryToCensoring(is.censored = is.na(Oc$Y1)),
      
      # Y1: L-node — observed focus-pop indicator; fill NA with 0
      # (value for censored rows doesn't matter; det_Q + C1 handle them)
      Y1 = ifelse(is.na(Oc$Y1), 0L, Oc$Y1),
      
      # L: time-dependent confounder — only observed when Delta=1
      # NA when Delta=0 (consistent with DGP)
      L  = Oc$L,
      
      # C2: censored when Delta=1 AND Y1=1 AND Delta2=0
      # (the only rows where Y2 is truly missing)
      C2 = BinaryToCensoring(
        is.censored = !is.na(Oc$Y1) & Oc$Y1 == 1 & is.na(Oc$Y2)
      ),
      
      # Y2: set to 0 when Y1=0 (screened, not at risk — deterministic);
      #     NA when Delta2=0 and Y1=1 (true missingness);
      #     Y2_star when Delta2=1 and Y1=1
      Y2 = ifelse(!is.na(Oc$Y1) & Oc$Y1 == 0, 0, Oc$Y2)
    )
    
    # ---- Deterministic Q function ----
    # Y1=0 and Delta=1 → Y2=0 deterministically
    # ltmle needs to know not to model Y2 for these rows
    det_Q_fun_num <- function(data_det, current.node, nodes, called.from.estimate.g) {
      is.deterministic <- current.node == "Y2" & data_det$Y1 == 0
      list(is.deterministic = is.deterministic, Q.value = 0)
    }
    
    # ---- G-forms ----
    gform_num <- c(
      "C1 ~ W1 + W2 + W3",          # P(Delta=1 | W) — screening mechanism
      "C2 ~ Y1 + L + W1 + W2 + W3"  # P(Delta2=1 | Y1, L, W) — follow-up mechanism
      # Y1 included: Delta2 depends on observed Y1
      # L included: TDC into follow-up measurement
    )
    
    # ---- ltmle call ----
    numerator_tmle <- tryCatch({
      fit <- suppressMessages(suppressWarnings(
        ltmle(
          data                     = Oc_num,
          Anodes                   = NULL,
          abar                     = NULL,
          Cnodes                   = c("C1", "C2"),
          Lnodes                   = c("Y1", "L"),   # both are intermediate L-nodes
          Ynodes                   = "Y2",
          SL.library               = SL.library,
          deterministic.Q.function = det_Q_fun_num,
          variance.method          = "ic",
          stratify                 = TRUE,
          estimate.time            = FALSE
        )
      ))
      v <- unname(fit$estimates["tmle"])
      if (is.na(v) || !is.finite(v)) numerator_naive else v
    }, error = function(e) numerator_naive)
    
 ###########################################   ################################   
    
    # ---- TMLE denominator: P(Y1*=1) — unchanged from original ----
    Oc_den <- Oc[, c(indv_cov, "Delta", "Y1")]
    Oc_den$Y1 <- ifelse(is.na(Oc_den$Y1), 0, Oc_den$Y1)
    den_tmle <- tryCatch({
      fit <- suppressMessages(suppressWarnings(
        ltmle(data           = Oc_den,
              Anodes         = "Delta",
              Ynodes         = "Y1",
              abar           = 1,
              SL.library     = SL.library,
              estimate.time  = FALSE,
              stratify       = TRUE)
      ))
      v <- unname(fit$estimates["tmle"])
      if (is.na(v) || !is.finite(v) || v < 1e-6) den_unadj else v
    }, error = function(e) den_unadj)
    
    # Final ratio: cluster-level estimate of E[Y2_star | Y1*=1, c]
    tmle_val <- if (den_tmle > 1e-6) numerator_tmle / den_tmle else NA_real_
    
    temp <- cbind(
      id = j, U = 1, alpha = 1,
      Oc[1, c(clust_cov, "n")],
      W1 = mean(Oc$W1), W2 = mean(Oc$W2), W3 = mean(Oc$W3),
      A  = Oc[1, "A"],
      screened  = screened,
      eligible  = eligible,
      unadj     = numerator_naive / den_unadj,
      tmle      = tmle_val,
      den_tmle  = den_tmle,
      numerator = numerator_tmle,
      numerator_naive=numerator_naive,
      den_unadj=den_unadj
    )
    stage1_out[j, ] <- temp
  }
  
  colnames(stage1_out) <- colnames(temp)
  stage1_out
}




