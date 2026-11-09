# Code by Joy Nakato & Laura B. Balzer
# jznakato@berkeley.edu & laura.balzer@berkeley.edu

##=======================================================================================================

# looper: function to implement all estimators
## Inputs:
# dgp: Data generating process
# J: Number of clusters
# N_mean: cluster size mean
# N_sd: cluster size sd
# SL.library: SuperLearner library, default (NULL) gives GLM
# verbose: indicator to print update
## Output: A list of estimates from all estimators with inference

looper<- function(dgp, J, N_mean, N_sd, psi=NULL, scale,level,
                  SL.library=NULL, verbose=F, m=5,maxit = 1, fcs = FALSE,b0 , bA ,passive=TRUE){
  

  maxit <- if (fcs) max(maxit, 10) else maxit
  # full data including counterfactuals
  
  O <- gen_data_until_valid(dgp=dgp, J=J, N_mean=N_mean, N_sd=N_sd, verbose=F,b0 = b0, bA = bA)
  
  indv_cov<-c("W1","W2","W3")  
  clust_cov<-c("E1","E2")
  
  screened <- eligible <- unadj <- tmle <- tmle_tmle <-  
    gee <- stmle_eligible <-   NULL
  
  #---------------------
  # TWO-STAGE ESTIMATORS
  #---------------------
  
  # do stage1
  t_tmle<- system.time({Yc <- getYc(O=O, indv_cov=indv_cov, clust_cov=clust_cov, SL.library=SL.library)
  
  
  cand_cov <- c('U','E1','E2')
  data.input <- cbind(Yc[,c('id','n','alpha','A','den_tmle', cand_cov)], Y=Yc$tmle)
  ## Add a line of getting appropriate weights
  
  # bookkeeping for the means' degrees of freedom
  n_clust <- length(unique(data.input$id))
  df_t    <- n_clust - 2
  
  
  
  ## INDIVIDUAL LEVEL EFFECTS
  data.input<- get.weights2(data=data.input, weighting='indv')
  
 Stage2_RD_ind<- Stage2(goal="RD", target='indv', data.input=data.input,
                            do.data.adapt=T, do.cv.variance=T,
                            cand.QAdj=cand_cov, cand.Qform='glm',
                            cand.gAdj=cand_cov, cand.gform='glm',
                            verbose=F, psi=psi[["RD_indiv"]])

  
Stage2_RR_ind<- Stage2(goal="aRR", target='indv', data.input=data.input,
                            do.data.adapt=T, do.cv.variance=T,
                            cand.QAdj=cand_cov, cand.Qform='glm',
                            cand.gAdj=cand_cov, cand.gform='glm',
                            verbose=F, psi=psi[["RR_indiv"]])
  
  tmle_tmle_RD_ind <- Stage2_RD_ind[,c("CV.est","CV.CI.lo","CV.CI.hi","CV.se", "CV.pval","CV.bias","CV.cover","CV.reject")]
  #tmle_tmle_RD_ind <- tmle_tmle_RD_ind[,c("est","CI.lo","CI.hi","se", "pval","bias","cover","reject")]
  tmle_tmle_RD_ind <-   tmle_tmle_RD_ind %>% dplyr::mutate(method = "tmle_tmle_RD_iATE", .before = 1)
  colnames(tmle_tmle_RD_ind) <- c("method",'est', 'CI.lo', 'CI.hi', 'se', 'pval',"bias","cover","reject")
  
  
  tmle_tmle_RR_ind <- Stage2_RR_ind[,c("CV.est","CV.CI.lo","CV.CI.hi","CV.se", "CV.pval","CV.bias","CV.cover","CV.reject")]
  
  # tmle_tmle_RR_ind <- tmle_tmle_RR_ind[,c("est","CI.lo","CI.hi","se", "pval","bias","cover","reject")]
  tmle_tmle_RR_ind <-   tmle_tmle_RR_ind %>% dplyr::mutate(method = "tmle_tmle_RR_iATE", .before = 1)
  colnames(tmle_tmle_RR_ind) <- c("method",'est', 'CI.lo', 'CI.hi', 'se', 'pval',"bias","cover","reject")
  
  # --- mean rows (standard / non-CV inference, pulled from RD run) ---
  tmle_tmle_mu1_iATE <- extract_tmle_mean( Stage2_RD_ind, "Txt",
                                          true_mean = psi[["mu1_indiv"]],
                                          df_for_t  = df_t) %>%
    dplyr::mutate(method = "tmle_tmle_mu1_iATE", .before = 1)
  
  tmle_tmle_mu0_iATE <- extract_tmle_mean(Stage2_RD_ind, "Con",
                                          true_mean = psi[["mu0_indiv"]],
                                          df_for_t  = df_t) %>%
    dplyr::mutate(method = "tmle_tmle_mu0_iATE", .before = 1)
  
  
  
  ## CLUSTER LEVEL EFFECTS
  data.input<- get.weights2(data=data.input, weighting='clust')
  
  stage2_RD_clust<- Stage2(goal="RD", target='clust', data.input=data.input,
                              do.data.adapt=T, do.cv.variance=T,
                              cand.QAdj=cand_cov, cand.Qform='glm',
                              cand.gAdj=cand_cov, cand.gform='glm',
                              verbose=F, psi= psi[["RD_clust"]])
  stage2_RR_clust<- Stage2(goal="aRR", target='clust', data.input=data.input,
                              do.data.adapt=T, do.cv.variance=T,
                              cand.QAdj=cand_cov, cand.Qform='glm',
                              cand.gAdj=cand_cov, cand.gform='glm',
                              verbose=F, psi=psi[["RR_clust"]])
  
  tmle_tmle_RD_clust <- stage2_RD_clust[,c("CV.est","CV.CI.lo","CV.CI.hi","CV.se", "CV.pval","CV.bias","CV.cover","CV.reject")]
  #tmle_tmle_RD_clust <- tmle_tmle_RD_clust[,c("est","CI.lo","CI.hi","se", "pval","bias","cover","reject")] 
  tmle_tmle_RD_clust <-   tmle_tmle_RD_clust %>% dplyr::mutate(method = "tmle_tmle_RD_cATE", .before = 1)
  colnames(tmle_tmle_RD_clust) <- c("method",'est', 'CI.lo', 'CI.hi', 'se', 'pval',"bias","cover","reject")
  
  
  tmle_tmle_RR_clust <- stage2_RR_clust[,c("CV.est","CV.CI.lo","CV.CI.hi","CV.se", "CV.pval","CV.bias","CV.cover","CV.reject")]
  # tmle_tmle_RR_clust <- tmle_tmle_RR_clust[,c("est","CI.lo","CI.hi","se", "pval","bias","cover","reject")]
  tmle_tmle_RR_clust <-   tmle_tmle_RR_clust %>% dplyr::mutate(method = "tmle_tmle_RR_cATE", .before = 1)
  colnames(tmle_tmle_RR_clust) <- c("method",'est', 'CI.lo', 'CI.hi', 'se', 'pval',"bias","cover","reject")
  
  
  # --- mean rows ---
  tmle_tmle_mu1_cATE <- extract_tmle_mean( stage2_RD_clust, "Txt",
                                          true_mean = psi[["mu1_clust"]],
                                          df_for_t  = df_t) %>%
    dplyr::mutate(method = "tmle_tmle_mu1_cATE", .before = 1)
  
  tmle_tmle_mu0_cATE <- extract_tmle_mean(stage2_RD_clust, "Con",
                                          true_mean = psi[["mu0_clust"]],
                                          df_for_t  = df_t) %>%
    dplyr::mutate(method = "tmle_tmle_mu0_cATE", .before = 1)}
  
)
  
  
  
  #---------------------
  # STANDARDIZED GEE FOR M IMPUTED DATASETS
  #---------------------
  
  t_fcs<-system.time({  
  # Create missing Y
  Obs1 <- O
 
  if (passive) {
    mi_spec <- make_mi_spec_passive(Obs1, fcs = fcs)
    vs <-  c("Y1", "L", "A_L", "Y2_star", "Y2")
    imp_datasets <- do.imputation.passive(data = mi_spec$dat, mi_spec = mi_spec,
                                          m = m, maxit = maxit, visitSequence = vs)
  } else {
    mi_spec <- make_mi_spec(Obs1, fcs = fcs, interaction = "jav")
    vs <- c("Y1", "L", "A_L", "Y2") 
    imp_datasets <- do.imputation(data = mi_spec$dat, mi_spec = mi_spec,
                                  m = m, maxit = maxit, visitSequence = vs)
  }
  
  imp_datasets <- lapply(imp_datasets, function(d) dplyr::filter(d, Y1 == 1))
  
  

  # one call per level, returns RD/RR/mu1/mu0 in one wide row
  res_imp_std_indiv <- dplyr::bind_rows(
    lapply(imp_datasets, run_gee_std, psi = psi, clust.lvl = FALSE), .id = "imp"
  )
  res_imp_std_clust <- dplyr::bind_rows(
    lapply(imp_datasets, run_gee_std, psi = psi, clust.lvl = TRUE),  .id = "imp"
  )
  
  
  pool_block <- function(res_imp, prefix, true_est, method_label) {
    cols <- grep(paste0("^", prefix, "_"), names(res_imp), value = TRUE)
    block <- res_imp[, c("imp", cols)]
    names(block) <- sub(paste0("^", prefix, "_"), "", names(block))   # strip prefix
    pool_rubin(block, true_est = true_est) %>%
      dplyr::transmute(method = method_label,
                       est, se, CI.lo, CI.hi, pval, bias, reject, cover)
  }
  
  # individual-level estimands
  gee_std_rd_ind  <- pool_block(res_imp_std_indiv, "RD",  psi[["RD_indiv"]],  "StdGEE_RD_iATE")
  gee_std_rr_ind  <- pool_block(res_imp_std_indiv, "RR",  psi[["RR_indiv"]],  "StdGEE_RR_iATE")
  gee_std_mu1_ind <- pool_block(res_imp_std_indiv, "mu1", psi[["mu1_indiv"]], "StdGEE_mu1_iATE")
  gee_std_mu0_ind <- pool_block(res_imp_std_indiv, "mu0", psi[["mu0_indiv"]], "StdGEE_mu0_iATE")
  
  # cluster-level estimands
  gee_std_rd_clust  <- pool_block(res_imp_std_clust, "RD",  psi[["RD_clust"]],  "StdGEE_RD_cATE")
  gee_std_rr_clust  <- pool_block(res_imp_std_clust, "RR",  psi[["RR_clust"]],  "StdGEE_RR_cATE")
  gee_std_mu1_clust <- pool_block(res_imp_std_clust, "mu1", psi[["mu1_clust"]], "StdGEE_mu1_cATE")
  gee_std_mu0_clust <- pool_block(res_imp_std_clust, "mu0", psi[["mu0_clust"]], "StdGEE_mu0_cATE")})
  
  # Combine everything
  
  res <- dplyr::bind_rows(
    # GEE side
    gee_std_rd_ind,    gee_std_rr_ind,
    gee_std_mu1_ind,   gee_std_mu0_ind,
    gee_std_rd_clust,  gee_std_rr_clust,
    gee_std_mu1_clust, gee_std_mu0_clust,
    # TMLE side
    tmle_tmle_RD_ind,  tmle_tmle_RR_ind,
    tmle_tmle_mu1_iATE, tmle_tmle_mu0_iATE,
    tmle_tmle_RD_clust, tmle_tmle_RR_clust,
    tmle_tmle_mu1_cATE, tmle_tmle_mu0_cATE
  )
  
  timings <- data.frame(
    method  = c("2stage_TMLE", "MI_GEE"),
    elapsed = c(t_tmle["elapsed"],t_fcs["elapsed"])
  )
  return(list(res = res, timings = timings))
  
}

# Pulls Txt (mu1) or Con (mu0) from a Stage2 output and returns a one-row
# data.frame matching the column structure of the contrast rows.
extract_tmle_mean <- function(stage2_out, side = c("Txt", "Con"),
                              true_mean, df_for_t, use_cv = TRUE) {
  side   <- match.arg(side)
  prefix <- if (use_cv) paste0(side, ".CV") else side
  
  est   <- stage2_out[[paste0(prefix, ".est")]]
  CI.lo <- stage2_out[[paste0(prefix, ".CI.lo")]]
  CI.hi <- stage2_out[[paste0(prefix, ".CI.hi")]]
  se    <- stage2_out[[paste0(prefix, ".se")]]
  
  tstat <- est / se
  pval  <- 2 * pt(abs(tstat), df = df_for_t, lower.tail = FALSE)
  
  data.frame(
    est    = est,
    CI.lo  = CI.lo,
    CI.hi  = CI.hi,
    se     = se,
    pval   = pval,
    bias   = est - true_mean,
    cover  = as.numeric(CI.lo <= true_mean & true_mean <= CI.hi),
    reject = as.numeric(pval < 0.05)
  )
}

