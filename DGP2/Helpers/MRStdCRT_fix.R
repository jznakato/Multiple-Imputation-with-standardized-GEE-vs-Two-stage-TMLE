##MRStdCRT_point: Modified function to include outputs of treatment specific means.



MRStdCRT_point <- function(formula, data, cluster, trt, trtprob,
                           family = gaussian(link = "identity"),
                           corstr, method = "GLM", scale){
    
    tryCatch(stopifnot(is.data.frame(data)),
             error = function(e) stop("Error: The provided object is not of class 'data.frame'."))
    tryCatch(stopifnot(inherits(family, "family")),
             error = function(e) stop("Error: The provided family is not a 'family' object."))
    stopifnot(is.character(cluster), is.character(trt))
    
    if (method == "GLMM") {
      fam <- family$family; lnk <- family$link
      is_allowed <- (fam == "gaussian" && lnk == "identity") ||
        (fam == "binomial" && lnk == "logit") ||
        (fam == "poisson"  && lnk == "log")
      if (!is_allowed)
        stop("`family` for method = \"GLMM\" must be one of:\n",
             "  - gaussian(link = \"identity\")\n",
             "  - binomial(link = \"logit\")\n",
             "  - poisson(link = \"log\")\n",
             "You supplied: ", fam, "(link = \"", lnk, "\")", call. = FALSE)
    }
    
    outcome      <- all.vars(formula[[2]])
    rhs_terms    <- terms(formula)
    term_labels  <- attr(rhs_terms, "term.labels")
    all_vars_in_formula <- setdiff(all.vars(formula), outcome)
    rhs_syms     <- setdiff(all_vars_in_formula, trt)
    
    trt_pat    <- paste0("(^|:)", trt, "(:|$)")
    int_labels <- term_labels[grepl(trt_pat, term_labels)]
    if (length(int_labels) > 0) {
      partners_raw <- unique(unlist(strsplit(int_labels, ":", fixed = TRUE)))
      int_partners <- intersect(setdiff(partners_raw, trt), rhs_syms)
    } else {
      int_partners <- character(0)
    }
    
    columns_to_check <- c(outcome, rhs_syms, cluster, trt)
    missing_columns  <- setdiff(columns_to_check, colnames(data))
    if (length(missing_columns) > 0)
      stop(paste("Error: The column(s)", paste(missing_columns, collapse = ", "),
                 "do not exist in the data."))
    
    data1 <- data %>%
      dplyr::select(dplyr::all_of(c(cluster, trt, outcome, rhs_syms))) %>%
      dplyr::rename(cluster = !!cluster, A = !!trt, Y = !!outcome) %>%
      dplyr::arrange(cluster)
    data1$prob <- trtprob
    
    if (length(rhs_syms) > 0) {
      adj_cov_names <- paste0("X", seq_along(rhs_syms))
      map_old2new   <- stats::setNames(adj_cov_names, rhs_syms)
      names(data1)[4:(3 + length(rhs_syms))] <- adj_cov_names
    } else {
      adj_cov_names <- NULL
      map_old2new   <- character(0)
    }
    
    data1 <- data1 %>% dplyr::group_by(cluster) %>%
      dplyr::mutate(N = dplyr::n()) %>% dplyr::ungroup()
    
    data_clus <- as.data.frame(
      data1 %>% dplyr::group_by(cluster) %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), mean))
    )
    
    cluster_vars <- adj_cov_names[sapply(adj_cov_names, function(nm) {
      all(tapply(data1[[nm]], data1$cluster,
                 function(v) length(unique(v[!is.na(v)])) == 1))
    })]
    ind_cov_names <- setdiff(adj_cov_names, cluster_vars)
    
    if (length(ind_cov_names) > 0) {
      data1 <- data1 %>% dplyr::group_by(cluster) %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(ind_cov_names),
                                    ~mean(.x), .names = "{.col}b")) %>%
        dplyr::mutate(dplyr::across(dplyr::all_of(ind_cov_names),
                                    ~.x - get(paste0(cur_column(), "b")),
                                    .names = "{.col}c")) %>%
        dplyr::ungroup()
    }
    
    partners_new <- unname(map_old2new[int_partners])
    partners_new <- partners_new[!is.na(partners_new)]
    
    main_terms_i <- c("A",
                      if (length(ind_cov_names) > 0) paste0(ind_cov_names, "c"),
                      if (length(ind_cov_names) > 0) paste0(ind_cov_names, "b"),
                      if (length(cluster_vars)  > 0) cluster_vars)
    ia_terms_i <- character(0)
    if (length(partners_new) > 0) {
      ia_terms_i <- c(
        if (length(intersect(partners_new, ind_cov_names)) > 0)
          c(paste0("A:", intersect(partners_new, ind_cov_names), "c"),
            paste0("A:", intersect(partners_new, ind_cov_names), "b")),
        if (length(intersect(partners_new, cluster_vars)) > 0)
          paste0("A:", intersect(partners_new, cluster_vars))
      )
    }
    formulai <- stats::as.formula(paste("Y ~", paste(c(main_terms_i, ia_terms_i),
                                                     collapse = " + ")))
    
    main_terms_c <- c("A",
                      if (length(ind_cov_names) > 0) ind_cov_names,
                      if (length(cluster_vars)  > 0) cluster_vars)
    ia_terms_c <- if (length(partners_new) > 0) paste0("A:", partners_new) else character(0)
    formulac <- stats::as.formula(paste("Y ~", paste(c(main_terms_c, ia_terms_c),
                                                     collapse = " + ")))
    
    model <- switch(method,
                    GLM  = try(stats::glm(formulac, data = data_clus, family = family), silent = TRUE),
                    LMM  = try(nlme::lme(as.formula(formulai), random = ~1 | cluster,
                                         data = data1), silent = TRUE),
                    GEE  = try(geepack::geeglm(as.formula(formulai), id = cluster,
                                               data = data1, corstr = corstr, family = family),
                               silent = TRUE),
                    GLMM = try(lme4::glmer(paste(formulai, "+ (1 | cluster)"),
                                           data = data1, family = family), silent = TRUE),
                    stop("Invalid method specified.")
    )
    
    if (inherits(model, "try-error")) {
      warning("All-0 or all-1 clusters yield infinite logit; falling back to non-parametric.")
      eta <- data_clus %>% dplyr::select(cluster) %>%
        dplyr::mutate(eta1 = 0, eta0 = 0) %>%
        dplyr::group_by(cluster) %>% as.data.frame()
    } else if (method == "GLM") {
      eta <- data_clus %>%
        dplyr::mutate(eta1 = stats::predict(model, newdata = dplyr::mutate(data_clus, A = 1),
                                            type = "response"),
                      eta0 = stats::predict(model, newdata = dplyr::mutate(data_clus, A = 0),
                                            type = "response")) %>%
        dplyr::select(cluster, eta1, eta0) %>%
        dplyr::group_by(cluster) %>% as.data.frame()
    } else if (method == "LMM") {
      nd1 <- data1; nd1$A <- 1
      nd0 <- data1; nd0$A <- 0
      beta <- lme4::fixef(model)
      mm1  <- model.matrix(formulai, data = nd1)
      mm0  <- model.matrix(formulai, data = nd0)
      p1 <- as.numeric(mm1 %*% beta)
      p0 <- as.numeric(mm0 %*% beta)
      eta <- dplyr::tibble(cluster = data1$cluster, eta1 = p1, eta0 = p0) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarise(eta1 = mean(eta1), eta0 = mean(eta0), .groups = "drop")
    } else if (method == "GEE") {
      eta <- as.data.frame(
        dplyr::summarise(
          dplyr::group_by(
            dplyr::select(
              dplyr::mutate(data1,
                            eta1 = stats::predict(model, newdata = dplyr::mutate(data1, A = 1),
                                                  type = "response"),
                            eta0 = stats::predict(model, newdata = dplyr::mutate(data1, A = 0),
                                                  type = "response")),
              cluster, eta1, eta0),
            cluster),
          eta1 = mean(eta1), eta0 = mean(eta0), .groups = "drop"))
    } else if (method == "GLMM") {
      pred1 <- predict(model, newdata = transform(data1, A = 1), re.form = NA, type = "link")
      pred0 <- predict(model, newdata = transform(data1, A = 0), re.form = NA, type = "link")
      if (family$link == "logit") {
        Pi <- pi
        s  <- sqrt(3 * as.numeric(lme4::VarCorr(model)) / Pi^2 + 1)
        eta <- data1 %>%
          dplyr::mutate(eta1 = plogis(pred1 / s),
                        eta0 = plogis(pred0 / s)) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(eta1 = mean(eta1), eta0 = mean(eta0)) %>%
          as.data.frame()
      } else if (family$link == "log") {
        v <- as.numeric(lme4::VarCorr(model))
        eta <- data1 %>%
          dplyr::mutate(eta1 = exp(pred1 + v / 2),
                        eta0 = exp(pred0 + v / 2)) %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarise(eta1 = mean(eta1), eta0 = mean(eta0)) %>%
          as.data.frame()
      }
    }
    
    mu_C1 <- data_clus$A / data_clus$prob       * (data_clus$Y - eta$eta1) + eta$eta1
    mu_C0 <- (1 - data_clus$A) / (1 - data_clus$prob) * (data_clus$Y - eta$eta0) + eta$eta0
    
    m1c <- mean(mu_C1)
    m0c <- mean(mu_C0)
    m1i <- mean(data_clus$N * mu_C1) / mean(data_clus$N)
    m0i <- mean(data_clus$N * mu_C0) / mean(data_clus$N)
    
    cate <- switch(scale,
                   RD = m1c - m0c,
                   RR = m1c / m0c,
                   OR = (m1c / (1 - m1c)) / (m0c / (1 - m0c)),
                   stop("Invalid scale specified."))
    
    iate <- switch(scale,
                   RD = m1i - m0i,
                   RR = m1i / m0i,
                   OR = (m1i / (1 - m1i)) / (m0i / (1 - m0i)),
                   stop("Invalid scale specified."))
    
    test_NICS <- cate - iate
    
    list(data1, data_clus,
         c(cate, iate, test_NICS, m1c, m0c, m1i, m0i))
  }
  



MRStdCRT_fit <- function(formula, data, cluster, trt,
                         trtprob = rep(0.5, nrow(data)),
                         method,
                         family = gaussian(link = "identity"),
                         corstr, scale, jack = 1, alpha = 0.05) {
  
  scale <- match.arg(scale, c("RD", "RR", "OR"))
  
  # ---- trtprob handling: unchanged ----
  if (is.null(trtprob)) {
    suppressWarnings({
      df_prob <- dplyr::mutate(dplyr::ungroup(dplyr::mutate(
        dplyr::group_by(data, .data[[cluster]]),
        first_trt = dplyr::first(.data[[trt]]))),
        prob_A_1 = mean(first_trt == 1, na.rm = TRUE),
        prob_A_0 = mean(first_trt == 0, na.rm = TRUE),
        assigned_value = ifelse(.data[[trt]] == 1, prob_A_1, prob_A_0))
    })
    trtprob <- df_prob$assigned_value
  } else {
    n_rows <- nrow(data); clus <- data[[cluster]]; uniq_c <- unique(clus)
    if (length(trtprob) == length(uniq_c)) {
      trtprob <- if (!is.null(names(trtprob)))
        trtprob[as.character(clus)] else trtprob[match(clus, uniq_c)]
    } else if (length(trtprob) != n_rows) {
      stop("`trtprob` must have length nrow(data) or number of clusters.")
    }
    if (any(!is.finite(trtprob)) || any(trtprob <= 0 | trtprob >= 1, na.rm = TRUE))
      warning("Some entries of `trtprob` are not in (0,1).")
  }
  
  # ---- point estimation: pes is now length 7 ----
  temp      <- MRStdCRT_point(formula, data, cluster, trt,
                                         trtprob, family, corstr, method, scale)
  data_clus <- temp[[2]]
  m         <- nrow(data_clus)
  pes       <- temp[[3]]            # c(cate, iate, test_NICS, m1c, m0c, m1i, m0i)
  k         <- length(pes)
  
  # ---- leave-one-cluster-out jackknife ----
  cluster_names  <- unique(data_clus$cluster)
  point_est_jack <- matrix(NA_real_, nrow = m, ncol = k)
  for (i in cluster_names) {
    keep      <- data[[cluster]] != i
    data_jack <- droplevels(data[keep, , drop = FALSE])
    point_est_jack[which(cluster_names == i), ] <-
     MRStdCRT_point(formula, data_jack, cluster, trt,
                                trtprob[keep], family, corstr, method, scale)[[3]]
  }
  
  jackse <- switch(as.character(jack),
                   "1" = sqrt(apply(point_est_jack, 2, var) * (m - 1)^2 / m),
                   "2" = {
                     rep_point <- matrix(pes, nrow = m, ncol = k, byrow = TRUE)
                     sqrt(colSums((point_est_jack - rep_point)^2) * (m - 1) / m)
                   },
                   "3" = {
                     rep_point <- matrix(pes, nrow = m, ncol = k, byrow = TRUE)
                     sqrt(colSums((point_est_jack - rep_point)^2) * (m - 1) / m) * sqrt(m / (m - 1))
                   },
                   stop("Invalid jackknife type specified.")
  )
  
  # ---- assemble table; drop row 3 (test_NICS — used only for ICS test) ----
  thres <- qt(1 - alpha / 2, df = m - 1)
  tab   <- cbind(Estimate    = pes,
                 `Std. Error` = jackse,
                 `CI lower`   = pes - thres * jackse,
                 `CI upper`   = pes + thres * jackse)[-3, , drop = FALSE]
  
  rownames(tab) <- c("cATE", "iATE",
                     "mu1_cluster", "mu0_cluster",
                     "mu1_indiv",   "mu0_indiv")
  
  tab <- as.data.frame(tab, check.names = FALSE)
  tab[["t"]]       <- tab[["Estimate"]] / tab[["Std. Error"]]
  tab[["p-value"]] <- 2 * pt(abs(tab[["t"]]), df = m - 1, lower.tail = FALSE)
  
  # ---- ICS test on element 3 of the original pes (NOT the trimmed table) ----
  test_sta <- pes[3] / jackse[3]
  p_val    <- 2 * min(pt(test_sta, df = m - 1), 1 - pt(test_sta, df = m - 1))
  ics_test <- c(test_sta, p_val)
  
  fit_list <- list(estimate = tab, m = m, N = temp[[2]]$N,
                   family = family, model = method,
                   ics = ics_test, alpha = alpha, scale = scale)
  class(fit_list) <- "MRS_obj"
  fit_list
}

