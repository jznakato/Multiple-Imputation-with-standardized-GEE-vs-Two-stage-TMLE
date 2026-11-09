rm(list = ls())
source("generate_data.R")

## ---- settings -------------------------------------------------------------
dgp     <- "main"
J       <- 100
N_mean  <- 200
N_sd    <- 10
b0      <- -1
bA      <- 1.5
R       <- 500
J_truth <- 5000
seed0   <- 423

 base_dir <- "/home/jznakato/MI_vs_TMLE_DGP2"
 scen     <- sprintf("dgp2_%s_J%d_b0%s_bA%s", dgp, J, b0, bA)
 data_dir <- file.path(base_dir, "sim_data", scen)
 dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
 ## ---- truth (once) ----------------------------------------------------------
 # set.seed(seed0)
 # truth <- get_truth(dgp = dgp, J = J_truth, N_mean = N_mean, N_sd = N_sd)
 # saveRDS(truth, file.path(data_dir, "truth.rds"))
 # print(truth)
# saveRDS(list(dgp = dgp, J = J, N_mean = N_mean, N_sd = N_sd, b0 = b0, bA = bA,
#              R = R, J_truth = J_truth, seed0 = seed0),
#         file.path(data_dir, "settings.rds"))


truth <- readRDS(file.path(data_dir, "truth.rds"))
print(truth)
## ---- generate R trials -----------------------------------------------------
manifest <- data.frame(rep = 1:R, seed = seed0 + 1:R, file = NA_character_, n_tries = NA_integer_)

for (r in 1:R) {
  set.seed(manifest$seed[r])                       # each trial reproducible on its own
  O <- gen_data_until_valid(dgp = dgp, J = J, N_mean = N_mean, N_sd = N_sd,
                            b0 = b0, bA = bA)
  f <- file.path(data_dir, sprintf("rep_%03d.rds", r))
  saveRDS(O, f)
  manifest$file[r]    <- f
  manifest$n_tries[r] <- attr(O, "n_tries")
  if (r %% 50 == 0) message("generated ", r, " / ", R)
}

write.csv(manifest, file.path(data_dir, "manifest.csv"), row.names = FALSE)
table(manifest$n_tries)        # how often regeneration was needed
