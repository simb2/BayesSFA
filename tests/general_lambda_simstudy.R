library(Rcpp)
library(RcppArmadillo)
library(BayesSFA)
library(MASS)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scoringRules)
library(furrr)
library(future)
source(here::here("tests", "sim_helpers.R"))
set.seed(20)

# Leave a couple of cores free; each worker holds its own copy of the fitted
# model's draws, so memory (not core count) is the binding constraint here.
future::plan(future::multisession, workers = min(6, future::availableCores() - 1))

# ---- Data simulation --------------------------------------------------------

sim_data_general <- function(N, V, q) {
  Lambda  <- matrix(rnorm(V * q), nrow = V, ncol = q)
  factors <- t(MASS::mvrnorm(N, rep(0, q), diag(q)))
  data    <- Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, V), diag(V)))
  list(data = data, factors = factors, Lambda = Lambda, Sigma = diag(V), N = N, V = V)
}

sim_data_3 <- function(N, Lambda, Sigma, factors) {
  Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, nrow(Lambda)), Sigma))
}

# ---- Goodness-of-fit helpers ------------------------------------------------

contributed_variance_noise <- function(L, y) {
  sum(diag(L %*% t(L))) / sum(diag(cov(t(y))))
}

avg_crps_list <- function(lambda_samps, true_val) {
  V <- nrow(true_val); q <- ncol(true_val); n <- length(lambda_samps)
  samps <- array(unlist(lambda_samps), dim = c(V, q, n))
  total <- 0
  for (i in seq_len(V))
    for (j in seq_len(q))
      total <- total + scoringRules::crps_sample(true_val[i, j], dat = samps[i, j, ])
  total / (V * q)
}

compute_mse <- function(N, Lambda, Sigma, factors, Lambda_est, sigma_est, factors_est) {
  mean(replicate(20, {
    samp  <- sim_data_3(N, Lambda,     Sigma,           factors)
    samp2 <- sim_data_3(N, Lambda_est, diag(sigma_est), factors_est)
    mean((samp - samp2)^2)
  }))
}

get_sparse_est <- function(fit, q) {
  est <- fit$estimates[[as.character(q)]]
  list(
    lambda_est   = est$lambda_est,
    sigma2_mean  = est$sigma2_mean,
    factors_est  = est$factors_est,
    lambda_draws = fit$draws$Lambda_test
  )
}

# ---- Simulation settings ----------------------------------------------------

N_TRUE <- 150
n_vars <- c(100, 150, 200)
Q_TRUE <- 10
N_REP  <- 50

settings <- tidyr::crossing(N = N_TRUE, V = n_vars, q = Q_TRUE, rep = seq_len(N_REP))
samples  <- purrr::pmap(settings, function(N, V, q, rep) sim_data_general(N, V, q))

saveRDS(samples, here::here("tests", "general_lambda_simstudy_samples.rds"))

# Fit one sparse model on one sample and immediately reduce it to the three
# scalar metrics. Keeping only these per row (instead of the full fit/draws
# for all ~150 rows at once) is what keeps memory bounded for this study.
fit_and_score <- function(samp, idx, constraint, label_suffix) {
  ensure_pkg_loaded()
  V <- nrow(samp$data)
  fit <- fitBSFA(
    y = samp$data, constraint = constraint, fixed = TRUE,
    q = Q_TRUE, n_runs = 10000,
    alpha = rep(1.5, V), beta = rep(1.5, V),
    theta.shape = 1.5, theta.rate = 1.5,
    hyperparams = list(aH = 2, bH = 2),
    thin = 5, burn = 5000
  )
  if (settings$rep[idx] == 1) {
    save_mcmc_diagnostics(fit$draws, paste0("general_lambda_V", V, "_", label_suffix))
  }
  est <- get_sparse_est(fit, q = Q_TRUE)
  list(
    mse   = compute_mse(settings$N[idx], samp$Lambda, samp$Sigma, samp$factors,
                         est$lambda_est, est$sigma2_mean, est$factors_est),
    crps  = avg_crps_list(est$lambda_draws, samp$Lambda),
    noise = contributed_variance_noise(est$lambda_est, samp$data)
  )
}

# ---- UGLT sparse model ------------------------------------------------------

uglt_results <- furrr::future_imap(samples, fit_and_score, constraint = "UGLT", label_suffix = "UGLT",
                                    .options = furrr::furrr_options(seed = TRUE))
uglt_mse   <- purrr::map_dbl(uglt_results, "mse")
uglt_crps  <- purrr::map_dbl(uglt_results, "crps")
uglt_noise <- purrr::map_dbl(uglt_results, "noise")

# ---- Sparse PLT model -------------------------------------------------------

splt_results <- furrr::future_imap(samples, fit_and_score, constraint = "PLT", label_suffix = "SparsePLT",
                                    .options = furrr::furrr_options(seed = TRUE))
splt_mse   <- purrr::map_dbl(splt_results, "mse")
splt_crps  <- purrr::map_dbl(splt_results, "crps")
splt_noise <- purrr::map_dbl(splt_results, "noise")

# ---- Results ----------------------------------------------------------------

results <- settings |> dplyr::mutate(
  uglt_mse  = uglt_mse,  uglt_crps = uglt_crps,  uglt_noise = uglt_noise,
  splt_mse  = splt_mse,  splt_crps = splt_crps,  splt_noise = splt_noise
)

print(results)
saveRDS(results, here::here("tests", "results_general_lambda.rds"))
