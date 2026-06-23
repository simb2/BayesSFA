devtools::load_all()
library(MASS)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scoringRules)
source(here::here("tests", "sim_helpers.R"))
set.seed(20)

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

settings <- tidyr::crossing(N = N_TRUE, V = n_vars, q = Q_TRUE)
samples  <- purrr::pmap(settings, sim_data_general)

# ---- UGLT sparse model ------------------------------------------------------

uglt_fits <- purrr::map(samples, function(samp) {
  V <- nrow(samp$data)
  fitBSFA(
    y = samp$data, constraint = "UGLT", fixed = TRUE,
    q = Q_TRUE, n_runs = 7000,
    alpha = rep(1.5, V), beta = rep(1.5, V),
    theta.shape = 1.5, theta.rate = 1.5,
    hyperparams = list(aH = 2, bH = 2),
    thin = 2, burn = 1000
  )
})

uglt_ests        <- purrr::map(uglt_fits, get_sparse_est, q = Q_TRUE)
uglt_lambda_est  <- purrr::map(uglt_ests, "lambda_est")
uglt_sigma_est   <- purrr::map(uglt_ests, "sigma2_mean")
uglt_factors_est <- purrr::map(uglt_ests, "factors_est")

uglt_mse <- purrr::pmap_dbl(
  tibble::tibble(
    N = settings$N, Lambda = map(samples, "Lambda"), Sigma = map(samples, "Sigma"),
    factors = map(samples, "factors"),
    Lambda_est = uglt_lambda_est, sigma_est = uglt_sigma_est, factors_est = uglt_factors_est
  ),
  compute_mse
)
uglt_crps  <- purrr::map2_dbl(uglt_ests, map(samples, "Lambda"),
  ~ avg_crps_list(.x$lambda_draws, .y))
uglt_noise <- purrr::map2_dbl(uglt_lambda_est, map(samples, "data"),
  contributed_variance_noise)

# ---- Sparse PLT model -------------------------------------------------------

splt_fits <- purrr::map(samples, function(samp) {
  V <- nrow(samp$data)
  fitBSFA(
    y = samp$data, constraint = "PLT", fixed = TRUE,
    q = Q_TRUE, n_runs = 7000,
    alpha = rep(1.5, V), beta = rep(1.5, V),
    theta.shape = 1.5, theta.rate = 1.5,
    hyperparams = list(aH = 2, bH = 2),
    thin = 2, burn = 1000
  )
})

splt_ests        <- purrr::map(splt_fits, get_sparse_est, q = Q_TRUE)
splt_lambda_est  <- purrr::map(splt_ests, "lambda_est")
splt_sigma_est   <- purrr::map(splt_ests, "sigma2_mean")
splt_factors_est <- purrr::map(splt_ests, "factors_est")

splt_mse <- purrr::pmap_dbl(
  tibble::tibble(
    N = settings$N, Lambda = map(samples, "Lambda"), Sigma = map(samples, "Sigma"),
    factors = map(samples, "factors"),
    Lambda_est = splt_lambda_est, sigma_est = splt_sigma_est, factors_est = splt_factors_est
  ),
  compute_mse
)
splt_crps  <- purrr::map2_dbl(splt_ests, map(samples, "Lambda"),
  ~ avg_crps_list(.x$lambda_draws, .y))
splt_noise <- purrr::map2_dbl(splt_lambda_est, map(samples, "data"),
  contributed_variance_noise)

# ---- Results ----------------------------------------------------------------

results <- settings |> dplyr::mutate(
  uglt_mse  = uglt_mse,  uglt_crps = uglt_crps,  uglt_noise = uglt_noise,
  splt_mse  = splt_mse,  splt_crps = splt_crps,  splt_noise = splt_noise
)

print(results)
saveRDS(results, here::here("tests", "results_general_lambda.rds"))
