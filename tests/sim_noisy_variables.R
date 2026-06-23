devtools::load_all()
library(MASS)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scoringRules)
source(here::here("tests", "sim_helpers.R"))
set.seed(8)

# ---- Data simulation ---------------------------------------------------------

sim_data_L <- function(N, V, q) {
  Lambda <- matrix(0, V, q)
  for (i in seq_len(V))
    for (j in seq_len(min(i, q)))
      Lambda[i, j] <- if (i == j) abs(rnorm(1, sd = sqrt(0.5))) else rnorm(1, sd = sqrt(0.5))
  factors <- t(MASS::mvrnorm(N, rep(0, q), diag(q)))
  data    <- Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, V), diag(V)))
  list(data = data, factors = factors, Lambda = Lambda, Sigma = diag(V), N = N, V = V)
}

sim_data_N <- function(N, V, q) {
  n_signal <- floor(V * 0.6)
  Lambda   <- matrix(0, V, q)
  for (i in seq_len(n_signal))
    for (j in seq_len(min(i, q)))
      Lambda[i, j] <- if (i == j) abs(rnorm(1, sd = sqrt(0.5))) else rnorm(1, sd = sqrt(0.5))
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

avg_crps_array <- function(lambda_samps, true_val) {
  V <- nrow(lambda_samps[1, , ]); q <- ncol(lambda_samps[1, , ])
  total <- 0
  for (i in seq_len(V))
    for (j in seq_len(q))
      total <- total + scoringRules::crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
  total / (V * q)
}

compute_mse <- function(N, Lambda, Sigma, factors, Lambda_est, sigma_est, factors_est) {
  mean(replicate(20, {
    samp  <- sim_data_3(N, Lambda,     Sigma,           factors)
    samp2 <- sim_data_3(N, Lambda_est, diag(sigma_est), factors_est)
    mean((samp - samp2)^2)
  }))
}

# Extract modal-r estimates from a fitBSFA result (fixed = TRUE means r = q always)
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

n_subj    <- c(300, 400)
n_vars    <- c(40, 50)
n_factors <- 2

settings <- tidyr::crossing(N = n_subj, V = n_vars, q = n_factors, noisy = c(FALSE, TRUE))
samples  <- purrr::pmap(settings, function(N, V, q, noisy)
  if (noisy) sim_data_N(N, V, q) else sim_data_L(N, V, q))

# ---- Non-sparse PLT model ---------------------------------------------------

plt_starts <- purrr::map(samples, function(samp) {
  centered <- samp$data - rowMeans(samp$data)
  sv       <- svd(centered)
  Lambda_0 <- sv$u[, seq_len(n_factors)]
  W_0      <- t(sv$v[, seq_len(n_factors)]) * sv$d[seq_len(n_factors)]
  sigma2_0 <- pmax(diag(cov(t(centered - Lambda_0 %*% W_0))), 1e-6)
  list(Lambda_0 = Lambda_0, sigma2_0 = sigma2_0)
})

plt_fits <- purrr::map2(samples, plt_starts, function(samp, start) {
  fitBayesPLT(
    Lambda_0 = start$Lambda_0, sigma2_0 = start$sigma2_0,
    n_runs = 5000, data = samp$data,
    q = n_factors, nu = 3, s2 = 0.5, c_0 = 1,
    thin = 2, burn = 500
  )
})

plt_lambda_est   <- purrr::map(plt_fits, ~ apply(.x$Lambda,    c(2, 3), mean))
plt_sigma_est    <- purrr::map(plt_fits, ~ apply(.x$Variances, 2,       mean))
plt_factors_est  <- purrr::map(plt_fits, ~ apply(.x$Factors,   c(2, 3), mean))

plt_mse <- purrr::pmap_dbl(
  tibble::tibble(
    N = settings$N, Lambda = map(samples, "Lambda"), Sigma = map(samples, "Sigma"),
    factors = map(samples, "factors"),
    Lambda_est = plt_lambda_est, sigma_est = plt_sigma_est, factors_est = plt_factors_est
  ),
  compute_mse
)
plt_crps  <- purrr::map2_dbl(plt_fits, map(samples, "Lambda"),
  ~ avg_crps_array(.x$Lambda, .y))
plt_noise <- purrr::map2_dbl(plt_lambda_est, map(samples, "data"),
  contributed_variance_noise)

# ---- UGLT sparse model ------------------------------------------------------

uglt_fits <- purrr::map(samples, function(samp) {
  V <- nrow(samp$data)
  fitBSFA(
    y = samp$data, constraint = "UGLT", fixed = TRUE,
    q = n_factors, n_runs = 5000,
    alpha = rep(1.5, V), beta = rep(1.5, V),
    theta.shape = 1.5, theta.rate = 1.5,
    hyperparams = list(aH = 2, bH = 2),
    thin = 2, burn = 500
  )
})

uglt_ests        <- purrr::map(uglt_fits, get_sparse_est, q = n_factors)
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
    q = n_factors, n_runs = 5000,
    alpha = rep(1.5, V), beta = rep(1.5, V),
    theta.shape = 1.5, theta.rate = 1.5,
    hyperparams = list(aH = 2, bH = 2),
    thin = 2, burn = 500
  )
})

splt_ests        <- purrr::map(splt_fits, get_sparse_est, q = n_factors)
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
  plt_mse   = plt_mse,   plt_crps  = plt_crps,   plt_noise  = plt_noise,
  uglt_mse  = uglt_mse,  uglt_crps = uglt_crps,  uglt_noise = uglt_noise,
  splt_mse  = splt_mse,  splt_crps = splt_crps,  splt_noise = splt_noise
)

print(results)
saveRDS(results, here::here("tests", "results_noisy_variables.rds"))
