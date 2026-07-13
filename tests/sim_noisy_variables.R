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
set.seed(8)

# Leave a couple of cores free; each worker holds its own copy of the fitted
# model's draws, so memory (not core count) is the binding constraint here.
future::plan(future::multisession, workers = min(6, future::availableCores() - 1))

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

# Build a descriptive diagnostics label for row `idx` of `settings`.
noisy_run_label <- function(settings, idx, model_name) {
  noisy_tag <- if (settings$noisy[idx]) "Noisy" else "NonNoisy"
  paste0("noisy_N", settings$N[idx], "_V", settings$V[idx], "_", noisy_tag, "_", model_name)
}

# ---- Simulation settings ----------------------------------------------------

n_subj    <- c(100, 150)
n_vars    <- c(220, 320)
n_factors <- 8
N_REP     <- 50

settings <- tidyr::crossing(N = n_subj, V = n_vars, q = n_factors,
                             noisy = c(FALSE, TRUE), rep = seq_len(N_REP))
samples  <- purrr::pmap(settings, function(N, V, q, noisy, rep)
  if (noisy) sim_data_N(N, V, q) else sim_data_L(N, V, q))

saveRDS(samples, here::here("tests", "sim_noisy_variables_samples.rds"))

# ---- Non-sparse PLT model ---------------------------------------------------

plt_starts <- purrr::map(samples, function(samp) {
  centered <- samp$data - rowMeans(samp$data)
  sv       <- svd(centered)
  Lambda_0 <- sv$u[, seq_len(n_factors)]
  W_0      <- t(sv$v[, seq_len(n_factors)]) * sv$d[seq_len(n_factors)]
  sigma2_0 <- pmax(diag(cov(t(centered - Lambda_0 %*% W_0))), 1e-6)
  list(Lambda_0 = Lambda_0, sigma2_0 = sigma2_0)
})

plt_results <- furrr::future_pmap(list(samples, plt_starts, seq_along(samples)), function(samp, start, idx) {
  ensure_pkg_loaded()
  fit <- fitBayesPLT(
    Lambda_0 = start$Lambda_0, sigma2_0 = start$sigma2_0,
    n_runs = 7000, data = samp$data,
    q = n_factors, nu = 3, s2 = 0.5, c_0 = 1,
    thin = 2, burn = 2000
  )
  if (settings$rep[idx] == 1) {
    save_mcmc_diagnostics(fit, noisy_run_label(settings, idx, "NonSparsePLT"))
  }
  lambda_est  <- apply(fit$Lambda,    c(2, 3), mean)
  sigma_est   <- apply(fit$Variances, 2,       mean)
  factors_est <- apply(fit$Factors,   c(2, 3), mean)
  list(
    mse   = compute_mse(settings$N[idx], samp$Lambda, samp$Sigma, samp$factors,
                         lambda_est, sigma_est, factors_est),
    crps  = avg_crps_array(fit$Lambda, samp$Lambda),
    noise = contributed_variance_noise(lambda_est, samp$data)
  )
}, .options = furrr::furrr_options(seed = TRUE))
plt_mse   <- purrr::map_dbl(plt_results, "mse")
plt_crps  <- purrr::map_dbl(plt_results, "crps")
plt_noise <- purrr::map_dbl(plt_results, "noise")

# Fit one sparse model on one sample and immediately reduce it to the three
# scalar metrics. Keeping only these per row (instead of the full fit/draws
# for all ~400 rows at once) is what keeps memory bounded for this study.
fit_and_score <- function(samp, idx, constraint, label_suffix) {
  ensure_pkg_loaded()
  V <- nrow(samp$data)
  fit <- fitBSFA(
    y = samp$data, constraint = constraint, fixed = TRUE,
    q = n_factors, n_runs = 7000,
    alpha = rep(1.5, V), beta = rep(1.5, V),
    theta.shape = 1.5, theta.rate = 1.5,
    hyperparams = list(aH = 2, bH = 2),
    thin = 2, burn = 2000
  )
  if (settings$rep[idx] == 1) {
    save_mcmc_diagnostics(fit$draws, noisy_run_label(settings, idx, label_suffix))
  }
  est <- get_sparse_est(fit, q = n_factors)
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
  plt_mse   = plt_mse,   plt_crps  = plt_crps,   plt_noise  = plt_noise,
  uglt_mse  = uglt_mse,  uglt_crps = uglt_crps,  uglt_noise = uglt_noise,
  splt_mse  = splt_mse,  splt_crps = splt_crps,  splt_noise = splt_noise
)

print(results)
saveRDS(results, here::here("tests", "results_noisy_variables.rds"))
