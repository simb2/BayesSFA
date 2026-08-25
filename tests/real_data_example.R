# Real-data comparison of the three models (Non-sparse PLT, Sparse PLT, UGLT)
# on the MOFA microbiome benchmark data (same source as
# BayesianWorkflow/Thesis_Code/real_data_application.qmd).
#
# Unlike the simulation studies in this directory, there is no ground-truth
# Lambda/Sigma/factors to compare against, so "Goodness of Fit" and "CRPS"
# here are out-of-sample predictive checks instead of parameter-recovery
# checks: for each of N_REP random train/test splits (by subject/column), we
# fit on the training subjects, draw fresh factors from the model's N(0, I)
# prior for the held-out subjects (their true factors are never observed),
# and score the resulting posterior-predictive replicates against the actual
# held-out data. "Noise Ratio" (tr(Lambda Lambda^T) / tr(cov(y))) needs no
# truth and is computed directly on the training fit.
#
# Runtime: this fits 3 models x N_REP splits, each an MCMC run on the full
# V x N_train microbiome matrix. Start with small N_REP/n_runs to sanity
# check before scaling up.

library(Rcpp)
library(RcppArmadillo)
library(BayesSFA)
library(data.table)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scoringRules)
library(furrr)
library(future)
source(here::here("tests", "sim_helpers.R"))
set.seed(1)

future::plan(future::multisession, workers = min(6, future::availableCores() - 1))

# ---- Load & prepare data -----------------------------------------------------
# Cached locally after the first run since the source is a remote FTP file.

data_cache <- here::here("tests", "microbiome_bacteria_data.rds")
if (file.exists(data_cache)) {
  dt_bact <- readRDS(data_cache)
} else {
  dt <- data.table::fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/microbiome/data.txt.gz")
  dt_bact <- dt[dt$view == "Bacteria", ] |>
    tidyr::pivot_wider(names_from = sample, values_from = "value") |>
    dplyr::select(-view)
  saveRDS(dt_bact, data_cache)
}

y_mat <- as.matrix(dplyr::select(dt_bact, -feature))

# The samplers assume a fully-observed V x N matrix; drop any taxa with gaps
# left by pivot_wider() rather than imputing them.
complete_rows <- rowSums(is.na(y_mat)) == 0
if (!all(complete_rows)) {
  message(sum(!complete_rows), " of ", nrow(y_mat),
          " features have missing values across samples; dropping them.")
  y_mat <- y_mat[complete_rows, ]
}

V <- nrow(y_mat)
N <- ncol(y_mat)
cat("Data matrix: V =", V, "features x N =", N, "samples\n")

Q_FIT <- 10

# ---- Full-data fits + diagnostics --------------------------------------------
# One fit per model on all N samples, for inspecting recovered factor
# structure directly (no held-out split, so no predictive metrics here).

y_full <- y_mat - rowMeans(y_mat)  # no intercept term in the model, so center

FULL_MCMC_ARGS <- list(
  n_runs = 14000, alpha = rep(1, V), beta = rep(1, V),
  theta.shape = 1.5, theta.rate = 1.5, hyperparams = list(aH = 2, bH = 2),
  thin = 7, burn = 2000, fixed = FALSE
)

fit_uglt_full <- do.call(fitBSFA, c(list(y = y_full, constraint = "UGLT", q = Q_FIT), FULL_MCMC_ARGS))
fit_splt_full <- do.call(fitBSFA, c(list(y = y_full, constraint = "PLT",  q = Q_FIT), FULL_MCMC_ARGS))

plt_full_start <- {
  sv       <- svd(y_full)
  Lambda_0 <- sv$u[, seq_len(Q_FIT)]
  W_0      <- t(sv$v[, seq_len(Q_FIT)]) * sv$d[seq_len(Q_FIT)]
  list(Lambda_0 = Lambda_0, sigma2_0 = pmax(diag(cov(t(y_full - Lambda_0 %*% W_0))), 1e-6))
}
fit_plt_full <- fitBayesPLT(
  Lambda_0 = plt_full_start$Lambda_0, sigma2_0 = plt_full_start$sigma2_0,
  n_runs = 14000, data = y_full, q = Q_FIT, nu = 3, s2 = 0.5, c_0 = 1,
  thin = 7, burn = 2000
)

save_mcmc_diagnostics(fit_uglt_full$draws, "real_data_UGLT")
save_mcmc_diagnostics(fit_splt_full$draws, "real_data_SparsePLT")
save_mcmc_diagnostics(list(T_stat = fit_plt_full$T_stat), "real_data_NonSparsePLT")

modal_r <- function(r_draws) as.numeric(names(which.max(table(r_draws))))
modal_r_uglt <- modal_r(fit_uglt_full$draws$r)
modal_r_splt <- modal_r(fit_splt_full$draws$r)

cat("Recovered factor count (r) on full data:\n")
cat("  UGLT:       mean =", mean(fit_uglt_full$draws$r), " modal =", modal_r_uglt, "\n")
cat("  Sparse PLT: mean =", mean(fit_splt_full$draws$r), " modal =", modal_r_splt, "\n")

noise_full <- c(
  "Non-sparse PLT" = contributed_variance_noise(apply(fit_plt_full$Lambda, c(2, 3), mean), y_full),
  "Sparse PLT"     = contributed_variance_noise(fit_splt_full$estimates[[as.character(modal_r_splt)]]$lambda_est, y_full),
  "UGLT"           = contributed_variance_noise(fit_uglt_full$estimates[[as.character(modal_r_uglt)]]$lambda_est, y_full)
)
cat("Noise ratio on full data:\n")
print(noise_full)

# ---- PSIS-LOO model comparison (entry-wise, single full-data fit) ------------
# See BayesianLOOCV-3.pdf SS3.6 (Vehtari et al. 2016): importance-weight each
# model's full-data posterior draws to approximate the leave-one-out
# predictive density per (feature, subject) entry, Pareto-smoothed via
# loo::loo() for stability. No refitting or held-out split needed here,
# unlike the split-CV comparison below.

log_lik_uglt <- pointwise_loglik_sparse(fit_uglt_full, y_full)
log_lik_splt <- pointwise_loglik_sparse(fit_splt_full, y_full)
log_lik_plt  <- pointwise_loglik_plt(fit_plt_full, y_full)

loo_uglt <- run_psis_loo(log_lik_uglt)
loo_splt <- run_psis_loo(log_lik_splt)
loo_plt  <- run_psis_loo(log_lik_plt)

loo_by_model <- list("Non-sparse PLT" = loo_plt, "Sparse PLT" = loo_splt, "UGLT" = loo_uglt)

cat("PSIS-LOO comparison (entry-wise elpd_loo, higher is better):\n")
print(loo::loo_compare(loo_by_model))

cat("Pareto k diagnostics (share of entries where the PSIS-LOO estimate is unreliable):\n")
print(purrr::map_dfr(loo_by_model, pareto_k_summary, .id = "model"))

saveRDS(loo_by_model, here::here("tests", "results_real_data_psis_loo.rds"))

# ---- Out-of-sample predictive helpers ----------------------------------------

# Held-out subjects' factors are unobserved (unlike in the simulation studies,
# where the generating factors are known), so replicate data is generated by
# drawing fresh factors from the N(0, I) prior.
predictive_mse <- function(y_test, lambda_est, sigma2_est, n_reps = 20) {
  q      <- ncol(lambda_est)
  n_test <- ncol(y_test)
  mean(replicate(n_reps, {
    factors_new <- matrix(rnorm(q * n_test), nrow = q, ncol = n_test)
    y_rep       <- sim_data_3(n_test, lambda_est, diag(sigma2_est), factors_new)
    mean((y_test - y_rep)^2)
  }))
}

# lambda_draws/sigma2_draws: lists of equal length (fitBSFA's fit$draws$Lambda_test
# / fit$draws$sigma_test). One posterior-predictive replicate per draw, using
# freshly-sampled factors for the held-out subjects.
predictive_crps_list <- function(lambda_draws, sigma2_draws, y_test) {
  n_test  <- ncol(y_test)
  n_draws <- length(lambda_draws)
  Vt      <- nrow(y_test)
  pred_samps <- array(NA_real_, dim = c(Vt, n_test, n_draws))
  for (m in seq_len(n_draws)) {
    q_m <- ncol(lambda_draws[[m]])
    factors_new <- matrix(rnorm(q_m * n_test), nrow = q_m, ncol = n_test)
    pred_samps[, , m] <- sim_data_3(n_test, lambda_draws[[m]], diag(sigma2_draws[[m]]), factors_new)
  }
  total <- 0
  for (i in seq_len(Vt))
    for (j in seq_len(n_test))
      total <- total + scoringRules::crps_sample(y_test[i, j], dat = pred_samps[i, j, ])
  total / (Vt * n_test)
}

# lambda_draws/sigma2_draws: arrays (fitBayesPLT's fit$Lambda / fit$Variances).
predictive_crps_array <- function(lambda_draws, sigma2_draws, y_test) {
  n_test  <- ncol(y_test)
  n_draws <- dim(lambda_draws)[1]
  Vt      <- dim(lambda_draws)[2]
  q       <- dim(lambda_draws)[3]
  pred_samps <- array(NA_real_, dim = c(Vt, n_test, n_draws))
  for (m in seq_len(n_draws)) {
    factors_new <- matrix(rnorm(q * n_test), nrow = q, ncol = n_test)
    pred_samps[, , m] <- sim_data_3(n_test, lambda_draws[m, , ], diag(sigma2_draws[m, ]), factors_new)
  }
  total <- 0
  for (i in seq_len(Vt))
    for (j in seq_len(n_test))
      total <- total + scoringRules::crps_sample(y_test[i, j], dat = pred_samps[i, j, ])
  total / (Vt * n_test)
}

# Split by subject (column); center by the training subjects' row means only,
# so no information from the held-out subjects leaks into the fit.
make_split <- function(y, prop_test = 0.2) {
  n        <- ncol(y)
  test_idx <- sample.int(n, size = round(prop_test * n))
  train    <- y[, -test_idx, drop = FALSE]
  test     <- y[,  test_idx, drop = FALSE]
  row_mean <- rowMeans(train)
  list(train = train - row_mean, test = test - row_mean)
}

# ---- Cross-validated comparison ----------------------------------------------

N_REP     <- 10
PROP_TEST <- 0.2

CV_MCMC_ARGS <- list(
  n_runs = 6000, alpha = rep(1.5, V), beta = rep(1.5, V),
  theta.shape = 1.5, theta.rate = 1.5, hyperparams = list(aH = 2, bH = 2),
  thin = 3, burn = 2000, fixed = TRUE
)
CV_PLT_ARGS <- list(n_runs = 6000, nu = 3, s2 = 0.5, c_0 = 1, thin = 3, burn = 2000)

splits <- purrr::map(seq_len(N_REP), ~ make_split(y_mat, PROP_TEST))

fit_and_score_sparse <- function(split, constraint) {
  ensure_pkg_loaded()
  fit <- do.call(fitBSFA, c(list(y = split$train, constraint = constraint, q = Q_FIT), CV_MCMC_ARGS))
  est <- get_sparse_est(fit, q = Q_FIT)
  list(
    mse   = predictive_mse(split$test, est$lambda_est, est$sigma2_mean),
    crps  = predictive_crps_list(fit$draws$Lambda_test, fit$draws$sigma_test, split$test),
    noise = contributed_variance_noise(est$lambda_est, split$train)
  )
}

fit_and_score_plt <- function(split) {
  ensure_pkg_loaded()
  sv       <- svd(split$train)
  Lambda_0 <- sv$u[, seq_len(Q_FIT)]
  W_0      <- t(sv$v[, seq_len(Q_FIT)]) * sv$d[seq_len(Q_FIT)]
  sigma2_0 <- pmax(diag(cov(t(split$train - Lambda_0 %*% W_0))), 1e-6)

  fit <- do.call(fitBayesPLT, c(
    list(Lambda_0 = Lambda_0, sigma2_0 = sigma2_0, data = split$train, q = Q_FIT),
    CV_PLT_ARGS
  ))
  lambda_est <- apply(fit$Lambda,    c(2, 3), mean)
  sigma_est  <- apply(fit$Variances, 2,       mean)
  list(
    mse   = predictive_mse(split$test, lambda_est, sigma_est),
    crps  = predictive_crps_array(fit$Lambda, fit$Variances, split$test),
    noise = contributed_variance_noise(lambda_est, split$train)
  )
}

plt_results <- furrr::future_map(splits, fit_and_score_plt,
                                  .options = furrr::furrr_options(seed = TRUE))
plt_mse   <- purrr::map_dbl(plt_results, "mse")
plt_crps  <- purrr::map_dbl(plt_results, "crps")
plt_noise <- purrr::map_dbl(plt_results, "noise")

uglt_results <- furrr::future_map(splits, fit_and_score_sparse, constraint = "UGLT",
                                   .options = furrr::furrr_options(seed = TRUE))
uglt_mse   <- purrr::map_dbl(uglt_results, "mse")
uglt_crps  <- purrr::map_dbl(uglt_results, "crps")
uglt_noise <- purrr::map_dbl(uglt_results, "noise")

splt_results <- furrr::future_map(splits, fit_and_score_sparse, constraint = "PLT",
                                   .options = furrr::furrr_options(seed = TRUE))
splt_mse   <- purrr::map_dbl(splt_results, "mse")
splt_crps  <- purrr::map_dbl(splt_results, "crps")
splt_noise <- purrr::map_dbl(splt_results, "noise")

results <- tibble::tibble(
  rep       = seq_len(N_REP),
  plt_mse   = plt_mse,   plt_crps  = plt_crps,   plt_noise  = plt_noise,
  uglt_mse  = uglt_mse,  uglt_crps = uglt_crps,  uglt_noise = uglt_noise,
  splt_mse  = splt_mse,  splt_crps = splt_crps,  splt_noise = splt_noise
)
print(results)
saveRDS(results, here::here("tests", "results_real_data_example.rds"))

# ---- Comparison boxplot -------------------------------------------------------

pivot_real_data_results <- function(results) {
  results |>
    tidyr::pivot_longer(
      cols          = c(plt_mse, plt_crps, plt_noise,
                        uglt_mse, uglt_crps, uglt_noise,
                        splt_mse, splt_crps, splt_noise),
      names_to      = c("model", "metric"),
      names_pattern = "^(plt|uglt|splt)_(.*)$",
      values_to     = "value"
    ) |>
    dplyr::mutate(
      model = dplyr::case_when(
        model == "plt"  ~ "Non-sparse PLT",
        model == "splt" ~ "Sparse PLT",
        model == "uglt" ~ "UGLT"
      ),
      metric = dplyr::case_when(
        metric == "mse"   ~ "Predictive MSE",
        metric == "crps"  ~ "Predictive CRPS",
        metric == "noise" ~ "Noise Ratio"
      ),
      model  = factor(model,  levels = names(MODEL_COLORS)),
      metric = factor(metric, levels = c("Predictive CRPS", "Predictive MSE", "Noise Ratio"))
    )
}

plot_real_data_comparison <- function(results) {
  long <- pivot_real_data_results(results)

  ggplot(long, aes(x = model, y = value, fill = model)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = MODEL_COLORS, guide = "none") +
    labs(x = NULL, y = NULL,
         title = "Out-of-sample model comparison on the microbiome data",
         subtitle = paste0(N_REP, " random ", round(PROP_TEST * 100), "% held-out splits")) +
    bsfa_boxplot_theme()
}

p <- plot_real_data_comparison(results)
p
save_bsfa_plot(p, "real_data_comparison_boxplot.png", width = 9, height = 6)
