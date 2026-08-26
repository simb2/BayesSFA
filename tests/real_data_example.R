# Bayesian LOO-CV comparison of the three models (Non-sparse PLT, Sparse PLT,
# UGLT) on the MOFA microbiome benchmark data (same source as
# BayesianWorkflow/Thesis_Code/real_data_application.qmd).
#
# Follows BayesianLOOCV-3.pdf SS3.6 (Vehtari, Mononen, Tolvanen, Sivula &
# Winther, 2016, JMLR): rather than refitting per held-out point, each
# model's full-data posterior draws are importance-weighted to approximate
# the leave-one-out predictive density per (feature, subject) entry, then
# Pareto-smoothed (loo::loo()) for stability -- no held-out split or refit
# needed, since the likelihood factorizes over entries given each draw's
# Lambda, F, sigma2. Applies here rather than the Laplace/EP-specific
# machinery in the paper because these models are fit by full MCMC, not a
# Gaussian/EP posterior approximation.
#
# Runtime: this fits 3 models once each on the full V x N microbiome matrix
# at n_runs = 14000. Reduce n_runs for a quick sanity check before scaling up.

library(Rcpp)
library(RcppArmadillo)
library(BayesSFA)
library(data.table)
library(purrr)
library(tidyr)
library(dplyr)
source(here::here("tests", "sim_helpers.R"))
set.seed(1)

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

# rowSums(is.na(y_mat)) - none of the entries are missing. 

V <- nrow(y_mat)
N <- ncol(y_mat)
cat("Data matrix: V =", V, "features x N =", N, "samples\n")

Q_FIT <- 10

# ---- Full-data fits -----------------------------------------------------------
# One fit per model on all N samples; these posterior draws are the sole
# input to the PSIS-LOO comparison below (no held-out split needed).

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

# ---- PSIS-LOO model comparison (entry-wise, single full-data fit) ------------

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

p_loo <- plot_elpd_comparison(loo_by_model, estimate = "elpd_loo",
                               title = "PSIS-LOO model comparison (microbiome data)")
save_bsfa_plot(p_loo, "real_data_psis_loo_elpd.png", width = 7, height = 5)

# ---- WAIC model comparison (cross-check against PSIS-LOO) -------------------
# Same log_lik matrices as above; see BayesianLOOCV-3.pdf SS3.8. Cheaper than
# PSIS-LOO (no importance weighting) but a Taylor-series approximation to LOO
# rather than an estimate of it, and has no per-point reliability diagnostic.

waic_uglt <- run_waic(log_lik_uglt)
waic_splt <- run_waic(log_lik_splt)
waic_plt  <- run_waic(log_lik_plt)

waic_by_model <- list("Non-sparse PLT" = waic_plt, "Sparse PLT" = waic_splt, "UGLT" = waic_uglt)

cat("WAIC comparison (elpd_waic, higher is better):\n")
print(loo::loo_compare(waic_by_model))

saveRDS(waic_by_model, here::here("tests", "results_real_data_waic.rds"))
