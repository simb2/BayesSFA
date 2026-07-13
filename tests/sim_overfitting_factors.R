library(Rcpp)
library(RcppArmadillo)
library(BayesSFA)
library(MASS)
library(purrr)
library(tidyr)
library(ggplot2)
library(furrr)
library(future)
source(here::here("tests", "sim_helpers.R"))
set.seed(8)

# Leave a couple of cores free; each worker holds its own copy of the fitted
# model's draws, so memory (not core count) is the binding constraint here.
future::plan(future::multisession, workers = min(6, future::availableCores() - 1))

N_TRUE <- 60
V_TRUE <- 180
Q_TRUE <- 10
Q_FIT  <- Q_TRUE + 3   # overfit by 3 extra factors
N_REP  <- 50

ALPHA <- rep(1.5, V_TRUE)
BETA <- rep(1.5, V_TRUE)
MCMC_ARGS <- list(
  n_runs = 10000,
  alpha = ALPHA,
  beta = BETA,
  theta.shape = 1.5,
  theta.rate  = 1.5,
  hyperparams = list(aH = 2, bH = 2),
  thin        = 5,
  burn        = 5000,
  fixed       = FALSE
)

# ---- One replicate: simulate both truths, fit both models on each -----------

modal_r <- function(draws) as.numeric(names(which.max(table(draws$r))))

run_replicate <- function(rep_id) {
  ensure_pkg_loaded()
  samp_uglt <- sim_data_UGLT(N_TRUE, V_TRUE, Q_TRUE, plt_structure = FALSE)
  samp_plt  <- sim_data_UGLT(N_TRUE, V_TRUE, Q_TRUE, plt_structure = TRUE)

  fit_uglt_on_uglt <- do.call(fitBSFA,
    c(list(y = samp_uglt$data, constraint = "UGLT", q = Q_FIT), MCMC_ARGS))
  fit_plt_on_uglt  <- do.call(fitBSFA,
    c(list(y = samp_uglt$data, constraint = "PLT",  q = Q_FIT), MCMC_ARGS))
  fit_uglt_on_plt  <- do.call(fitBSFA,
    c(list(y = samp_plt$data, constraint = "UGLT", q = Q_FIT), MCMC_ARGS))
  fit_plt_on_plt   <- do.call(fitBSFA,
    c(list(y = samp_plt$data, constraint = "PLT",  q = Q_FIT), MCMC_ARGS))

  if (rep_id == 1) {
    save_mcmc_diagnostics(fit_uglt_on_uglt$draws, "overfit_UGLTdata_UGLTmodel")
    save_mcmc_diagnostics(fit_plt_on_uglt$draws,  "overfit_UGLTdata_SparsePLTmodel")
    save_mcmc_diagnostics(fit_uglt_on_plt$draws,  "overfit_PLTdata_UGLTmodel")
    save_mcmc_diagnostics(fit_plt_on_plt$draws,   "overfit_PLTdata_SparsePLTmodel")
  }

  list(
    truth = list(uglt = samp_uglt, plt = samp_plt),
    r_summary = tibble::tibble(
      rep       = rep_id,
      truth     = rep(c("UGLT data", "PLT data"), each = 2),
      fit_model = rep(c("UGLT", "Sparse PLT"), times = 2),
      mean_r    = c(
        mean(fit_uglt_on_uglt$draws$r), mean(fit_plt_on_uglt$draws$r),
        mean(fit_uglt_on_plt$draws$r),  mean(fit_plt_on_plt$draws$r)
      ),
      modal_r   = c(
        modal_r(fit_uglt_on_uglt$draws), modal_r(fit_plt_on_uglt$draws),
        modal_r(fit_uglt_on_plt$draws),  modal_r(fit_plt_on_plt$draws)
      )
    )
  )
}

# ---- Run all replicates ------------------------------------------------------

replicates <- furrr::future_map(seq_len(N_REP), run_replicate,
                                 .options = furrr::furrr_options(seed = TRUE))

overfit_truth_samples <- purrr::map(replicates, "truth")
overfit_r_results     <- purrr::map_dfr(replicates, "r_summary")

saveRDS(overfit_truth_samples, here::here("tests", "sim_overfitting_factors_samples.rds"))
saveRDS(overfit_r_results,     here::here("tests", "sim_overfitting_factors_results.rds"))

print(overfit_r_results |> dplyr::group_by(truth, fit_model) |>
        dplyr::summarise(mean_r = mean(mean_r), .groups = "drop"))
