devtools::load_all()
library(MASS)
library(purrr)
library(tidyr)
library(ggplot2)
source(here::here("tests", "sim_helpers.R"))
set.seed(8)

N_TRUE <- 2000
V_TRUE <- 30
Q_TRUE <- 4
Q_FIT  <- Q_TRUE + 3   # overfit by 3 extra factors

ALPHA  <- rep(1.5, V_TRUE)
BETA   <- rep(1.5, V_TRUE)
MCMC_ARGS <- list(
  n_runs      = 6000,
  alpha       = ALPHA,
  beta        = BETA,
  theta.shape = 1.5,
  theta.rate  = 1.5,
  hyperparams = list(aH = 2, bH = 2),
  thin        = 2,
  burn        = 1000,
  fixed       = FALSE
)

# ---- UGLT-structured truth --------------------------------------------------

samp_uglt <- sim_data_UGLT(N_TRUE, V_TRUE, Q_TRUE, plt_structure = FALSE)

fit_uglt_on_uglt <- do.call(fitBSFA,
  c(list(y = samp_uglt$data, constraint = "UGLT", q = Q_FIT), MCMC_ARGS))
fit_plt_on_uglt  <- do.call(fitBSFA,
  c(list(y = samp_uglt$data, constraint = "PLT",  q = Q_FIT), MCMC_ARGS))

saveRDS(fit_uglt_on_uglt$draws, here::here("tests", "overfit_uglt_struct_draws_UGLT.rds"))
saveRDS(fit_plt_on_uglt$draws,  here::here("tests", "overfit_uglt_struct_draws_PLT.rds"))
saveRDS(samp_uglt,              here::here("tests", "overfit_uglt_struct_truth.rds"))

# ---- PLT-structured truth ---------------------------------------------------

samp_plt <- sim_data_UGLT(N_TRUE, V_TRUE, Q_TRUE, plt_structure = TRUE)

fit_uglt_on_plt  <- do.call(fitBSFA,
  c(list(y = samp_plt$data, constraint = "UGLT", q = Q_FIT), MCMC_ARGS))
fit_plt_on_plt   <- do.call(fitBSFA,
  c(list(y = samp_plt$data, constraint = "PLT",  q = Q_FIT), MCMC_ARGS))

saveRDS(fit_uglt_on_plt$draws, here::here("tests", "overfit_plt_struct_draws_UGLT.rds"))
saveRDS(fit_plt_on_plt$draws,  here::here("tests", "overfit_plt_struct_draws_PLT.rds"))
saveRDS(samp_plt,              here::here("tests", "overfit_plt_struct_truth.rds"))

# ---- Summary: effective factor count per model ------------------------------

r_table <- function(draws, label) {
  cat(label, "— r distribution:\n")
  print(table(draws$r))
  cat("\n")
}

r_table(fit_uglt_on_uglt$draws, "UGLT data / UGLT model")
r_table(fit_plt_on_uglt$draws,  "UGLT data / sparse PLT model")
r_table(fit_uglt_on_plt$draws,  "PLT data  / UGLT model")
r_table(fit_plt_on_plt$draws,   "PLT data  / sparse PLT model")


overfit_res_list <- list(r_table(fit_uglt_on_uglt$draws, "UGLT data / UGLT model"),
                         r_table(fit_plt_on_uglt$draws,  "UGLT data / sparse PLT model"),
                         r_table(fit_uglt_on_plt$draws,  "PLT data  / UGLT model"),
                         r_table(fit_plt_on_plt$draws,   "PLT data  / sparse PLT model"))
saveRDS(overfit_res_list, file = 'overfit_res_list.rds')
