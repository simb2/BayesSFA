# Shared helpers for simulation studies. Source this file at the top of each
# sim study script.
#
# Dependencies: MASS, ggplot2, reshape2, scoringRules, tibble, furrr, future

# ---- Parallel worker setup ---------------------------------------------------

# `future`'s multisession workers are fresh R processes. `library(BayesSFA)`
# in the main session normally gets auto-detected and re-attached in each
# worker by `future` itself, but call this at the top of any function passed
# to future_map()/future_pmap()/future_imap() as a cheap, explicit guarantee
# that the package (and its compiled Rcpp/RcppArmadillo code) is loaded.
ensure_pkg_loaded <- function() {
  if (!isTRUE(getOption("bsfa_worker_loaded", FALSE))) {
    # Rcpp/RcppArmadillo must be attached *before* BayesSFA, or BayesSFA's
    # compiled code fails at runtime with "function 'enterRNGScope' not
    # provided by package 'Rcpp'" - a namespace-load-order quirk, not a
    # build issue (reproduces even from a clean R CMD INSTALL).
    suppressPackageStartupMessages({
      library(Rcpp)
      library(RcppArmadillo)
      library(BayesSFA)
    })
    options(bsfa_worker_loaded = TRUE)
  }
}

# ---- Data simulators --------------------------------------------------------

sim_data_3 <- function(N, Lambda, Sigma, factors) {
  Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, nrow(Lambda)), Sigma))
}

sim_data_general <- function(N, V, q) {
  Lambda  <- matrix(rnorm(V * q), nrow = V, ncol = q)
  factors <- t(MASS::mvrnorm(N, rep(0, q), diag(q)))
  data    <- Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, V), diag(V)))
  list(data = data, factors = factors, Lambda = Lambda, Sigma = diag(V), N = N, V = V)
}

sim_data_UGLT <- function(N, V, q, plt_structure = FALSE) {
  sigma2 <- rep(1, V)
  theta  <- rep(1, q)
  tau    <- seq_len(q) / (2 * q + 1)
  pivots <- if (plt_structure) seq_len(q) else sort(sample(seq_len(V - 3), q))

  delta <- matrix(0L, V, q)
  for (j in seq_len(q)) {
    delta[pivots[j], j] <- 1L
    below <- seq.int(pivots[j] + 1L, V)
    if (length(below))
      delta[below, j] <- rbinom(length(below), 1, 1 - tau[j])
  }
  Lambda <- matrix(0, V, q)
  for (i in seq_len(V)) {
    active <- which(delta[i, ] == 1L)
    if (length(active))
      Lambda[i, active] <- MASS::mvrnorm(
        1, rep(0, length(active)),
        diag(theta[active], length(active)) * sigma2[i]
      )
  }
  Lambda  <- Lambda %*% diag(diag(sign(Lambda[pivots, ])), q)
  factors <- t(MASS::mvrnorm(N, rep(0, q), diag(q)))
  data    <- Lambda %*% factors + t(MASS::mvrnorm(N, rep(0, V), diag(sigma2)))
  list(data = data, factors = factors, Lambda = Lambda, Sigma = diag(sigma2),
       N = N, V = V, delta = delta, pivots = pivots, theta = theta)
}

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

# ---- CRPS -------------------------------------------------------------------

# For fitBSFA draws: lambda_samps is a list of V x q matrices.
avg_crps <- function(lambda_samps, true_val) {
  V <- nrow(true_val); q <- ncol(true_val); n <- length(lambda_samps)
  samps <- array(unlist(lambda_samps), dim = c(V, q, n))
  total <- 0
  for (i in seq_len(V))
    for (j in seq_len(q))
      total <- total + scoringRules::crps_sample(true_val[i, j], dat = samps[i, j, ])
  total / (V * q)
}

# For fitBayesPLT draws: lambda_samps is a (kept x V x q) array.
avg_crps_array <- function(lambda_samps, true_val) {
  V <- nrow(lambda_samps[1, , ]); q <- ncol(lambda_samps[1, , ])
  total <- 0
  for (i in seq_len(V))
    for (j in seq_len(q))
      total <- total + scoringRules::crps_sample(true_val[i, j], dat = lambda_samps[, i, j])
  total / (V * q)
}

# ---- Goodness-of-fit --------------------------------------------------------

contributed_variance_noise <- function(L, y) {
  sum(diag(L %*% t(L))) / sum(diag(cov(t(y))))
}

# True signal-to-total-variance ratio tr(Lambda Lambda^T) / tr(Omega), with
# Omega = Lambda Lambda^T + Sigma the true generative total covariance.
# Paired with the `*_noise` columns already saved by the sim studies (which
# are tr(Lambda_est Lambda_est^T) / tr(cov(y)), i.e. the same ratio computed
# from the fit, using tr(cov(y)) as the practical stand-in for tr(Omega_hat)
# since Sigma_hat was never saved to disk) to get the Noise Error metric
# |true_ratio - est_ratio|.
true_noise_ratio <- function(Lambda, Sigma) {
  signal <- sum(diag(Lambda %*% t(Lambda)))
  signal / (signal + sum(diag(Sigma)))
}

# ... absorbs extra pmap columns (V, q, noisy, etc.)
compute_goodness_of_fit <- function(N, Lambda, Sigma, factors,
                                    Lambda_est, sigma_est, factors_est, ...) {
  mse_vec <- replicate(20, {
    samp  <- sim_data_3(N, Lambda,     Sigma,           factors)
    samp2 <- sim_data_3(N, Lambda_est, diag(sigma_est), factors_est)
    mean((samp - samp2)^2)
  })
  list(mean = mean(mse_vec), draws = mse_vec)
}

# ---- Posterior extraction ---------------------------------------------------

# Extract estimates from a fitBSFA result with fixed = TRUE (r = q always).
get_sparse_est <- function(fit, q) {
  est <- fit$estimates[[as.character(q)]]
  list(
    lambda_est   = est$lambda_est,
    sigma2_mean  = est$sigma2_mean,
    factors_est  = est$factors_est,
    lambda_draws = fit$draws$Lambda_test
  )
}

# Extract posterior means from a fitBayesPLT result.
get_plt_est <- function(fit) {
  list(
    lambda_est   = apply(fit$Lambda,    c(2, 3), mean),
    sigma2_mean  = apply(fit$Variances, 2,       mean),
    factors_est  = apply(fit$Factors,   c(2, 3), mean)
  )
}

# ---- PSIS-LOO (entry-wise) ---------------------------------------------------
#

# fitBSFA draws (one list element per retained iteration): Lambda_test[[i]] is
# V x r_i, W[[i]] is r_i x N (r_i can vary by draw after spurious-factor
# filtering -- fine here since only the V x N product Lambda %*% F is used),
# sigma_test[[i]] is length V. Returns an (n_draws x V*N) pointwise log-lik
# matrix, column-major over (feature, subject) to match as.vector(y).
pointwise_loglik_sparse <- function(fit, y) {
  N <- ncol(y)
  draws <- fit$draws
  n_draws <- nrow(draws)
  log_lik <- matrix(NA_real_, nrow = n_draws, ncol = length(y))
  for (i in seq_len(n_draws)) {
    pred <- draws$Lambda_test[[i]] %*% draws$W[[i]]
    sd_i <- rep(sqrt(draws$sigma_test[[i]]), times = N)
    log_lik[i, ] <- dnorm(as.vector(y), mean = as.vector(pred), sd = sd_i, log = TRUE)
  }
  log_lik
}

# fitBayesPLT draws (arrays): Lambda is (n_draws x V x q), Factors is
# (n_draws x q x N), Variances is (n_draws x V).
pointwise_loglik_plt <- function(fit, y) {
  N <- ncol(y)
  n_draws <- dim(fit$Lambda)[1]
  log_lik <- matrix(NA_real_, nrow = n_draws, ncol = length(y))
  for (i in seq_len(n_draws)) {
    pred <- fit$Lambda[i, , ] %*% fit$Factors[i, , ]
    sd_i <- rep(sqrt(fit$Variances[i, ]), times = N)
    log_lik[i, ] <- dnorm(as.vector(y), mean = as.vector(pred), sd = sd_i, log = TRUE)
  }
  log_lik
}

# Thin wrapper so callers don't need to know the loo:: call signature.
run_psis_loo <- function(log_lik) {
  loo::loo(log_lik)
}

# WAIC (BayesianLOOCV-3.pdf SS3.8): a cheaper, non-importance-sampling
# alternative on the same log-lik matrix -- a Taylor-series approximation to
# LOO with no refitting or reweighting, but also no per-point reliability
# diagnostic like PSIS-LOO's Pareto k, and known to be biased under weak
# priors / hierarchical structure (which the sparse shrinkage prior on Lambda
# arguably is). Useful as a cross-check against the PSIS-LOO ranking, not a
# replacement for it.
run_waic <- function(log_lik) {
  loo::waic(log_lik)
}

# Share of entries with Pareto k > 0.7 (loo's threshold above which the PSIS
# estimate for that entry is considered unreliable; see BayesianLOOCV-3.pdf
# SS3.6.3), plus the max k, per model.
pareto_k_summary <- function(loo_obj) {
  k <- loo_obj$diagnostics$pareto_k
  c(n_high_k = sum(k > 0.7), max_k = max(k), prop_high_k = mean(k > 0.7))
}

# ---- Plots ------------------------------------------------------------------

# Shared model color/theme so every sim-study boxplot looks consistent.
MODEL_COLORS <- c(
  "Non-sparse PLT" = "#283593",
  "Sparse PLT"     = "#2e7d32",
  "UGLT"           = "#4db6ac"
)

bsfa_boxplot_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background   = ggplot2::element_blank(),
      strip.text.y.right = ggplot2::element_text(angle = 0, hjust = 0),
      strip.text.x       = ggplot2::element_text(face = "bold"),
      legend.position     = "bottom",
      panel.grid.major.x  = ggplot2::element_blank(),
      panel.spacing       = ggplot2::unit(0.5, "lines")
    )
}

# Point-range plot of an elpd estimate +/- 1 SE per model, from a named list
# of loo:: results as returned by run_psis_loo()/run_waic() (e.g.
# list("Non-sparse PLT" = ..., "Sparse PLT" = ..., "UGLT" = ...)). `estimate`
# selects the row of loo_obj$estimates to plot -- "elpd_loo" for PSIS-LOO
# results, "elpd_waic" for WAIC results.
plot_elpd_comparison <- function(loo_by_model, estimate = "elpd_loo",
                                  title = "PSIS-LOO model comparison") {
  df <- purrr::imap_dfr(loo_by_model, function(loo_obj, model) {
    est <- loo_obj$estimates[estimate, ]
    tibble::tibble(model = model, elpd = est[["Estimate"]], se = est[["SE"]])
  })
  df$model <- factor(df$model, levels = names(MODEL_COLORS))

  ggplot2::ggplot(df, ggplot2::aes(x = model, y = elpd, color = model)) +
    ggplot2::geom_pointrange(
      ggplot2::aes(ymin = elpd - se, ymax = elpd + se),
      linewidth = 0.8, size = 0.6
    ) +
    ggplot2::scale_color_manual(values = MODEL_COLORS, guide = "none") +
    ggplot2::labs(x = NULL, y = estimate, title = title,
                  subtitle = "Point = estimate, bars = ± 1 SE (higher is better)") +
    bsfa_boxplot_theme()
}

# Save a plot into tests/figures/ (created if missing), as PNG.
save_bsfa_plot <- function(plot, filename, width = 8, height = 6) {
  fig_dir <- here::here("tests", "figures")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(fig_dir, filename), plot = plot, width = width, height = height, dpi = 300)
}

# Save MCMC diagnostic plots (trace, ACF, and factor-count trace if present)
# for one fit's draws into tests/figures/diagnostics/.
# draws: a fitBSFA draws tibble (has $T_stat, $r) or a fitBayesPLT-style list
# (has $T_stat only).
save_mcmc_diagnostics <- function(draws, label) {
  diag_dir <- here::here("tests", "figures", "diagnostics")
  dir.create(diag_dir, showWarnings = FALSE, recursive = TRUE)

  T_stat <- draws$T_stat
  trace_df <- tibble::tibble(iteration = seq_along(T_stat), T_stat = T_stat)

  trace_plot <- ggplot2::ggplot(trace_df, ggplot2::aes(x = iteration, y = T_stat)) +
    ggplot2::geom_line(linewidth = 0.4, alpha = 0.7) +
    ggplot2::labs(x = "Iteration (post burn-in/thin)", y = "T statistic",
                  title = paste0(label, " — trace")) +
    ggplot2::theme_minimal()

  has_r <- !is.null(draws$r)
  if (has_r) {
    r_df <- tibble::tibble(iteration = seq_along(draws$r), r = draws$r)
    r_plot <- ggplot2::ggplot(r_df, ggplot2::aes(x = iteration, y = r)) +
      ggplot2::geom_step(linewidth = 0.4) +
      ggplot2::labs(x = "Iteration (post burn-in/thin)", y = "Number of active factors (r)",
                    title = paste0(label, " — factor count")) +
      ggplot2::theme_minimal()
  }

  if (has_r && requireNamespace("patchwork", quietly = TRUE)) {
    combined <- trace_plot / r_plot
    ggplot2::ggsave(file.path(diag_dir, paste0(label, "_trace.png")),
                     combined, width = 8, height = 6, dpi = 300)
  } else {
    ggplot2::ggsave(file.path(diag_dir, paste0(label, "_trace.png")),
                     trace_plot, width = 8, height = 4, dpi = 300)
    if (has_r) {
      ggplot2::ggsave(file.path(diag_dir, paste0(label, "_r_trace.png")),
                       r_plot, width = 8, height = 4, dpi = 300)
    }
  }

  png(file.path(diag_dir, paste0(label, "_acf.png")), width = 800, height = 600, res = 150)
  acf(T_stat, main = paste0(label, " — ACF of T statistic"))
  dev.off()

  invisible(NULL)
}

make_heat_map <- function(Lambda) {
  Lambda_df       <- reshape2::melt(Lambda)
  Lambda_df$Var1  <- factor(Lambda_df$Var1, levels = rev(unique(Lambda_df$Var1)))
  Lambda_df$Var2  <- factor(Lambda_df$Var2)
  ggplot2::ggplot(Lambda_df, ggplot2::aes(y = Var1, x = Var2, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "j", y = "i") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
                   legend.position = "none")
}

# post_draws: fit$draws for fitBSFA, or list(T_stat = fit$T_stat) for fitBayesPLT.
make_trace_acf_plot <- function(post_draws) {
  T_stat_df <- tibble::tibble(
    Value     = post_draws$T_stat,
    Iteration = seq_along(post_draws$T_stat)
  )
  list(
    acf   = acf(post_draws$T_stat, plot = FALSE),
    trace = ggplot2::ggplot(T_stat_df, ggplot2::aes(x = Iteration, y = Value)) +
      ggplot2::geom_line(linewidth = 0.4, alpha = 0.5) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Iteration", y = "T statistic")
  )
}
