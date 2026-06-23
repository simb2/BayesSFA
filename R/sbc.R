# Simulation-Based Calibration (SBC) for BayesSFA models.
# Reference: Talts et al. (2018) "Validating Bayesian Inference Algorithms
# with Simulation-Based Calibration", arXiv:1804.06788.
#
# Procedure per simulation s = 1, ..., S:
#   1. Draw (theta_s, Y_s) from the joint prior-predictive p(theta) p(Y | theta).
#   2. Run the MCMC sampler to obtain L posterior draws theta_s^{(1:L)}.
#   3. For each scalar summary g, record rank r_s = #{l : g(theta_s^{(l)}) < g(theta_s)}.
# Under a correct sampler r_s ~ Uniform{0, 1, ..., L}, so rank histograms
# should be flat.

# ---- Prior samplers --------------------------------------------------------

#' Sample from the prior of the non-sparse PLT factor model.
#'
#' Draws Lambda, sigma2, and Factors jointly from the prior and returns
#' simulated data Y = Lambda F + epsilon.
#'
#' @param V Number of observed variables.
#' @param q Number of factors.
#' @param N Number of observations.
#' @param nu Degrees-of-freedom hyperparameter; sigma2_k ~ Inv-Gamma(nu/2, nu*s2/2).
#' @param s2 Scale hyperparameter for sigma2.
#' @param c_0 Prior variance for loading entries.
#' @return List with Lambda (V x q), sigma2 (length V), Factors (q x N), Y (V x N).
#' @keywords internal
sim_prior_PLT <- function(V, q, N, nu, s2, c_0) {
  sigma2 <- 1 / rgamma(V, shape = nu / 2, rate = nu * s2 / 2) # correct

  # PLT structure: lower-triangular with positive diagonal
  Lambda <- matrix(0, V, q) # correct
  for (j in seq_len(q)) {
    Lambda[j, j] <- abs(rnorm(1, 0, sqrt(c_0)))  
    if (j < V) {
      rows_below <- seq(j + 1, V)
      Lambda[rows_below, j] <- rnorm(length(rows_below), 0, sqrt(c_0))
    }
  }

  Factors <- matrix(rnorm(q * N), q, N)
  eps <- matrix(rnorm(V * N) * sqrt(rep(sigma2, N)), V, N)
  Y <- Lambda %*% Factors + eps

  list(Lambda = Lambda, sigma2 = sigma2, Factors = Factors, Y = Y)
} # so far so good. 

#' Sample from the prior of the sparse UGLT/PLT factor model.
#'
#' @param V Number of observed variables.
#' @param q Upper bound on factors.
#' @param N Number of observations.
#' @param alpha Length-V shape parameters; sigma2_k ~ Inv-Gamma(alpha_k, beta_k).
#' @param beta Length-V rate parameters for sigma2.
#' @param theta.shape Shape for theta_j ~ Inv-Gamma(theta.shape, theta.rate).
#' @param theta.rate Rate for theta_j.
#' @param hyperparams List with aH and bH for tau_j ~ Beta(aH, bH).
#' @return List with Lambda (V x q), delta (V x q), sigma2, tau, theta, Factors (q x N), Y (V x N).
#' @keywords internal
sim_prior_sparse <- function(V, q, N, alpha, beta, theta.shape, theta.rate,
                             hyperparams, is_uglt = FALSE) {

  repeat {
    sigma2 <- 1 / rgamma(V, shape = alpha, rate = beta)
    theta <- 1 / rgamma(q, shape = theta.shape, rate = theta.rate) # okay so far but let's keep in mind shape vs scale
    # PLT pivot structure: pivot for column j is row j
    # pivots <- seq_len(q)
    pivots <- if (is_uglt) as.integer(sample(1:V, q, replace = FALSE)) else seq_len(q)
    tau    <- rbeta(q, hyperparams$aH, hyperparams$bH)
    delta  <- matrix(0L, V, q)
    Lambda <- matrix(0, V, q)
    for (j in seq_len(q)) {
      piv <- pivots[j]
      delta[piv, j] <- 1L  # pivot row always active
      # Rows below pivot: each independently Bernoulli(tau_j)
      if (piv < V) {
        rows_below <- seq(piv + 1, V)
        delta[rows_below, j] <- rbinom(length(rows_below), 1L, tau[j])
      }

      active <- which(delta[, j] == 1L)
      if (length(active) > 0) {
        sd_vec <- sqrt(theta[j] * sigma2[active])
        vals <- rnorm(length(active), 0, sd_vec)
        # Pivot entry must be positive (identifiability constraint)
        vals[active == piv] <- abs(vals[active == piv])
        Lambda[active, j] <- vals
      }
    }

    if (!is_uglt) break
    delta_id <- delta[rowSums(delta) > 0, , drop = FALSE] + 0
    if (sparvaride::counting_rule_holds(delta_id)) break
  }
  Factors <- matrix(rnorm(q * N), q, N)
  eps <- matrix(rnorm(V * N) * sqrt(rep(sigma2, N)), V, N)
  Y <- Lambda %*% Factors + eps

  list(Lambda = Lambda, delta = delta, sigma2 = sigma2, tau = tau,
       theta = theta, Factors = Factors, Y = Y)
}

# ---- SBC runner: PLT (non-sparse) ------------------------------------------

#' Simulation-based calibration for the non-sparse PLT model
#'
#' Runs S SBC trials for \code{fitBayesPLT}. Each trial draws ground-truth
#' parameters from the prior, simulates data, runs the sampler, and records
#' the rank of the true T-statistic among the \code{n_runs} posterior draws.
#' Under a correctly implemented sampler the ranks are Uniform\{0, ..., L\}.
#'
#' @param S Number of SBC simulations.
#' @param n_runs Number of MCMC iterations per simulation (= L posterior draws
#'   when \code{burn = 1, thin = 1}).
#' @param V Number of observed variables.
#' @param q Number of latent factors.
#' @param N Number of observations per simulation.
#' @param nu Degrees-of-freedom hyperparameter for \eqn{\sigma^2}.
#' @param s2 Scale hyperparameter for \eqn{\sigma^2}.
#' @param c_0 Prior variance for loading entries.
#' @param burn First iteration to retain. Default 1.
#' @param thin Thinning interval. Default 1.
#' @return An \code{sbc_result} object.
#' @export
sbc_PLT <- function(S, n_runs, V, q, N, nu, s2, c_0, burn = 1, thin = 1) {
  ranks  <- numeric(S)
  L_kept <- NULL
  cli::cli_progress_bar("SBC (PLT)", total = S)

  for (s in seq_len(S)) {
    prior  <- sim_prior_PLT(V, q, N, nu, s2, c_0)
    T_true <- sum(diag(prior$Lambda %*% t(prior$Lambda) + diag(prior$sigma2)))

    init <- sim_prior_PLT(V, q, N, nu, s2, c_0)
    fit  <- fitBayesPLT(
      Lambda_0 = init$Lambda, sigma2_0 = init$sigma2,
      n_runs = n_runs, data = prior$Y,
      q = q, nu = nu, s2 = s2, c_0 = c_0,
      thin = thin, burn = burn
    )

    if (is.null(L_kept)) L_kept <- length(fit$T_stat)
    ranks[s] <- sum(fit$T_stat < T_true)
    cli::cli_progress_update()
  }

  structure(list(ranks = ranks, L = L_kept, S = S, model = "PLT", T_stat = fit$T_stat),
            class = "sbc_result")
}

# ---- SBC runner: sparse PLT / UGLT -----------------------------------------

#' Simulation-based calibration for the sparse factor model
#'
#' Runs S SBC trials for \code{fitBSFA} with \code{fixed = TRUE} (number of
#' factors held at q). Each trial draws ground-truth parameters from the prior,
#' simulates data, runs the sampler, and records the rank of the true
#' T-statistic among the posterior draws.
#'
#' @param S Number of SBC simulations.
#' @param n_runs Number of MCMC iterations per simulation.
#' @param V Number of observed variables.
#' @param q Number of latent factors.
#' @param N Number of observations per simulation.
#' @param alpha Length-V shape parameters for \eqn{\sigma^2} prior.
#' @param beta Length-V rate parameters for \eqn{\sigma^2} prior.
#' @param theta.shape Shape hyperparameter for \eqn{\theta} prior.
#' @param theta.rate Rate hyperparameter for \eqn{\theta} prior.
#' @param hyperparams List with \code{aH} and \code{bH} for the \eqn{\tau} prior.
#' @param constraint Either \code{"PLT"} or \code{"UGLT"}.
#' @param burn First iteration to retain. Default 1.
#' @param thin Thinning interval. Default 1.
#' @return An \code{sbc_result} object.
#' @export
sbc_sparse <- function(S, n_runs, V, q, N, alpha, beta, theta.shape, theta.rate,
                       hyperparams, constraint = "PLT", burn = 1, thin = 1) {
  is_uglt <- constraint == "UGLT"
  ranks       <- numeric(S)
  T_true_vec  <- numeric(S)
  T_stat_list <- vector("list", S)
  n_skipped   <- 0L
  L_kept      <- NULL
  cli::cli_progress_bar(paste0("SBC (sparse ", constraint, ")"), total = S)

  for (s in seq_len(S)) {
    prior  <- sim_prior_sparse(V, q, N, alpha, beta, theta.shape, theta.rate,
                               hyperparams, is_uglt = is_uglt)
    T_true <- sum(diag(0.5 * (prior$Lambda %*% t(prior$Lambda) + diag(prior$sigma2))))

    fit <- fitBSFA(
      y = prior$Y, constraint = constraint, fixed = TRUE,
      q = q, n_runs = n_runs,
      alpha = alpha, beta = beta,
      theta.shape = theta.shape, theta.rate = theta.rate,
      hyperparams = hyperparams,
      thin = thin, burn = burn
    )

    T_post <- fit$draws$T_stat
    if (length(T_post) == 0L) {
      ranks[s]        <- NA_real_
      T_true_vec[s]   <- NA_real_
      T_stat_list[[s]] <- numeric(0)
      n_skipped       <- n_skipped + 1L
    } else {
      if (is.null(L_kept)) L_kept <- length(T_post)
      ranks[s]        <- sum(T_post < T_true)
      T_true_vec[s]   <- T_true
      T_stat_list[[s]] <- T_post
    }
    cli::cli_progress_update()
  }

  valid        <- !is.na(ranks)
  structure(list(ranks   = ranks[valid], L = L_kept, S = sum(valid),
                 S_total = S, n_skipped = n_skipped,
                 T_true  = T_true_vec[valid], T_stat = T_stat_list[valid],
                 model   = paste0("sparse ", constraint)),
            class = "sbc_result")
}

# ---- S3 class: sbc_result --------------------------------------------------

#' Print an sbc_result object
#' @param x An \code{sbc_result} object.
#' @param ... Ignored.
#' @export
print.sbc_result <- function(x, ...) {
  cat("SBC result (", x$model, ")\n", sep = "")
  cat("  S =", x$S, "valid simulations")
  if (!is.null(x$n_skipped) && x$n_skipped > 0L)
    cat("  (", x$n_skipped, "skipped — no identifiable draws)")
  cat(",  L =", x$L, "posterior draws\n")
  cat("  Rank range: [", min(x$ranks), ",", max(x$ranks), "]\n")
  invisible(x)
}

#' Plot rank histogram for an sbc_result object
#'
#' Draws a rank histogram for the T-statistic. The red dashed line shows the
#' expected count per bin under a uniform rank distribution.
#'
#' @param x An \code{sbc_result} object.
#' @param ... Passed to \code{hist()}.
#' @export
plot.sbc_result <- function(x, ...) {
  L_vec <- lengths(x$T_stat)
  L <- max(x$ranks, na.rm = TRUE)
  S <- sum(L_vec > 0)
  p <- 1 / (L + 1)
  expected <- S * p
  lo <- qbinom(0.005, S, p)
  hi <- qbinom(0.995, S, p)
  hist(x$ranks,
       breaks = seq(-0.5, L + 0.5, by = 1),
       xlab = "rank", ylab = "count",
       main = paste("SBC rank histogram -", x$model),
       freq = TRUE, ...)
  abline(h = expected, col = "red", lty = 2)
  abline(h = lo, col = "red", lty = 3)
  abline(h = hi, col = "red", lty = 3)
}

