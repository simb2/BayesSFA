#' Run MCMC sampler for the non-sparse PLT factor model
#'
#' Draws posterior samples of the loadings matrix \eqn{\Lambda}, latent factors
#' \eqn{F}, and idiosyncratic variances \eqn{\sigma^2} under the PLT (lower-triangular)
#' factor model. Each iteration applies a parameter-expansion boost step to
#' improve mixing, then sign-corrects draws for identifiability.
#'
#' @param Lambda_0 \eqn{v \times q} matrix of initial loadings.
#' @param sigma2_0 Length-v vector of initial idiosyncratic variances.
#' @param n_runs Number of MCMC iterations. Default 1000.
#' @param data \eqn{v \times N} data matrix \eqn{Y}.
#' @param q Number of latent factors.
#' @param nu Degrees-of-freedom parameter for the inverse-gamma prior on \eqn{\sigma^2}.
#' @param s2 Scale parameter for the inverse-gamma prior on \eqn{\sigma^2}.
#' @param c_0 Prior variance on the free loadings.
#' @param thin Thinning interval applied after burn-in. Default 1 (keep every draw).
#' @param burn Number of initial iterations to discard. Default 1.
#'
#' @return A list with elements:
#' \describe{
#'   \item{Factors}{Array of dimension \code{(kept, q, N)} of factor draws.}
#'   \item{Lambda}{Array of dimension \code{(kept, v, q)} of loadings draws.}
#'   \item{Variances}{Matrix of dimension \code{(kept, v)} of variance draws.}
#'   \item{T_stat}{Numeric vector of the trace statistic at each kept iteration.}
#' }
#' @export
fitBayesPLT <- function(Lambda_0, sigma2_0, n_runs = 1000, data, q, nu, s2, c_0, thin = 1, burn = 1) {
  cli::cli_progress_bar("Sampling from Posterior", total = n_runs)
  W <- array(data = NA, dim = c(n_runs, q, dim(data)[2]))
  L <- array(data = NA, dim = c(n_runs, dim(data)[1], q))
  L[1, , ] <- Lambda_0
  S <- array(data = NA, dim = c(n_runs, dim(data)[1]))
  S[1, ] <- sigma2_0
  T_stat <- numeric(n_runs)

  W[1, , ] <- sample_factors(Lambda_0, sigma2_0, data, q)
  L[1, , ] <- sample_loadings(data, sigma2_0, W[1, , ], c_0)
  S[1, ] <- sample_variances(nu, s2, data, W[1, , ], L[1, , ])
  boost_1 <- boost(Lambda = L[1, , ], factors = W[1, , ], c_0)
  L[1, , ] <- boost_1$Lambda_new
  W[1, , ] <- boost_1$factors_new
  T_stat[1] <- sum(diag(L[1, , ] %*% t(L[1, , ]) + diag(S[1, ])))
  cli::cli_progress_update()

  for (c in 2:n_runs) {
    W[c, , ] <- sample_factors(L[c - 1, , ], S[c - 1, ], data, q)
    L[c, , ] <- sample_loadings(data, S[c - 1, ], W[c, , ], c_0)
    S[c, ] <- sample_variances(nu, s2, data, W[c, , ], L[c, , ])
    boost <- boost(Lambda = L[c, , ], factors = W[c, , ], c_0)
    L[c, , ] <- boost$Lambda_new
    W[c, , ] <- boost$factors_new
    T_stat[c] <- sum(diag(L[c, , ] %*% t(L[c, , ]) + diag(S[c, ])))
    cli::cli_progress_update()
  }
  for (l in 1:n_runs) {
    for (i in 1:dim(data)[2]) {
      W[l, , i] <- diag(sign(diag(L[l, , ]))) %*% W[l, , i]
    }
    L[l, , ] <- L[l, , ] %*% diag(sign(diag(L[l, , ])))
  }
  thin_burn <- seq(burn, n_runs, by = thin)
  return(list(
    Factors = W[thin_burn, , ], Lambda = L[thin_burn, , ],
    Variances = S[thin_burn, ], T_stat = T_stat[thin_burn]
  ))
}
