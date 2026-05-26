#' Fits a Bayesian Sparse Factor models under either the UGLT or PLT constraint.
#'
#' Implements a partially collapsed gibbs sampler. Each iteration Each iteration samples \eqn{\tau}, \eqn{\delta},
#' \eqn{l}, \eqn{\Lambda}, \eqn{\sigma^2}, \eqn{F}, and \eqn{\theta} in sequence, followed by
#' a boosting step step. Applies post-processing: spurious factor removal, sign fixing, orders columns, and filters out draws that
#' fail the 3579 counting rule and filters out unidentified draws.

#' @param y \eqn{v \times N} data matrix \eqn{Y}.
#' @param constraint Either 'UGLT' or 'PLT'.
#' @param fixed Either TRUE or FALSE. Determines whether to return estimates with spurious factors removed.
#' @param q Number of factors \eqn{q}.
#' @param n_runs Total number of MCMC iterations.
#' @param alpha Length-v shape parameters \eqn{\alpha} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param beta Length-v rate parameters \eqn{\beta} for the \eqn{G^{-1}} prior on \eqn{\sigma^2}.
#' @param theta.shape Scalar shape hyperparameter \eqn{a_\theta} for the \eqn{G^{-1}} prior on \eqn{\theta}.
#' @param theta.rate Scalar rate hyperparameter \eqn{b_\theta} for the \eqn{G^{-1}} prior on \eqn{\theta}.
#' @param hyperparams List with \code{aH} and \code{bH} for the Beta prior on \eqn{\tau}.
#' @param thin Thinning interval (retain every \code{thin}-th draw after burn-in).
#' @param burn Index of the first draw to retain.
#' @return List with \code{estimates} (posterior summaries by factor dimension r)
#'   and \code{draws} (tibble of retained MCMC draws).
#' @export
fitBSFA <- function(y, constraint, fixed, q, n_runs, alpha, beta, theta.shape, theta.rate, hyperparams, thin = 1, burn = 1) {
  switch(constraint,
    'UGLT' = run_mcmc_UGLT(q, n_runs, alpha, beta, theta.shape, theta.rate, hyperparams, y, thin = thin, burn = burn, fixed),
    'PLT' = run_mcmc_sparse_PLT(q, n_runs, alpha, beta, theta.shape, theta.rate, hyperparams, y, thin = thin, burn = burn, fixed)
  )
}

