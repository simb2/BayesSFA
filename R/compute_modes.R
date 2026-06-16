#' @keywords internal
#' Compute the posterior mode of \eqn{\delta}
#'
#' Returns the most frequently visited sparsity pattern across MCMC draws.
#'
#' @param x Tibble of draws (filtered to a fixed r), with a \code{delta_test} list-column.
#' @return \eqn{v \times r} integer matrix of the modal sparsity pattern.
post_mode_delta <- function(x) {
  delta_arr <- simplify2array(x$delta_test)  # V x r x n_draws
  patterns <- apply(delta_arr, 3, function(d) paste(d, collapse = ","))
  mode_pattern <- names(which.max(table(patterns)))
  delta_arr[, , match(mode_pattern, patterns)]
}

#' Compute the posterior mean of \eqn{\Lambda} under the modal \eqn{\delta}
#'
#' Averages loading matrices over draws where the sparsity pattern matches
#' \code{post_mode_delta} exactly.
#'
#' @param x Tibble of draws with \code{Lambda_test} and \code{delta_test} list-columns.
#' @param post_mode_delta \eqn{v \times r} modal sparsity matrix from \code{post_mode_delta()}.
#' @return \eqn{v \times r} posterior mean loading matrix.
post_est_lambda <- function(x, post_mode_delta) {
  lambda_arr <- simplify2array(x$Lambda_test)  # V x r x n_draws
  delta_arr <- simplify2array(x$delta_test)    # V x r x n_draws
  r <- ncol(post_mode_delta)
  result <- matrix(0, nrow(post_mode_delta), r)
  for (j in seq_len(r)) {
    col_match <- apply(delta_arr[, j, , drop = FALSE], 3,
                       function(d) all(d == post_mode_delta[, j]))
    if (any(col_match)) {
      result[, j] <- rowMeans(lambda_arr[, j, col_match, drop = FALSE])
    } else {
      # fallback: marginal conditional mean for this column
      d_sum <- rowSums(delta_arr[, j, ])
      l_sum <- rowSums(lambda_arr[, j, ] * delta_arr[, j, ])
      result[, j] <- (l_sum / pmax(d_sum, 1)) * post_mode_delta[, j]
    }
  }
  result
}

#' Compute posterior means of \eqn{\sigma^2}, \eqn{\tau}, \eqn{\theta}, and \eqn{F}
#'
#' Averages \eqn{\sigma^2}, \eqn{\tau}, \eqn{\theta}, and \eqn{F} across all draws in \code{x}.
#'
#' @param x Tibble of draws with list-columns \code{sigma_test}, \code{tau_test},
#'   \code{theta_test}, and \code{W}.
#' @return List with \code{sigma2_mean}, \code{tau_mean}, \code{theta_mean}, \code{factors_est}.
compute_post_means <- function(x) {
  list(
    sigma2_mean = rowMeans(do.call(cbind, x$sigma_test)),
    tau_mean = rowMeans(do.call(cbind, x$tau_test)),
    theta_mean = rowMeans(do.call(cbind, x$theta_test)),
    factors_est = apply(simplify2array(x$W), c(1, 2), mean)
  )
}
