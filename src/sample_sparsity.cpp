// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::export]]
arma::umat sample_sparsity(arma::mat y, arma::mat factors, arma::vec tau, arma::vec theta,
                           arma::umat delta, arma::uvec pivots, arma::vec alpha, arma::vec beta,
                           arma::vec inner_prod_y) {

  vec inner_prod_f = sum(square(factors), 1);
  vec log_prior_odds = log(tau) - log(1 - tau);
  uvec all_cols = regspace<uvec>(0, delta.n_cols - 1);

  for (uword j = 0; j < delta.n_cols; j++) {
    if (pivots(j) >= delta.n_rows) {
      continue;
    }

    uvec I_j = regspace<uvec>(pivots(j), delta.n_rows - 1);

    umat delta_j = delta;
    delta_j.shed_col(j);
    uvec other_active = sum(delta_j.rows(I_j), 1);
    uvec zero_mask = (other_active == 0);

    vec O_ij(I_j.n_elem);
    vec posterior_shape = alpha(I_j) + 0.5 * y.n_cols;

    if (any(zero_mask)) {
      // fix: use find() to convert boolean uvec to index uvec
      uvec zr = I_j(find(zero_mask));

      double B_iT = theta(j) / (1 + theta(j) * inner_prod_f(j));
      double D_0 = (-0.5) * std::log1p(theta(j) * inner_prod_f(j));

      vec m_iT = y.rows(zr) * factors.row(j).t();
      vec posterior_rate_null = beta(zr) + 0.5 * inner_prod_y(zr);
      vec posterior_rate_1 = beta(zr) + 0.5 * (inner_prod_y(zr) - square(m_iT) * B_iT);

      uvec zero_idx = find(zero_mask);
      O_ij(zero_idx) = posterior_shape(zero_idx) % (log(posterior_rate_null) - log(posterior_rate_1)) + D_0;
    }

    uvec non_zero_ix = find(zero_mask == 0);

    for (uword k : non_zero_ix) {
      uword i = I_j(k);
      uvec cols_no_j = all_cols;
      cols_no_j.shed_row(j);

      uvec delta_row_i = conv_to<uvec>::from(delta.row(i));
      uvec j_singleton(1);
      j_singleton(0) = j;
      uvec ord_cols = join_vert(cols_no_j(find(delta_row_i(cols_no_j) == 1)), j_singleton);
      uword q_i = ord_cols.n_elem;

      mat X_i = factors.rows(ord_cols).t();

      mat P_i = X_i.t() * X_i + diagmat(1 / theta(ord_cols));
      vec m_i = X_i.t() * y.row(i).t();

      mat L_i = chol(P_i).t();
      vec x_i = solve(arma::trimatl(L_i), m_i);

      double posterior_rate_1_i = beta(i) + 0.5 * (inner_prod_y(i) - sum(square(x_i)));

      double x_star = x_i(q_i - 1);
      double C_0_i = posterior_rate_1_i + 0.5 * x_star * x_star;

      double D_ij = -log(L_i(q_i - 1, q_i - 1)) - 0.5 * std::log(theta(j));

      O_ij(k) = posterior_shape(k) * log(C_0_i / posterior_rate_1_i) + D_ij;
    }

    vec O_post = O_ij + log_prior_odds(j);
    vec log_u = log(randu<vec>(I_j.n_elem));
    uvec cur = conv_to<uvec>::from(delta.col(j))(I_j);

    uvec to_1 = (cur == 0) % (log_u <= O_post);
    uvec to_0 = (cur == 1) % (log_u <= -O_post);

    // fix: directly index into delta rather than filling a temporary
    uvec idx_to_1 = I_j(find(to_1));
    uvec idx_to_0 = I_j(find(to_0));
    for (uword idx : idx_to_1) delta(idx, j) = 1;
    for (uword idx : idx_to_0) delta(idx, j) = 0;
  }

  return delta;
}
