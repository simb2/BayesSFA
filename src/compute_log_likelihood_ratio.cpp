// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::export]]
double compute_log_likelihood_ratio(umat delta, uword l_new, uword j, mat factors, mat y,
  vec inner_prod_y, vec alpha, vec beta, vec theta) {

  uword N = y.n_cols;

  // R passes 1-based indices; convert to 0-based for Armadillo
  l_new -= 1;
  j -= 1;

  delta(l_new, j) = 1;

  if (sum(delta.row(l_new)) == 1) {
    // --- Null case: j is the only active column in row l_new ---
    // Compares delta[l_new,j]=1 (only j active) vs delta[l_new,j]=0 (empty model).
    vec X_i_1 = factors.row(j).t();                        // N-vector
    double L_1_inv = dot(X_i_1, X_i_1) + 1.0 / theta(j);
    double L_1 = 1.0 / L_1_inv;
    double m_i_1 = as_scalar(y.row(l_new) * X_i_1);          // scalar inner product

    return 0.5 * std::log(L_1 / theta(j)) -
      (N / 2.0 + alpha(l_new)) * (
        std::log(beta(l_new) + 0.5 * (inner_prod_y(l_new) - m_i_1 * m_i_1 * L_1)) -
        std::log(beta(l_new) + 0.5 * inner_prod_y(l_new)));

  } else {
    // --- Non-null case: row l_new has other active columns besides j ---
    // Compares full model (active set includes j) vs reduced model (j excluded).

    // delta[l_new,j]=1 quantities  (j already set to 1 above)
    uvec active_1 = find(delta.row(l_new).t() == 1);
    mat X_i_1 = factors.rows(active_1).t();                 // N x q_1
    mat L_1_inv = X_i_1.t() * X_i_1 + diagmat(1.0 / theta(active_1));
    mat L_1 = inv(L_1_inv);
    vec m_i_1 = X_i_1.t() * y.row(l_new).t();              // q_1-vector

    // delta[l_new,j]=0 quantities
    delta(l_new, j) = 0;
    uvec active_0 = find(delta.row(l_new).t() == 1);
    mat X_i_0 = factors.rows(active_0).t();                 // N x q_0
    mat L_0_inv = X_i_0.t() * X_i_0 + diagmat(1.0 / theta(active_0));
    mat L_0 = inv(L_0_inv);
    vec m_i_0 = X_i_0.t() * y.row(l_new).t();              // q_0-vector

    // Schur complement of L_0'_inv in L_1_inv (block matrix det lemma):
    // det(L_1) / det(L_0') = 1 / schur, so log-det contribution is -log(schur) - log(theta_j)
    double schur = as_scalar(
      factors.row(j) * factors.row(j).t() +
      1.0 / theta(j) -
      factors.row(j) * X_i_0 * L_0 * X_i_0.t() * factors.row(j).t()
    );

    double term_2 = (N / 2.0 + alpha(l_new)) * (
      std::log(beta(l_new) + 0.5 * (inner_prod_y(l_new) - as_scalar(m_i_1.t() * L_1 * m_i_1))) -
      std::log(beta(l_new) + 0.5 * (inner_prod_y(l_new) - as_scalar(m_i_0.t() * L_0 * m_i_0))));

    return 0.5 * (-std::log(schur) - std::log(theta(j))) - term_2;
  }
}
