// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;

// [[Rcpp::export]]
Rcpp::List sample_loadings_variances(mat y, mat factors, mat delta, vec theta, 
  vec alpha, vec beta, vec inner_prod_y) {
  mat Lambda_new = mat(delta.n_rows, delta.n_cols, fill::zeros);
  Rcpp::NumericVector sigma2_new(delta.n_rows);
  
  for (uword i = 0; i < delta.n_rows; i ++) {
    if (sum(delta.row(i) != 0)) {
      uvec filter = find(delta.row(i) == 1);
      vec theta_a = theta.rows(filter);
      mat L_0_inv = diagmat(1/theta_a);
      mat X_i_delta = factors.rows(filter).t();
      mat P_i = X_i_delta.t() * X_i_delta + L_0_inv;

      vec m_iN = X_i_delta.t() * y.row(i).t();
      mat L_chol = chol(P_i).t();

      vec x = solve(trimatl(L_chol), m_iN);
      sigma2_new(i) = 1.0 / R::rgamma(alpha(i) + y.n_cols/2.0,
        1.0 / (beta(i) + 0.5*(inner_prod_y(i) - sum(square(x)))));

      vec mi = solve(trimatl(L_chol).t(), x);
      vec sample = mi + sqrt(sigma2_new(i)) * solve(trimatl(L_chol).t(), randn<vec>(filter.n_elem));
      Lambda_new.submat(uvec{i}, filter) = sample.t();
    } else {
      sigma2_new(i) = 1.0 / R::rgamma(alpha(i) + y.n_cols/2.0, 1.0 / (beta(i) + 0.5 * inner_prod_y(i)));
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("Lambda_new") = Lambda_new,
    Rcpp::Named("sigma2_new") = sigma2_new
  );
}