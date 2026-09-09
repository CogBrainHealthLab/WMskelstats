#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List fast_rint_reg_multi_y_cpp(const arma::mat& Y, 
                                     const arma::mat& X, 
                                     const arma::uvec& group_offsets, 
                                     const arma::uvec& group_sizes,
                                     double gamma = 0.5) {
  int N = Y.n_rows;
  int M = Y.n_cols;
  int K = group_sizes.n_elem;
  int p = X.n_cols;
  
  // Precompute weights per group
  arma::vec w(K);
  for (int i = 0; i < K; ++i) {
    w(i) = gamma / (1.0 + group_sizes[i] * gamma);
  }
  
  // Build and invert (X^T V^-1 X) ONCE for all M outcomes
  arma::mat XtVinvX = X.t() * X;
  for (int i = 0; i < K; ++i) {
    int start = group_offsets[i];
    int n_i = group_sizes[i];
    arma::rowvec Sx_i = arma::sum(X.rows(start, start + n_i - 1), 0);
    XtVinvX -= w(i) * (Sx_i.t() * Sx_i);
  }
  
  arma::mat XtVinvX_inv = arma::inv_sympd(XtVinvX);
  arma::vec xtvx_inv_diag = XtVinvX_inv.diag();
  
  // Group sums for Y across all M outcomes
  arma::mat Sy(K, M, arma::fill::zeros);
  for (int i = 0; i < K; ++i) {
    int start = group_offsets[i];
    int n_i = group_sizes[i];
    Sy.row(i) = arma::sum(Y.rows(start, start + n_i - 1), 0);
  }
  
  // Compute X^T V^-1 Y (p x M matrix)
  arma::mat XtVinvY = X.t() * Y;
  for (int i = 0; i < K; ++i) {
    int start = group_offsets[i];
    int n_i = group_sizes[i];
    arma::rowvec Sx_i = arma::sum(X.rows(start, start + n_i - 1), 0);
    XtVinvY -= w(i) * (Sx_i.t() * Sy.row(i));
  }
  
  // Compute Coefficients B (p x M matrix)
  arma::mat Beta = XtVinvX_inv * XtVinvY;
  
  // Compute Residual Variances (1 x M vector)
  arma::rowvec YtVinvY = arma::sum(Y % Y, 0);
  arma::rowvec w_Sy2 = arma::sum(Sy % Sy % arma::repmat(w, 1, M), 0);
  YtVinvY -= w_Sy2;
  
  arma::rowvec beta_XtVinvY = arma::sum(Beta % XtVinvY, 0);
  arma::rowvec sigma2_e = (YtVinvY - beta_XtVinvY) / (N - p);
  
  // Compute t-statistics (p x M matrix)
  arma::mat t_stats(p, M);
  for (int j = 0; j < M; ++j) {
    double s2 = sigma2_e(j);
    for (int i = 0; i < p; ++i) {
      double se = std::sqrt(s2 * xtvx_inv_diag(i));
      t_stats(i, j) = Beta(i, j) / se;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("t_stat") = t_stats,
    Rcpp::Named("coefficients") = Beta
  );
}