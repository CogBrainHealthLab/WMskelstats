#' Fast Linear Mixed Effects Model with Random Intercepts for Multiple Outcomes
#'
#' Fits random-intercept linear mixed models across multiple outcome variables 
#' (columns of Y) simultaneously against a single design matrix X. Highly efficient 
#' for high-throughput screening (e.g., GWAS, EWAS, transcriptomics).
#'
#' @param Y A numeric matrix or data frame of dimensions \code{N x M}, where N is 
#'   the number of observations and M is the number of outcome variables.
#' @param X A numeric matrix or data frame of dimensions \code{N x p} containing 
#'   the predictor variables.
#' @param id A vector of length N indicating group/cluster membership for each observation.
#' @param add_intercept Logical. If \code{TRUE} (default), checks for an existing 
#'   intercept column of 1s in \code{X} and prepends one if missing.
#' @param gamma Double. Variance ratio parameter (\eqn{\sigma^2_u / \sigma^2_e}). 
#'   Defaults to 0.5.
#'
#' @return A list containing:
#' \item{t_stat}{A matrix of dimensions \code{p x M} containing t-statistics for each predictor and outcome.}
#' \item{coefficients}{A matrix of dimensions \code{p x M} containing estimated fixed-effect coefficients (\eqn{\beta}).}
#'
#' @importFrom Rcpp cppFunction
#' @importFrom stats aggregate
#' @export
lme_fast <- function(Y, X, id, add_intercept = TRUE, gamma = 0.5) {
  # 1. Compile C++ function inline if not already compiled in current session
  if (!exists("fast_rint_reg_multi_y_cpp", mode = "function")) {
    Rcpp::cppFunction(
      code = '
      Rcpp::List fast_rint_reg_multi_y_cpp(const arma::mat& Y, const arma::mat& X, 
                                           const arma::uvec& group_offsets, const arma::uvec& group_sizes,
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
      }',
      depends = "RcppArmadillo",
      env = globalenv()
    )
  }
  
  # 2. Type validation and matrix conversions
  if (missing(Y) || missing(X) || missing(id)) {
    stop("Arguments 'Y', 'X', and 'id' must all be provided.")
  }
  
  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)
  
  if (nrow(X_mat) != nrow(Y_mat) || length(id) != nrow(Y_mat)) {
    stop("Row dimensions of 'X', 'Y', and length of 'id' must match.")
  }
  
  # 3. Automatically handle intercept column
  if (add_intercept) {
    has_intercept <- any(apply(X_mat, 2, function(col) all(col == 1)))
    if (!has_intercept) {
      X_mat <- cbind("(Intercept)" = 1, X_mat)
    }
  }
  
  # 4. Sort inputs sequentially by cluster ID
  ord <- order(id)
  Y_sorted <- Y_mat[ord, , drop = FALSE]
  X_sorted <- X_mat[ord, , drop = FALSE]
  id_sorted <- id[ord]
  
  # 5. Compute C++ group offsets and sizes
  sizes <- as.numeric(table(id_sorted))
  offsets <- c(0, cumsum(sizes)[-length(sizes)])
  
  # 6. Execute compiled C++ routine
  res <- fast_rint_reg_multi_y_cpp(
    Y = Y_sorted,
    X = X_sorted,
    group_offsets = offsets,
    group_sizes = sizes,
    gamma = gamma
  )
  
  # 7. Assign predictor and outcome names to output matrices
  x_names <- colnames(X_mat)
  if (is.null(x_names)) {
    x_names <- c("(Intercept)", paste0("X", seq_len(ncol(X_mat) - 1)))
  }
  
  y_names <- colnames(Y_mat)
  if (is.null(y_names)) {
    y_names <- paste0("Y", seq_len(ncol(Y_mat)))
  }
  
  dimnames(res$coefficients) <- list(x_names, y_names)
  dimnames(res$t_stat)       <- list(x_names, y_names)
  
  return(res)
}
