library(Rcpp)
library(RcppArmadillo)

# Single-threaded C++ engine retaining algebraic O(N)-free EM iterations
Rcpp::sourceCpp(code = '
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List fast_rint_reg_multi(const arma::mat& Y, const arma::mat& X, 
                               const arma::uvec& group_offsets, const arma::uvec& group_sizes,
                               bool reml = true, double tol = 1e-6, int maxiters = 100) {
    int N = Y.n_rows;
    int m = Y.n_cols;
    int K = group_sizes.n_elem;
    int p = X.n_cols;
    
    double df_e = reml ? (N - p) : N;
    double df_b = reml ? std::max(1.0, (double)(K - 1)) : (double)K;
    arma::vec group_sizes_d = arma::conv_to<arma::vec>::from(group_sizes);
    
    // --- 1. Global Pre-computations (O(N) operations performed ONCE) ---
    arma::mat XtX = X.t() * X;
    arma::mat XtY = X.t() * Y;                             // p x m
    arma::vec Yty = arma::sum(arma::square(Y), 0).t();      // m x 1 (y_j^T * y_j)
    
    // Group sums of X (K x p)
    arma::mat Sx(K, p, arma::fill::zeros);
    for (int i = 0; i < K; ++i) {
        int start = group_offsets[i];
        int n_i = group_sizes[i];
        Sx.row(i) = arma::sum(X.rows(start, start + n_i - 1), 0);
    }
    
    // Group sums of Y (K x m)
    arma::mat Sy(K, m, arma::fill::zeros);
    for (int i = 0; i < K; ++i) {
        int start = group_offsets[i];
        int n_i = group_sizes[i];
        Sy.row(i) = arma::sum(Y.rows(start, start + n_i - 1), 0);
    }

    // Output Storage Containers
    arma::mat beta_mat(p, m, arma::fill::zeros);
    arma::mat se_mat(p, m, arma::fill::zeros);
    arma::mat tstat_mat(p, m, arma::fill::zeros);
    arma::vec sigma2_e_vec(m);
    arma::vec sigma2_b_vec(m);

    // --- 2. Sequential Single-Threaded Outcome Loop ---
    for (int j = 0; j < m; ++j) {
        arma::vec XtY_j = XtY.col(j);
        arma::vec Sy_j = Sy.col(j);
        double yty_j = Yty(j);
        
        double sigma2_e = yty_j / N;
        double sigma2_b = sigma2_e * 0.5;
        
        arma::vec beta(p, arma::fill::zeros);
        arma::mat XtVinvX(p, p);
        
        for (int iter = 0; iter < maxiters; ++iter) {
            double gamma = sigma2_b / sigma2_e;
            
            // Fast Matrix-Vector Operations O(K * p + p^2)
            arma::vec w = gamma / (1.0 + group_sizes_d * gamma);
            XtVinvX = XtX - Sx.t() * (Sx.each_col() % w);
            arma::vec XtVinvY = XtY_j - Sx.t() * (w % Sy_j);
            
            arma::vec beta_new = arma::solve(XtVinvX, XtVinvY, arma::solve_opts::fast);
            
            // Algebraic Quadratic Form for ||y - X*beta||^2 (O(p^2), independent of N)
            double res_sq = yty_j - 2.0 * arma::dot(beta_new, XtY_j) + arma::as_scalar(beta_new.t() * XtX * beta_new);
            arma::vec s_e = Sy_j - Sx * beta_new;
            
            double rss = res_sq - arma::dot(w, arma::square(s_e));
            arma::vec u_hat = (sigma2_b / (sigma2_e + group_sizes_d * sigma2_b)) % s_e;
            double sum_post_var = arma::sum(sigma2_b / (1.0 + group_sizes_d * gamma));
            
            double new_sigma2_e = std::max(1e-8, rss / df_e);
            double new_sigma2_b = std::max(1e-8, (arma::dot(u_hat, u_hat) + sum_post_var) / df_b);
            
            if (arma::norm(beta_new - beta, "inf") < tol && 
                std::abs(new_sigma2_e - sigma2_e) < tol) {
                beta = beta_new;
                sigma2_e = new_sigma2_e;
                sigma2_b = new_sigma2_b;
                break;
            }
            
            beta = beta_new;
            sigma2_e = new_sigma2_e;
            sigma2_b = new_sigma2_b;
        }
        
        // Exact REML Covariance Matrix Calculation
        arma::mat inv_XtVinvX = arma::inv_sympd(XtVinvX);
        arma::vec se = arma::sqrt(sigma2_e * inv_XtVinvX.diag());
        arma::vec tstats = beta / se;
        
        beta_mat.col(j) = beta;
        se_mat.col(j) = se;
        tstat_mat.col(j) = tstats;
        sigma2_e_vec(j) = sigma2_e;
        sigma2_b_vec(j) = sigma2_b;
    }

    return Rcpp::List::create(
        Rcpp::Named("coefficients") = beta_mat,
        Rcpp::Named("std_errors")   = se_mat,
        Rcpp::Named("tstats")       = tstat_mat,
        Rcpp::Named("sigma2_e")     = sigma2_e_vec,
        Rcpp::Named("sigma2_b")     = sigma2_b_vec
    );
}
')

# R Wrapper Interface
lme_fast <- function(Y, X, id, intercept = TRUE, reml = TRUE, tol = 1e-6, maxiters = 100) {
  Y_mat <- as.matrix(Y)
  X_mat <- as.matrix(X)
  
  if (intercept) {
    has_intercept <- all(X_mat[, 1] == 1)
    if (!has_intercept) {
      X_mat <- cbind("(Intercept)" = 1, X_mat)
    } else if (is.null(colnames(X_mat)) || colnames(X_mat)[1] == "") {
      colnames(X_mat)[1] <- "(Intercept)"
    }
  }
  
  ord <- order(id)
  Y_sorted <- Y_mat[ord, , drop = FALSE]
  X_sorted <- X_mat[ord, , drop = FALSE]
  id_sorted <- id[ord]
  
  sizes <- as.numeric(table(id_sorted))
  offsets <- c(0, cumsum(sizes)[-length(sizes)])
  
  res <- fast_rint_reg_multi(
    Y = Y_sorted,
    X = X_sorted,
    group_offsets = offsets,
    group_sizes = sizes,
    reml = reml,
    tol = tol,
    maxiters = maxiters
  )
  
  pred_names <- colnames(X_mat)
  if (is.null(pred_names)) pred_names <- paste0("X", 1:ncol(X_mat))
  
  resp_names <- colnames(Y_mat)
  if (is.null(resp_names)) resp_names <- paste0("Y", 1:ncol(Y_mat))
  
  dimnames(res$coefficients) <- list(pred_names, resp_names)
  dimnames(res$std_errors)   <- list(pred_names, resp_names)
  dimnames(res$tstats)       <- list(pred_names, resp_names)
  names(res$sigma2_e)        <- resp_names
  names(res$sigma2_b)        <- resp_names
  
  return(res)
}