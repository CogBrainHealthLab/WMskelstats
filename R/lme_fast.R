#' Fast Linear Mixed Effects Model with Random Intercepts for Multiple Outcomes
#'
#' Fits random-intercept linear mixed models across multiple outcome variables 
#' (columns of Y) simultaneously against a single design matrix X.
#'
#' @param Y A numeric matrix or data frame of dimensions N x M.
#' @param X A numeric matrix or data frame of dimensions N x p.
#' @param id A vector of length N indicating group/cluster membership.
#' @param add_intercept Logical. Checks for an existing intercept column and prepends one if missing.
#' @param gamma Double. Variance ratio parameter. Defaults to 0.5.
#'
#' @return A list containing matrices for t_stat and coefficients.
#' @useDynLib WMskelstats, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
lme_fast <- function(Y, X, id, add_intercept = TRUE, gamma = 0.5) {
  if (missing(Y) || missing(X) || missing(id)) {
    stop("Arguments 'Y', 'X', and 'id' must all be provided.")
  }
  
  X_mat <- as.matrix(X)
  Y_mat <- as.matrix(Y)
  
  if (nrow(X_mat) != nrow(Y_mat) || length(id) != nrow(Y_mat)) {
    stop("Row dimensions of 'X', 'Y', and length of 'id' must match.")
  }
  
  if (add_intercept) {
    has_intercept <- any(apply(X_mat, 2, function(col) all(col == 1)))
    if (!has_intercept) {
      X_mat <- cbind("(Intercept)" = 1, X_mat)
    }
  }
  
  ord <- order(id)
  Y_sorted <- Y_mat[ord, , drop = FALSE]
  X_sorted <- X_mat[ord, , drop = FALSE]
  id_sorted <- id[ord]
  
  sizes <- as.numeric(table(id_sorted))
  offsets <- c(0, cumsum(sizes)[-length(sizes)])
  
  # Execute pre-compiled package C++ routine
  res <- fast_rint_reg_multi_y_cpp(
    Y = Y_sorted,
    X = X_sorted,
    group_offsets = offsets,
    group_sizes = sizes,
    gamma = gamma
  )
  
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