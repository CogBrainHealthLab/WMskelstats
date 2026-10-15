#' Multiple-outcome linear regression t statistics
#'
#' Fits ordinary least squares models with one common design matrix using QR.
#' @param X Numeric design matrix with observations in rows and predictors in
#'   columns. Include the intercept explicitly, for example with model.matrix().
#' @param Y Numeric outcome matrix with observations in rows and outcomes in
#'   columns, or a numeric vector for a single outcome. Rows must align with X.
#' @param block_size Positive integer specifying outcomes per processing block.
#' @param tol Reciprocal-condition-number cutoff for the column-scaled design.
#'   This is not the same rank tolerance as in .lm.fit().
#' @param return_coefficients Logical. If TRUE, also return fitted coefficients
#'   in the original units of X. Defaults to FALSE.
#' @return By default, a numeric matrix of conventional coefficient t statistics
#'   for a null value of zero. If return_coefficients = TRUE, a list with matrices
#'   named t and coefficients. All matrices have predictors in rows and outcomes
#'   in columns, with names inherited from X and Y.
#' @details Requires finite data, full column rank and nrow(X) > ncol(X).
#'   Residual degrees of freedom are nrow(X) - ncol(X). No missing observations
#'   are removed. These are ordinary, not robust or mixed-model, standard errors.
#'   Exact or nearly perfect fits produce undefined or unstable t statistics.
#'   The QR decomposition is shared across outcomes within each call.
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' set.seed(1)
#' X <- cbind(Intercept = 1, age = rnorm(100))
#' Y <- matrix(rnorm(100 * 5), nrow = 100)
#' t_values <- lm_fast(X, Y)
#' t_values["age", ]
#' fit <- lm_fast(X, Y, return_coefficients = TRUE)
#' fit$coefficients["age", ]
lm_fast <- function(X, Y, block_size = 2048L, tol = 1e-12,
                    return_coefficients = FALSE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (!is.numeric(X) || !is.numeric(Y) || is.complex(X) || is.complex(Y)) {
    stop("X and Y must contain real numeric values.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"
  if (ncol(X) < 1L || nrow(X) <= ncol(X) ||
      nrow(Y) != nrow(X) || ncol(Y) < 1L) {
    stop("Require nrow(X) == nrow(Y), nrow(X) > ncol(X), and nonempty columns.",
         call. = FALSE)
  }
  if (!is.numeric(block_size) || is.complex(block_size) ||
      length(block_size) != 1L || !is.finite(block_size) ||
      block_size < 1 || block_size > .Machine$integer.max ||
      block_size != floor(block_size)) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || is.complex(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("Require 0 < tol < 1.", call. = FALSE)
  }
  if (!is.logical(return_coefficients) || length(return_coefficients) != 1L ||
      is.na(return_coefficients)) {
    stop("return_coefficients must be TRUE or FALSE.", call. = FALSE)
  }
  # lm_multi_t_cpp is generated in R/RcppExports.R by compileAttributes().
  ans <- lm_multi_t_cpp(X, Y, as.integer(block_size), as.double(tol),
                        return_coefficients)
  dn <- list(colnames(X), colnames(Y))
  if (return_coefficients) {
    dimnames(ans$t) <- dn
    dimnames(ans$coefficients) <- dn
  } else {
    dimnames(ans) <- dn
  }
  ans
}


###############################################################################
#################################INVERSE TOOLS#################################


#' Prepare voxel data for the inverse linear model
#'
#' Removes the nuisance covariates from every voxel once, so each permutation
#' only needs to process the contrast variable. Voxels that are fully explained
#' by the covariates (nothing left after adjustment) are set to zero, which
#' gives them a t of 0.
#'
#' @param vox_data N x V numeric matrix of voxel-wise data (one row per
#'   observation, one column per voxel).
#' @param Z N x p numeric nuisance design including the intercept, i.e. the
#'   model without the contrast of interest.
#' @param tol Tolerance passed to qr() when checking that Z has full column rank.
#' @return A list with Vunit (voxels adjusted for Z and scaled to unit length),
#'   qr (the QR decomposition of Z), and df (residual degrees of freedom of the
#'   full model, N minus p minus 1).
#' @noRd
lm_inverse_prep <- function(vox_data, Z, tol = 1e-12) {
  Z <- as.matrix(Z); vox_data <- as.matrix(vox_data)
  n <- nrow(Z)
  qrZ <- qr(Z, tol = tol)
  if (qrZ$rank < ncol(Z)) stop("Nuisance design is rank deficient.", call. = FALSE)
  df <- n - ncol(Z) - 1L
  if (df < 1L) stop("Not enough observations for the inverse model.", call. = FALSE)
  
  Vres <- qr.resid(qrZ, vox_data)
  nrm  <- sqrt(colSums(Vres^2))
  scl  <- pmax(1, sqrt(colSums(vox_data^2)))
  bad  <- nrm <= sqrt(.Machine$double.eps) * scl   # voxel fully explained by Z
  nrm[bad] <- 1
  Vunit <- Vres / rep(nrm, each = n)
  Vunit[, bad] <- 0
  colnames(Vunit) <- colnames(vox_data)
  list(Vunit = Vunit, qr = qrZ, df = df)
}


#' t statistics for voxel -> y, adjusting for Z (inverse linear model)
#'
#' Computes the t statistic of each voxel as a predictor of y, given the
#' covariates in Z. This equals the partial correlation between y and the voxel
#' converted to a t value, which is the same number an ordinary regression of
#' y on Z and that voxel would give.
#'
#' @param prep Output of lm_inverse_prep().
#' @param y Numeric vector of length N: the contrast variable, acting as the
#'   outcome (the observed values, or a permuted version of them).
#' @return Numeric vector of V t statistics, one per voxel, named after the
#'   columns of vox_data. All zeros if y is fully explained by Z.
#' @importFrom stats setNames
#' @noRd
lm_inverse_t <- function(prep, y) {
  yr <- qr.resid(prep$qr, as.numeric(y))
  ny <- sqrt(sum(yr^2))
  V  <- ncol(prep$Vunit)
  if (ny < sqrt(.Machine$double.eps)) return(stats::setNames(numeric(V), colnames(prep$Vunit)))
  r <- drop(crossprod(prep$Vunit, yr)) / ny
  r <- pmax(pmin(r, 1), -1)
  r * sqrt(prep$df / pmax(1 - r^2, .Machine$double.eps))
}