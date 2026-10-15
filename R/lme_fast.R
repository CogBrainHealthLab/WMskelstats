#' Voxel-wise REML random-intercept models
#'
#' @param Y Numeric N x M matrix of outcomes.
#' @param X Numeric N x P design matrix, optionally without an intercept.
#' @param id Participant/group identifiers of length N.
#' @param add_intercept Add an intercept unless one already exists.
#' @param gamma_max Upper bound for the random/residual variance ratio.
#' @param tol Optimization tolerance on log(1 + gamma).
#' @param max_iter Maximum iterations per optimization bracket.
#' @param grid_size Number of initial variance-ratio search points.
#'
#' @returns List containing coefficients, standard errors, t-statistics,
#'   variance estimates, and fitting status for each outcome.
#'
#' @examples
#' set.seed(1)
#' X <- cbind(Intercept = 1, age = rnorm(100))
#' Y <- matrix(rnorm(100 * 5), nrow = 100)
#' id <- paste0(rep(x = "sub-", 100), floor(runif(100, min=0, max=100)))
#' fit <- suppressWarnings(lme_fast(X, Y, id))
#' 
#' @importFrom utils head
#' @export
lme_fast <- function(
    Y, X, id,
    add_intercept = TRUE,
    gamma_max = 1e8,
    tol = 1e-5,
    max_iter = 100L,
    grid_size = 9L) {
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  
  if (!is.numeric(X) || !is.numeric(Y)) {
    stop("X and Y must be numeric; dummy-code categorical predictors.")
  }
  
  if (nrow(X) != nrow(Y) || length(id) != nrow(Y)) {
    stop("Rows of X and Y and length(id) must match.")
  }
  
  if (ncol(X) < 1L || ncol(Y) < 1L) {
    stop("X and Y must each contain at least one column.")
  }
  
  if (anyNA(id) || any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("Missing/nonfinite values are not supported. Prepare a common complete dataset first.")
  }
  
  if (!is.logical(add_intercept) ||
      length(add_intercept) != 1L ||
      is.na(add_intercept)) {
    stop("add_intercept must be TRUE or FALSE.")
  }
  
  if (length(gamma_max) != 1L ||
      !is.finite(gamma_max) || gamma_max <= 0 ||
      length(tol) != 1L || !is.finite(tol) || tol <= 0) {
    stop("gamma_max and tol must be positive finite scalars.")
  }
  
  if (length(max_iter) != 1L ||
      !is.finite(max_iter) || max_iter < 1 ||
      max_iter != floor(max_iter) ||
      length(grid_size) != 1L ||
      !is.finite(grid_size) || grid_size < 5 ||
      grid_size != floor(grid_size)) {
    stop("max_iter and grid_size must be integers, with grid_size >= 5.")
  }
  
  if (is.null(colnames(X))) {
    colnames(X) <- paste0("X", seq_len(ncol(X)))
  }
  
  if (add_intercept &&
      !any(vapply(seq_len(ncol(X)),
                  function(k) all(X[, k] == 1),
                  logical(1)))) {
    X <- cbind("(Intercept)" = 1, X)
  }
  
  if (nrow(X) <= ncol(X)) {
    stop("Insufficient residual degrees of freedom.")
  }
  
  # Scale design columns for numerical stability.
  # This does not change the fitted model.
  x_scale <- sqrt(colMeans(X^2))
  
  if (any(!is.finite(x_scale)) || any(x_scale == 0)) {
    stop("X contains a zero or numerically invalid column.")
  }
  
  X_scaled <- sweep(X, 2L, x_scale, "/")
  
  if (qr(X_scaled)$rank < ncol(X_scaled)) {
    stop("X is rank deficient; remove redundant predictors.")
  }
  
  # Integer codes avoid ordering problems with numeric versus text IDs.
  id_code <- match(id, unique(id))
  ord <- order(id_code)
  sizes <- tabulate(id_code)
  
  if (length(sizes) < 2L || !any(sizes > 1L)) {
    stop("Variance-component estimation requires multiple participants and repeated observations.")
  }
  
  offsets <- c(0, utils::head(cumsum(sizes), -1L))
  
  res <- fast_rint_reg_multi_y_cpp(
    Y = Y[ord, , drop = FALSE],
    X = X_scaled[ord, , drop = FALSE],
    group_offsets = offsets,
    group_sizes = sizes,
    gamma_max = gamma_max,
    tol = tol,
    max_iter = as.integer(max_iter),
    grid_size = as.integer(grid_size)
  )
  
  # Convert coefficients and SEs back to original predictor units.
  res$coefficients <- sweep(res$coefficients, 1L, x_scale, "/")
  res$std_error <- sweep(res$std_error, 1L, x_scale, "/")
  
  y_names <- colnames(Y)
  if (is.null(y_names)) {
    y_names <- paste0("Y", seq_len(ncol(Y)))
  }
  
  for (nm in c("coefficients", "std_error", "t_stat")) {
    dimnames(res[[nm]]) <- list(colnames(X), y_names)
  }
  
  for (nm in c("var_random", "var_residual", "gamma", "status")) {
    names(res[[nm]]) <- y_names
  }
  
  if (any(res$status != 0L)) {
    warning(
      "Some voxels require inspection: status 1 = upper bound; ",
      "2 = failed/degenerate/unidentified; 3 = iteration limit."
    )
  }
  
  res
}



###############################################################################
#################################INVERSE TOOLS#################################


#' Prepare voxel data for the inverse mixed-effects model
#'
#' Precomputes everything that does not change between permutations: the voxel
#' data, its per-subject means, the nuisance design and the subject coding.
#'
#' @param vox_data N x V numeric matrix of voxel-wise data (one row per
#'   observation, one column per voxel).
#' @param Z N x p numeric nuisance design including the intercept, i.e. the
#'   model without the contrast of interest.
#' @param id Vector of length N identifying the subject for each row. Rows with
#'   the same id are repeated measures of the same subject.
#' @return A list used by lme_inverse_t().
#' @noRd
#' 
lme_inverse_prep <- function(vox_data, Z, id) {
  vox_data <- as.matrix(vox_data); Z <- as.matrix(Z)
  n <- nrow(Z)
  if (length(id) != n) stop("length(id) must equal nrow(Z).", call. = FALSE)
  gid <- match(id, unique(id))
  ng  <- tabulate(gid)
  df  <- n - ncol(Z) - 1L
  if (df < 1L) stop("Not enough observations for the inverse model.", call. = FALSE)
  list(V = vox_data,
       Vmu = rowsum(vox_data, gid, reorder = TRUE) / ng,   # G x V group means
       Z = Z, gid = gid, ng = ng, df = df)
}


#' t statistics for voxel -> y, adjusting for Z, with a random subject intercept
#'
#' @param prep Output of lme_inverse_prep().
#' @param y Numeric vector of length N: the contrast variable, acting as the outcome.
#' @param gamma Non-negative number: the ratio of between-subject variance to
#'   residual variance for y, estimated beforehand with lme_fast(). 0 means
#'   independent observations; larger values mean stronger within-subject similarity.
#' @param block_size Number of voxels processed at a time (limits memory use).
#' @return Numeric vector of V t statistics, one per voxel.
#' @importFrom stats setNames
#' @noRd
lme_inverse_t <- function(prep, y, gamma, block_size = 5000L) {
  gid <- prep$gid; ng <- prep$ng
  kr  <- (1 - 1 / sqrt(1 + ng * gamma))[gid]      # per observation
  
  # T a = a - k_g * mean_g(a); T'T is proportional to V^-1 for V = I + gamma*K
  tw <- function(a) {
    a  <- as.matrix(a)
    mu <- rowsum(a, gid, reorder = TRUE) / ng
    a - kr * mu[gid, , drop = FALSE]
  }
  Zw  <- tw(prep$Z)
  yw  <- tw(y)[, 1L]
  qrZ <- qr(Zw)
  if (qrZ$rank < ncol(Zw)) stop("Whitened nuisance design is rank deficient.", call. = FALSE)
  yr <- qr.resid(qrZ, yw)
  ny <- sqrt(sum(yr^2))
  
  V   <- ncol(prep$V)
  out <- numeric(V)
  if (ny >= sqrt(.Machine$double.eps)) {
    for (s in seq.int(1L, V, by = block_size)) {
      cols <- s:min(V, s + block_size - 1L)
      Vw   <- prep$V[, cols, drop = FALSE] - kr * prep$Vmu[gid, cols, drop = FALSE]
      Vr   <- qr.resid(qrZ, Vw)
      nv   <- sqrt(colSums(Vr^2))
      bad  <- nv <= sqrt(.Machine$double.eps) * pmax(1, sqrt(colSums(Vw^2)))
      r    <- drop(crossprod(Vr, yr)) / (ifelse(bad, 1, nv) * ny)
      r[bad] <- 0
      r    <- pmax(pmin(r, 1), -1)
      out[cols] <- r * sqrt(prep$df / pmax(1 - r^2, .Machine$double.eps))
    }
  }
  stats::setNames(out, colnames(prep$V))
}