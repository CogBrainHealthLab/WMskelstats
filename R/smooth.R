#' 3D Spatial Smoothing for Voxel Matrices
#'
#' @param data_mat Numeric matrix of size N x V (N subjects, V voxels).
#' @param coords Numeric or integer matrix of size V x 3 containing X, Y, Z voxel grid coordinates.
#' @param fwhm Full Width at Half Maximum (scalar or length-3 vector in mm).
#' @param sigma Gaussian standard deviation (scalar or length-3 vector in mm). Ignored if `fwhm` is given.
#' @param voxel_size Voxel dimension in mm (scalar or length-3 vector). Default is c(1, 1, 1).
#'
#' @return A smoothed N x V matrix.
#' @useDynLib WMskelstats, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
smooth_vox <- function(data_mat, 
                                coords, 
                                fwhm = NULL, 
                                sigma = NULL, 
                                voxel_size = c(1, 1, 1)) {
  
  # Validate matrix input dimensions
  data_mat <- as.matrix(data_mat)
  coords <- as.matrix(coords)
  
  if (ncol(data_mat) != nrow(coords)) {
    stop("Number of columns in data_mat (V) must match number of rows in coords.")
  }
  if (ncol(coords) != 3) {
    stop("coords must have exactly 3 columns (X, Y, Z).")
  }
  
  # Convert FWHM to Sigma if provided
  if (!is.null(fwhm)) {
    sigma <- fwhm / (2 * sqrt(2 * log(2)))
  } else if (is.null(sigma)) {
    stop("Must specify either 'fwhm' or 'sigma'.")
  }
  
  # Ensure length-3 vectors for anisotropic support
  if (length(sigma) == 1) sigma <- rep(sigma, 3)
  if (length(voxel_size) == 1) voxel_size <- rep(voxel_size, 3)
  
  # Scale sigma from mm to voxel units
  sigma_voxels <- sigma / voxel_size
  
  # Execute single-threaded C++ routine
  smoothed_mat <- cpp_smooth_voxel_matrix(
    data_mat     = data_mat,
    coords       = coords,
    sigma_voxels = sigma_voxels
  )
  
  return(smoothed_mat)
}