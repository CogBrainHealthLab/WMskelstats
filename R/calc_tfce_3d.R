#' Threshold-Free Cluster Enhancement (TFCE) in 3D
#'
#' @param t_stat 3D array of t-statistics.
#' @param E Extent exponent (default = 0.5).
#' @param H Height exponent (default = 2.0).
#' @param dh Step size (default = 0.1).
#' @param connectivity Neighborhood connectivity: 6, 18, or 26 (default = 26).
#'
#' @return 3D array of TFCE statistics.
#' @export
calc_tfce_3d <- function(t_stat, E = 0.5, H = 2.0, dh = 0.1, connectivity = 26) {
  if (!is.array(t_stat) || length(dim(t_stat)) != 3) {
    stop("Input 't_stat' must be a 3D array.")
  }
  if (!connectivity %in% c(6, 18, 26)) {
    stop("Connectivity must be 6, 18, or 26.")
  }
  
  calc_tfce_cpp(
    t_stat = t_stat,
    dims = dim(t_stat),
    E = E,
    H = H,
    dh = dh,
    connectivity = connectivity
  )
}