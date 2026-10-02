# Optional parity check, run BEFORE uninstalling rpyANTs/Python.
# Rscript tests/compare_with_rpyants.R fixed.nii.gz moving.nii.gz transform.h5
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args) == 3L, requireNamespace("rpyANTs", quietly=TRUE),
          requireNamespace("Rcpp", quietly=TRUE))
source("R/qsi_native.R")
options(WMskelstats.qsi_native_dir=normalizePath("inst/qsi-native"))
legacy <- rpyANTs::ants_apply_transforms(fixed=args[1], moving=args[2],
              transformlist=list(args[3]), imagetype=0L, interpolator="linear")
native <- qsi_apply_transforms(fixed=args[1], moving=args[2],
              transformlist=list(args[3]), imagetype=0L, interpolator="linear")
x <- as.array(legacy[])
y <- as.array(native)
stopifnot(identical(dim(x), dim(y)), identical(is.finite(x), is.finite(y)))
keep <- is.finite(x)
if (!any(keep)) stop("No finite voxels to compare")
error <- abs(x[keep]-y[keep])
print(c(max_absolute_error=max(error), mean_absolute_error=mean(error),
        rmse=sqrt(mean(error^2))))
tolerance <- 1e-5 + 1e-5*max(abs(x[keep]))
if (max(error) > tolerance)
  stop("Parity check exceeded tolerance; investigate orientation, precision and ITK versions")
message("Scalar image parity check passed at tolerance ", tolerance)
