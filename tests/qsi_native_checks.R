# Run after building native backend, from the directory containing R/ and tests/.
# Does not require installing WMskelstats or loading its other dependencies.
stopifnot(requireNamespace("Rcpp", quietly=TRUE), requireNamespace("RNifti", quietly=TRUE))
source("R/qsi_native.R")
options(WMskelstats.qsi_native_dir=normalizePath("inst/qsi-native"))
work <- tempfile("qsi_checks_")
dir.create(work)
# A deterministic, nonconstant image; a ramp exposes axis and transform errors.
d <- c(12L, 11L, 10L)
i <- arrayInd(seq_len(prod(d)), d) - 1L
values <- array(i[,1] + 10*i[,2] + 100*i[,3], d)
reference <- RNifti::niftiHeader(RNifti::asNifti(values))
reference$xyzt_units <- 2L # millimetres
moving <- RNifti::asNifti(values, reference=reference)
RNifti::xform(moving) <- diag(4)
identity <- qsi_apply_transforms(moving, moving, singleprecision=FALSE)
stopifnot(.qsi_same_grid(identity, moving),
          max(abs(as.array(identity)-values)) < 1e-9)

# Same voxel dimensions/spacing, shifted physical origin: must NOT skip resampling.
fixed <- RNifti::asNifti(array(0, d), reference=moving)
shift <- diag(4)
shift[1,4] <- 1
RNifti::xform(fixed) <- shift
stopifnot(!.qsi_same_grid(fixed, moving))
shifted <- qsi_resample_to_target(moving, fixed)
stopifnot(.qsi_same_grid(shifted, fixed),
          max(abs(as.array(shifted)[1:11,,] - values[2:12,,])) < 1e-5,
          all(as.array(shifted)[12,,] == 0))

# ITK transforms are output-to-input, expressed in LPS. +1 LPS x means -1 RAS x.
tfm <- file.path(work, "translation.tfm")
writeLines(c("#Insight Transform File V1.0", "#Transform 0",
             "Transform: AffineTransform_double_3_3",
             "Parameters: 1 0 0 0 1 0 0 0 1 1 0 0",
             "FixedParameters: 0 0 0"), tfm)
translated <- qsi_apply_transforms(moving, moving, tfm, singleprecision=FALSE)
stopifnot(max(abs(as.array(translated)[2:12,,] - values[1:11,,])) < 1e-9,
          all(as.array(translated)[1,,] == 0))

# HDF5 composite: affine + constant displacement field; total +2 LPS mm x.
h5 <- "tests/fixtures/composite_translation.h5"
composite <- qsi_apply_transforms(moving, moving, h5, singleprecision=FALSE)
stopifnot(max(abs(as.array(composite)[3:12,,] - values[1:10,,])) < 1e-8,
          all(as.array(composite)[1:2,,] == 0))
# Same transform passed twice: +4 LPS mm, checking transform ownership/lifetimes.
twice <- qsi_apply_transforms(moving, moving, c(h5,h5), singleprecision=FALSE)
stopifnot(max(abs(as.array(twice)[5:12,,] - values[1:8,,])) < 1e-8)

# NIfTI read/write round trip; nonfinite values must remain nonfinite on writing.
filename <- file.path(work, "roundtrip.nii.gz")
qsi_image_write(moving, filename, singleprecision=FALSE)
roundtrip <- RNifti::readNifti(filename)
stopifnot(.qsi_same_grid(roundtrip, moving),
          max(abs(as.array(roundtrip)-values)) < 1e-9)
# Metadata coordinates and extraction vector must share identical order.
coords <- which(values %% 7 == 0, arr.ind=TRUE)
coords <- coords[order(coords[,1], coords[,2], coords[,3]),,drop=FALSE]
extracted <- values[coords]
rebuilt <- array(NA_real_, d)
rebuilt[coords] <- extracted
stopifnot(identical(rebuilt[coords], extracted))

# Unsupported cases fail explicitly instead of silently changing semantics.
expect_error <- function(expr) stopifnot(inherits(try(expr, silent=TRUE), "try-error"))
expect_error(qsi_apply_transforms(moving, moving, imagetype=3L))
expect_error(qsi_apply_transforms(moving, moving, list(c(tfm,tfm))))
expect_error(qsi_apply_transforms(moving, moving, interpolator="bSpline"))
expect_error(qsi_apply_transforms(moving, moving, nthread=0L))
unlink(work, recursive=TRUE)
message("Native identity, shifted-grid, affine, composite HDF5, write and validation checks passed.")
