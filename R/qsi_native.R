# Private DLL, separate from WMskelstats's existing Rcpp/OpenMP DLL.
.qsi_native_state <- new.env(parent = emptyenv())

.qsi_native_call <- function(symbol, ...) {
  if (is.null(.qsi_native_state$dll)) {
    root <- getOption("WMskelstats.qsi_native_dir",
                      system.file("qsi-native", package = "WMskelstats"))
    dll <- file.path(root, paste0("qsiNative", .Platform$dynlib.ext))
    if (!nzchar(root) || !file.exists(dll))
      stop("Native QSI backend is missing. Run tools/build_qsi_native.R in the package source, then reinstall WMskelstats.")
    .qsi_native_state$dll <- dyn.load(normalizePath(dll, mustWork = TRUE), local = TRUE)
  }
  address <- getNativeSymbolInfo(paste0("_qsiNative_", symbol),
                                PACKAGE = .qsi_native_state$dll)$address
  .Call(address, ...)
}

.qsi_image_input <- function(x, directory, name) {
  if (is.character(x)) {
    if (length(x) != 1L || is.na(x) || !file.exists(x))
      stop(name, " must name exactly one existing NIfTI file")
    path <- normalizePath(x, winslash = "/", mustWork = TRUE)
    image <- RNifti::readNifti(path, internal = TRUE)
  } else {
    if (!inherits(x, "niftiImage"))
      stop(name, " must be a NIfTI filename or RNifti niftiImage with spatial metadata")
    image <- x
    path <- file.path(directory, paste0(name, ".nii"))
  }
  if (length(dim(image)) != 3L)
    stop(name, " must be a scalar 3-D image; 4-D/vector/tensor images are unsupported")
  if (!is.character(x)) RNifti::writeNifti(image, path, datatype = "float64")
  path
}

#' Apply ITK transforms to a scalar 3-D NIfTI image
#'
#' Uses the fixed image's complete spatial grid. An empty transform list performs
#' identity resampling in physical coordinates. Composite QSIprep .h5 transforms
#' are read by ITK, including their embedded displacement fields. Transform order
#' is ANTs/ITK order (last listed applied first); no transforms are inverted.
#' @param fixed,moving NIfTI filename or RNifti niftiImage.
#' @param transformlist Character vector or list of transform filenames.
#' @param interpolator "linear" or "nearestNeighbor".
#' @param imagetype Must be 0 (scalar 3-D).
#' @param defaultvalue Value outside the moving image; default 0.
#' @param nthread Positive number of ITK work units.
#' @param output Optional output NIfTI path.
#' @param singleprecision Write float32 output when TRUE; interpolation remains double precision.
#' @return RNifti niftiImage on the fixed grid.
#' @export
qsi_apply_transforms <- function(fixed, moving, transformlist = character(),
                                 interpolator = "linear", imagetype = 0L,
                                 defaultvalue = 0, nthread = 1L,
                                 output = NULL, singleprecision = TRUE) {
  if (length(imagetype) != 1L || is.na(imagetype) || imagetype != 0L)
    stop("Only imagetype = 0 (scalar 3-D) is supported")
  interpolator <- match.arg(interpolator, c("linear", "nearestNeighbor"))
  if (length(nthread) != 1L || !is.finite(nthread) || nthread < 1 ||
      nthread > .Machine$integer.max || nthread != as.integer(nthread))
    stop("nthread must be a positive integer")
  if (length(defaultvalue) != 1L || !is.finite(defaultvalue)) stop("Invalid defaultvalue")
  if (!is.logical(singleprecision) || length(singleprecision) != 1L || is.na(singleprecision))
    stop("singleprecision must be TRUE or FALSE")
  if (is.list(transformlist)) {
    if (any(lengths(transformlist) != 1L))
      stop("Each transformlist element must contain exactly one filename")
    transformlist <- unlist(transformlist, use.names = FALSE)
  }
  if (is.null(transformlist)) transformlist <- character()
  if (!is.character(transformlist) || anyNA(transformlist) || any(!file.exists(transformlist)))
    stop("transformlist must contain existing transform filenames")
  if (any(!grepl("\\.(h5|hdf5|mat|tfm|txt)$", transformlist, ignore.case = TRUE)))
    stop("Use ITK .h5/.hdf5/.mat/.tfm/.txt transforms; standalone warp NIfTI files are unsupported")
  transformlist <- vapply(transformlist, normalizePath, character(1),
                          winslash = "/", mustWork = TRUE, USE.NAMES = FALSE)
  work <- tempfile("qsi_native_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  fixed_path <- .qsi_image_input(fixed, work, "fixed")
  moving_path <- .qsi_image_input(moving, work, "moving")
  result <- file.path(work, "result.nii")
  .qsi_native_call("qsi_apply_transforms_cpp", fixed_path, moving_path,
                   transformlist, result, interpolator, as.double(defaultvalue),
                   as.integer(nthread), singleprecision)
  image <- RNifti::readNifti(result, internal = FALSE)
  if (!is.null(output)) qsi_image_write(image, output, singleprecision)
  image
}

#' Resample a scalar image to a reference grid
#' @param moving,fixed NIfTI filename or RNifti niftiImage.
#' @param interpolator "linear" or "nearestNeighbor".
#' @param nthread Positive number of ITK work units.
#' @return RNifti niftiImage.
#' @export
qsi_resample_to_target <- function(moving, fixed, interpolator = "linear", nthread = 1L) {
  qsi_apply_transforms(fixed, moving, character(), interpolator, nthread = nthread)
}

#' Write a scalar 3-D NIfTI image using ITK
#' @param image NIfTI filename or RNifti niftiImage.
#' @param filename Output .nii or .nii.gz filename.
#' @param singleprecision Write float32 rather than float64.
#' @return Output filename, invisibly.
#' @export
qsi_image_write <- function(image, filename, singleprecision = TRUE) {
  if (!is.character(filename) || length(filename) != 1L || is.na(filename) ||
      !grepl("\\.nii(\\.gz)?$", filename, ignore.case = TRUE))
    stop("filename must be one .nii or .nii.gz path")
  if (!is.logical(singleprecision) || length(singleprecision) != 1L || is.na(singleprecision))
    stop("singleprecision must be TRUE or FALSE")
  if (!dir.exists(dirname(filename))) stop("Output directory does not exist")
  filename <- file.path(normalizePath(dirname(filename), winslash = "/", mustWork = TRUE),
                        basename(filename))
  work <- tempfile("qsi_write_")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- .qsi_image_input(image, work, "input")
  # Always stage output separately, permitting input and output to be the same file.
  result <- file.path(work, if (grepl("\\.gz$", filename)) "result.nii.gz" else "result.nii")
  .qsi_native_call("qsi_image_write_cpp", input, result, singleprecision)
  if (!file.copy(result, filename, overwrite = TRUE)) stop("Failed to write ", filename)
  invisible(filename)
}

.qsi_same_grid <- function(x, y, tolerance = 1e-5) {
  identical(dim(x), dim(y)) &&
    isTRUE(all.equal(unname(RNifti::xform(x)), unname(RNifti::xform(y)),
                     tolerance = tolerance, check.attributes = FALSE)) &&
    identical(bitwAnd(as.integer(RNifti::niftiHeader(x)$xyzt_units), 7L),
              bitwAnd(as.integer(RNifti::niftiHeader(y)$xyzt_units), 7L))
}
