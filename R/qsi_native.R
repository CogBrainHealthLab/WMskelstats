# QSI backend using an existing ANTs executable and RNifti.
# Replaces the entire earlier R/qsi_native.R; no custom QSI DLL is loaded.
.qsi_ants_executable <- function() {
  executable <- getOption("WMskelstats.antsApplyTransforms", "antsApplyTransforms")
  if (!is.character(executable) || length(executable) != 1L ||
      is.na(executable) || !nzchar(executable))
    stop("Set options(WMskelstats.antsApplyTransforms='/full/path/to/antsApplyTransforms')")
  if (!file.exists(executable)) executable <- unname(Sys.which(executable))
  if (!nzchar(executable) || !file.exists(executable))
    stop("antsApplyTransforms was not found. Add its bin directory to PATH or set options(WMskelstats.antsApplyTransforms='/full/path/to/antsApplyTransforms').")
  normalizePath(executable, winslash = "/", mustWork = TRUE)
}

.qsi_run_ants <- function(fixed, moving, transforms, output,
                          interpolator, defaultvalue, nthread, log_file) {
  executable <- .qsi_ants_executable()
  old_threads <- Sys.getenv("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(old_threads)) Sys.unsetenv("ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS") else
      Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS = old_threads)
  }, add = TRUE)
  Sys.setenv(ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS = as.character(nthread))
  args <- c("-d", "3", "-e", "0", "-i", moving, "-r", fixed,
            "-o", output, "-n", if (interpolator == "linear") "Linear" else "NearestNeighbor",
            "-f", format(defaultvalue, digits = 17, trim = TRUE, decimal.mark = "."),
            "--float", "0")
  if (!length(transforms)) transforms <- "identity"
  for (transform in transforms) args <- c(args, "-t", transform)
  status <- system2(executable, args = vapply(args, shQuote, character(1)),
                    stdout = log_file, stderr = log_file)
  if (!identical(as.integer(status), 0L) || !file.exists(output)) {
    details <- if (file.exists(log_file)) paste(readLines(log_file, warn = FALSE), collapse = "\n") else ""
    stop("antsApplyTransforms failed (exit status ", status, ").\n", details, call. = FALSE)
  }
  invisible(output)
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

#' Apply transforms using the ANTs command-line executable
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
#' @param nthread Positive number of threads requested for the ANTs subprocess.
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
  .qsi_run_ants(fixed_path, moving_path, transformlist, result,
                 interpolator, defaultvalue, as.integer(nthread),
                 file.path(work, "ants.log"))
  # Match the old wrapper's output-storage option; computation is double precision.
  image <- RNifti::asNifti(RNifti::readNifti(result, internal = FALSE),
                           datatype = if (singleprecision) "float32" else "float64",
                           internal = FALSE)
  if (!is.null(output)) qsi_image_write(image, output, singleprecision)
  image
}

#' Resample a scalar image to a reference grid
#' @param moving,fixed NIfTI filename or RNifti niftiImage.
#' @param interpolator "linear" or "nearestNeighbor".
#' @param nthread Positive number of threads requested for the ANTs subprocess.
#' @return RNifti niftiImage.
#' @export
qsi_resample_to_target <- function(moving, fixed, interpolator = "linear", nthread = 1L) {
  qsi_apply_transforms(fixed, moving, character(), interpolator, nthread = nthread)
}

#' Write a scalar 3-D NIfTI image using RNifti
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
  RNifti::writeNifti(RNifti::readNifti(input, internal = FALSE), result,
                     datatype = if (singleprecision) "float32" else "float64")
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
