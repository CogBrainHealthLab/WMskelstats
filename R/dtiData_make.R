# Replacement for the original dtiData_make(); same existing positional arguments.
# dwi_dir is optional, but qsi_extract supplies the active session's dwi directory.
dtiData_make <- function(sub_s, subfiles, silent = FALSE, dwi_dir = NULL) {
  if (!is.null(dwi_dir)) {
    if (length(dwi_dir) != 1L || is.na(dwi_dir)) stop("Invalid dwi_dir")
    subfiles <- list.files(dwi_dir, recursive = FALSE, full.names = TRUE)
  }
  # Exact dwi-directory membership avoids finding gradients in another modality.
  paths <- gsub("\\", "/", subfiles, fixed = TRUE)
  files <- subfiles[basename(dirname(paths)) == "dwi"]
  # Keep exact subject/session identity, but allow acq-*, run-*, dir-* etc.
  files <- files[startsWith(basename(files), paste0(sub_s, "_"))]
  files <- unique(files[file.exists(files)])
  fail <- function(message) {
    if (!silent) warning(sub_s, ": ", message, call. = FALSE)
    list(NA, NA)
  }

  names <- basename(files)
  dwi <- files[grepl("_space-ACPC_", names) &
                 grepl("_desc-preproc_dwi\\.nii(\\.gz)?$", names)]
  if (!length(dwi)) return(fail("No ACPC preprocessed DWI image found in the active dwi directory."))
  if (length(dwi) > 1L) {
    stop(sub_s, ": multiple ACPC DWI images found. Select one acquisition/run before extraction:\n",
         paste(dwi, collapse = "\n"), call. = FALSE)
  }

  # Derive sidecars from the chosen image, rather than collecting all gradients.
  stem <- sub("\\.nii(\\.gz)?$", "", dwi)
  bval_file <- paste0(stem, ".bval")
  bvec_file <- paste0(stem, ".bvec")
  btable_file <- paste0(stem, ".b_table")
  if (!file.exists(bvec_file)) return(fail(paste("Missing matching bvec:", bvec_file)))
  if (!file.exists(bval_file) && !file.exists(btable_file)) {
    return(fail(paste("Missing matching bval or b_table for", basename(dwi))))
  }

  mask_stem <- sub("_desc-preproc_dwi$", "_desc-brain_mask", stem)
  mask <- c(paste0(mask_stem, ".nii"), paste0(mask_stem, ".nii.gz"))
  mask <- mask[file.exists(mask)]
  if (!length(mask)) {
    # QSIprep may provide one shared ACPC mask without an acquisition/run entity.
    mask <- files[grepl("_space-ACPC_", names) &
                    grepl("_desc-brain_mask\\.nii(\\.gz)?$", names)]
  }
  if (!length(mask)) return(fail("No matching ACPC brain mask found."))
  if (length(mask) > 1L) stop(sub_s, ": multiple possible brain masks; cannot choose safely.")

  bvec <- as.matrix(read.table(bvec_file))
  if (file.exists(bval_file)) {
    bval <- scan(bval_file, quiet = TRUE)
  } else {
    bval <- as.numeric(read.table(btable_file)[[1L]])
  }
  if (!silent) {
    message("  DWI:  ", dwi)
    message("  bvec: ", bvec_file)
    message("  bval: ", if (file.exists(bval_file)) bval_file else btable_file)
  }
  dtiDataobj <- dti::readDWIdata(
    gradient = bvec, bvalue = bval, dirlist = dwi, format = "NIFTI"
  )
  dtiDataobj <- dti::setmask(dtiDataobj, mask)
  list(dtiDataobj, RNifti::readNifti(dwi))
}
