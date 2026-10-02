# Look only in this session, then in the shared subject/anat directory.
# Never use the inverse transform or a transform from a different session.
.qsi_find_transform <- function(root, subid, session = "") {
  if (is.null(root)) return(character())
  suffix <- "_from-ACPC_to-MNI152NLin2009cAsym_mode-image_xfm.h5"
  subject_dir <- file.path(root, subid)
  pick <- function(files, prefix, exact = FALSE) {
    names <- basename(files)
    keep <- if (exact) names == paste0(prefix, suffix) else
      startsWith(names, paste0(prefix, "_")) & endsWith(names, suffix)
    candidates <- unique(files[keep])
    if (length(candidates) > 1L)
      stop("Ambiguous ACPC-to-MNI152 transforms for ", prefix, ":\n",
           paste(candidates, collapse = "\n"), call. = FALSE)
    candidates
  }
  if (nzchar(session)) {
    files <- list.files(file.path(subject_dir, session), recursive = TRUE,
                        full.names = TRUE)
    found <- pick(files, paste(subid, session, sep = "_"))
    if (length(found)) return(found)
  }
  files <- list.files(file.path(subject_dir, "anat"), recursive = TRUE,
                      full.names = TRUE)
  pick(files, subid, exact = TRUE)
}
