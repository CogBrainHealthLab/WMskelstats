# Run from the WMskelstats source root:
# Rscript tools/build_qsi_native.R . /path/to/ITK/lib/cmake/ITK-5.x
# Rebuild on each target OS, R architecture and compiler toolchain.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1L]] else ".", winslash = "/", mustWork = TRUE)
itk_dir <- if (length(args) >= 2L) args[[2L]] else Sys.getenv("ITK_DIR")
if (!nzchar(itk_dir) || !file.exists(file.path(itk_dir, "ITKConfig.cmake")))
  stop("Supply the directory containing ITKConfig.cmake as argument 2 or ITK_DIR")
itk_dir <- normalizePath(itk_dir, winslash = "/", mustWork = TRUE)
if (!requireNamespace("Rcpp", quietly = TRUE)) stop("Install Rcpp first")
if (!requireNamespace("RNifti", quietly = TRUE)) stop("Install RNifti first")
cmake <- Sys.which("cmake")
if (!nzchar(cmake)) stop("Install CMake >= 3.18 and add it to PATH")
source_dir <- file.path(root, "tools", "qsi_native")
cpp_source <- file.path(root, "src", "qsi_native", "qsi_native.cpp")
if (!file.exists(cpp_source)) stop("Missing native source")
work <- file.path(root, "tools", "qsi_native_build")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
stage <- file.path(work, "qsiNative")
dir.create(file.path(stage, "src"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(stage, "R"), showWarnings = FALSE)
writeLines(c("Package: qsiNative", "Version: 0.0.1", "Title: Native QSI Build Staging",
             "Description: Temporary package used to generate registered Rcpp bindings.",
             "License: GPL-3", "Imports: Rcpp", "LinkingTo: Rcpp"), file.path(stage, "DESCRIPTION"))
writeLines(c("useDynLib(qsiNative, .registration=TRUE)", "importFrom(Rcpp, evalCpp)"),
           file.path(stage, "NAMESPACE"))
stopifnot(file.copy(cpp_source, file.path(stage, "src"), overwrite = TRUE))
Rcpp::compileAttributes(stage, verbose = TRUE)
# Link against this R installation, never a different R or Python's ITK.
lib_candidates <- c(file.path(R.home("lib"), c("libR.so", "libR.dylib", "libR.dll.a", "x64/libR.dll.a")),
                    file.path(R.home(), c("bin/x64/R.dll", "bin/R.dll")))
r_library <- lib_candidates[file.exists(lib_candidates)][1L]
if (is.na(r_library)) stop("Cannot locate this R installation's link library")
output <- file.path(root, "inst", "qsi-native")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
build <- file.path(work, "build")
cmake_args <- c("-S", source_dir, "-B", build,
                paste0("-DITK_DIR=", itk_dir), paste0("-DQSI_NATIVE_SOURCE=", cpp_source), "-DCMAKE_BUILD_TYPE=Release",
                paste0("-DR_INCLUDE_DIR=", R.home("include")),
                paste0("-DRCPP_INCLUDE_DIR=", system.file("include", package = "Rcpp")),
                paste0("-DRCPP_EXPORTS=", file.path(stage, "src", "RcppExports.cpp")),
                paste0("-DR_LIBRARY=", r_library),
                paste0("-DR_DLL_SUFFIX=", .Platform$dynlib.ext),
                paste0("-DCMAKE_INSTALL_PREFIX=", output))
if (.Platform$OS.type == "windows") {
  # ITK must also be built with this Rtools toolchain. MSVC-built ITK is incompatible.
  compiler <- Sys.which("g++")
  make <- Sys.which("make")
  c_compiler <- Sys.which("gcc")
  if (any(!nzchar(c(compiler, make, c_compiler)))) stop("Put the matching Rtools compiler and make on PATH")
  cmake_args <- c(cmake_args, "-G", "MinGW Makefiles",
                  paste0("-DCMAKE_CXX_COMPILER=", compiler),
                  paste0("-DCMAKE_C_COMPILER=", c_compiler),
                  paste0("-DCMAKE_MAKE_PROGRAM=", make))
}
run <- function(arguments) {
  status <- system2(cmake, vapply(arguments, shQuote, character(1)))
  if (status != 0L) stop("Native build failed (exit status ", status, ")")
}
run(cmake_args)
run(c("--build", build, "--config", "Release", "--parallel", "2"))
run(c("--install", build, "--config", "Release"))
dll <- file.path(output, paste0("qsiNative", .Platform$dynlib.ext))
stopifnot(file.exists(dll))
# Test the compiled exports now. Rcpp must already be loaded for its callable API.
loaded <- dyn.load(dll, local = TRUE)
for (symbol in c("qsi_apply_transforms_cpp", "qsi_image_write_cpp"))
  stopifnot(inherits(getNativeSymbolInfo(paste0("_qsiNative_", symbol), PACKAGE = loaded), "NativeSymbolInfo"))
dyn.unload(dll)
message("Native backend built and load-checked: ", dll)
message("Now install WMskelstats with R CMD INSTALL . and run tests/qsi_native_checks.R")
