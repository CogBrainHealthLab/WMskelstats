#' @title Demo data downloader
#'
#' @description Downloads all demo data used in the package's examples and vignettes at once.
#' @param path Path to which the demo data folder will be downloaded. Default is tempdir().
#' @param quiet A boolean stating whether to print the progress and outcome of the downloading. Default is FALSE.
#' @return A string containing the path to the downloaded demo_data folder
#' @examples
#' dl_demo()
#' @importFrom utils menu untar unzip vignette
#' @export

dl_demo=function(path=tempdir(), quiet=FALSE){
  # Folder name
  out_dir <- paste0(path,"/demo_data")
  
  # Create demo_data directory if it doesn't exist
  if (!dir.exists(out_dir)) {dir.create(out_dir)}
  
  # GitHub repository URL
  base_url <- "https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/refs/heads/main/inst/demo_data/"
  
  files_to_download <- c(
    "ds007090_behdata.tsv",
    "skel_fa_FMRIB58_skeleton_2mm_t0.2.rds",
    "skel_md_FMRIB58_skeleton_2mm_t0.2.rds")
  
  # Download each file
  for (file in files_to_download) {
    if(!file.exists(paste0(out_dir,'/',file)))
    {
      # Construct raw GitHub URL (different from the tree URL)
      raw_url <- paste0(base_url, file)
      
      # Destination path
      dest_path <- file.path(out_dir, file)
      
      tryCatch({
        # Download the file
        if (quiet==FALSE)
        {
          download.file(url = raw_url, destfile = dest_path, mode = "wb")
          cat("Downloaded:", file, "\n")
        } else
        {
          download.file(url = raw_url, destfile = dest_path, mode = "wb",
                        quiet = TRUE)
        }
        
        # Check if file is a zip file and extract it
        if (grepl("\\.zip$", file, ignore.case = TRUE)) {
          #unzip or untar as unzip struggles with long paths
          zip_result <- suppressWarnings(try(unzip(zipfile = dest_path, exdir = out_dir), silent = TRUE))
          if (is.null(zip_result) | inherits(zip_result, "try-error")==TRUE) {untar(tarfile = dest_path, exdir = out_dir)}
          
          if (quiet==FALSE){cat("Unzipped:", file, "\n")}
          unlink(dest_path) #remove zip file once unzipped
        }
      }, error = function(e) 
      {
        if (quiet==FALSE){
          cat("Failed to download:", file, "\n");
          cat("Error:", e$message, "\n")}
      })
    }
  }
  
  if (quiet==FALSE){
    cat("\nDownload and extraction process completed.\n")
    cat("\nDemo files were saved in:", normalizePath(out_dir), "\n")
  }
  
  return(normalizePath(out_dir))
}