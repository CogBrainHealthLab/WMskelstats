#Functions to download, load and report streamline endpoints using reference from the HCP full connectome tractography (Elias et al. (2024) <doi:10.1038/s41597-024-03197-0>)

###############################################################################
###############################################################################
###############################################################################

#' @title Database check for WMskelstats HCP connectome streamline data
#'
#' @description Function to check whether streamline data library is available. The data is downloaded and stored in the path indicated by system.file('extdata', package='VertexWiseR'). 
#' @examples
#' if (interactive()) {
#' streamline_data_check()
#' streamline_data=load_streamline_data(file.path(system.file('extdata', package='WMskelstats'), "/templates/dTOR-985_full-connectome_MNI152-2mm_R.h5"))
#' }
#' @importFrom utils menu
#' @export
streamline_data_check=function()
{
  #Check if data is present or prompt for download
  extdatadir=system.file('extdata', package='WMskelstats')
  if (!dir.exists(paste0(extdatadir,'/templates/'))) {dir.create(paste0(extdatadir,'/templates/'))}
  dest_path=paste0(extdatadir,'/templates/dTOR-985_full-connectome_MNI152-2mm_R.h5')
  #prompt
  if (!file.exists(dest_path))
  {
    prompt = utils::menu(c("Yes", "No"), title=paste0(
      "\nThe HCP connectome streamline data was not detected inside WMskelstats's package directory (",dest_path,"). It is needed to identify streamlines crossing voxels and localise their endpoints.\n\nDo you want the streamline data (~442.3 MB) to be downloaded now?"))
    if (prompt==1) {
      
      #function to check if url exists
      #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
      #https://stackoverflow.com/a/60627969
      valid_url <- function(url_in, t = 2) {
        con <- url(url_in)
        check <- suppressWarnings(try(open.connection(con, open = "rb", timeout = t), silent = TRUE)[1])
        suppressWarnings(try(close.connection(con), silent = TRUE))
        is.null(check)}
      
      #downloader 
      #Check if URL works and avoid returning error but only print message as requested by CRAN:
      SERVER_URL="https://researchdata.ntu.edu.sg/"
      PERSISTENT_ID="doi:10.21979/N9/QIJFOZ"
      version="2.0"
      file_id="571906"
      url=paste0(SERVER_URL, "/api/access/datafile/", file_id)
      if (valid_url(url)) {
        #extend time out of 60s temporarily
        old <- options(timeout = max(3600, getOption("timeout"))) 
        on.exit(options(old), add = TRUE)
        
        #keep only a tmp file at first in case it fails
        tmp <- paste0(dest_path, ".part") #prevent partial file at dest_path
        on.exit(unlink(tmp), add = TRUE) # no-op after a successful rename
        successfully_dled <- tryCatch({
          download.file(url, 
                        destfile = tmp, 
                        mode = "wb"); 
          TRUE}, error = function(e) FALSE)
        if (!successfully_dled) stop(paste0("dTOR-985_full-connectome_MNI152-2mm_R.h5 failed to be downloaded from the dataverse repository. Please check your internet connection. Alternatively, you may visit https://doi.org/10.21979/N9/QIJFOZ and download the object manually."))
        
        #if part is entirely downloaded successfully, can adopt true name
        invisible(file.rename(tmp, dest_path))
        
      } else { 
        stop(paste0("dTOR-985_full-connectome_MNI152-2mm_R.h5 failed to be downloaded from the dataverse repository. Please check your internet connection. Alternatively, you may visit https://doi.org/10.21979/N9/QIJFOZ and download the object manually.")) 
      }
      #if user refuses, stops if required, just returns a message if optional at this stage
    } else {
      stop("\nStreamline endpoints can only be identified with the HCP connectome streamline data.\n") 
    } 
  } 
}

###############################################################################
###############################################################################
###############################################################################

#' @title Loader for WMskelstats HCP connectome streamline data
#'
#' @description This script reads the WMskelstats HCP connectome streamline data, which is a lookup table inside a .h5 file, containing the voxel coordinates of every streamline in the HCP full connectome tract file (Elias et al. (2024) <doi:10.1038/s41597-024-03197-0>), transposed on WMskelstats' MNI 152 template.
#' @returns A large list containing (see details):
#' - A 1x902,630 numeric vector corresponding to vox_ptr
#' - A 1x11,820,000 numeric vector corresponding to start_r
#' - A 1x11,820,000 numeric vector corresponding to end_r
#' - A 1x2,119,927,747 numeric vector corresponding to sids
#' @example db=load_streamline_data()
#' @details The library is a voxel-to-streamline lookup table for the dTOR-985 Human Connectome Project connectome (Elias et al. (2024) <doi:10.1038/s41597-024-03197-0>) with 11,820,000 streamlines. For each voxel of the MNI152_T1_2mm template (FSL mask), it identifies crossing streamlines, and gives their endpoints' coordinates. It was built from the dTOR-985 structural connectome full tractogram .trk file (<doi:10.6084/m9.figshare.25209947>) on the 91×109×91 MNI152 2mm grid used by WMskelstats.
#'Content: 
#'- sids: vector listing streamline IDs crossing (1-based, same numbering as the .trk), voxel per voxel 
#'- vox_ptr: Voxel index indicating for each voxel, where it starts and stops appearing in sids (helps indexing all streamlines that cross each voxel) 
#'- start_r: voxel (linear index) of each streamline’s first point. 0 means that end fell outside the grid. 
#'- end_r: The same for the last point of each streamline.
#'- grid_shape: shape of the grid for indexing reference 91 109 91.
#' @importFrom hdf5r H5File h5attr
#' @noRd

load_streamline_data <- function() 
{
  streamline_data_check()
  path=file.path(system.file('extdata', package='WMskelstats'), "/templates/dTOR-985_full-connectome_MNI152-2mm_R.h5")
  
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all())
  list(ptr   = as.numeric(f[["vox_ptr"]]$read()), #voxel indices to read sids
       start = f[["start_r"]]$read(), #for each streamline, first voxel end
       end   = f[["end_r"]]$read(), #for each streamline, last voxel end
       sids  = f[["sids"]]$read(), #streamline IDs, ordered by voxel they are crossing
       dims  = as.integer(hdf5r::h5attr(f, "grid_shape"))) #grid shape 91 109 91
}


###############################################################################
###############################################################################
###############################################################################
#' @title Streamline endpoints identifier
#' @description
#' This function allows users to enter a mask of skeleton clusters or t-stat data (e.g. in a map outputted by `TFCE_threshold()`), and to obtain a map of voxels in the same template space where streamlines crossing through the cluster end in both directions. The streamlines used as reference come from the dTOR-985 Human Connectome Project connectome (Elias et al. (2024) <doi:10.1038/s41597-024-03197-0>) with 11,820,000 streamlines, mapped on WMskelstats' MNI 152 2mm template space.
#' @param cluster_mask A binary numeric or boolean vector aligned with skeleton coordinates (e.g. model$pos_mask), which indicates the voxels of interest, through which the streamlines are crossing.
#' @param skel_coords A V x 3 matrix, containing the X, Y and Z coordinates of the skeleton data matching the `cluster_mask`. Alternatively, a list object containing $skel_coords. They will be used to cross-reference the template coordinates.
#' @param min_size An integer stating the minimum number of voxels composing the cluster(s) from `cluster_mask` for crossing streamlines to be included. Default is 20.
#' @param split_clusters A boolean object stating whether to return coordinates and endpoints maps separately for each cluster in `cluster_mask`, each in a list called cluster_N, with a coordinates matrix and a vector of endpoints values. If FALSE, the function only returns total_coords and total_values maps, which contain endpoints of all streamlines crossing any voxels in `cluster_mask`. Default is FALSE. Note that generating maps for each cluster increases memory usage (potentially a few hundreds MB when there are many clusters).
#' @param weighted A boolean object stating whether to return endpoints voxel values as the number of their streamlines crossing the cluster_map (FALSE), or as the number of streamlines multiplied by the sum of values in the cluster-voxels they crossed (TRUE). Unless `stats_map` is given, each cluster voxel value = 1. Default is FALSE.
#' @param stats_map An optional numeric vector of statistics, matching  `skel_coords` (e.g. a `thresholded_tstat_map` returned by `TFCE_threshold()`). It is used to report clusters' peak values, and as weight for the endpoints map (except if `weighted` is FALSE)
#' @param keep_skel A boolean object stating whether to keep the skeleton voxels coordinates from `skel_coords` in the output maps, so they can be plotted together with the endpoints (with NA values). Endpoint voxels that overlap with the skeleton are not affected. Default is TRUE.
#' @returns A list containing:
#' 
#' - clusters: a data frame with one row per cluster (`cluster`, `cluster_size`, `n_streamlines`, `n_endpoint_voxels` and, if `stats_map` is given, the peak's `x`, `y`, `z` and `max_value`).
#' - total_coords: a V x 3 integer matrix, in the same layout as `skel_coords`, of the voxels where the streamlines crossing the cluster(s) end.
#' - total_values: for each row of total_coords, the number of streamline ends in that voxel (NA for the skeleton voxels appended by `keep_skel`). If `weighted`==TRUE, the streamline count weighted by sum of values in crossed cluster voxels. 
#' - grid: dimensions of the template grid.
#' @importFrom mmand shapeKernel components
#' @export

get_endpoints <- function(cluster_mask, 
                          skel_coords, 
                          min_size = 20, 
                          split_clusters = FALSE,
                          weighted = FALSE,
                          stats_map = NULL,
                          keep_skel = TRUE) 
{
  #Load the template lookup table with streamline IDs and endpoint coordinates matching every voxel of the MNI 152 2mm template 
  db=load_streamline_data()
  dims <- db$dims #91 109 91
  
  #Skeleton coordinates that the clusters should correspond to
  coords <- if (is.list(skel_coords)) skel_coords$skel_coords else skel_coords
  coords <- as.matrix(coords)[, 1:3, drop = FALSE]
  storage.mode(coords) <- "integer"
  if (length(cluster_mask) != nrow(coords))
    stop("cluster_mask has ", length(cluster_mask), " values but skel_coords has ", nrow(coords), " rows.")
  if (!is.null(stats_map) && length(stats_map) != nrow(coords)) stop("stats_map must match skel_coords's", nrow(coords)," length.")
  
  #Check that the skeleton's coordinates fall within the streamline data template's own coordinates
  if (anyNA(coords) || any(coords < 1L) || any(coords[, 1] > dims[1]) ||
      any(coords[, 2] > dims[2]) || any(coords[, 3] > dims[3]))
    stop("The skeleton's coordinates fall outside the lookup grid (", paste(dims, collapse = "x"), "). The model may have been run on a different template.")
  
  #Prints the skeleton voxels that do match a streamline in the dTOR HCP full connectome (about 6.5% for our FMRIB58 2mm FA skeleton do not)
  #cov_vox <- diff(db$ptr) > 0
  #lin_sk  <- coords[, 1] + dims[1] * (coords[, 2] - 1) + dims[1] * dims[2] * (coords[, 3] - 1)
  #message(sprintf("%.1f%% of model skeleton voxels have streamlines in the lookup;.", 100 * mean(cov_vox[lin_sk])))
  
  #Turn cluster vector to 3D array, using the coordinates
  img <- array(0, dim = dims)
  img[coords] <- ifelse(is.na(cluster_mask), 0, as.numeric(cluster_mask))
  
  #Weight map: the stats_map if given, otherwise just the same cluster_mask binarised (will be sums of 1s)
  w_src  <- if (!is.null(stats_map)) stats_map else cluster_mask
  wt_img <- array(0, dim = dims)
  wt_img[coords] <- ifelse(is.na(w_src), 0, as.numeric(w_src))
  
  #Identify clusters (connectivity hard coded as 26 for now)
  connectivity=26
  if (!connectivity %in% c(6, 26)) {warning("Connectivity must be 6 or 26. Set to 26."); connectivity <- 26}
  kernel <- if (connectivity == 6) mmand::shapeKernel(c(3, 3, 3), type = "diamond") else mmand::shapeKernel(c(3, 3, 3), type = "box") 
  #assign ID to clusters based on kernel neighbours
  components <- mmand::components(img != 0, kernel = kernel) 
  tab <- table(components[components > 0]) #count voxels for each cluster ID
  ids <- as.integer(names(tab)[tab >= min_size]) #clusters above min_size
  if (!length(ids)) stop("No clusters >= min_size.")
  
  #Placeholder 3D arrays where endpoints will be written
  nv <- prod(dims)
  total <- array(0L, dim = dims); 
  total_w <- array(0, dim = dims)
  info <- list()
  by_cluster <- list() #used if split_clusters is TRUE
  
  #For each cluster kept,  
  for (cl in ids) {
    vox <- which(components == cl) #get voxels from that cluster 
    n <- db$ptr[vox + 1] - db$ptr[vox] #from the h5: use ptr (voxel's first and last index of matching streamline IDs) to count streamlines that cross through each of those voxel
    
    #prepare 3D arrays
    m  <- array(0L, dim = dims)
    mw <- array(0,  dim = dims)
    s  <- integer(0)
    
    #per-streamline weight: sum of the crossed voxels stat (from stats_map), or crossed voxels count
    if (sum(n) > 0) {
      idx <- rep(db$ptr[vox], n) + sequence(n) #sids matching the cluster voxels
      sid <- db$sids[idx] #extract database info for the matching streamlines
      wrs <- rowsum(rep(wt_img[vox], n), sid) #weighted map voxels' summed up for each relevant streamline
      
      #New indices for matching streamlines
      s   <- as.integer(rownames(wrs)); #streamline indices
      w <- wrs[, 1] #weights
      
      #get endpoints
      startpt <- db$start[s]; 
      endpt <- db$end[s]
      #drops the ~59 streamlines whose both ends fell outside the MNI152 mask
      ok <- startpt > 0 
      ends_idx  <- c(startpt[ok], endpt[ok]) #voxel indices of both ends
      ends_weights <- c(w[ok],  w[ok])  #weight at each end
      
      if (length(ends_idx)) {
        m[] <- tabulate(ends_idx, nbins = nv) #streamline count map
        rs2 <- rowsum(ends_weights, ends_idx) #weighted count (sum of voxel values)
        mw[as.integer(rownames(rs2))] <- rs2[, 1] #weighted map
      }
    }
    
    #append to all clusters maps
    total   <- total + m
    total_w <- total_w + mw
    
    #If split_clusters is TRUE, will keep cluster's voxels and counts here
    if (split_clusters) {
      cl_idx <- which(m > 0, arr.ind = TRUE)[, 1:3, drop = FALSE]
      by_cluster[[as.character(cl)]] <- list(idx = cl_idx, 
                                             counts = m[cl_idx], 
                                             counts_weighted = mw[cl_idx])
    }
    
    #Write a cluster summary similar to model output
    row <- data.frame(cluster = cl, 
                      cluster_size = length(vox),
                      n_streamlines = length(s), 
                      n_endpoint_voxels = sum(m > 0))
    if (!is.null(stats_map)) {  #peak matching cluster summary
      peakv <- vox[which.max(abs(wt_img[vox]))]
      #convert linear index back into xyz columns for skel_viewer
      lin2ijk <- function(v, dims) {
        v0 <- v - 1
        cbind(i = v0 %% dims[1] + 1,
              j = (v0 %/% dims[1]) %% dims[2] + 1,
              k = v0 %/% (dims[1] * dims[2]) + 1)}
      row[c("x", "y", "z")] <- as.list(lin2ijk(peakv, dims))
      row$max_value <- round(wt_img[peakv], 4)
    }
    info[[length(info) + 1]] <- row
  }
  
  #Lays out endpoint voxels like skel_coords (V x 3, named x, y, z), with one value per row
  layout_endpoints <- function(idx, counts) {
    coords_out <- idx
    dimnames(coords_out) <- list(NULL, c("x","y","z"))
    
    #Optionally append the skeleton voxels that are not endpoints (NA values), without touching the endpoint voxels that overlap the skeleton
  if (keep_skel) {
      lin <- function(a) a[, 1] + dims[1] * (a[, 2] - 1) + dims[1] * dims[2] * (a[, 3] - 1)
      skel_only <- coords[!(lin(coords) %in% lin(coords_out)), , drop = FALSE]
      colnames(skel_only) <- c("x","y","z")
      coords_out <- rbind(coords_out, skel_only)
      counts <- c(counts, rep(NA, nrow(skel_only)))
    }
    list(coords = coords_out, vals = counts)
  }
  
  #All clusters together: endpoint voxels are defined by the streamline count, the value reported is the count or the weighted count
  tot_idx <- which(total > 0, arr.ind = TRUE)[, 1:3, drop = FALSE]
  tot <- layout_endpoints(tot_idx, if (weighted) total_w[tot_idx] else total[tot_idx])
  
  endpoints_obj <- list(clusters     = do.call(rbind, info),
                        total_coords = tot$coords,
                        total_values = tot$vals,
                        grid         = dims)
  
  #split_clusters: for every cluster, append its own coords and (weighted) counts, named after its cluster ID
  if (split_clusters) {
    for (cl in names(by_cluster)) {
      b <- by_cluster[[cl]]
      endpoints_obj[[paste0("cluster_", cl)]] <- layout_endpoints(b$idx, if (weighted) b$counts_weighted else b$counts)
    }
  }
  
return(endpoints_obj)
}