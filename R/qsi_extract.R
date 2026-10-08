#' @title QSI diffusion-weighted imaging metrics extractor
#'
#' @description Extracts diffusion weighted imaging based metrics across a whole cohort datasets from QSIprep or QSIrecon pipeline outputs, masked into a common skeleton template, and merging them into single RDS files for each metric of interest.
#' @details For QSIprep output, the function makes use of tools from the `dti` R package to build a diffusion-weighted map and estimate DTI or DKI tensors. This requires bvec, bval and the dwi images to all be present in the subject directory. The map estimated from it is then coregistered to MNI 152 2mm space using the native ITK C++ backend and the ACPC-to-MNI152 transforms that QSIprep generates. 
#' For QSIrecon outputs, maps are already present in MNI 152 space, and the qsi_extract() function regrids MNI152 maps whenever their physical grid differs from the skeleton template. 
#' The FA skeleton is based on FSL's FMRIB58_FA-skeleton_1mm downsampled to 2mm.
#' @param inputdir A string object containing the path to the QSIprep or QSIrecon output dataset. For QSIrecon, specify "derivatives/qsirecon-*/" instead of the parent directory, as some files have identical suffixes and cannot be disentangled across reconstructions.
#' @param outputdir A string object containing the path of the directory where the final cohort-wise RDS will be stored for each metric (as well as metrics maps if `keep_maps` is set as TRUE). Default is 'cohort_skeletons' in the R temporary directory (tempdir()).
#'@param metrics A string object or vector of string objects containing the name (s) of the metric(s) to be estimated, in lower case. For QSIprep outputs, only the dti package's "dtiIndices/dkiIndices" metrics apply ("fa","ga","md","k1","k2","k3","mk","mk2","kaxial","kradial","fak"); for QSIrecon outputs, it can be any reconstruction output that has a "*param-\*_dwimap.nii.gz" suffix in MNI152 space (e.g., 'icvf', 'od' etc.). Default is c('fa', 'md').
#'@param skeleton_fathreshold A numerical object with the Fractional Anisotropy (FA) threshold value to apply on the template FA skeleton (FMRIB58 2mm). Default is 0.2.
#'@param dti_tensor A string object stating the tensor to be used for applicable metrics estimation ('dtiTensor' or 'dkiTensor'). Default is dtiTensor. Argument ignored for QSIrecon output.
#'@param dti_method A string object containing the method to be used for tensor-based estimations. If `dti_tensor` is 'dtiTensor', options include "nonlinear", "linear" (default), "quasi-likelihood"; if `dti_tensor` is 'dkiTensor', options include "CLLS-QP" (default), "CLLS-H", "ULLS", "QL", "NLR". Argument ignored for QSIrecon output.
#'@param dti_sigma An integer specifying the sigma value (scale parameter of the signal's distribution) to be used as part of the tensor estimation. Default is NULL. Argument ignored for QSIrecon output. 
#'@param dti_L An integer specifying the effective degrees of freedom for the tensor estimation. Default is 1.  Argument ignored for QSIrecon output.
#'@param nthread Number of CPU threads for tensor estimation and ITK work units for image resampling.
#'@param keep_maps A logical object to determine whether files such as estimated tensor maps and coregistered maps are to be written in the `outputdir`. Default is FALSE.
#'@param qsiprep_path A string containing the path to the QSIprep output (Optional). Ignored if inputdir and qsiprep_path are the same. Its purpose is for QSIrecon processing to retrieve ACPC-to-MNI152 transformation matrices if MNI152 coregistration was not done by QSIrecon. The QSIprep folder must have the same subjects and sessions as the the QSIrecon output.
#'@param silent A logical object to determine whether messages will be silenced. Default is FALSE.
#'
#' @returns A list of 2D matrices, each matrix corresponding to one metric from 
#' `metrics`. Each element (skel_matrices$fa, skel_matrices$md, skel_matrices$ga, ...) is its own separate matrix: rows = subjects (and sessions), columns = voxels. Additionally, the list contains the coordinates of the skeleton voxels (skel_coords matrix), the skeleton template they are based on, and the FA threshold selected. 
#' 
#' @importFrom dti readDWIdata dtiTensor dkiTensor dtiIndices dkiIndices setmask
#' @importFrom RNifti asNifti readNifti writeNifti pixdim
#' @export 

qsi_extract=function(inputdir,
                     outputdir, 
                     metrics=c('fa', 'md'),
                     skeleton_fathreshold=0.2,
                     dti_tensor='dtiTensor', 
                     dti_method, 
                     dti_sigma=NULL,
                     dti_L=1, 
                     nthread=4, 
                     keep_maps=FALSE,
                     qsiprep_path=NULL,
                     silent=FALSE){
  
  #if silent is TRUE: will silence all dti package functions/system prints
  if (silent) {
    shush <- file(nullfile(), open = "wb")
    sink(shush, type = "output")
    #if function breaks, disable
    on.exit({ sink(type = "output"); close(shush) }, add = TRUE)
  }
  
  #Output directory
  if (missing("outputdir")) {
    warning(paste0('No outputdir argument was given. The matrix objects will be saved in a directory named "cohort_skeletons" inside the R temporary directory (tempdir()).\n'))
    outputdir=file.path(tempdir(), 'cohort_skeletons')
  } else {
    dir.create(outputdir, showWarnings=FALSE, recursive=TRUE)
  }
  
  dir.create(outputdir, showWarnings=FALSE, recursive=TRUE)
  if (!dir.exists(outputdir)) stop("Cannot create outputdir")

  #Preload skeleton template
  template='FMRIB58_FA-skeleton_2mm'
  skeleton_template=RNifti::readNifti(paste0(system.file('extdata',package='WMskelstats'),'/templates/', template))
  #Premake thresholded skeleton mask 
  skeleton_mask = skeleton_masker(skeleton_template=skeleton_template, 
                                  skeleton_fathreshold=skeleton_fathreshold)
  #reorder coordinates for later use (will also be reordered for subject data)
  skel_coords=skeleton_mask[[2]] 
  skel_coords=skel_coords[order(skel_coords[,'x'], skel_coords[,'y'], skel_coords[,'z']), , drop=FALSE]
  skeleton_mask[[2]]=skel_coords 
  
  #The native ITK DLL is loaded lazily by qsi_apply_transforms().
  
  #save metadata including template, threshold, skeleton mask coordinates in a list to be appended for later rebuild
  metadata=list(skel_coords,template, skeleton_fathreshold)
  names(metadata)=c('skel_coords','skel_template','skel_threshold')
  
  #clear dtiDataobj in case these are also in the environment
  if (exists('dtioutput', inherits = FALSE)) { remove(dtioutput, dtiDataobj) }
  if (exists('dtiTensorobj', inherits = FALSE)) { remove(dtiTensorobj) }
  if (exists('dkiTensorobj', inherits = FALSE)) { remove(dkiTensorobj) }
  
  #prepare grand skeleton list (will contain cohort matrices for each metrics) 
  skel_list <- setNames(vector("list", length(metrics)),
                          paste0("skel_", metrics))
    
  # Discover subject directories, then keep every session's files separate.
  subject_dirs <- list.dirs(inputdir, recursive = FALSE, full.names = TRUE)
  subject_dirs <- subject_dirs[grepl("^sub-", basename(subject_dirs))]
  sublist <- basename(subject_dirs)

  if (!is.null(qsiprep_path) &&
      (length(qsiprep_path) != 1L || is.na(qsiprep_path) || !dir.exists(qsiprep_path))) {
    stop("qsiprep_path must be NULL or one existing directory")
  }

  for (subid in sublist)
  {
    subject_dir <- file.path(inputdir, subid)
    session_dirs <- list.dirs(subject_dir, recursive = FALSE, full.names = TRUE)
    session_dirs <- session_dirs[grepl("^ses-", basename(session_dirs))]
    if (!length(session_dirs)) session_dirs <- subject_dir

    for (session_dir in session_dirs)
    {
      session_name <- if (identical(session_dir, subject_dir)) "" else basename(session_dir)
      sub_s <- if (nzchar(session_name)) paste(subid, session_name, sep = "_") else subid
      subfiles <- list.files(session_dir, recursive = TRUE, full.names = TRUE)
      if (!length(subfiles)) {
        warning("No files found for ", sub_s, ". Skipping.")
        next
      }

      subfiles_qsiprep <- character()
      if (!is.null(qsiprep_path)) {
        prep_dir <- file.path(qsiprep_path, subid)
        if (nzchar(session_name)) prep_dir <- file.path(prep_dir, session_name)
        subfiles_qsiprep <- list.files(prep_dir, recursive = TRUE, full.names = TRUE)
      }
      if(!silent){message("\nProcessing ", sub_s,"...")}
      
      for (m in metrics)
      {
        if(!silent){message(paste0( " Metric: ", m))}
        #############################
        #create map from DWI file for QSIprep outputs
        #If map not already computed (QSIPREP), compute if applicable
        metric_map=startsWith(basename(subfiles), paste0(sub_s, "_")) & grepl(paste0("_space-[^_]+_model-.*_param-", m,"_dwimap\\.nii(\\.gz)?$"),subfiles)
        if (length(which(metric_map)) == 0)
        {
          if(!exists('dtiDataobj')){if(!silent){message(paste0("  => No preexisting map found, trying to build a dti object..."))}}
          
          #if metric is not computable by dti package, skip
          if (! m %in% c("fa","ga","md","k1","k2","k3","mk","mk2","kaxial","kradial","fak")){
            if(!silent){message(paste0(
            "  /!\\ ", m," cannot be computed here. Options are: 
    - fa,ga,md (dtiIndices/dkiIndices), 
    - k1,k2,k3,mk,mk2,kaxial,kradial,fak (dkiIndices). 
    Skipping"))}
            next
          }
          
          #############################
          #making the dti object (dtiData_make function)
          #if dtiDataobj has been created already, skip: will be reused across metrics
          #and cleared before next subject
          if(!exists('dtioutput')){
            if(!silent){message("  => Fetching individual DWI data ...")}
            dtioutput=dtiData_make(sub_s, subfiles, silent,
                                      dwi_dir = file.path(session_dir, "dwi"))
            dtiDataobj=dtioutput[[1]]
            dwivol=dtioutput[[2]] #will be reused later for coreg
          }
          #if dtiDataobj has been attempted to be created, but failed, skip
          if (!inherits(dtiDataobj, "dtiData")) {
            if(!silent){warning(paste0('No DTI/DKI map could be computed for',sub_s,', as either bval, bvec or brain_mask files are missing'))}
            next
          } 
          
          #############################
          #computing metrics
          #dti has two algorithm for metrics computation
          #if already computed, reuse the tensor (all relevant metrics are available)
          
          if (dti_tensor=='dtiTensor' & !exists('dtiTensorobj')) {
            #DTI
            if(!silent){message(paste0("  => Computing diffusion tensor using ", 
                                       dti_tensor, "..."))}
            
            if (missing(dti_method)){dti_method=c("linear")}
            dtiTensorobj  <- dti::dtiTensor(dtiDataobj, method=dti_method, 
                                            L=dti_L, sigma=dti_sigma, 
                                            mc.cores = nthread)
            Indicesobj <- dti::dtiIndices(dtiTensorobj, 
                                          mc.cores = nthread) 
          } else if (dti_tensor=='dkiTensor' & !exists('dkiTensorobj')) {
            #DKI
            if(!silent){message(paste0("  => Computing diffusion kurtosis tensor (and diffusion tensor)  using ", dti_tensor, "..."))}
            
            if (missing(dti_method)){dti_method=c("CLLS-QP")}
            dkiTensorobj  <- dti::dkiTensor(dtiDataobj, method=dti_method, 
                                            L=dti_L, sigma=dti_sigma, 
                                            mc.cores = nthread) 
            Indicesobj <- dti::dkiIndices(dkiTensorobj, 
                                          mc.cores = nthread)
            
          } else if (dti_tensor!='dtiTensor' & dti_tensor!='dkiTensor') 
          {stop('The dti_tensor argument must either be dtiTensor or dkiTensor')
          } else if (exists('dtiTensorobj') | exists('dkiTensorobj')){
            if(!silent){message(paste0("  => Reusing previously computed tensor..."))}
          }
          
          if (! m %in% slotNames(Indicesobj))
          {warning(paste0('  => ', m,' did not get outputted in the Indices. It may be an issue with the dti package.
                            Skipping')); next}  
          
          #############################
          #create map
          metricmap <- slot(Indicesobj, m) #3D array
          niivol <- RNifti::asNifti(metricmap, reference = dwivol)
          
          #save map if wanted
          if(keep_maps){
            mapdir=paste0(outputdir,'/',m,'_maps')
            dir.create(mapdir, showWarnings=FALSE)
            if(!silent){message(paste0("  => Writing metrics map to ",mapdir))}
            mapfile=file.path(mapdir, paste0(sub_s,"_",m,"_map.nii.gz"))
            RNifti::writeNifti(niivol, mapfile)
          } 
          mapfile=niivol
          
          #############################
          #coregister map to MNI 152 using QSIprep's transforms
          if(!silent){message("  => Coregistering metrics map to MNI152NLin2009cAsym...")}
          #looking for transformation matrix, either in default path 
          #or in QSIprep path if specified
          transform_path <- .qsi_find_transform(inputdir, subid, session_name)
          if (!length(transform_path) && !is.null(qsiprep_path)) {
            transform_path <- .qsi_find_transform(qsiprep_path, subid, session_name)
          }
          #if still not found, skip
          if(length(transform_path)==0)
          { if (!silent){message(paste0("  No valid transformation matrix found for", sub_s, ", ('*_from-ACPC_to-MNI152NLin2009cAsym_mode-image_xfm.h5')."))}
            break
          }
          
          #coregister with the transform_path found
          finalmap=ACPC_to_MNI152(mapfile, transform_path, keep_maps=keep_maps,
                                   fixed=skeleton_template, outputdir=outputdir,
                                   m=m, sub_s=sub_s, nthread=nthread)
        
        } else {
          
          #check if QSIrecon map exists in MNI152
          metric_map_MNI152=startsWith(basename(subfiles), paste0(sub_s, "_")) & grepl(paste0("_space-MNI152NLin2009cAsym_model-.*_param-", m,"_dwimap\\.nii(\\.gz)?$"),subfiles)
          
          #if only ACPC, coregister
          if(length(which(metric_map_MNI152))==0)
          {
            #############################
            #coregister ACPC maps to MNI 152 using QSIprep's transforms
            if(!silent){message("  => Coregistering metrics map to MNI152NLin2009cAsym...")}
            #looking for transformation matrix, either in default path 
            #or in QSIprep path if specified
            transform_root <- if (is.null(qsiprep_path)) inputdir else qsiprep_path
            transform_path <- .qsi_find_transform(transform_root, subid, session_name)
            #if not found, skip
            if(length(transform_path)==0)
            {
              if (!silent){ 
              {message(paste0("  No valid transformation matrix found for", sub_s, ", ('*_from-ACPC_to-MNI152NLin2009cAsym_mode-image_xfm.h5'), even in the given qsiprep_path."))}}
              break
            }
              
            #coregister with the transform_path found
            finalmap=ACPC_to_MNI152(subfiles[metric_map & grepl("_space-ACPC_", subfiles)], transform_path,
                                    keep_maps=keep_maps, fixed=skeleton_template,
                                    outputdir=outputdir, m=m, sub_s=sub_s, nthread=nthread)
            
          } else {
            ####################################
            #already in MNI152 for QSIrecon maps
            
            if(!silent){message("  => Using preexisting map:");
                        message(paste0("    ", basename(subfiles[which(metric_map_MNI152==TRUE)])))}
            mapfile_coreg=subfiles[which(metric_map_MNI152==TRUE)]
            if (length(mapfile_coreg) != 1L) stop("Multiple MNI maps found for ", sub_s, " / ", m)
            orig_vol <- RNifti::readNifti(mapfile_coreg)
            
            
            #Resample whenever the complete physical grid differs, including shifted 2mm grids
            if (!.qsi_same_grid(orig_vol, skeleton_template)) {
              if(!silent){message("  => Resampling metrics map to the skeleton grid...")}
                resampled_vol <- qsi_apply_transforms(
                fixed = skeleton_template, #exact target grid
                moving = orig_vol,
                interpolator = "linear",
                transformlist = list(),  #identity transform in physical coordinates
                nthread = nthread
              )
              
              #save to dedicated folder if needed
              if(keep_maps){
                mapdirmni152=file.path(outputdir, paste0(m,'_maps_MNI152'))
                dir.create(mapdirmni152, showWarnings=FALSE)
                mapfile_coreg=file.path(mapdirmni152,paste0(sub_s,"_",m,"_map_MNI152.nii.gz"))
                qsi_image_write(resampled_vol,  file.path(mapdirmni152,paste0(sub_s,"_",m,"_map_MNI152.nii.gz")))
              }
              finalmap=resampled_vol #either way
    
            } else {
              #If already on the same physical grid, use file directly
              #save to dedicated folder if needed
              if(keep_maps){
                mapdirmni152=file.path(outputdir, paste0(m,'_maps_MNI152'))
                dir.create(mapdirmni152, showWarnings=FALSE)
                if(!silent){message(paste0("  => Copying map to ", mapdirmni152))}
                RNifti::writeNifti(orig_vol, file.path(mapdirmni152,paste0(sub_s,"_",m,"_map_MNI152.nii.gz")))
              }
              finalmap=orig_vol
              
            }
          }
        }
        
        ####################################
        if(!silent){message("  => Extracting values using the FMRIB58 FA 2mm template skeleton...")}
        #Extract skeleton of the map for each metric separately
        #safeguard
        if (!.qsi_same_grid(finalmap, skeleton_template)) stop("Final map and skeleton have different physical grids")
        metrics_array=as.array(finalmap)
        if(!identical(dim(metrics_array), dim(skeleton_mask[[1]]))){
        stop("The FA skeleton template does not share the subject's map dimensions. The downsampling to 2mm may have failed.")}
        #Get subject values in the template skeleton mask
        #vectorise values and give it the name of subject/ses
        subj_skeleton <- matrix(metrics_array[skeleton_mask[[2]]], nrow = 1,
                                dimnames = list(sub_s, NULL))
        skel_list[[paste0('skel_',m)]][[sub_s]] = subj_skeleton
      }
      
      #clear dtiDataobj of that subject if created
      if (exists('dtioutput', inherits = FALSE)) { remove(dtioutput, dtiDataobj) }
      if (exists('dtiTensorobj', inherits = FALSE)) { remove(dtiTensorobj) }
      if (exists('dkiTensorobj', inherits = FALSE)) { remove(dkiTensorobj) }
      
    }
  }
  ####################################
  #set to individual matrices, return the grand list
  #empty metrics matrices will be removed
  if(!silent){
    message("Metrics without applicable subjects removed: ",
            paste(names(skel_list)[sapply(skel_list, is.null)], collapse = ", "))
  }
  skel_list <- skel_list[!sapply(skel_list, is.null)]
  #turn to matrices
  skel_matrices = lapply(names(skel_list), function(m) do.call(rbind, skel_list[[m]])) 
  names(skel_matrices) = names(skel_list)
  
  #save to outputdir separately per metric
  skeldir=paste0(outputdir,'/cohort_metrics_skeletons')
  dir.create(skeldir, showWarnings=FALSE)
  for (skel_m in names(skel_matrices)){
    #define file name
    file_path=paste0(skeldir,'/',skel_m,'_FMRIB58_skeleton_2mm_t', skeleton_fathreshold,'.rds')
    #append metadata
    rds_file=append(list(skel_matrices[[skel_m]]),metadata)
    names(rds_file)[1]=skel_m
    #Save
    saveRDS(object=rds_file, 
            file = file_path)
    if(!silent){message(paste0("\u2713 Final ", skel_m, " cohort matrix saved to", file_path))}
  }
  return(append(skel_matrices,metadata))
}

#################################################################################
#################################################################################
#################################################################################

#' @title dtiData maker
#'
#' @description Function to create a dtiData class object from the dti package with the right files (automatically fetched), which can then be used to compute tensor or DKI metrics if needed.
#' @param sub_s A string indicating the subject ID and their session if applicable (e.g. "sub-0001_ses-1"). Will be used to read inside BIDS-formatted file names.
#' @param subfiles A list of files, with full path names, from a QSIprep output directory.
#' @param silent Whether to print messages and warnings or not. Default is FALSE.
#' @importFrom dti readDWIdata setmask
#' @noRd

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
  btable_candidates <- paste0(stem, c(".b_table", ".b_table.txt"))
  btable_existing <- btable_candidates[file.exists(btable_candidates)]
  btable_file <- if (length(btable_existing)) btable_existing[1L] else btable_candidates[1L]
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

#################################################################################
#################################################################################
#################################################################################

#' @title Fractional anisotropy skeleton masker
#'
#' @description Function to get mask of a template fractional anisotropy (FA) skeleton at a desired threshold. The skeleton is based on FSL's FMRIB58_FA-skeleton_1mm and downsampled to 2mm.
#' @param sub_s A string indicating the subject ID and their session if applicable (e.g. "sub-0001_ses-1"). Will be used to read inside BIDS-formatted file names.
#' @param metrics_map A metrics map in MNI 152 2mm space
#' @param skeleton_fathreshold A numerical object with the (Fractional Anisotropy) threshold value with which to apply the template skeleton. Default is 0.2.
#' @param skeleton_template A string object naming the template FA skeleton (currently, only 'FMRIB58_FA-skeleton_2mm') or a niftiImage object of the template loaded with RNifti. It is optional as qsi_extract already preloads it when running skeleton_masker, but the latter function can load it by on its own if needed.
#' @param silent Whether to print messages and warnings or not. Default is FALSE.
#' @importFrom RNifti readNifti
#' @examples
#' skeleton_mask=skeleton_masker(skeleton_fathreshold=0.2) 
#' @export 

skeleton_masker=function(skeleton_template, skeleton_fathreshold=0.2){
  
  #load template skeleton if not provided
  if(missing(skeleton_template)){
    skeleton_template=RNifti::readNifti(paste0(system.file('extdata',package='WMskelstats'),'/templates/FMRIB58_FA-skeleton_2mm.nii'))
  } else if (!inherits(skeleton_template,'niftiImage'))
  {
    skeleton_template=RNifti::readNifti(paste0(system.file('extdata',package='WMskelstats'),'/templates/',skeleton_template,'.nii'))
  } 
  
  #make binary mask out of skeleton and threshold based on FA values: 
  skeleton_bin=array(0L, dim=dim(skeleton_template))
  thresh=skeleton_fathreshold*10000
  skeleton_bin[skeleton_template>=thresh]=1L
  
  #keep voxel coordinates of the mask for later rebuild
  skeleton_bin_coords <- which(skeleton_bin == 1, arr.ind = TRUE)
  #LAS coordinates for FMRIB (see RNifti::orientation(skeleton_template):
  colnames(skeleton_bin_coords)=c('x','y','z')
  
  return(list(skeleton_bin,skeleton_bin_coords))
}


#################################################################################
#################################################################################
#################################################################################

# Applies a precomputed spatial transform; this does not estimate registration.
ACPC_to_MNI152=function(mapfile, transform_path, qsiprep_path=NULL,
                        keep_maps=FALSE, fixed=NULL, outputdir=NULL,
                        m=NULL, sub_s=NULL, nthread=1L){
  if (length(transform_path) != 1L || !file.exists(transform_path))
    stop("Expected exactly one existing ACPC-to-MNI152 composite transform")
  if (is.null(fixed)) {
    fixed=system.file("extdata", "templates", "FMRIB58_FA-skeleton_2mm.nii",
                      package="WMskelstats")
    if (!nzchar(fixed)) stop("FMRIB58 skeleton reference is missing")
  }
  warped_vol=qsi_apply_transforms(fixed=fixed, moving=mapfile,
                                  imagetype=0L, interpolator="linear",
                                  transformlist=list(transform_path), nthread=nthread)
  if (keep_maps) {
    if (is.null(outputdir) || is.null(m) || is.null(sub_s))
      stop("outputdir, m and sub_s are required when keep_maps=TRUE")
    mapdir=file.path(outputdir, paste0(m, "_maps_MNI152"))
    dir.create(mapdir, showWarnings=FALSE, recursive=TRUE)
    qsi_image_write(warped_vol, file.path(mapdir, paste0(sub_s,"_",m,"_map_MNI152.nii.gz")))
  }
  warped_vol
}

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
