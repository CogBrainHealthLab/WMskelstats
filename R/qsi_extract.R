#' @title QSI diffusion-weighted imaging metrics extractor
#'
#' @description Extracts diffusion weighted imaging based metrics across a whole cohort datasets from QSIprep or QSIrecon pipeline outputs, masked into a common skeleton template, and merging them into single RDS files for each metric of interest.
#' @details For QSIprep output, the function makes use of tools from the `dti` R package to build a diffusion-weighted map and estimate DTI or DKI tensors. This requires bvec, bval and the dwi images to all be present in the subject directory. The map estimated from it is then coregistered to MNI 152 2mm space using the `rpyANTs` package and the ACPC-to-MNI152 transforms that QSIprep generates. 
#' For QSIrecon outputs, maps are already present in MNI 152 space, and the qsi_extract() function will only regrid the MNI152 maps to 2 mm if resolution is higher. 
#' The FA skeleton is based on FSL's FMRIB58_FA-skeleton_1mm downsampled to 2mm.
#' @param inputdir A string object containing the path to the QSIprep or QSIrecon output dataset. For QSIrecon, specify "derivatives/qsirecon-*/" instead of the parent directory, as some files have identical suffixes and cannot be disentangled across reconstructions.
#' @param outputdir A string object containing the path of the directory where the final cohort-wise RDS will be stored for each metric (as well as metrics maps if `keep_maps` is set as TRUE). Default is 'cohort_skeletons' in the R temporary directory (tempdir()).
#'@param metrics A string object or vector of string objects containing the name (s) of the metric(s) to be estimated, in lower case. For QSIprep outputs, only the dti package's "dtiIndices/dkiIndices" metrics apply ("fa","ga","md","k1","k2","k3","mk","mk2","kaxial","kradial","fak"); for QSIrecon outputs, it can be any reconstruction output that has a "*param-\*_dwimap.nii.gz" suffix in MNI152 space (e.g., 'icvf', 'od' etc.). Default is c('fa', 'md').
#'@param skeleton_fathreshold A numerical object with the Fractional Anisotropy (FA) threshold value to apply on the template FA skeleton (FMRIB58 2mm). Default is 0.2.
#'@param dti_tensor A string object stating the tensor to be used for applicable metrics estimation ('dtiTensor' or 'dkiTensor'). Default is dtiTensor. Argument ignored for QSIrecon output.
#'@param dti_method A string object containing the method to be used for tensor-based estimations. If `dti_tensor` is 'dtiTensor', pptions include "nonlinear" (default), "linear", "quasi-likelihood"; if `dti_tensor` is 'dkiTensor', options include "CLLS-QP" (default), "CLLS-H", "ULLS", "QL", "NLR". Argument ignored for QSIrecon output.
#'@param dti_sigma An integer specifying the sigma value (scale parameter of the signal's distribution) to be used as part of the tensor estimation. Default is NULL. Argument ignored for QSIrecon output. 
#'@param dti_L An integer specifying the effective degrees of freedom for the tensor estimation. Default is 1.  Argument ignored for QSIrecon output.
#'@param nthread Number of CPU threads for the dti package to use when estimating the tensor and metrics from the data. Argument ignored for QSIrecon output.
#'@param keep_maps A logical object to determine whether files such as estimated tensor maps and coregistered maps are to be written in the `outputdir`. Default is FALSE.
#'@param silent A logical object to determine whether messages will be silenced. Default is FALSE.
#'
#' @returns A list of 2D matrices, each matrix corresponding to one metric from 
#' `metrics`. Each element (skel_matrices$fa, skel_matrices$md, skel_matrices$ga, ...) is its own separate matrix: rows = sub_ses, columns = voxels. Additionally, the list contains the coordinates of the skeleton voxels (skel_coords matrix), the skeleton template they are based on, and the FA threshold selected. 
#' 
#' @examples
#' SCMvextract(sdirpath = "subcortexmesh_output_metrics", 
#' outputdir=paste0(tempdir(), "\\subcortices"), template='fsaverage', measure="surfarea") 
#' @importFrom dti readDWIdata dtiTensor dkiTensor dtiIndices dkiIndices setmask
#' @importFrom rpyANTs load_ants ants_apply_transforms
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
    outputdir=paste0(tempdir(),'\\cohort_skeletons')
  } else {
    dir.create(paste0(outputdir), showWarnings=FALSE)
  }
  
  #Preload skeleton template
  template='FMRIB58_FA-skeleton_2mm'
  skeleton_template=RNifti::readNifti(paste0(system.file('extdata',package='WMskelstats'),'/templates/', template))
  #Premake thresholded skeleton mask 
  skeleton_mask = skeleton_masker(skeleton_template=skeleton_template, 
                                  skeleton_fathreshold=skeleton_fathreshold)
  #save metadata including template, threshold, skeleton mask coordinates in a list to be appended for later rebuild
  metadata=list(skeleton_mask[[2]],template, skeleton_fathreshold)
  names(metadata)=c('skel_coords','skel_template','skel_threshold')
  
  #prepare grand skeleton list (will contain cohort matrices for each metrics) 
  skel_list <- setNames(vector("list", length(metrics)),
                          paste0("skel_", metrics))
    
  #subject list
  sublist=list.files(path = inputdir, recursive = F)
  sublist=unique(stringr::str_extract(sublist, "sub-[^/]+"))
  sublist=sublist[!is.na(sublist)]
  
  for (subid in sublist)
  {
    #get all files for that subject
    subdirs=list.dirs(path=paste0(inputdir,'/',subid), recursive = FALSE,
                      full.names = TRUE)
    subfiles= list.files(path = subdirs, recursive = TRUE, full.names = TRUE)
    if (length(subfiles)==0) {warning(paste0('No files found for ',subid,'. Skipping')); next}
    
    #if ses- directories present, compute both ses separately
    if (length(which(grepl('ses-', basename(subdirs), ignore.case = TRUE)))>0)
    {
      subdirs=subdirs[grepl('ses-', basename(subdirs), ignore.case = TRUE)]
      subses=paste0(subid,'_',basename(subdirs))
    } else
    {
      subses=subid  
    }
    
    for (sub_s in subses)
    {
      if(!silent){message("Processing ", sub_s,"...")}
      
      for (m in metrics)
      {
        if(!silent){message(paste0( " Metric: ", m))}
        #############################
        #create map from DWI file for QSIprep outputs
        #If map not already computed (QSIPREP), compute if applicable
        metric_map=grepl(paste0(sub_s,"_space-MNI152NLin2009cAsym_model-.*_param-", m,"_dwimap\\.nii(\\.gz)?$"),subfiles)
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
            dtioutput=dtiData_make(sub_s, subfiles, silent)
            dtiDataobj=dtioutput[[1]]
            dwivol=dtioutput[[2]] #will be reused later for coreg
          }
          #if dtiDataobj has been attempted to be created, but failed, skip
          if (!inherits(dtiDataobj, "dtiData")) {next} 
          
          #############################
          #computing metrics
          #dti has two algorithm for metrics computation
          
          if (dti_tensor=='dtiTensor') {
            #DTI
            if(!silent){message(paste0("  => Computing diffusion tensor using ", 
                                       dti_tensor, "..."))}
            
            if (missing(dti_method)){dti_method=c("nonlinear")}
            dtiTensorobj  <- dti::dtiTensor(dtiDataobj, method=dti_method, 
                                            L=dti_L, sigma=dti_sigma, 
                                            mc.cores = setCores(nthread,reprt = FALSE))
            Indicesobj <- dti::dtiIndices(dtiTensorobj, 
                                          mc.cores = setCores(nthread,reprt = FALSE)) 
          } else if (dti_tensor=='dkiTensor') {
            #DKI
            if(!silent){message(paste0("  => Computing diffusion kurtosis tensor (and diffusion tensor)  using ", dti_tensor, "..."))}
            
            if (missing(dti_method)){dti_method=c("CLLS-QP")}
            dkiTensorobj  <- dti::dkiTensor(dtiDataobj, method=dti_method, 
                                            L=dti_L, sigma=dti_sigma, 
                                            mc.cores = setCores(nthread,reprt = FALSE)) 
            Indicesobj <- dti::dkiIndices(dkiTensorobj, 
                                          mc.cores = setCores(nthread,reprt = FALSE))
            
          } else {stop('The dti_tensor argument must either be dtiTensor or 
                       dkiTensor')}
          
          if ((! m %in% slotNames(Indicesobj)) & silent==FALSE)
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
            mapfile=paste0(mapdir,"\\",sub_s,"_",m,"_map.nii.gz")
            RNifti::writeNifti(niivol, mapfile)
          } 
          mapfile=niivol
          
          #############################
          #coregister map to MNI 152 using QSIprep's transforms
          if(!silent){message("  => Coregistering metrics map to MNI152NLin2009cAsym...")}
          #Check if MNI152 is available
          #$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz but downloadable from git
          mni152_2mmvol= paste0(system.file('extdata',package='WMskelstats'),'/templates/MNI152_T1_2mm.nii.gz')
          if (!file.exists(mni152_2mmvol)){
            prompt = utils::menu(c("Yes", "No"), title=paste0(
              "\nThe MNI 152 template (2mm) was not detected inside the package directory (", paste0(system.file('extdata',package='WMskelstats'),'/MNI152_T1_2mm.nii.gz'), "). It is needed for coregistration. It can be downloaded from github.\n\nDo you want the template (~1.34 MB) to be downloaded now?"))
            if (prompt==1) {
              #function to check if url exists
              #courtesy of Schwarz, March 11, 2020, CC BY-SA 4.0:
              #https://stackoverflow.com/a/60627969
              valid_url <- function(url_in,t=2){
                con <- url(url_in)
                check <- suppressWarnings(try(open.connection(con,open="rt",timeout=t),silent=TRUE)[1])
                suppressWarnings(try(close.connection(con),silent=TRUE))
                ifelse(is.null(check),TRUE,FALSE)}
              
              #Check if URL works and avoid returning error but only print message as requested by CRAN:
              url="https://raw.githubusercontent.com/CogBrainHealthLab/WMskelstats/main/inst/extdata/templates/MNI152_T1_2mm.nii.gz"
              if(valid_url(url)) {
                download.file(url=url,destfile = paste0(system.file(package='VertexWiseR'),'/extdata/templates/MNI152_T1_2mm.nii.gz'))
              }  else { 
                stop("The template failed to be downloaded from the github repository. Please check your internet connection. Alternatively, you may visit https://github.com/CogBrainHealthLab/WMskelstats/tree/main/inst/extdata/templates and download the file manually.") #ends function
              }
            } else {
            stop("Coregistration cannot be done without the MNI 152 template.\n")}
          }
          
          #Uses rpyANTs (python version of ANTs read via reticulate in R, instead of the R version that requires ITK compiling, which can fail)
          transform_path=grep(paste0(sub_s,"_from-ACPC_to-MNI152NLin2009cAsym_mode-image_xfm.h5"), subfiles, value = TRUE)
          warped_vol <- rpyANTs::ants_apply_transforms(
            fixed = mni152_2mmvol, 
            moving = mapfile,
            imagetype = 0,
            transformlist = list(transform_path)
          )
          ants <- rpyANTs::load_ants()
          
          #save coregistered map if needed
          if(keep_maps){
            mapdirmni152=paste0(outputdir,'\\',m,'_maps_MNI152')
            dir.create(mapdirmni152, showWarnings=FALSE)
            mapfile_coreg=paste0(mapdirmni152,"\\",sub_s,"_",m,"_map_MNI152.nii.gz")
            coreg=ants$image_write(warped_vol, mapfile_coreg)
          } 
          finalmap=warped_vol
        
        } else {
          
          ####################################
          #already in MNI152 for QSIrecon maps
          
          if(!silent){message("  => Using preexisting map:");
                      message(paste0("    ", basename(subfiles[which(metric_map==TRUE)])))}
          mapfile_coreg=subfiles[which(metric_map==TRUE)]
          orig_vol <- RNifti::readNifti(mapfile_coreg)
          
          #downsample it to 2mm, if higher resolution
          if (any(RNifti::pixdim(orig_vol) < 2)) {
            if(!silent){message("  => Downsampling metrics map to 2 mm...")}
              resampled_vol <- rpyANTs::ants_apply_transforms(
              fixed = skeleton_template, #same grid as MNI 152 so no dl needed
              moving = orig_vol,
              interpolator = "linear",
              transformlist = list()   #no transform needed as same grid
            )
            ants <- rpyANTs::load_ants()
            
            #save to dedicated folder if needed
            if(keep_maps){
              mapdirmni152=paste0(outputdir,'\\',m,'_maps_MNI152')
              dir.create(mapdirmni152, showWarnings=FALSE)
              mapfile_coreg=paste0(mapdirmni152,"\\",sub_s,"_",m,"_map_MNI152.nii.gz")
              finalmap=ants$image_write(resampled_vol,  paste0(mapdirmni152,"\\",sub_s,"_", m,"_map_MNI152.nii.gz"))
            }
            finalmap=resampled_vol #either way
  
          } else {
            #If already 2 mm, use file directly
            #save to dedicated folder if needed
            if(keep_maps){
              mapdirmni152=paste0(outputdir,'\\',m,'_maps_MNI152')
              dir.create(mapdirmni152, showWarnings=FALSE)
              if(!silent){message(paste0("  => Copying map to ", mapdirmni152))}
              RNifti::writeNifti(orig_vol, paste0(mapdirmni152,"\\",sub_s,"_", m,"_map_MNI152.nii.gz"))
            }
            finalmap=orig_vol
            
          }
        }
        
        ####################################
        if(!silent){message("  => Extracting values using the FMRIB58 FA 2mm template skeleton...")}
        #Extract skeleton of the map for each metric separately
        #safeguard
        metrics_array=finalmap[]
        if(!identical(dim(metrics_array), dim(skeleton_mask[[1]]))){
        stop("The FA skeleton template does not share the subject's map dimensions. The downsampling to 2mm may have failed.")}
        #Get subject values in the template skeleton mask
        #vectorise values and give it the name of subject/ses
        subj_skeleton=matrix(metrics_array[skeleton_mask[[1]]==1], 
                           nrow=1, dimnames=list(sub_s, NULL))
        skel_list[[paste0('skel_',m)]][[sub_s]] = subj_skeleton
      }
      
      #clear dtiDataobj of that subject if created
      if (exists('dtioutput', inherits = FALSE)) { remove(dtioutput, dtiDataobj) }
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

dtiData_make=function(sub_s, 
                      subfiles, 
                      silent=FALSE){
  
  #If enough data to compute DTI/DKI map, do it
  if (length(grep(paste0(sub_s,"_space-ACPC_desc-preproc_dwi.nii*"),
                  subfiles, value=TRUE))>0 & 
      length(grep(paste0(sub_s,"_space-ACPC_desc-preproc_dwi.bval"),
                  subfiles, value=TRUE))>0 & 
      length(grep(paste0(sub_s,"_space-ACPC_desc-preproc_dwi.bvec"),
                  subfiles, value=TRUE))>0 &
      length(grep(paste0("(?=.*/dwi/)(?=.*",sub_s,"_space-ACPC_desc-brain_mask\\.nii(\\.gz)?)"), subfiles, perl = TRUE,value = TRUE))>0
  )
  {
    #define DWI volume and associated bvals and bvec
    bvec <- as.matrix(read.table(grep(paste0(sub_s,"_space-ACPC_desc-preproc_dwi.bvec"), subfiles, value = TRUE)))
    bval <- scan(grep(paste0(sub_s,"_space-ACPC_desc-preproc_dwi.bval"), subfiles,value = TRUE), quiet=TRUE)
    dwivol <- grep(paste0(sub_s,"_space-ACPC_desc-preproc_dwi.nii*"), subfiles, value = TRUE)
    #create dti package base object
    dtiDataobj <- dti::readDWIdata(
      gradient = bvec,
      bvalue   = bval,
      dirlist  = dwivol,
      format   = "NIFTI")
    #mask out DWI data using the brain mask in output
    dtiDataobj <- dti::setmask(dtiDataobj, grep(paste0("(?=.*/dwi/)(?=.*",sub_s,"_space-ACPC_desc-brain_mask\\.nii(\\.gz)?)"), subfiles, perl = TRUE, value = TRUE))
    return(list(dtiDataobj, RNifti::readNifti(dwivol)))
  } else {
    dtiDataobj=NA
    if (!silent){
      warning(paste0('No DTI/DKI map could be computed for',sub_s,', as either bval, bvec or brain_mask files are missing'))
      return(list(NA,NA))
    }
  }
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
#' @param skeleton_template The template FA skeleton (currently, only FMRIB58 2mm). It is optional as qsi_extract already preloads it when running skeleton_masker, but the latter function can load it by on its own if needed.
#' @param silent Whether to print messages and warnings or not. Default is FALSE.
#' @importFrom RNifti readNifti
#' @examples
#' skeleton_mask=skeleton_masker(skeleton_fathreshold=0.2) 
#' @export 

skeleton_masker=function(skeleton_template, skeleton_fathreshold=0.2){
  
  #load template skeleton if not provided
  if(missing(skeleton_template)){
    skeleton_template=RNifti::readNifti(paste0(system.file('extdata',package='WMskelstats'),'/templates/FMRIB58_FA-skeleton_2mm.nii'))
  }
  
  #make binary mask out of skeleton and threshold based on FA values: 
  skeleton_bin=array(0L, dim=dim(skeleton_template))
  thresh=skeleton_fathreshold*10000
  skeleton_bin[skeleton_template>=thresh]=1L
  
  #keep voxel coordinates of the mask for later rebuild
  skeleton_bin_coords <- which(skeleton_bin == 1, arr.ind = TRUE)
  
  return(list(skeleton_bin,skeleton_bin_coords))
}