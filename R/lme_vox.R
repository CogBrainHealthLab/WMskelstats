############################################################################################################################
############################################################################################################################
#' @title Voxel-wise analysis with threshold-free cluster enhancement (mixed effect)
#'
#' @description Fits a linear mixed effects model with the voxel-wise data as the predicted outcome, and returns t-stat and threshold-free cluster enhancement (TFCE) statistical maps for the selected contrast.
#' 
#' @details This TFCE method is adapted from the \href{https://github.com/nilearn/nilearn/blob/main/nilearn/mass_univariate/_utils.py#L7C8-L7C8}{'Nilearn' Python library}. 
#' 
#' @param model An N X P data.frame object containing N rows for each subject and P columns for each predictor included in the model.This data.frame should not include the random effects variable.
#' @param contrast A N x 1 numeric vector or object containing the values of the predictor of interest. Its length should equal the number of subjects in model (and can be a single column from model). The t-stat and TFCE maps will be estimated only for this predictor.
#' @param random A N x 1 numeric vector or object containing the values of the random variable (optional). Its length should be equal to the number of subjects in model (it should NOT be inside the model data.frame).
#' @param formula An optional string or formula object describing the predictors to be fitted against the surface data, replacing the model, contrast, or random arguments. If this argument is used, the formula_dataset argument must also be provided.
#' - The dependent variable (DV) is not needed, and the formula will start with ~. The DV will be the surface data value by default, but it can be swapped with contrast as IV via the "inverse" argument.
#' - The first independent variable in the formula will always be interpreted as the contrast of interest for which to estimate cluster-thresholded t-stat maps. 
#' - Only one random regressor can be given and must be indicated as '(1|variable_name)'.
#' @param formula_dataset An optional data.frame object containing the independent variables to be used with the formula (the IV names in the formula must match their column names in the dataset).
#' @param smooth_FWHM A numeric value specifying the desired smoothing width in mm. It should not be specified if the vox_data has been smoothed previously with smooth_vox(), because this results in vox_data being smoothed twice.
#' @param vox_data A N x V matrix object containing the voxel-wse data (N row for each subject, V for each voxel)
#' @param coords A V x 3 matrix, containing the X, Y and Z coordinates of the voxels used in `vox_data`.
#' @param nperm An integer specifying the number of permutations generated for the subsequent thresholding procedures (default = 100)
#' @param tail An integer specifying whether to test a one-sided positive (1), one-sided negative (-1) or two-sided (2) hypothesis
#' @param nthread An integer specifying the number of CPU threads to allocate 
#' @param perm_type A string object specifying whether to permute the rows ("row"), between subjects ("between"), within subjects ("within") or between and within subjects ("within_between") for random subject effects. Default is "row". 
#'
#' @returns A list object containing the t-stat and the TFCE statistical maps which can then be subsequently thresholded using TFCE_threshold()
#' 
#' @useDynLib WMskelstats, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom foreach foreach 
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom doSNOW registerDoSNOW
#' @export


##Main function

LME_TFCE=function(model,contrast, random, formula, formula_dataset, vox_data, coords, smooth_FWHM, nperm=1000, tail=2, nthread=4,  perm_type="row")
{
  #if the user chooses to use a formula, run the formula reader
  #and output appropriate objects
  if (!missing(formula) & !missing(formula_dataset))
  {
    formula_model=model_formula_reader(formula, formula_dataset) 
    model=formula_model$model
    contrast=formula_model$contrast
    if (!is.null(formula_model$random)) 
    {random=formula_model$random}
  } else if ((missing(formula) & !missing(formula_dataset)) | (!missing(formula) & missing(formula_dataset)))
  {stop('The formula and the formula_dataset arguments must both be provided to work.')}
  
  #run all checks for correct structure, recode variables when needed with
  #model_check()
  if (missing(random)) {stop('The random argument must be provided')};
  if (missing(smooth_FWHM)) {smooth_FWHM=NULL};
  model_summary=model_check(model=model, contrast=contrast, smooth_FWHM=smooth_FWHM,coords=coords,
                            random=random, vox_data=vox_data)
  model=model_summary$model
  contrast=model_summary$contrast
  vox_data=model_summary$vox_data
  colno=model_summary$colno

  if (!is.null(model_summary$random)) {random=model_summary$random}
  
  ##unpermuted TFCE
  start=Sys.time()
  message("Estimating unpermuted TFCE image...")
  mod=lme_fast(Y = vox_data,X = model,id=random)
  TFCE.orig=calc_tfce_3d(WMskelstats:::df_to_vol(coords,mod$t_stat[colno+1,]), tail=tail)
  end=Sys.time()
  
  message(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))
  
  ##permuted model
  #generating permutation sequences  
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  
  if(perm_type=="within_between") {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(random)}} 
  else if(perm_type=="within") {for (perm in 1:nperm)  {permseq[,perm]=perm_within(random)}} 
  else if(perm_type=="between") {for (perm in 1:nperm)  {permseq[,perm]=perm_between(random)}} 
  else if(perm_type=="row") {for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}}
  
  #activate parallel processing
  unregister_dopar = function() {
    .foreachGlobals <- utils::getFromNamespace(".foreachGlobals", "foreach"); 
    env =  .foreachGlobals;
    #rm(list=ls(name=env), pos=env) #handled by foreach::registerDoSEQ()
  }
  unregister_dopar()
  
  cl=parallel::makeCluster(nthread)
  doParallel::registerDoParallel(cl)
  #preload variables for the cluster workers
  parallel::clusterExport(cl, c("vox_data","permseq"), envir=environment())
  `%dopar%` = foreach::`%dopar%`
  
  #progress bar
  doSNOW::registerDoSNOW(cl)
  pb=txtProgressBar(max = nperm, style = 3)
  progress=function(n) setTxtProgressBar(pb, n)
  opts=list(progress = progress)
  
  #fitting permuted model and extracting max-TFCE values in parallel streams
  start=Sys.time()
  TFCE.max=foreach::foreach(perm=1:nperm, .combine="c",.packages = "WMskelstats", .options.snow = opts)  %dopar%
    {
      mod.perm=lme_fast(Y = vox_data[permseq[, perm]],X = model,id=random)
      return(max(abs(calc_tfce_3d(WMskelstats:::df_to_vol(coords,mod.perm$t_stat[colno+1,]), tail=tail)),na.rm = T))
    }
  end=Sys.time()
  message(paste("\nCompleted in ",round(difftime(end, start, units='mins'),1)," minutes \n",sep=""))
  parallel::stopCluster(cl)
  unregister_dopar()

  
  ##saving list objects
  returnobj=list(mod$t_stat[colno+1,],vol_to_df(TFCE.orig), TFCE.max,tail)
  names(returnobj)=c("t_stat","TFCE.orig","TFCE.max","tail")
  
  return(returnobj)
}