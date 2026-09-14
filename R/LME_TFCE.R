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
#' @importFrom mori share
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
  ## Freedman–Lane preparation: fit the reduced nuisance model once
  model = as.matrix(model)
  
  # Explicit intercept, including when the reduced model is intercept-only.
  # As in your original function, model should not already contain an intercept.
  X_full = cbind("(Intercept)" = 1, model)
  X_null = X_full[, -(colno + 1L), drop = FALSE]
  
  mod.null = lme_fast(
    Y = vox_data,
    X = X_null,
    id = random,
    add_intercept = FALSE
  )
  
  # Fixed-effects predictions ONLY: do not add fitted random intercepts.
  fitted_null = X_null %*% mod.null$coefficients
  
  # Marginal residuals retain participant-level variation.
  residuals_null = vox_data - fitted_null
  
  
  mod = lme_fast(Y = vox_data, X = model, id = random)
  
  t_orig = mod$t_stat[colno + 1, ]
  
  # Calculate one height step for the entire analysis.
  t_finite = t_orig[is.finite(t_orig)]
  
  if (!length(t_finite)) {
    stop("The observed t-statistic map contains no finite values.")
  }
  
  peak = switch(
    as.character(tail),
    "1"  = max(0, t_finite),
    "-1" = max(0, -t_finite),
    "2"  = max(abs(t_finite)),
    stop("tail must be 1, -1, or 2.")
  )
  
  # A positive fallback is needed if the relevant observed map is all zero.
  dh_fixed = if (peak > 0) peak / 100 else 0.1
  
  TFCE.orig = calc_tfce_3d(
    WMskelstats:::df_to_vol(coords, t_orig),
    dh = dh_fixed,
    tail = tail
  )
  end=Sys.time()
  
  message(paste("Completed in",round(difftime(end,start, units="secs"),1),"secs\nEstimating permuted TFCE images...\n",sep=" "))
  
  ##permuted model
  #generating permutation sequences  
  permseq=matrix(NA, nrow=NROW(model), ncol=nperm)
  
  if(perm_type=="within_between") {for (perm in 1:nperm)  {permseq[,perm]=perm_within_between(random)}} 
  else if(perm_type=="within") {for (perm in 1:nperm)  {permseq[,perm]=perm_within(random)}} 
  else if(perm_type=="between") {for (perm in 1:nperm)  {permseq[,perm]=perm_between(random)}} 
  else if(perm_type=="row") {for (perm in 1:nperm)  {permseq[,perm]=sample.int(NROW(model))}}
  
  # 1. Convert large datasets to shared memory objects
  fitted_null_shm = mori::share(fitted_null)
  residuals_null_shm = mori::share(residuals_null)
  permseq_shm = mori::share(permseq)
  
  rm(fitted_null, residuals_null, mod.null)
  
  # 2. Initialize cluster and attach single backend (doSNOW)
  cl=parallel::makeCluster(nthread)
  `%dopar%` = foreach::`%dopar%`
  doSNOW::registerDoSNOW(cl)
  
  # 3. Export shared object pointers and lightweight metadata once
  parallel::clusterEvalQ(cl, {
    loadNamespace("mori")
    loadNamespace("WMskelstats")
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = c(
      "fitted_null_shm",
      "residuals_null_shm",
      "permseq_shm",
      "model", "random", "coords", "colno", "tail"
    ),
    envir = environment()
  )
  
  # Progress bar setup
  pb=txtProgressBar(max = nperm, style = 3)
  opts=list(progress = function(n) setTxtProgressBar(pb, n))
  
  # 4. Execute parallel loop using shared memory handles
  start=Sys.time()
  TFCE.max = foreach::foreach(
    perm = seq_len(nperm),
    .combine = "c",
    .packages = c("WMskelstats", "mori"),
    .noexport = c("vox_data", "permseq","fitted_null_shm", "residuals_null_shm", "permseq_shm","model", "random", "coords", "colno", "tail"),
    .options.snow = opts,
    .errorhandling = "stop"
  ) %dopar% {
    
    perm_idx = permseq_shm[, perm]
    
    # Freedman–Lane: permute reduced-model marginal residuals,
    # then restore the original fixed nuisance effects.
    vox_data_perm = fitted_null_shm +residuals_null_shm[perm_idx, , drop = FALSE]
    
    mod.perm = lme_fast(
      Y = vox_data_perm,
      X = model,
      id = random
    )
    
    # For the updated fitter, stop on unsuccessful optimization.
    if (!is.null(mod.perm$status) &&
        any(is.na(mod.perm$status) | mod.perm$status != 0L)) {
      stop("Model fitting failed in permutation ", perm)
    }
    
    t_perm = mod.perm$t_stat[colno + 1L, ]
    
    if (any(!is.finite(t_perm))) {
      stop("Nonfinite t-statistics in permutation ", perm)
    }
    
    vol_stat = WMskelstats:::df_to_vol(coords, t_perm)
    
    max(
      abs(calc_tfce_3d(vol_stat, dh = dh_fixed, tail = tail)),
      na.rm = TRUE
    )
  }
  
  close(pb)
  end=Sys.time()
  
  message(sprintf("\nCompleted in %.1f minutes \n", difftime(end, start, units = "mins")))
  
  parallel::stopCluster(cl)

  
  ##saving list objects
  returnobj=list(mod$t_stat[colno+1,],vol_to_df(TFCE.orig), TFCE.max,tail)
  names(returnobj)=c("t_stat","TFCE.orig","TFCE.max","tail")
  
  return(returnobj)
}