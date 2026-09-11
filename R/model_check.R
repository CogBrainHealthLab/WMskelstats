
#' @title Model structure check
#' @description Ensures the voxel data, contrast, model and/or random objects in the analyses have the appropriate structures to enable vertex analyses to be run properly on the given variables.
#' @returns An error message if an issue or discrepancy is found.
#' @noRd

model_check=function(contrast, model, random, vox_data,smooth_FWHM,coords)
{
  #check if coords structure is correct
  if(!sum(dim(coords)==c(NROW(vox_data),3))!=2) {stop("coords must be provided in the form of a V x 3 data.matrix object")}
  
  #If the contrast/model is a tibble (e.g., taken from a read_csv output)
  #converts the columns to regular data.frame column types
  if ('tbl_df' %in% class(contrast) == TRUE) {
    if (inherits(contrast[[1]],"character")==TRUE) {contrast = contrast[[1]]
    } else {contrast = as.numeric(contrast[[1]])}
  } 
  if ('tbl_df' %in% class(model) == TRUE) {
    model=as.data.frame(model)
    if (NCOL(model)==1) {model = model[[1]]
    } else { for (c in 1:NCOL(model)) { 
      if(inherits(model[,c],"double")==TRUE) {model[,c] = as.numeric(model[,c])}
    }  }
  }
  
  #if contrast or model is a data.frame with 1 column, the variable needs to be flattened for the inherits() checks to work properly
  if(inherits(contrast,"data.frame")==TRUE & NCOL(contrast)==1){contrast=contrast[[1]]}
  if(inherits(model,"data.frame")==TRUE & NCOL(model)==1){model=model[[1]]}
  
  #numerise contrast
  if(inherits(contrast,"integer")==TRUE) {contrast=as.numeric(contrast)}
  
  #check if nrow is consistent for model and vox_data
  if(NROW(vox_data)!=NROW(model))  {stop(paste("The number of rows for vox_data (",NROW(vox_data),") and model (",NROW(model),") are not the same",sep=""))}
  
  #recode random variable to numeric
  if(!is.null(random)) { random=match(random,unique(random)) }
  
  ##checks
  #check contrast for consistency with the model data.frame
  if(NCOL(model)>1)
  {
    for(colno in 1:(NCOL(model)+1))
    {
      if(colno==(NCOL(model)+1))  {warning("contrast is not contained within model")}
      
      if(inherits(contrast,"character")==TRUE) 
      {
        if(identical(contrast,model[,colno]))  {break} 
      } else 
      {
        if(identical(suppressWarnings(as.numeric(contrast)),suppressWarnings(as.numeric(model[,colno]))))  {break}
      }
    }
  }  else
  {
    if(inherits(contrast,"character")==TRUE) 
    {
      if(identical(contrast,model))  {colno=1} 
      else  {stop("contrast is not contained within model")}
    } else
    {
      if(identical(as.numeric(contrast),as.numeric(model)))  {colno=1}
      else  {stop("contrast is not contained within model")}
    }
  }
  
  #incomplete data check
  idxF=which(complete.cases(model)==FALSE)
  if(length(idxF)>0)
  {
    message(paste("The model contains",length(idxF),"subjects with incomplete data. Subjects with incomplete data will be excluded from the current analysis\n"))
    model=model[-idxF,]
    contrast=contrast[-idxF]
    vox_data=vox_data[-idxF,]
    if(!is.null(random)) {random=random[-idxF]}
  }
  
  #check categorical and recode variable
  if(NCOL(model)>1)
  {
    for (column in 1:NCOL(model))
    {
      if(inherits(model[,column],"character")==TRUE | inherits(model[,column],"factor")==TRUE)
      {
        if(length(unique(model[,column]))==2)
        {
          message(paste("The binary variable '",colnames(model)[column],"' will be recoded with ",unique(model[,column])[1],"=0 and ",unique(model[,column])[2],"=1 for the analysis\n",sep=""))
          
          recode=rep(0,NROW(model))
          recode[model[,column]==unique(model[,column])[2]]=1
          model[,column]=recode
          contrast=model[,colno]
        } else if(length(unique(model[,column]))>2)    {stop(paste("The categorical variable '",colnames(model)[column],"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
      }      
    }
  } else
  {
    
    if(inherits(model,"character")==TRUE | inherits(model,"factor")==TRUE)
    {
      if(length(unique(model))==2)
      {
        message(paste("The model variable is binary and will be recoded such that ",unique(model)[1],"=0 and ",unique(model)[2],"=1 for the analysis\n",sep=""))
        
        recode=rep(0,NROW(model))
        recode[model==unique(model)[2]]=1
        model=recode
        contrast=model
      } else if(length(unique(model))>2)    {stop(paste("The categorical variable '",colnames(model),"' contains more than 2 levels, please code it into binarized dummy variables",sep=""))}
    }      
  }
  
  
  #check if vox_data is a multiple-rows matrix and NOT a vector
  if (is.null(nrow(vox_data)) | nrow(vox_data)==1)
  {stop("The voxel data must be a matrix containing multiple participants (rows).")}
  
  
  ##smoothing
  if(is.null(smooth_FWHM))
  {
    message("smooth_FWHM argument was not given. vox_data will not be smoothed here.\n")
  } else if(smooth_FWHM==0) 
  {
    message("smooth_FWHM set to 0: vox_data will not be smoothed here.\n")
  } else if(smooth_FWHM>0) 
  {
    message(paste("vox_data will be smoothed using a ",smooth_FWHM,"mm FWHM kernel", sep=""))
    vox_data=smooth_vox(data_mat,coords, fwhm = smooth_FWHM, sigma = NULL, voxel_size = c(1, 1, 1))
  }
  
  
  ##########################################
  #Output the right elements to be analysed
  model_summary=list(model=model, contrast=contrast,
                     random=random, vox_data=vox_data,
                     colno=colno)
  
  return(model_summary) 
}