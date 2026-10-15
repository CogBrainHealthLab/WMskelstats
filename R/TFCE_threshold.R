############################################################################################################################
#' @title Thresholding TFCE output
#'
#' @description Threshold TFCE maps from the TFCE_vertex_analysis() output and identifies significant clusters at the desired threshold. 
#' 
#' @param TFCEoutput An object containing the output from TFCE_vertex_analysis()
#' @param p A numeric object specifying the p-value to threshold the results (Default is 0.05)
#' @param k Cluster-forming threshold (Default is 20)
#'
#' @returns A list object containing the cluster level results, unthresholded t-stat map, thresholded t-stat map, and positive, negative and bidirectional cluster maps.
#' 
#' @examples
#' \donttest{
#'  demodata=dl_demo(quiet=TRUE)
#'  #voxel-wise DTI data
#'  dti_data=readRDS(paste0(demodata,"/skel_fa_FMRIB58_skeleton_2mm_t0.2.rds"))
#'  skeleton=dti_data$skel_fa
#'  
#'  #behavioural data to match the dti data
#'  dat_beh=read.csv(paste0(demodata,"/ds007090_behdata.csv"))
#'  dat_beh$id_ses=paste0(dat_beh$participant_id, '_', dat_beh$session_id)
#'  dat_beh=dat_beh[which(dat_beh$id_ses %in% row.names(skeleton)),]
#'  dat_beh=dat_beh[order(dat_beh$participant_id),]
#'  
#'  model=LM_TFCE(model=dat_beh[,c('age_months','sex')],
#'                contrast=dat_beh[,'age_months'], 
#'                vox_data=skeleton[row.names(skeleton),], 
#'                coords=dti_data$skel_coords, 
#'                smooth_FWHM=5, 
#'                nperm=100, 
#'                tail=2, 
#'                nthread=2)
#'  model_thresh=TFCE_threshold(model, p=0.05,  k=20)
#'  model_thresh$cluster_level_results
#' }
#' 
#' @export

TFCE_threshold=function(TFCEoutput, p=0.05,  k=20)
{
  nperm=length(TFCEoutput$TFCE.max)
  n_vox=length(TFCEoutput$t_stat)
  
  #check if number of permutations is adequate
  if(nperm<1/p)  {warning(paste("Not enough permutations were carried out to estimate the p<",p," threshold precisely\nConsider setting an nperm to at least ",ceiling(1/p),sep=""))}
  
  ##generating p map
  tfce.p=rep(NA,n_vox)
  
  TFCEoutput$t_stat[is.na(TFCEoutput$t_stat)]=0
  for (vox in 1:n_vox)  {tfce.p[vox]= (1 + sum(TFCEoutput$TFCE.max >=abs(TFCEoutput$TFCE.orig$value[vox]))) / (nperm + 1)}
  TFCEoutput$t_stat[is.na(TFCEoutput$t_stat)]=0
  
  t_stat.thresholdedP=TFCEoutput$t_stat
  t_stat.thresholdedP[tfce.p>p]=0
  ##Cluster level results
  ##positive cluster
  if(TFCEoutput$tail==1 |TFCEoutput$tail==2)
  {
    #zeroing out all negative voxels
    pos.t_stat.thresholdedP=t_stat.thresholdedP
    pos.t_stat.thresholdedP[pos.t_stat.thresholdedP<0]=0
    pos.t_stat.thresholdedP.vol=df_to_vol(coords = TFCEoutput$coords, data=pos.t_stat.thresholdedP)   
    pos.clust.results=get_clusters(pos.t_stat.thresholdedP.vol,min_size = k)
    if(NROW(pos.clust.results)==0)
    {
      pos.clust.results="No significant clusters"
    }
  }
  if(TFCEoutput$tail==-1 |TFCEoutput$tail==2)
  {
    #zeroing out all positive voxels
    neg.t_stat.thresholdedP=t_stat.thresholdedP
    neg.t_stat.thresholdedP[neg.t_stat.thresholdedP>0]=0
    neg.t_stat.thresholdedP.vol=df_to_vol(coords = TFCEoutput$coords, data=neg.t_stat.thresholdedP)   
    neg.clust.results=get_clusters(neg.t_stat.thresholdedP.vol,min_size = k)
    if(NROW(neg.clust.results)==0)
    {
    neg.clust.results="No significant clusters"
    }
  }
  if(TFCEoutput$tail==2)
  {
    clust.results=list(pos.clust.results,neg.clust.results)
    t_stat.thresholded.return=t_stat.thresholdedP
    pos.mask=pos.t_stat.thresholdedP
    pos.mask[pos.mask>0]=1
    neg.mask=neg.t_stat.thresholdedP
    neg.mask[neg.mask<0]=1
    
  } else if(TFCEoutput$tail==1)
  {
    clust.results=list(pos.clust.results,"Negative contrast not analyzed, only positive one-tailed TFCE statistics were estimated")
    t_stat.thresholded.return=pos.t_stat.thresholdedP
    pos.mask=pos.t_stat.thresholdedP
    pos.mask[pos.mask>0]=1
    neg.mask=NULL
  } else if(TFCEoutput$tail==-1)
  {
    clust.results=list("Positive contrast not analyzed, only negative one-tailed TFCE statistics were estimated",neg.clust.results)
    t_stat.thresholded.return=neg.t_stat.thresholdedP
    pos.mask=NULL
    neg.mask=neg.t_stat.thresholdedP
    neg.mask[neg.mask<0]=1
  }
  
  names(clust.results)=c("Positive contrast", "Negative contrast")
  t_stat.thresholded.return[t_stat.thresholded.return==0]=NA
  returnobj=list(clust.results,TFCEoutput$t_stat,t_stat.thresholded.return,pos.mask,neg.mask)
  names(returnobj)=c("cluster_level_results",
                     "unthresholded_tstat_map",
                     "thresholded_tstat_map",
                     "pos_mask","neg_mask")
  return(returnobj)
}

