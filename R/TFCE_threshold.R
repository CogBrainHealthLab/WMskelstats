############################################################################################################################
#' @title Thresholding TFCE output
#'
#' @description Threshold TFCE maps from the TFCE_vertex_analysis() output and identifies significant clusters at the desired threshold. 
#' 
#' @param TFCEoutput An object containing the output from TFCE_vertex_analysis()
#' @param p A numeric object specifying the p-value to threshold the results (Default is 0.05)
#' @param atlas An integer corresponding to the atlas of interest.  1=Desikan, 2=Destrieux-148, 3=Glasser-360, 4=Schaefer-100, 5=Schaefer-200, 6=Schaefer-400. Set to `1` by default. This argument is ignored for hippocampal surfaces. For applicable SubCortexMesh surfaces, 1=default base ROI, 2=anatomical atlas.
#' @param k Cluster-forming threshold (Default is 20)
#'
#' @returns A list object containing the cluster level results, unthresholded t-stat map, thresholded t-stat map, and positive, negative and bidirectional cluster maps.

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
  for (vox in 1:n_vox)  {tfce.p[vox]=length(which(TFCEoutput$TFCE.max>abs(TFCEoutput$TFCE.orig$value[vox])))/nperm}
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
    pos.t_stat.thresholdedP.vol=WMskelstats:::df_to_vol(coords = data.matrix(check[,c("x","y","z")]), data=pos.t_stat.thresholdedP)   
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
    neg.t_stat.thresholdedP.vol=WMskelstats:::df_to_vol(coords = data.matrix(check[,c("x","y","z")]), data=neg.t_stat.thresholdedP)   
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
    clust.results=list(pos.clust.results,"Negative contrast not analyzed, only negative one-tailed TFCE statistics were estimated")
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
  
  t_stat.thresholded.return[t_stat.thresholded.return==0]=NA
  returnobj=list(clust.results,TFCEoutput$t_stat,t_stat.thresholded.return,pos.mask,neg.mask)
  names(returnobj)=c("cluster_level_results",
                     "unthresholded_tstat_map",
                     "thresholded_tstat_map",
                     "pos_mask","neg_mask")
  return(returnobj)
}

