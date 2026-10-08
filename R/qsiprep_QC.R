#' @title qsiprep_QC
#'
#' @description collating all qsiprep QC data into a single .csv file, identifies subjects who fail QC, and if there are QC failures, a bash script (`del.sh`) will be generated to delete these subject folders
#'
#' @details This function collates all qsiprep QC data into a single .csv file and identifies subjects whose `raw_neighbor_corr` and/or `mean_FD` violates a preset criterion (i.e., number of SDs from the mean).Since there is a lack of consensus on the absolute cutoffs for quality controlling qsiprep output, a relative criteria was used instead. 
#'
#' @param filename Filename of the .csv file. Set to `QSIPREP_QC.csv` by default
#' @param N_SD_FD The relative cutoffs for `mean_FD` in terms of the number of SDs from the mean. Subjects whose `mean_FD` > (mean(mean_FD)+N_SD_FD*sd(mean_FD)) will be identified as QC failure cases. Set to `2` by default
#' @param N_SD_neigbour The relative cutoffs for `raw_neighbor_corr` in terms of the number of SDs from the mean. Subjects whose `raw_neighbor_corr` < (mean(raw_neighbor_corr)-N_SD_neighbor*sd(raw_neighbor_corr)) will be identified as QC failure cases. Set to `2` by default
#' @param absthresh_FD The absolute cutoffs for `mean_FD` . Subjects whose `mean_FD` > `absthresh_FD` will be identified as QC failure cases. Not specified by default. `N_SD_FD` will be ignored if specified
#' @param absthresh_neighbor The absolute cutoffs for `raw_neighbor_corr` . Subjects whose `raw_neighbor_corr` < `absthresh_neighbor` will be identified as QC failure cases. Not specified by default. `N_SD_neighbor` will be ignored if specified

#' @returns outputs a .csv file containing the aslprep QC data and the bash script `del.sh` if there are scans that fail the QC 
#'
#' @examples
#' \dontrun{
#' qsiprep_QC()
#' }
#' @export

########################################################################################################
########################################################################################################


qsiprep_QC=function(filename="QSIPREP_QC.csv", N_SD_FD=2, N_SD_neighbor=2, absthresh_FD, absthresh_neighbor)
{
  #modified rbind function 
  rbind_fill <- function(df1, df2) {
    all_cols <- union(names(df1), names(df2))
    for (col in setdiff(all_cols, names(df1))) df1[[col]] <- NA
    for (col in setdiff(all_cols, names(df2))) df2[[col]] <- NA
    rbind(df1, df2[, names(df1)])  # reorder df2 to match df1
  }
  
  #compile QC data from all subjects and sessions
  filelist=list.files(path = Sys.glob("sub-*"),pattern = "_desc-image_qc.tsv", recursive = T,full.names = T)
  subses=basename(gsub(pattern = "_desc-image_qc.tsv","",filelist))  
  for(sub in 1:length(filelist))
  {
    if(sub==1)
    {
      dat=read.table(filelist[sub],header = T, sep="\t") 
    } else
    {
      dat.temp=read.table(filelist[sub],header = T, sep="\t")
      dat=rbind_fill(dat,dat.temp)
      remove(dat.temp)
    }
  }
  #identify scans that fail QC
  if(missing(absthresh_FD)) {FD_crit=mean(dat$mean_fd)+(N_SD_FD*sd(dat$mean_fd))} 
  else  {FD_crit=absthresh_FD}
  
  if(missing(absthresh_neighbor)) {neighbor_crit=mean(dat$raw_neighbor_corr)-(N_SD_neighbor*sd(dat$raw_neighbor_corr))}
  else  {neighbor_crit=absthresh_neighbor}
  
  idx=unique(c(which(dat$raw_neighbor_corr < neighbor_crit),
               which(dat$mean_FD > FD_crit)))
  if(length(idx)>0)
  {
    ##if runs are detected
    if(!all(is.na(dat$run_id)))
    {
      if(length(unique(dat$run))==1)
      {
        del.sh=paste0("rm -rf ",dat$subject_id[idx],"/",dat$ses[idx])    
      } else
      {
        del.sh=paste0("rm -rf ",dat$subject_id,"/",dat$ses[idx],"/",dat$subject_id[idx],"_",dat$session_id[idx],"_run-",dat$run[idx],"*")
        del.sh=gsub("run-NA","",del.sh)
      }
    } else
    {
      del.sh=paste0("rm -rf ",dat$subject_id[idx],"/",dat$session_id[idx])    
    }

    cat(paste0(length(del.sh), " DWI scans failed QC. QC criteria:mean_FD > ",round(FD_crit,2), " and raw_neighbor_corr < ",round(neighbor_crit,2)))
    write.table(del.sh,row.names=F, col.names=F,quote=F, file="del.sh")
  } else { cat("all DWI scans passed QC")}
  write.table(dat, file=filename, row.names = F, sep=",")
}
