get_clusters <- function(img_3d, min_size = 2, connectivity = 26) {
  # Validate input
  if (!is.array(img_3d) || length(dim(img_3d)) != 3) {
    stop("Input must be a 3D array.")
  }
  
  # Define neighborhood kernel (6-face or 26-corner/edge/face connected)
  if (connectivity == 6) {
    kernel <- mmand::shapeKernel(c(3, 3, 3), type = "diamond")
  } else if (connectivity == 26) {
    kernel <- mmand::shapeKernel(c(3, 3, 3), type = "box")
  } else {
    stop("Connectivity must be either 6 or 26.")
  }
  
  # Step 1: Binarize array and label contiguous non-zero clusters
  binary_map <- img_3d != 0
  labeled_map <- mmand::components(binary_map, kernel = kernel)
  
  # Step 2: Extract non-zero coordinates, labels, and values
  nz_coords <- which(labeled_map > 0, arr.ind = TRUE)
  
  # Empty result dataframe structure
  empty_df <- data.frame(
    cluster_id   = integer(0),
    cluster_size = integer(0),
    max_value    = numeric(0),
    x            = integer(0),
    y            = integer(0),
    z            = integer(0)
  )
  
  if (nrow(nz_coords) == 0) return(empty_df)
  
  df <- data.frame(
    x       = nz_coords[, 1],
    y       = nz_coords[, 2],
    z       = nz_coords[, 3],
    val     = img_3d[nz_coords],
    cluster = labeled_map[nz_coords]
  )
  
  # Step 3: Filter clusters meeting the minimum voxel threshold
  cluster_counts <- table(df$cluster)
  valid_clusters <- as.numeric(names(cluster_counts[cluster_counts >= min_size]))
  
  if (length(valid_clusters) == 0) return(empty_df)
  
  df_valid <- df[df$cluster %in% valid_clusters, ]
  
  # Step 4: Summarize cluster sizes, maxima, and peak coordinates
  cluster_list <- split(df_valid, df_valid$cluster)
  
  results <- lapply(cluster_list, function(sub_df) {
    max_idx <- which.max(abs(sub_df$val))
    peak <- sub_df[max_idx, ]
    
    data.frame(
      cluster_id   = peak$cluster,
      cluster_size = nrow(sub_df),
      max_value    = peak$val,
      x            = peak$x,
      y            = peak$y,
      z            = peak$z
    )
  })
  
  # Combine, sort by cluster size, and return
  out_df <- do.call(rbind, results)
  out_df <- out_df[order(-out_df$cluster_size), ]
  rownames(out_df) <- NULL
  
  return(out_df)
}



##converting between dataframe and 3D array formats

vol_to_df <- function(data) {
  #coords=which(data!=min(data), arr.ind = T)
  coords=which(!is.na(data), arr.ind = T)
  # Convert 1D indices to 3D matrix coordinates
  
  res <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    z = coords[, 3],
    value = data[coords]
  )
  res <- res[order(res$x, res$y, res$z), ]
  
  return(res)
}
df_to_vol=function(coords,data)
{
  
  img_array=array(NA, dim = c(max(coords[,1]),max(coords[,2]),max(coords[,3])))
  img_array[coords]=data
  return(img_array)
}


## permutation functions for random subject effects
## Paired/grouped data points are first shuffled within subjects, then these pairs/groups are shuffled between subjects
perm_within_between=function(random)
{
  ##for groups of 2 or more (subjects with 2 or more measurements)
  perm.idx=rep(NA, length(random))
  for(count in 2:max(table(random)))
  {
    if(length(which(table(random)==count))>0)
    {
      sub.id=as.numeric(which(table(random)==count))
      if(length(sub.id)>1)
      {
        ##between group shuffling
        recode.vec=sample(sub.id)
        vec.idx=1
        for(sub in sub.id)
        {
          perm.idx[which(random==sub)]=sample(which(random==recode.vec[vec.idx])) ##sample— within subject shuffling
          vec.idx=vec.idx+1
        }   
        remove(vec.idx,recode.vec)  
      } else 
      {
        ##if only one subject has a certain count, between subject shuffling will not be possible, only within-subject shuffling will be carried out
        perm.idx[which(random==sub.id)]=sample(which(random==sub.id)) ##sample— within subject shuffling
      }
    }
  }
  ##for subjects with a single measurement
  sub.idx=which(is.na(perm.idx))
  if(length(sub.idx)>1)
  {
    perm.idx[sub.idx]=sample(sub.idx)  
  } else 
  {
    perm.idx[sub.idx]=sub.idx
  }
  return(perm.idx)
}

## Paired/grouped data points are shuffled within subjects, order of subjects in the dataset remains unchanged
perm_within=function(random)
{
  ##for groups of 2 or more (subjects with 2 or more measurements)
  perm.idx=rep(NA, length(random))
  
  for(count in 2:max(table(random)))
  {
    if(length(which(table(random)==count)>0))
    {
      sub.id=as.numeric(which(table(random)==count))
      for(sub in sub.id)
      {
        perm.idx[which(random==sub)]=sample(which(random==sub))
      }  
    }
  }
  return(perm.idx)
}

## Paired/grouped data points are shuffled between subjects, order of data points within subjects remains unchanged.
perm_between=function(random)
{
  ##for groups of 2 or more (subjects with 2 or more measurements)
  perm.idx=rep(NA, length(random))
  for(count in 2:max(table(random)))
  {
    if(length(which(table(random)==count))>0)
    {
      sub.id=as.numeric(which(table(random)==count))
      if(length(sub.id)>1)
      {
        ##between group shuffling
        recode.vec=sample(sub.id)
        vec.idx=1
        for(sub in sub.id)
        {
          perm.idx[which(random==sub)]=which(random==recode.vec[vec.idx])
          vec.idx=vec.idx+1
        }   
        remove(vec.idx,recode.vec)  
      }
    }
  }
  ##for subjects with a single measurement
  sub.idx=which(is.na(perm.idx))
  if(length(sub.idx)>1)
  {
    perm.idx[sub.idx]=sample(sub.idx)  
  } else 
  {
    perm.idx[sub.idx]=sub.idx
  }
  return(perm.idx)
}