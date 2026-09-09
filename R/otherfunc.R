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
