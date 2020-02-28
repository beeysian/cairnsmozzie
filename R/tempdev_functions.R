#' Given a latitude and longitude, determines the corresponding grid ID.
#' Assumes grid.df is a global variable.
#' This is to find the grid ID of an agent. assign_grid() is for microclimates.
#' @param testlat Latitude being considered.
#' @param testlong Longitude being considered.
#' @return Corresponding grid ID number.
get_gridID <- function(testlat, testlong){
  #closelong <- which(abs(grid.df$long - testlong) == min(abs(grid.df$long - testlong)))
  #Changed to the below on 21/1/20:
  closelong <- which(abs(grid.df$V2 - testlong) == min(abs(grid.df$V2 - testlong)))
  if(length(closelong) > 1){
    # The below used to have grid.df$lat but looks like we need to use V1, V2
    closelat <- closelong[which(abs(grid.df$V1[closelong] - testlat) == min(abs(grid.df$V1[closelong] - testlat)))]
  }else{closelat <- closelong }

#  else if(length(closelong) == 0){
#    closelat <- closelong
#  }
  if(length(closelat) > 1){
    closelat <- closelat[1]
  }
return(closelat)
}
