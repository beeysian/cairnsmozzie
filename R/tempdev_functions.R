#' Given a latitude and longitude, determines the corresponding grid ID.
#' Assumes grid.df is a global variable.
#' @param testlat Latitude being considered.
#' @param testlong Longitude being considered.
#' @return Corresponding grid ID number.
get_gridID <- function(testlat, testlong){
  closelong <- which(abs(grid.df$long - testlong) == min(abs(grid.df$long - testlong)))
  if(length(closelong) > 1){
    closelat <- closelong[which(abs(grid.df$lat[closelong] - testlat) == min(abs(grid.df$lat[closelong] - testlat)))]
  }
  else if(length(closelong) == 0){
    closelat <- closelong
  }
return(closelat)
}
