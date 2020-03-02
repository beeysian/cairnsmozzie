#' Generates random 'starting' positions for adult agents.
#' Uses boundary data for simulation region.
#' Draws a polygon from the boundary points given in boundaryDat
#'     and picks N randomly spaced points within it.
#' Includes a boundary check:
#'   will keep drawing points until there are N in the region/polygon.
#' Requires package 'sp' for Polygon object.
#' @param boundaryDat Geographic boundary of simulation region.
#' @param N Number of agents.
#' @return A data frame of lat and long positions for each adult agent.
init_position <- function(boundaryDat, N){
  i <- -1
  while(i == -1){
    polyg     <- Polygon(boundaryDat)
    samplePts <- spsample(polyg, N, type="random") #see documentation for other options for 'type'
    #check if points are in bounds:
    if(sum(point.in.polygon(samplePts$x,samplePts$y,boundaryDat$Long,boundaryDat$Lat)) >= N){
      i      <- 0
      posdf  <- as.data.frame(cbind(long=samplePts$x, lat=samplePts$y))
      griddf <- mapply(FUN = get_gridID, testlat = posdf$lat, testlong = posdf$long)
      posdf  <- cbind(posdf, gridID = griddf)
      return(posdf)
    }
  }
}

#' Roughly determines the stage of a juvenile clutch given numerical 'age'.
#' Needs checking with literature, or replace with better method.
#'
#' @param age 'Age' of clutch in days.
#' @return Juvenile stage: 0 for egg, 1 for larval, 3 for pupal.
init_juv_stage <- function(age){
  if(age >= 0 & age <= 5){
    stage <- 1
  }
  else if(age > 5 & age <= 11){
    stage <- 2
  }
  else{
    stage <- 3
  }
  return(stage)
}

#' Grabs the land type of a mosquito given its grid ID.
#' Used for microclimates work.
#'
#' @param gridID Grid ID of an agent
#' @return Land type
grid_to_land_type <- function(gridID){
  landtype <- landtype.df$LandType[gridID]
  return(landtype)
}
