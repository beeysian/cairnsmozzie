#' Converts a polygon boundary into a grid mesh, but in dataframe form
#'
#' @param bdary List of points along the boundary of the simulation area (in lat/long)
#' @param gridSize Grid size. Assumes square grid. In decimal form (eg 0.00025 = 25m^2 grid)
#' @return A dataframe of lat/long coordinates representing the region described by bdary
#'         with space sizes gridSize.
make_grid <- function(bdary, gridSize){
  latval  <- seq(min(bdary$lat), max(bdary$lat), gridSize)
  longval <- seq(min(bdary$long), max(bdary$long), gridSize)
  mesh    <- expand.grid(latval,longval) #turns list of lat and long into grid
  points  <- point.in.polygon(mesh$Var1, mesh$Var2, bdary$lat, bdary$long) #identifies grid points that lie within our bounday: NB does not include those on the boundary
  grid    <- cbind(mesh$Var1[which(points==1)],mesh$Var2[which(points==1)]) #puts these grid spaces into a df
  grid    <- as.data.frame(grid)
  return(grid)
}

#' Determines the grid ID number for each grid space within a bounding box.
#' Useful for, eg, large areas of greenery that we can easily draw a boundary
#'    around, instead of grabbing the lat/long of every point within the region.
#' Returns a dataframe of grid ID numbers that will be associated with that land type.
#'
#' @param bb Bounding box of region, in lat/long
#' @param gridSize Grid size. Assumes square grid. In decimal form (eg 0.00025 = 25m^2 grid)
#' @param grid Dataframe of area grid, as per what's created in make_grid
#' @return A dataframe of grid IDs associated with that land type.
bb_to_matchedIDs <- function(bb, gridSize, grid){
  # Passing "grid" might be slow, but let's just do it for now so I can get this code done
  mesh    <- expand.grid(seq(min(bb$lat), max(bb$lat), gridSize), seq(min(bb$long), max(bb$long), gridSize))
  points  <- point.in.polygon(mesh$Var1, mesh$Var2, bb$lat, bb$long) #identifies grid points that lie within our bounday: NB does not include those on the boundary

  bb.grid <- cbind(mesh$Var1[which(points > 0)], mesh$Var2[which(points > 0)]) #puts these grid spaces into a df
  bb.grid <- unique(as.data.frame(bb.grid)) # Get rid of duplicates

  matched.IDs <- lapply(1:nrow(bb.grid), function(x){assign_grid(bb.grid[x, ], grid)}) #the 'grid' is our simulation region 'grid'
  grid.pos    <- unlist(matched.IDs)
  matched.IDs <- as.data.frame(grid.pos)
  return(matched.IDs)
}

#' Determines the grid ID number corresponding to each lat/long point in the input.
#' Returns a dataframe of grid ID numbers that will be associated with that land type.
#'
#' @param points Lat/long points to be matched to a grid ID.
#' @param grid Dataframe of area grid, as per what's created in make_grid
#' @return A dataframe of grid IDs associated with that land type.
points_to_matchedIDs <- function(points, grid){
  IDs         <- lapply(1:nrow(points), function(x){assign_grid(points[x, ], grid)})
  grid.pos    <- unlist(IDs)
  matched.IDs <- as.data.frame(grid.pos)
  return(grid.pos)
}

#' Assigns a point to a grid ID.
#' Used for functions points_to_matchedIDs and bb_to_matchedIDs.
#' @param pos A point in lat/long coordinates
#' @param grid Dataframe of area grid, as per what's created in make_grid
#' @return A grid ID that has been assigned to that point. Should be 1 integer.
assign_grid <- function(pos, grid.df){
  close <- which((abs(pos$V2 - grid.df$V2) < gridSize) & abs(pos$V1 - grid.df$V1 < gridSize))
  grid.match <- -1
  if(length(close > 0)){
    close.dist <- list(abs(pos$V2 - grid.df$V2[close]), abs(pos$V1 - grid.df$V1[close]))
    closest    <- which(close.dist[[1]] == min(close.dist[[1]]) & close.dist[[2]] == min(close.dist[[2]]))
    grid.match <- close[closest]
  }
  else if(length(close == 0)){
    grid.match <- -1
  }
  return(grid.match)
}

#' Constructs the EKM chart with microclimates.
#' Uses data about the area of the region we're simulating.
#' If you're using this function for a different area, you will need to change
#'  all of the data that's read in in the first part of the function.
#' The lookup table is deterministic, and is made of nested lists.
#' The order will be: [stage][land type][day]
#' So when we update the EKS for each agent, we add the value that is at df[[stage]][land type, ][day]
#' Eg: for an egg that is in land type 1 at time step t=2, we add the value at dt[0][1][2] to the agents' EKS
#' To use this list, indexing works like: EKM_chart[[stage]][land_type, ][day]
#' Land types (as of 16/1/19) are:
#' \enumerate{
#' \item{Uncovered greenery}
#' \item{Vegetation}
#' \item{Roads}
#' \item{Water bodies}
#' \item{Carpark/uncovered concrete}
#' \item{Large buildings}
#' \item{Houses (tempdev informed by house type)}
#' }
#' @param bdary Boundary of simulation region.
#' @param tempData Vector of average temperatures over the simulation period.
#' @param typeTempDevs Temperature deviates for each land type.
#' @param noDays Number of days in the simulation period.
#' @param noLandTypes Number of different land types we consider.
#' @param gridSize Grid size, as per make_grid and bb_to_matchedIDs functions
#' @param const Constants for EKM calculation.
#' @return A list of 4 dataframes for each development stage, of the EKS for each day of the simulation per each land type
initialise_enzyme_wmicroclim <- function(bdary, tempData, typeTempDevs, noDays, noLandTypes, gridSize, const){
  # ---- Read in data
  houses <- read.table("inst/exdata/pp_housing_data.txt", header=TRUE) #Location of houses

  # The data taken from Google Maps ----
  greenery.bdary <- read.table("inst/exdata/greenery_bdary.txt", header=TRUE) # Bounding box of greenery area
  martyn.bdary   <- read.table("inst/exdata/martyn_park_bdary.txt", header=TRUE) # Martyn sports park
  tc.bdary       <- read.table("inst/exdata/martin_tennis_court.txt", sep=",", header=TRUE) #Tennis court next to Martyn sports park
  #severin.st     <- read.table("inst/exdata/severinst.txt", sep=",", header=FALSE) # Now included in microclim/roads.txt
  #roads          <- read.table("inst/exdata/roads.txt", sep=",", header=FALSE) # Now included in microclim/roads.txt

  # The data taken from OSM ----
  newroads <- read.table("inst/exdata/microclim/roads.txt", sep=",", header=TRUE)
  river    <- read.table("inst/exdata/microclim/river.txt", sep=",", header=TRUE)
  ogreen   <- read.table("inst/exdata/microclim/ogreen.txt", sep=",", header=TRUE)
  oconc    <- read.table("inst/exdata/microclim/concrete.txt", sep=",", header=TRUE)
  trees    <- read.table("inst/exdata/microclim/trees.txt", sep=",", header = TRUE)
  builds   <- read.table("inst/exdata/microclim/buildings.txt", sep=",", header = TRUE)
  # ---- Construct grid of spacing gridSize
  grid <- make_grid(bdary, gridSize)

  # ---- Determine location of houses
  houses.sorted <- houses[with(houses, order(long)), ]
  # This process is much the same as constructing 'grid'
  ## but we don't use the functions in case we want to vary the temp devs of
  ## individual houses.
  houses.unique     <- unique(houses)
  houses.in.region  <- point.in.polygon(houses.unique$long, houses.unique$lat, bdary$long, bdary$lat)
  final.houses      <- cbind(houses.unique$lat[which(houses.in.region ==1)], houses.unique$long[which(houses.in.region == 1)]) #puts these grid spaces into a df
  final.houses      <- as.data.frame(final.houses)
  houses.sorted     <- unique(final.houses[order(final.houses$V2, final.houses$V1), ])
  matched.IDs       <- lapply(1:nrow(houses.sorted), function(x){assign_grid(houses.sorted[x, ], grid)})
  grid.pos          <- unlist(matched.IDs)
  houses.matchedIDs <- as.data.frame(grid.pos)

  # ---- Determine grid spaces from data sets that are bounding boxes
  greenery.matchedIDs <- bb_to_matchedIDs(greenery.bdary, gridSize, grid)
  martyn.matchedIDs   <- bb_to_matchedIDs(martyn.bdary, gridSize, grid)
  tc.matchedIDs       <- bb_to_matchedIDs(tc.bdary, gridSize, grid)

  # ---- Determine grid spaces from data sets that are lists of points
  #severin.matchedIDs <- points_to_matchedIDs(severin.st, grid) # This dataset is deprecated
  #roads.matchedIDs   <- points_to_matchedIDs(roads, grid) # This dataset is deprecated

  # ---- Determine grid spaces from data sets that are lists of points: from OSM
  # The following code could be made into a function... one day.
  osmroads.dt         <- as.data.frame(cbind(newroads$X_lat, newroads$X_lon))
  osmroads.matchedIDs <- as.data.frame(points_to_matchedIDs(osmroads.dt, grid))
  osmroads.matchedIDs <- osmroads.matchedIDs[which(osmroads.matchedIDs != -1),] # Remove any -1 entries

  river.dt            <- as.data.frame(cbind(river$X_lat, river$X_lon))
  river.matchedIDs    <- as.data.frame(points_to_matchedIDs(river.dt, grid))
  river.matchedIDs    <- river.matchedIDs[which(river.matchedIDs != -1),] # Remove any -1 entries

  ogreen.dt           <- as.data.frame(cbind(ogreen$X_lat, ogreen$X_lon))
  ogreen.matchedIDs   <- as.data.frame(points_to_matchedIDs(ogreen.dt, grid))
  ogreen.matchedIDs   <- ogreen.matchedIDs[which(ogreen.matchedIDs != -1),] # Remove any -1 entries

  oconc.dt            <- as.data.frame(cbind(oconc$X_lat, oconc$X_lon))
  oconc.matchedIDs    <- as.data.frame(points_to_matchedIDs(oconc.dt, grid))
  oconc.matchedIDs    <- oconc.matchedIDs[which(oconc.matchedIDs != -1),] # Remove any -1 entries

  trees.dt            <- as.data.frame(cbind(trees$X_lat, trees$X_lon))
  trees.matchedIDs    <- as.data.frame(points_to_matchedIDs(trees.dt, grid))
  trees.matchedIDs    <- trees.matchedIDs[which(trees.matchedIDs != -1),] # Remove any -1 entries

  builds.dt           <- as.data.frame(cbind(builds$X_lat, builds$X_lon))
  builds.matchedIDs   <- as.data.frame(points_to_matchedIDs(builds.dt, grid))
  builds.matchedIDs   <- builds.matchedIDs[which(builds.matchedIDs != -1),] # Remove any -1 entries

  # ---- Construct a dataframe of each grid space and land type
  tempdata <- setNames(data.frame(matrix(ncol=4, nrow = nrow(grid))),c("Longitude","Latitude","LandType", "TemperatureDeviate"))
  tempdata$Latitude  <- grid$V1
  tempdata$Longitude <- grid$V2
  tempdata$TemperatureDeviate <- 0
  tempdata$LandType <- 1

  # ---- Set the different land types. There should be noLandTypes of land types
  tempdata$LandType[unlist(ogreen.matchedIDs)]   <- 1
  tempdata$LandType[unlist(martyn.matchedIDs)]   <- 1
  tempdata$LandType[unlist(greenery.matchedIDs)] <- 2
  tempdata$LandType[unlist(trees.matchedIDs)]    <- 2
  tempdata$LandType[unlist(tc.matchedIDs)]       <- 3
  tempdata$LandType[unlist(osmroads.matchedIDs)] <- 3
  tempdata$LandType[unlist(oconc.matchedIDs)]    <- 3
  tempdata$LandType[unlist(river.matchedIDs)]    <- 4
  tempdata$LandType[unlist(builds.matchedIDs)]   <- 5
  tempdata$LandType[unlist(houses.matchedIDs)]   <- 6

  # ---- Set the temperature deviates
  for (i in 1:noLandTypes){
    tempdata$TemperatureDeviate[which(tempdata$LandType == i)] <- typeTempDevs[i]
  }

  # ---- Construct the EKM chart
  # Apply temperature deviates to daily temperatures for each land type
  # Temperature at a particular land type = average daily temperature + its associated temperature deviate
  landTypeTemps <- vector(mode = "list", length = noLandTypes) # Initialising lists
  kelvinTemps <- vector(mode = "list", length = noLandTypes) # Initialising lists
  for(i in 1:noLandTypes){
    landTypeTemps[[i]] <- typeTempDevs[i] + tempData
  }
  # Convert each temperature from Celsius to Kelvin
  for(i in 1:noLandTypes){
    kelvinTemps[[i]] <- typeTempDevs[i] + tempData + const$KELV_CONV
  }
  # Converting kelvinTemps to a df for easy calculation
  temp.df <- data.frame(matrix(unlist(kelvinTemps), nrow=length(kelvinTemps), byrow=T))

  # ---- Calculating growth charts for each life stage
  EKM_egg    <- 24*((const$RHO_E*(temp.df/298)*exp((const$HA_E/const$R)*((1/298) - (1/temp.df))))/(1 + exp((const$HH_E/const$R)*(((1/const$THALF_E)-(1/temp.df))))))
  EKM_larval <- 24*((const$RHO_L*(temp.df/298)*exp((const$HA_L/const$R)*((1/298) - (1/temp.df))))/(1 + exp((const$HH_L/const$R)*(((1/const$THALF_L)-(1/temp.df))))))
  EKM_pupal  <- 24*((const$RHO_P*(temp.df/298)*exp((const$HA_P/const$R)*((1/298) - (1/temp.df))))/(1 + exp((const$HH_P/const$R)*(((1/const$THALF_P)-(1/temp.df))))))
  EKM_gono   <- 24*((const$RHO_G*(temp.df/298)*exp((const$HA_G/const$R)*((1/298) - (1/temp.df))))/(1 + exp((const$HH_G/const$R)*(((1/const$THALF_G)-(1/temp.df))))))

  # To use this list, indexing works like: EKM_chart[[stage]][land_type, ][day]
  EKM_chart  <- list(EKM_egg, EKM_larval, EKM_pupal, EKM_gono)

  return(c(EKM_chart, tempdata))
}
