#' Clean release data.
#' Since the method for estimating mosquito numbers is inprecise,
#' there are non-integer numbers of mosquitoes in each column.
#' Furthermore, this data set has no dates, rather the "day" as per Carla's work,
#' so I'm adding the dates back in.
#' @param releaseDat Mosquito release data.
#' @return Dataframe of cleaned release data.
release_clean <- function(releaseDat){
  #datevec <- c("10/1/2013", "17/1/2013", "25/01/2013", "31/1/2013", "7/2/2013", "14/2/2013", "21/2/2013", "28/2/2013", "7/3/2013", "14/3/2013", "21/3/2013", "28/3/2013", "4/4/2013","11/4/2013", "18/4/2013")
  colnames(releaseDat) <- c("GID", "Lat", "Long", "ReleaseSequence", "NoReleaseCups",
                            "NoMozzie", "NoFemale", "NoMale", "Date") #cohesive column names
  # Update 25/10/20: we use the real full data now so no need for manual dates B)
  #datevec <- c("10/1/2013", "17/1/2013", "31/1/2013", "7/2/2013", "13/2/2013", "21/2/2013", "28/2/2013", "7/3/2013", "14/3/2013", "21/3/2013", "25/3/2013", "28/3/2013")
  #release_dates <- as.Date(datevec, format = "%d/%m/%Y") # Same date format as trapping data, in case we want to plot

  #datevec <- releaseDat$Date
  #release_dates <- as.Date(datevec, format = "%d/%m/%Y")
  release_dates <- dmy(releaseDat$Date)


  releaseDat$NoMozzie <- round(releaseDat$NoMozzie) # Integer # of mozzies released
  releaseDat$NoFemale <- round(releaseDat$NoFemale) # Integer # of females released
  releaseDat$NoMale   <- releaseDat$NoMozzie - releaseDat$NoFemale # So noFemale + noMales= noReleased

  #' Convert Dates to simulation "Days" (1 to 110)
  #'  we use 110 because that's the number of days in the simulation.... currently
  daysvec <- seq(from = 1, to = 110, by = 1)
  datesvec <- seq(from = as.Date("2013-01-03"), to = as.Date("2013-04-21"), "days")

  uniquedates   <- unique(release_dates) # The dates when releases are made
  release_days <- which(datesvec %in% uniquedates) # The days in the simulation these correspond to
  Day <- rep(0, times = 110) #initialising vector

  for(i in 1:length(uniquedates)){
    Day[which(release_dates == uniquedates[i])] <- release_days[i]
  }
  releaseDat <- as.data.frame(cbind(releaseDat, Day))

  # below is old stuff ----

  #Adding date vector
  #Date  <- releaseDat$Day
  #dayunique <- unique(Date)
  #dayunique <- dayunique[order(Date)]
  #for(i in 1:length(dayunique)){
  #  Date[which(Date == dayunique[i])] <- datevec[i]
  #}
  #releaseDat <- as.data.frame(cbind(releaseDat, Date))

  return(releaseDat)
}
