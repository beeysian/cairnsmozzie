#' Clean release data.
#' Since the method for estimating mosquito numbers is inprecise,
#' there are non-integer numbers of mosquitoes in each column.
#' Furthermore, this data set has no dates, rather the "day" as per Carla's work,
#' so I'm adding the dates back in.
#' @param releaseDat Mosquito release data.
#' @return Dataframe of cleaned release data.
release_clean <- function(releaseDat){
  #datevec <- c("10/1/2013", "17/1/2013", "25/01/2013", "31/1/2013", "7/2/2013", "14/2/2013", "21/2/2013", "28/2/2013", "7/3/2013", "14/3/2013", "21/3/2013", "28/3/2013", "4/4/2013","11/4/2013", "18/4/2013")
  datevec <- c("10/1/2013", "17/1/2013", "31/1/2013", "7/2/2013", "13/2/2013", "21/2/2013", "28/2/2013", "7/3/2013", "14/3/2013", "21/3/2013", "25/3/2013", "28/3/2013")
  release_dates <- as.Date(datevec, format = "%d/%m/%Y") # Same date format as trapping data, in case we want to plot
  releaseDat$NoMozzie <- round(releaseDat$NoMozzie) # Integer # of mozzies released
  releaseDat$NoFemale <- round(releaseDat$NoFemale) # Integer # of females released
  releaseDat$NoMale   <- releaseDat$NoMozzie - releaseDat$NoFemale # So noFemale + noMales= noReleased

  #Adding date vector
  Date  <- releaseDat$Day
  dayunique <- unique(Date)
  dayunique <- dayunique[order(Date)]
  for(i in 1:length(dayunique)){
    Date[which(Date == dayunique[i])] <- datevec[i]
  }
  releaseDat <- as.data.frame(cbind(releaseDat, Date))

  return(releaseDat)
}
