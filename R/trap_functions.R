#' Clean trapping data.
#' Adds 3 columns at the end of the data CURRENTLY UNIMPLEMENTED:
#' \describe{
#'  \item{FWProp}{Proportion of trapped females that carry Wolbachia}
#'  \item{MWProp}{Proportion of trapped males that carry Wolbachia}
#'  \item{WProp}{Proportion of trapped mosquitoes that carry Wolbachia}
#' }
#' Adapted from code I wrote to create plots.
#' @param trapped Input trapping data, as a dataframe
#' @return Dataframe of cleaned trapping data.
trap_clean <- function(trapped){
  trapped$TrapStatus <- trapped$TrapStatus == "OK" # Change OK/Fail column into boolean

  # Now we need to remove the 0:00:00 timestamps at the end of each date format
  # for reference: regex bit is ctest <- sub("\\ .*", "", c),
  ## where c <- as.character(trapped$DateCollected[1])

  # Remove entries corresponding to failed traps
  trapped <- trapped[(which(trapped$TrapStatus == TRUE)),]
  trapped <- trapped[-(which(is.na(trapped$FemaleAeaeg))), ] #remove entries of NA female counts
  trapped <- trapped[-(which(is.na(trapped$MaleAeaeg))), ] #remove entries fo NA male counts
  #Maybe fix this later- just because there are NA female there may still be a real number of males?
  #Update 17/1/2020: Well, checking >which(is.na(trapped$FemaleAeaeg)) gives 4 entries, all have
  ## NA # of females and NA # of males except the last entry, where males = 0. So should be fine

  # Removing the timestamp from the dates since there is no discrepancy between entries
  trapped$DateCollected <- sub("\\ .*", "", trapped$DateCollected)
  trapped$DateDeployed <- sub("\\ .*", "", trapped$DateDeployed)

  # Ensuring date is in a format that R plays nicely with
  trapped$DateCollected <- as.Date(trapped$DateCollected, format = "%d/%m/%Y")
  trapped$DateDeployed <- as.Date(trapped$DateDeployed, format = "%d/%m/%Y")

  # Get the number of females & males trapped over time
  # Not needed atm...
  #femalesTrapped <- aggregate(trapped$FemaleAeaeg,
                              #by = list(Category = trapped$DateCollected), FUN = sum)

  #malesTrapped <- aggregate(trapped$MaleAeaeg,
                            #by = list(Category = trapped$DateCollected), FUN = sum)

  return(trapped)
}

#' Sets up traps:
#' Determines where unique traps positions are.
#' @param trapped Dataframe of trap data as per function trap_clean
#' @return Dataframe of trap locations in lat/long and its GID.
trap_setup <- function(trapped){
  traploc <- cbind(trapped$Lat, trapped$Long)
  traploc <- cbind(traploc, trapped$GID)
  traploc <- cbind(traploc, trapped$ToBufDist)
  traploc <- as.data.frame(traploc)
  traploc <- unique(traploc)
  colnames(traploc) <- c("Lat", "Long", "GID", "BufferZone")
  return(traploc)
}

#' Set up output file for trapping
#' It's going to be in longform
#' Variables: Lat, Long, Sex, Day, Wolbachia status, Number trapped.
#' @param trapped Cleaned trapping data
#' @param traploc Vector of trap locations
#' @param number_of_traps Number of traps.
#' @param noTimeSteps Number of time steps.
#' @return Longform dataframe to record trapping data.
trap_data_setup <- function(trapped, traploc, number_of_traps, noTimeSteps){
  N <- number_of_traps * noTimeSteps * 2 * 2 #*2 for Male and Female, *2 for Wolbachia and no Wolbachia

  #trapoutput <- setNames(data.frame(matrix(ncol = 5, nrow = N)),
                         #c("lat", "long","sex", "day","wolbachia"))

  return(trapoutput)
}

#' Finds any possible adult mosquitoes that are trapped by a particular trap
#' for a particular timestep.
#' If mosquitoes are within distance phi of a trap, there is a 0.5 probability
#' that a mosquito is trapped.
#' Only adults can be trapped.
#' @param traplat Lat of a particular trap.
#' @param traplong Long of a particular trap.
#' @param phi Distance that a mosquito is susceptible to a trap.
#' @param mozzie.dt Data.table of adult mosquitoes.
#' @param prob Probability of being trapped when within trap radius.
#' @return Vector of mosquito IDs to be trapped. Expect to be length 0 for the first few days. If length 0 return -1.
find_trapped <- function(traplat, traplong, phi, mozzie.dt, prob){
  bb         <- c(traplat - phi, traplat + phi, traplong - phi, traplong + phi) # Draw a box around trap of distance phi

  #to.trap    <- which(mozzie.dt$lat >= bb[1] && mozzie.dt$lat <= bb[2] && mozzie.dt$long <= bb[3] && mozzie.dt$long >= bb[4])
  to.trap    <- which(mozzie.dt$lat >= bb[1] & mozzie.dt$lat <= bb[2] & mozzie.dt$long >= bb[3] & mozzie.dt$long <= bb[4] & mozzie.dt$typeDeath == -1)

  no.to.trap <- length(to.trap) # To calculate random vector
  #Change below to a bernoulli?
  rand       <- runif(no.to.trap, min = 0, max = 1) #Calculates a Uniform[0,1] number for each mozzie in trap radius
  to.trap    <- to.trap[which(rand <= prob)]

  trapped.mozz <- -1 #Error handling
  if(length(to.trap) != 0){
    trapped.mozz <- mozzie.dt$ID[to.trap]
  }
  return(trapped.mozz)
}

#' Gets Wolbachia infection proportion of traps at appropriate time intervals.
#' In reality, traps were delpoyed and then emptied on certain days.
#' So we can't just look at the proportion of Wolbachia-carriers over time for the
#' entire graveyard
#' @param trapClean Cleaned trapping data as per function trap_clean.
trap_empty <- function(trapClean){
# Indices of traps that were functional during the simulation period
#timelyTrapped <- trapClean[which(trapClean$DateCollected > as.Date("2013/01/01") & trapClean$DateCollected < as.Date("2013/04/01")), ]
# Edited to account for 110 days:
timelyTrapped <- trapClean[which(trapClean$DateCollected > as.Date("2013/01/01") & trapClean$DateCollected < as.Date("2013/04/20")), ]
noEmpties     <- length(unique(timelyTrapped$DateCollected)) # Number of days that traps are emptied
emptyDates <- unique(timelyTrapped$DateCollected) # The dates that traps are emptied
diff       <- as.numeric(abs(diff(as.Date(emptyDates)))) # This gives us a vector of the difference in days between each trap empty
diff <- c(0, diff)
diff <- cumsum(diff) # Entry i is the difference in days between that date and the last trap empty date

DayDeployed   <- as.numeric(as.Date(timelyTrapped$DateDeployed) - as.Date("2013/01/01"))
timelyTrapped <- as.data.frame(cbind(timelyTrapped, DayDeployed))
DayCollected  <- as.numeric(as.Date(timelyTrapped$DateCollected) - as.Date("2013/01/01"))
timelyTrapped <- as.data.frame(cbind(timelyTrapped, DayCollected))
DaysActive    <- as.numeric(as.Date(timelyTrapped$DateCollected) - as.Date(timelyTrapped$DateDeployed))
timelyTrapped <- as.data.frame(cbind(timelyTrapped, DaysActive))

#date.df <- seq(from = as.Date("2013/01/02"), to = as.Date("2013/04/01"), by = "day")
#days.vec <- seq(from = 1, to = 90)
#date.df <- as.data.frame(cbind(as.Date(date.df), days.vec))
#colnames(date.df) <- c("Date", "Day")

 return(timelyTrapped)
}

#' Returns a vector of the dates where we're "cleaning" out the traps
#' i.e. when we actually "observe" the data in the graveyard
#' @param trapClean Cleaned trapping data as per function trap_clean.
#' @return Vector of days where we empty the traps.
when_trap_empty <- function(trapClean){
  # Indices of traps that were functional during the simulation period
  # Account for 110 days: edit this code to not use magic numbers
  timelyTrapped <- trapClean[which(trapClean$DateCollected > as.Date("2013/01/01") & trapClean$DateCollected < as.Date("2013/04/20")), ]
  #timelyTrapped <- trapClean[which(trapClean$DateCollected > as.Date("2013/01/01") & trapClean$DateCollected < as.Date("2013/04/01")), ]
  noEmpties     <- length(unique(timelyTrapped$DateCollected)) # Number of days that traps are emptied
  emptyDates <- unique(timelyTrapped$DateCollected) # The dates that traps are emptied
  diff       <- as.numeric(abs(diff(as.Date(emptyDates)))) # This gives us a vector of the difference in days between each trap empty
  diff <- c(0, diff)
  diff <- cumsum(diff) # Entry i is the difference in days between that date and the last trap empty date
  return(diff)
}

#' Not all traps are active during the entire simulation.
#' This function determines if a trap is active at time t before it attempts to trap anything.
#' @param currentTraps Cleaned trapping data with active dates converted to numerical days
#' @return Boolean TRUE if active and FALSE if not active.
when_trap_active <- function(currentTraps){


}

