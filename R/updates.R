#' Random dispersal of a single agent according to distance formulae.
#'
#' Method we use: distance and bearing, Haversine method
#' Note we work in units of radians and kilometres.
#' \enumerate{
#'  \item{Calculate 'd', the distance travelled}
#'  \item{Calculate 'theta', the direction or bearing travelled}
#'  \item{Convert to lat and long.}
#' }
#' @param lat Latitude or 'y' coordinate of agent.
#' @param long Longitude or 'x' coordinate of agent.
#' @param lambda Shape parameter for distance calculation. ABC parameter.
#' @return Updated \code{lat} and \code{long} of agent, as a list.
random_dispersal <- function(lat, long, lambda){
  d        <- rexp(1,1/lambda) #formula given in documentation
  #d <- rexp(1,lambda)
  theta.d  <- runif(1,0,2*pi) #formula given in documentation - distance in radians
  # We need to get theta.d (which is in radians) and convert it into
  # "clockwise from north" bearing
 # degs <- theta.d*(180/pi)

  #ratio    <- d/(1000*6371.01) #ratio between distance and the radius of the earth
  ratio    <- d/(6371.01) # ratio between distance and the radius of the earth in km

  latrad   <- (lat*pi)/180 #converting from degrees to radians
  longrad  <- (long*pi)/180
  #thetarad <- (theta.d*pi)/180

  latconv  <- asin(sin(latrad)*cos(ratio)+cos(latrad)*sin(ratio)*cos(theta.d))
  #longconv <- longrad+atan2(sin(thetarad)*cos(latrad)*sin(ratio), cos(ratio)-sin(latrad)*sin(longrad))
  longconv <- longrad+atan2(sin(theta.d)*cos(latrad)*sin(ratio), cos(ratio)-sin(latrad)*sin(latconv))

  # for testing
  #radlat <- (latconv*180)/pi
 # radlong <- (longconv*180)/pi

  #'below uses DEGREES in angle for calculation and not radians
#  latconv  <- asin(sin(latrad)*cos(ratio)+cos(latrad)*sin(ratio)*cos(degs))
#  longconv <- longrad+atan2(sin(degs)*cos(latrad)*sin(ratio), cos(ratio)-sin(latrad)*sin(latconv))


# ""When converting radians back to degrees (deg = rad⋅180/π), West is negative if using signed decimal degrees"
  #' do we need to do anything about that???
  newlat   <- (latconv*180)/pi #converting from radians to degrees
  newlong  <- (longconv*180)/pi

  # for testing
 # test <- as.data.frame(c(lat, radlat, newlat))
 # test <- cbind(test, as.data.frame(c(long, radlong, newlong)))
 # colnames(test) <- c("lat", "long")


  new.pos  <- list(lat=newlat, long=newlong)
  #COMMENT THIS OUT BELOW AFTER TESTING
  #new.pos <- list(lat = newlat, long = newlong, d = d, theta = theta.d)
  #CHECK boundaries
  return(new.pos)
}

#' Trapping of mosquitoes.
#'
#' This algorithm cycles through every trap in the region that is
#' active at time t. Every adult mosquito within 'radius' phi/2 of the
#' trap has a 50% chance of being trapped.
#' ('Radius' here referring to a square of lenth phi.)
#' Returns list of mosquito IDs that are to be trapped.
#' @param trapID unique ID of trap
#' @param phi phi/2 is the radius of influence of the trap.
#' @param mozzie.dt data.table of adult agents
#' @return List of mosquito IDs that are to be trapped.
#' @export
mosquito_trapping <- function(trapID, phi, mozzie.dt){


}

#' Gets the amount by which the Enzyme Kinetic Score for an agent will increase.
#' This value is added to their current EKS.
#' This function works as a 'lookup' for the EKS data, which is all calculated
#'     at the beginning of the simulation- this function just grabs it
#'
#' @param stage Stage of agent. 1: egg, 2: larvae, 3: pupae, 4: adult.
#' @param timestep Current timestep.
#' @return Amount by which an agent's EKS will increase for that timestep.
update_enzyme <- function(stage, timestep){
  if(stage == 1){
    EKSUpdate <- EKMChart[[1]][timestep]
  }
  else if(stage == 2){
    EKSUpdate <- EKMChart[[2]][timestep]
  }

  else if(stage == 3){
    EKSUpdate <- EKMChart[[3]][timestep]
  } else{
    EKSUpdate <- EKMChart[[4]][timestep]
  }
  return(EKSUpdate)
}

#' Gets the amount by which the Enzyme Kinetic Score for an agent will increase.
#' This value is added to their current EKS.
#' This fuction is used when microclimates are enabled.
#' This function works as a 'lookup' for the EKS data, which is all calculated
#'     at the beginning of the simulation- this function just grabs it
#'
#' @param stage Stage of agent. 1: egg, 2: larvae, 3: pupae, 4: adult.
#' @param timestep Current timestep.
#' @param landType Land Lype that the agent is on.
#' @return Amount by which an agent's EKS will increase for that timestep.
update_enzyme_microclim <- function(stage, timestep, landType){
  if(stage == 1){
    EKSUpdate <- EKMChart[[1]][landType, ][timestep]
  }
  else if(stage == 2){
    EKSUpdate <- EKMChart[[2]][landType, ][timestep]
  }

  else if(stage == 3){
    EKSUpdate <- EKMChart[[3]][landType, ][timestep]
  } else{
    EKSUpdate <- EKMChart[[4]][landType, ][timestep] # Adults
  }
  return(as.numeric(EKSUpdate))
}

#' Finds a mate for a female mosquitoes of breeding age.
#'
#' We draw a box of 'radius' \code{k} around the agent.
#' \code{k} is an ABC parameter and usually between 11m to 100m.
#' Note that here we mainly deal with ID number, an actual data.table variable,
#'  and not the agent's index in the main data.table.
#'
#' @section Algorithm:
#' Algorithm for finding a mate:
#' \itemize{
#'  \item{If there are no males within box, return -1 for no mate.}
#'  \item{If there is 1 male, that is the new mate.}
#'  \item{If there is >1 male, randomly choose one.}
#' }
#' Note that there is a max. number of agents that a male can mate with in a day,
#'  given by \code{max_daily_mates}.
#' It is important in testing to ensure that the \code{gonoCycle} variable for
#'  each male mate is properly updated on a global scope.
#' We assume that females only mate once.
#'
#' @section Who can mate?:
#' Recall that in order to mate, a female mosquito must have:
#' \itemize{
#'  \item{\code{enzyme} between 1.2 and 1.8 inclusive}
#'  \item{\code{gonoCycle} 1, 2 or 3}
#' }
#'
#' @param femID ID of female looking for mate.
#' @param k Distance within which female mosquito 'searches' for mate. ABC parameter.
#' @param max_daily_mates Number of times a male can mate in a day.
#' @param mozzie.dt Data.table of adult mosquitoes.
#' @return \code{mateID} of agent, or -1 if agent does not find a mate.
find_mate <- function(femID, k, max_daily_mates, mozzie.dt){
  new.mate <- -1 #Function returns -1 if it does not find a mate
  print("stop 1")
  print(mozzie.dt)
  #print(mozzie.dt[["lat"]][femID])
  position <- as.numeric(c("lat"  = mozzie.dt$lat[femID],
                           "long" = mozzie.dt$long[femID]))
print("stop 1.5")
  #CHANGE: magic number to k
  bdary    <- c("latmin" = position[1]-0.00011, "latmax" = position[1]+0.00011,
                "longmin" = position[2]-0.00011 , "longmax" = position[2]+0.00011)
  #CHECK EKS of males that can mate?
  possible.mates <- which(mozzie.dt$lat >= bdary["latmin"] & mozzie.dt$lat <= bdary["latmax"] & mozzie.dt$long >= bdary["longmin"] & mozzie.dt$long <= bdary["longmax"] & mozzie.dt$gender == 0)
  no.bachelors   <- length(possible.mates)
 print("stop 2")
   #print(paste0("no.bachelors: ",no.bachelors))
  if(no.bachelors == 0){
    new.mate <- -1
    print(paste0(femID, " no mate found"))
  }else if(no.bachelors == 1){
    if(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates){
      new.mate <- -1 #even though there is only one male, he has mated too much today
      print(paste0(femID, " 1 mate found, unsuitable"))
    }else{
      new.mate <- as.integer(possible.mates)
      mozzie.dt$gonoCycle[new.mate] <<- mozzie.dt$gonoCycle[new.mate] + 1
      #return(mozziedf$ID[possible.mates]) #New mate is just the single male they found
      print(paste0(femID, " 1 mate found"))
    }
  }else{
    if(length(which(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates)) == 0){
      #Randomly permutes the list of possible mates and then picks the one at the top of the pile
      new.mate <- sample(possible.mates)[1]
      #Increment number of mates of male by 1
      mozzie.dt$gonoCycle[new.mate] <<- mozzie.dt$gonoCycle[new.mate] + 1
      print(paste0(femID, " multiple mates, all ok"))
    }else{
      #Remove males who have mated more than "max_daily_mates" number of times in a day
      possible.mates <- possible.mates[-(which(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates))]
      #Now we have one more condition to check for:
      ##If we drop some males from "possible.mates" we might end up dropping them all,
      ##so if we dropped them all we return -1
      if(length(possible.mates) == 0){
        #print(paste0(femID, " multiple mates, had to drop some AND then ended up with none"))
        new.mate <- 1
      }else{
        #Randomly permutes the list of possible mates and then picks the one at the top of the pile
        new.mate <- sample(possible.mates)[1]
        mozzie.dt$gonoCycle[new.mate] <<- mozzie.dt$gonoCycle[new.mate] + 1 #increment number of mates by 1
        print(paste0(femID, " multiple mates, had to drop some"))
      }
    }
  }
  return(new.mate)
}

#' Second try function for finding a mate, for pararellisation
#' Takes in a female mosquito's position, then finds a single mate for them.
#' Returns the ID of the male mate.
#' NB the input of this function is the row position of a female in the current
#' data.frame of agents. The output is the ID number (a column/variable) of a
#' mate, and not its row position.
#' @param femID ID of female mosquito to find a mate for.
#' @param k 2*k is the width of the bounding box drawn around the mosquito.
#' This function will attempt to find a mate within this box.
#' @param max_daily_mates Maximum number of times a male mosquito can mate in one day.
#' @return The ID number of the male mate.
find_a_mate <- function(femID, k, max_daily_mates, mozzie.dt){
  new.mate <- -1 #Function returns -1 if it does not find a mate

  # Grab position of female
  to.mate.pos    <- as.numeric(c("lat"  = mozzie.dt$lat[femID],
                                 "long" = mozzie.dt$long[femID]))
  # Draw boundary around female
  to.mate.bdary  <- c("latmin" = to.mate.pos[1] - k, "latmax" = to.mate.pos[1] + k,
                      "longmin" = to.mate.pos[2] - k , "longmax" = to.mate.pos[2] + k)
  #' Find males inside that boundary
  #' I find it easier to do a dplyr::filter and then change to ID numbers etc,
  #' system.time of the below on nrow(mozzie.dt) = 91795 gives:
  #'   user  system elapsed
  #' 0.008   0.008   0.018
  possible.mates <- filter(mozzie.dt, gender == 0, gonoCycle <= max_daily_mates,
                           lat >= to.mate.bdary["latmin"],
                           lat <= to.mate.bdary["latmax"],
                           long >= to.mate.bdary["longmin"],
                           long <= to.mate.bdary["longmax"])
  no_bachelors   <- nrow(possible.mates)

  if(no_bachelors == 0){
    #' no_bachelors == 0 means that there are no suitable mates.
    mateID <- -1
  }else if(no_bachelors == 1){
    #' no_bachelors == 1 means there is one suitable bachelor
    mateID <- as.integer(possible.mates$ID)
  }else{
    #' no_bachelors > 1 means there are multiple suitable bachelors
    #' We randomly permute them and then take the top of the pile.
    mateID <- as.integer(sample(possible.mates$ID)[1])
  }

  #' Return the ID :)
  return(mateID)
}



#' Takes a juvenile agent and splits it into adult ages
#' This is to simulate 'emergence' of juveniles form aquatic stage.
#' Many variables stay the same, like motherID, but new ones are added or changed
#' Agents also undergo one motility event when they emerge.
#' This becomes their new lat/long.
#'
#' @section Data.table variables and initialisation:
#' ID:     Unique ID number of agent.
#' gender: Male is 0, female is 1. Sampled by 1 Binomial trial as opposed to a
#'         Bernoulli trial as Bernoulli requires another package.
#' mateID: Unique ID of their mate. Since no new adults will have a mate
#'         yet, it is initialised as -1.
#'         Males will always have mateID as -1 since they can have multiple mates.
#' enzyme: Enzyme Kinetic Score. See \code{init} for explanation. Initialised
#'         at 0 since they have just moved up from previous stage.
#' age:    Age in days. Given from their juvenile agent entry.
#' gonoCycle: Gonotrophic cycle. Means something different for males and females.
#'            males: number of times they've mated in a day, to be reset daily
#'            females: how many times they've laid a clutch of eggs
#'            Starts at 0 for everyone since they've just emerged.
#' timeDeath: Timestep they died: initialised as -1 as they are alive.
#' typeDeath: Random mortality/trapped death/death due to old age: which type?
#' whereTrapped: in the event of trapped death, where did they die? -1 otherwise.
#' motherID: Unique ID of mother. Taken from their juvenile agent entry.
#' fatherID: Unique ID of father. Taken from their juvenile agent entry.
#' infStatus: 1 if they carry Wolbachia, 0 if no Wolbachia, -1 for CI
#'            There should not be any with -1 because they should not have hatched.
#' lat:      Initial north/south or 'y' coordinate of agent.
#'            Should start with 145. Determined by motility event.
#' long:     Initial east/west of 'x' coordinate of agent.
#'            Should start with -16. Determined by motility event.
#' @param juvID The index of the agent in the juvenile data.table to be converted to adult agents.
#' @param idStart Where to start 'counting' the ID numbers of new agents from.
#' @param pmale Probability of being male.
#' @param lambda Shape parameter for distance calculation. ABC parameter.
#' @return A data.table of \code{N} adult agents.
juv_to_adult <- function(juvID, idStart, pmale, lambda, juv.dt, grid.df){

  new.dt <- data.table(ID=idStart:(idStart+juv.dt$clutchSize[[juvID]]-1),
                       gender=numeric(juv.dt$clutchSize[[juvID]]),
                       lat=numeric(juv.dt$clutchSize[[juvID]]),
                       long=numeric(juv.dt$clutchSize[[juvID]]),
                       mateID=numeric(juv.dt$clutchSize[[juvID]]),
                       enzyme=numeric(juv.dt$clutchSize[[juvID]]),
                       age=numeric(juv.dt$clutchSize[[juvID]]),
                       gonoCycle=numeric(juv.dt$clutchSize[[juvID]]),
                       timeDeath=numeric(juv.dt$clutchSize[[juvID]]),
                       typeDeath=numeric(juv.dt$clutchSize[[juvID]]),
                       whereTrapped=numeric(juv.dt$clutchSize[[juvID]]),
                       motherID=numeric(juv.dt$clutchSize[[juvID]]),
                       fatherID=numeric(juv.dt$clutchSize[[juvID]]),
                       infStatus=numeric(juv.dt$clutchSize[[juvID]]),
                       releaseLoc=numeric(juv.dt$clutchSize[[juvID]]),
                       gridID = numeric(juv.dt$clutchSize[[juvID]]))
  #print(length(juv.dt))
#print(juv.dt$clutchSize)
  # new.dt <- data.table(ID=idStart:(idStart+juv.dt$clutchSize[juvID]-1),
  #                      gender=numeric(juv.dt$clutchSize[juvID]),
  #                      lat=numeric(juv.dt$clutchSize[juvID]),
  #                      long=numeric(juv.dt$clutchSize[juvID]),
  #                      mateID=numeric(juv.dt$clutchSize[juvID]),
  #                      enzyme=numeric(juv.dt$clutchSize[juvID]),
  #                      age=numeric(juv.dt$clutchSize[juvID]),
  #                      gonoCycle=numeric(juv.dt$clutchSize[juvID]),
  #                      timeDeath=numeric(juv.dt$clutchSize[juvID]),
  #                      typeDeath=numeric(juv.dt$clutchSize[juvID]),
  #                      whereTrapped=numeric(juv.dt$clutchSize[juvID]),
  #                      motherID=numeric(juv.dt$clutchSize[juvID]),
  #                      fatherID=numeric(juv.dt$clutchSize[juvID]),
  #                      infStatus=numeric(juv.dt$clutchSize[juvID]),
  #                      releaseLoc=numeric(juv.dt$clutchSize[juvID]),
  #                      gridID = numeric(juv.dt$clutchSize[juvID]))

  new.dt$gender <- as.integer(lapply(new.dt$gender, function(x) x<- rbinom(1,1,1-pmale))) #probability of male is calculated above. since female mozzies are represented by 1 (a success) we have 1-pmale

  #We assume that mozzies disperse a bit from their original position when they hatch
  #CHECK boundaries
  positions   <- lapply(1:juv.dt$clutchSize[[juvID]], function(x) random_dispersal(juv.dt$lat[[juvID]],juv.dt$long[[juvID]], lambda)) #creates a list of lats and longs
  positions   <- do.call(rbind,positions)
  new.dt$lat  <- as.numeric(positions[,1])
  new.dt$long <- as.numeric(positions[,2])
  #new.dt$gridID <- mapply(FUN = get_gridID, new.dt$lat, new.dt$long)
  new.dt$gridID <- mapply(FUN = get_gridID, new.dt$lat, new.dt$long, MoreArgs = list(gridlong = grid.df$V2, gridlat = grid.df$V1))

  new.dt$mateID      <- new.dt$mateID[new.dt$gender == 0] <- -1 #males don't have mates
  new.dt$enzyme      <- 0 #enzyme resets because they're moving to the next stage

  new.dt$age <- juv.dt$age[[juvID]] #CHANGE
  #new.dt$age <- 14 + t
  new.dt$gonoCycle    <- 0 #gonoCycle will start at 0 for everyone when they become an adult
  new.dt$whereTrapped <- -1

  #new.dt$motherID <- -1
  #new.dt$fatherID <- -1

  new.dt$motherID     <- juv.dt$mother[[juvID]]
  new.dt$fatherID     <- juv.dt$father[[juvID]]

  new.dt$infStatus    <- as.numeric(lapply(new.dt$infStatus, function(x) x<-rbinom(1,1,juv.dt$infProb[juvID])))
  #new.dt$infStatus <- juv.dt$infProb[[juvID]] #FIX
  new.dt$releaseLoc   <- -1
  new.dt$timeDeath    <- -1
  new.dt$typeDeath    <- -1

  return(new.dt)
}

#' Updates the development stage of juvenile agents.
#'
#' If a juvenile agent's Enzyme Kinetic Score is greater than 0.95,
#' they age up to the next development stage:
#' \itemize{
#' \item{Egg to larvae OR}
#' \item{Larvae to pupae.}
#' }
#' Pupae development into adult agents should have already been handled.
#' This function runs on the entire list of juvenile agents at once.
#' @param stage List of all juvenile stages from juv.dt.
#' @param enzyme List of all juvenile agent enzyme score from juv.dt.
#' @return Updated \code{stage} and \code{enzyme} of agents, as a list.
#' @export
update_juv_stage <- function(stage, enzyme){

  stage[which(enzyme >= 0.95)] <- stage[which(enzyme >= 0.95)] +1

}

#' Applies daily natural death rate to a juvenile clutch.
#'
#' Uses the exponential model of mortality (constant death rate at each timestep)
#' @param clutch.size Size of clutch.
#' @param alpha Natural death rate.
#' @return Updated \code{clutch.size}.
#' @export
resize_clutch <- function(clutch.size, alpha){
  new.size <- rbinom(1, size = clutch.size, prob = (1 - alpha))
  return(new.size)
}
