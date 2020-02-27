#' Creates an initial data.table of adult agents.
#' Since this is an initial spread, some data is made up/dummy:
#' for example, we don't know the ID of the parents.
#' Latitude/Longitude of agents is randomly determined through
#' another function \code{init_position} and added to the data.table
#' after construction.
#'
#' @section Data.table variables and initialisation:
#' \itemize{
#' \item{ID:     Unique ID number of agent.}}
#' gender: Male is 0, female is 1. Sampled by 1 Binomial trial as opposed to a
#'         Bernoulli trial as Bernoulli requires another package.
#' mateID: Unique ID of their mate. Since no initial mosquitoes will have a mate
#'         yet, it is initialised as -1.
#'         Males will always have mateID as -1 since they can have multiple mates.
#' enzyme: Enzyme Kinetic Score. See \code{init} for explanation. Initialised
#'         uniform randomly.
#' age:    Age in days. Initialised uniform randomly around what we would expect
#'         young adults to be.
#' gonoCycle: Gonotrophic cycle. Means something different for males and females.
#'            males: number of times they've mated in a day, to be reset daily
#'            females: how many times they've laid a clutch of eggs
#'            we roughly estimate the gonoCycle of females based on \code{age}.
#' timeDeath: Timestep they died: initialised as -1 as they are alive.
#' typeDeath: Random mortality/trapped death/death due to old age: which type?
#' whereTrapped: in the event of trapped death, where did they die? -1 otherwise.
#' motherID: Unique ID of mother. -1 for initial adult data.table
#' fatherID: Unique ID of father. -1 for initial adult data.table.
#' infStatus: 1 if they carry Wolbachia, 0 if no Wolbachia, -1 for CI
#'            for initial wild type we assume they all start with 0.
#' lat:      Initial north/south or 'y' coordinate of agent. Should start with -16.
#' long:     Initial east/west of 'x' coordinate of agent. Should start with 145.
#' @param N The number of initial adult agents.
#' @param pmale Probability of being male.
#' @return A data.table of \code{N} adult agents.
initialise_adults <- function(N,pmale){
  noVariables <- 13
  init.dt <- setNames(data.frame(matrix(ncol = noVariables, nrow = N)),
                 c("ID", "gender","mateID", "enzyme","age","gonoCycle","timeDeath","typeDeath","whereTrapped","motherID","fatherID","infStatus","releaseLoc"))

  init.dt$gender       <- as.integer(lapply(init.dt$gender, function(x) x<- rbinom(1,1,1-pmale))) #Need to convert to not be a list
  init.dt$ID           <- 1:N
  init.dt$age          <- as.integer(lapply(init.dt$age, function (x) x<-round(truncnorm::rtruncnorm(1,mean=20,sd=2,a=14,b=30))))
  init.dt$motherID     <- -1
  init.dt$fatherID     <- -1
  init.dt$releaseLoc   <- -1
  init.dt$timeDeath    <- -1
  init.dt$typeDeath    <- -1
  init.dt$infStatus    <- 0
  init.dt$whereTrapped <- -1
  init.dt$mateID       <- -1
  init.dt$enzyme       <- as.double(lapply(init.dt$enzyme, function(x) x <- runif(1, min=0, max=1)))  #CHANGE: want to use initial age distribution to get spread of initial enzymes

  #FIX this shit and make it work better :(
  init.dt$gonoCycle[which(init.dt$gender == 0)] <- 0 #males start off at 0 because 'gonoCycle' tracks number of mating events in a day for males
  init.dt$gonoCycle[which(init.dt$gender == 1 & init.dt$age >=14 & init.dt$age < 20)] <- 0
  init.dt$gonoCycle[which(init.dt$gender == 1 & init.dt$age >=20 & init.dt$age < 26)] <- 1
  init.dt$gonoCycle[which(init.dt$gender == 1 & init.dt$age >=26)] <- 2

  #initial position of each mosquito, then combine with rest of dataframe
  position.dt <- init_position(boundaryDat,N)
  init.dt     <- cbind(init.dt,position.dt)

  return(init.dt)
}

#' Creates an initial data.table of juvenile agents.
#' Since this is an initial spread, some data is made up/dummy:
#' for example, we don't know the ID of the parents.
#' Latitude/Longitude of agents is randomly determined and added to the
#'  data.table after construction.
#'
#' @section Data.table variables and initialisation:
#' \describe{
#' \item{motherID:}{ID of mother.}
#' \item{fatherID:}{ID of father.}
#' \item{age:}{Age in days. Initialised uniform randomly around what we would expect
#'              young adults to be.}
#' \item{stage:}{Development stage of clutch. 1: egg, 2: larvae, 3: pupae.}
#' \item{infProb:}{Probability of carrying Wolbachia. 0: no Wolbachia, -1: Cytoplasmic
#'               Incompatability, else infProb is nonzero.}
#' \item{lat:}{Initial north/south or 'y' coordinate of agent. Should start with -16.}
#' \item{long:}{Initial east/west of 'x' coordinate of agent. Should start with 145.}
#' \item{clutchSize:}{Number of juveniles in the clutch.}
#' \item{enzyme:}{Enzyme Kinetic Score of agents. Initialised uniform randomly. FIX}
#' \item{pDeath:}{Probability of death.}
#' }
#' @param Njuv The number of initial juvenile agents.
#' @return A data.table of \code{Njuv} juvenile agents.
initialise_juveniles <- function(Njuv){
  juv.dt <- setNames(data.frame(matrix(ncol = 8, nrow = Njuv)),
          c("motherID","fatherID","age","stage" ,"infProb", "clutchSize","enzyme","pDeath"))
  juv.dt$motherID <- -1
  juv.dt$fatherID <- -1
  juv.dt$age      <- as.integer(lapply(juv.dt$age, function(x) x <- round(runif(1, min=0, max=14),0))) #Uniform between 0 and 14 as a quick fix- FIX
  juv.dt$stage    <- mapply(FUN = init_juv_stage, juv.dt$age)

  #this can be replaced with a binary for inf/noninf since fringe cases lead to CI
  juv.dt$infProb  <- 0 #dummy setup, this will be sampled in RABC (also need correct value from Carla)

  juv.dt$clutchSize <- param$eta_1 #since these will always be wild type
  juv.dt$enzyme   <- as.double(lapply(juv.dt$enzyme, function(x) x <- runif(1, min=0, max=0.95))) #CHANGE: this should not be uniform
  juv.dt$pDeath   <- param$alpha_j #fix this. see removeNatDeath in Carla's code

  ##BUT for the moment we're just using the same method as the adult lat/long to get working prototype
  #juv.dt$lat  <- lapply(juv.dt$lat, function (x) x<-round(runif(1, min = min(boundaryDat$Lat), max = max(boundaryDat$Lat)),6))
  #juv.dt$long <- lapply(juv.dt$long, function (x) x<-round(runif(1, min = min(boundaryDat$Long), max = max(boundaryDat$Long)),6))

  #The following two lines calculate lat, long and gridID
  position.dt <- init_position(boundaryDat, Njuv)
  juv.dt      <- cbind(juv.dt,position.dt)


  return(juv.dt)
}

#' Creates a data.table of newly laid eggs to add to the juvenile data.table.
#'
#' Input is the indices of mothers attempting to lay eggs.
#' Much of this code is regarding the handling of CI:
#' \itemize{
#'  \item{Neither mother nor father carry Wolbachia: no Wolbachia in offspring}
#'  \item{Mother carries Wolbachia: offspring have nonzero probability of carrying}
#'  \item{Mother does not carry, father does: offspring suffer CI and won't hatch}
#' }
#' Latitude/Longitude of agents is determined as the same position that the
#'  mother is currently at.
#' Clutch sizes, as per literature, differ depending on Wolbachia status of mother.
#' This is handled by the ABC parameter eta_1 and eta_2.
#'
#' @section data.table variables and initialisation:
#' \itemize{
#' \item{motherID:   ID of mother.}
#' \item{fatherID:   ID of father.}
#' \item{age:        Age in days. Initialised at 0 since they're new!}
#' \item{stage:      Development stage of clutch. Should only be 1 since these are eggs.}
#' \item{infProb:    Probability of carrying Wolbachia. 0: no Wolbachia, -1: Cytoplasmic
#'                   Incompatability, else infProb is nonzero.}
#' \item{lat:        Initial north/south or 'y' coordinate of agent.
#'                   Should start with -16. Should be same as mother.}
#' \item{long:       Initial east/west of 'x' coordinate of agent.
#'                   Should start with 145. Should be same as mother.}
#' \item{clutchSize: Number of juveniles in the clutch.}
#' \item{pDeath:     Probability of death.}
#' }
#' @param toLay List of the indices of mothers attempting to lay a clutch.
#' @return A data.table of juvenile agents in stage 1 corresponding to each mother in \code{toLay}.
initialise_eggs <- function(toLay){
  noMothers <- length(toLay)
  eggs.dt <- data.table(motherID = numeric(noMothers), fatherID = numeric(noMothers), motherStatus = numeric(noMothers), fatherStatus = numeric(noMothers), age = numeric(noMothers), stage = numeric(noMothers), infProb = numeric(noMothers) , lat = double(noMothers), long = double(noMothers), clutchSize =  numeric(noMothers), enzyme = double(noMothers), pDeath = double(noMothers), gridID = integer(noMothers))

  eggs.dt$motherID <- mozziedf$ID[toLay]
  eggs.dt$fatherID <- mozziedf$mateID[toLay]
  eggs.dt$age    <- 0 #because they're new!
  eggs.dt$stage  <- 1 #because they're eggs!
  eggs.dt$lat    <- mozziedf$lat[toLay] #mother's position
  eggs.dt$long   <- mozziedf$long[toLay] #mother's position
  eggs.dt$gridID <- mapply(FUN = get_gridID, eggs.dt$lat, eggs.dt$long)

  # This code needs to be removed- fix code so I don't have to do this --------
  eggs.dt[eggs.dt=="NULL"] <- 0
  eggs.dt <- na.omit(eggs.dt) #CHANGE THIS it shouldnt be happening
  #eggs.dt <- eggs.dt[complete.cases(eggs.dt$infProb),]
  eggs.dt$clutchSize[which(mozziedf$infStatus[toLay] == 1)] <- eta_2 #finish
  eggs.dt$clutchSize[which(mozziedf$infStatus[toLay] == 0)] <- eta_1 #finish

  eggs.dt$enzyme <- 0 #Since they are new eggs, they haven't had the "chance" to accumulate enzyme yet

  # CI ----
  eggs.dt$fatherStatus <- mozziedf$infStatus[which(mozziedf$ID %in% mozziedf$mateID[toLay])] #should give a list of 1s and 0s corresponding to their father's Wolbachia status
  eggs.dt$motherStatus <- mozziedf$infStatus[toLay] #same but for mother

  #If mother is a Wolbachia carrier then offspring should be (with probability p_1)
  eggs.dt$infProb[which(eggs.dt$motherStatus == 1)] <- (1-p_1)
  eggs.dt$pDeath[which(eggs.dt$motherStatus == 1)]  <- 0.01 #natural death rate

  #FIX : HACKY FOR LOOP for ANZIAM 2/2/19

  motheruninf <- which(eggs.dt$motherStatus == 0)

  for(i in 1:length(motheruninf)){
    if(eggs.dt$fatherStatus[motheruninf[i]] == 0){
      eggs.dt$infProb[i] <- 0
      eggs.dt$pDeath <- alpha_j
    }
    else if(eggs.dt$fatherStatus[motheruninf[i]] == 1){
      eggs.dt$infProb[i] <- -1
      eggs.dt$pDeath[i] <- -1
    }
    else{
      #error handling case; offspring are just wild CHECK
      eggs.dt$infProb[i] <- 0
      eggs.dt$pDeath <- alpha_j
    }
  }
  #Get rid of "extra information" in data table that we no longer need
  eggs.dt <- eggs.dt[,c("motherStatus","fatherStatus") := NULL]

  return(eggs.dt)
}

#' Simulates a release of Wolbachia-carrying mosquitoes.
#' Releases in this experiment are always adults.
#' Much of this is functionally similar to initialise_adults.
#' Since this is a release, some data is made up/dummy:
#' for example, we don't know the ID of the parents.
#' Latitude/Longitude of agents is given by data supplied by the WMP.
#'
#' @section Data.table variables and initialisation:
#' \describe{
#' \item{ID:}{Unique ID number of agent.}
#' \item{gender:}{Male is 0, female is 1. We know proportions of female/males released
#'          from data.}
#' \item{mateID:}{Unique ID of their mate. Since no initial mosquitoes will have a mate
#'         yet, it is initialised as -1.
#'         Males will always have mateID as -1 since they can have multiple mates.}
#' \item{enzyme:}{Enzyme Kinetic Score. See \code{init} for explanation. Initialised
#'         uniform randomly.}
#' \item{age:}{Age in days. Initialised uniform randomly around what we would expect
#'         young adults to be.}
#' \item{gonoCycle:}{Gonotrophic cycle. Means something different for males and females.
#'            males: number of times they've mated in a day, to be reset daily
#'            females: how many times they've laid a clutch of eggs
#'            we roughly estimate the gonoCycle of females based on \code{age}.}
#' \item{timeDeath:}{Timestep they died: initialised as -1 as they are alive.}
#' \item{typeDeath:}{Random mortality/trapped death/death due to old age: which type?}
#' \item{whereTrapped:}{In the event of trapped death, where did they die? -1 otherwise.}
#' \item{motherID:}{Unique ID of mother. -1 since wild release.}
#' \item{fatherID:}{Unique ID of father. -1 since wild release.}
#' \item{infStatus:}{1 if they carry Wolbachia, 0 if no Wolbachia, -1 for CI.
#'            They should mostly be 1 with a few 0 due to incomplete transmission.}
#' \item{lat:}{Initial north/south or 'y' coordinate of agent. Should start with -16.}
#' \item{long:}{Initial east/west of 'x' coordinate of agent. Should start with 145.}
#' }
#' @param noReleased The number of initial adult agents.
#' @param noMale Number of males.
#' @param noFemale Number of females.
#' @param lat Latitude of release site.
#' @param long Longitude of release site.
#' @return A data.table of \code{N} adult agents.
initialise_release <- function(noReleased, noMale, noFemale, lat, long, idStart, GID){
  propInfRelease <- 1 #proportion of released mosquitoes that carry Wolbachia (CHECK: get actual value from WMP)
  noVariables <- 16 #number of variables used to track state of mozzie minus 2 (lat & long). this is so I make sure to remember to change it as necessary
  #18/7: change release.dt so we just have a numeric age for mosquitoes rather than timeBirth/timeAdult/timeDeath
  release.dt <- setNames(data.frame(matrix(ncol = noVariables, nrow = noReleased)), c("ID", "gender", "lat","long","mateID", "enzyme","age","gonoCycle","timeDeath","typeDeath","whereTrapped","motherID","fatherID","infStatus","releaseLoc","gridID"))

  release.dt$gender <- 0 #Initialise all as male
  release.dt$gender[1:noFemale] <- 1 #Make #noFemale entries female
  #release.dt$gender <- lapply(release.dt$gender, function(x) x <- rbinom(1,1,1 - pmale)) #probability of male is calculated above. since female mozzies are represented by 1 (a success) we have 1-pmale
  #release.dt$ID     <- (max(mozziedf$ID) + 1):((max(mozziedf$ID)) + noReleased) #unique ID for each mosquito #CHECK
  release.dt$ID <- idStart:(idStart+noReleased-1) #Unique ID for each mosquito

  #release.dt$age <- lapply(release.dt$age, function (x) x<-round(rtruncnorm(1,mean=20,sd=2,a=14,b=30)))
  release.dt$age        <- as.integer(lapply(release.dt$age, function(x) x <- round(runif(1,min=18,max=21),0)))
  release.dt$motherID   <- -1
  release.dt$fatherID   <- -1
  release.dt$releaseLoc <- GID
  #release.dt$long <- 145.760
  #release.dt$lat <- -16.918
  release.dt$mateID     <- -1
  release.dt$long       <- long
  release.dt$lat        <- lat
  releaseGrid           <- get_gridID(testlat = lat, testlong = long) #Every mozzie released will be at the same grid
  release.dt$gridID     <- releaseGrid
  release.dt$timeDeath  <- -1
  release.dt$typeDeath  <- -1
  #release.dt$infStatus    <- lapply(release.dt$infStatus, function(x) x <-rbinom(1,1,propInfRelease))
  release.dt$infStatus  <- 1
  #infStatus: 1 for wolbachia, 0 for no wolbachia, -1 for CI
  release.dt$whereTrapped <- -1

  release.dt$enzyme <- as.double(lapply(release.dt$enzyme, function(x) x <- runif(1, min=0, max=1)))  #CHANGE: want to use initial age distribution to get spread of initial enzymes

  release.dt$gonoCycle[which(release.dt$gender == 0)] <- 0 #males start off at 0 because 'gonoCycle' tracks number of mating events in a day for males
  release.dt$gonoCycle[which(release.dt$gender == 1 & release.dt$age >=14 & release.dt$age < 20)] <- 0
  release.dt$gonoCycle[which(release.dt$gender == 1 & release.dt$age >=20 & release.dt$age < 26)] <- 1
  release.dt$gonoCycle[which(release.dt$gender == 1 & release.dt$age >=26)] <- 2

  return(release.dt)

}

#' Initialises graveyard data.table
#' The graveyard is where the entries of all dead adult mozzies go
#' @return An empty data.table that's ready to be rbind()-ed
initialise_graveyard <- function(){

  graveyard <- data.table(
      ID           = integer(),
      gender       = integer(),
      mateID       = double(),
      enzyme       = double(),
      age          = integer(),
      gonoCycle    = numeric(),
      timeDeath    = numeric(),
      typeDeath    = numeric(),
      whereTrapped = double(),
      motherID     = numeric(),
      fatherID     = numeric(),
      infStatus    = numeric(),
      releaseLoc   = numeric(),
      long         = double(),
      lat          = double(),
      gridID       = numeric()
  )
return(graveyard)
}

#' Initialises juvenile graveyard data.table
#' The entries here are for egg clutches that won't hatch due to CI
#' This does not include juveniles that died of natural causes
#'   (since that's just a number in the juv.dt data.table)
#'  @return Initialised data.table of the juvenile graveyard.
initialise_juv_graveyard <- function(){

 juv.graveyard <- data.table(
   motherID   = numeric(),
   fatherID   = numeric(),
   age        = integer(),
   stage      = numeric(),
   infProb    = numeric(),
   clutchSize = numeric(),
   enzyme     = double(),
   pDeath     = numeric(),
   long       = double(),
   lat        = double(),
   gridID     = integer()
   )

  return(juv.graveyard)
}


