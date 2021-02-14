#' Creates an initial data.table of adult agents.
#' Since this is an initial spread, some data is made up/dummy:
#' for example, we don't know the ID of the parents.
#' Latitude/Longitude of agents is randomly determined through
#' another function init_position and added to the data.table
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
#' enzyme: Enzyme Kinetic Score. See init for explanation. Initialised
#'         uniform randomly.
#' age:    Age in days. Initialised uniform randomly around what we would expect
#'         young adults to be.
#' gonoCycle: Gonotrophic cycle. Means something different for males and females.
#'            males: number of times they've mated in a day, to be reset daily
#'            females: how many times they've laid a clutch of eggs
#'            we roughly estimate the gonoCycle of females based on age.
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
#' @return A data.table of N adult agents.
initialise_adults <- function(N, pmale, boundaryDat, grid.df){
  noVariables <- 13
  init.dt <- setNames(data.frame(matrix(ncol = noVariables, nrow = N)),
                 c("ID", "gender","mateID", "enzyme","age","gonoCycle","timeDeath","typeDeath","whereTrapped","motherID","fatherID","infStatus","releaseLoc"))

  init.dt$gender       <- as.integer(lapply(init.dt$gender, function(x) x<- rbinom(1,1,1-pmale))) #Need to convert to not be a list
  init.dt$ID           <- 1:N
  init.dt$age          <- as.integer(lapply(init.dt$age, function (x) x<-round(truncnorm::rtruncnorm(1,mean=20,sd=2,a=15,b=30))))
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
  #print(head(grid.df))
  position.dt <- init_position(boundaryDat, N, grid.df)
  init.dt     <- cbind(init.dt, position.dt)

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
#' @return A data.table of Njuv juvenile agents.
initialise_juveniles <- function(Njuv, param, grid.df, boundaryDat){
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
  position.dt <- init_position(boundaryDat, Njuv, grid.df)
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
#' \describe{
#' \item{motherID:}   {ID of mother.}
#' \item{fatherID:}   {ID of father.}
#' \item{age:}        {Age in days. Initialised at 0 since they're new!}
#' \item{stage:}      {Development stage of clutch. Should only be 1 since these are eggs.}
#' \item{infProb:}    {Probability of carrying Wolbachia. 0: no Wolbachia, -1: Cytoplasmic
#'                   Incompatability, else infProb is nonzero.}
#' \item{lat:}        {Initial north/south or 'y' coordinate of agent.
#'                   Should start with -16. Should be same as mother.}
#' \item{long:}       {Initial east/west of 'x' coordinate of agent.
#'                   Should start with 145. Should be same as mother.}
#' \item{clutchSize:} {Number of juveniles in the clutch.}
#' \item{pDeath:}     {Probability of death.}
#' }
#' @param toLay List of the indices of mothers attempting to lay a clutch.
#' @param eta_1 Clutch size of non-Wolbachia carrying mothers.
#' @param eta_2 Clutch size Wolbachia carrying mothers.
#' @param p_1   Probability of complete maternal Wolbachia transmission.
#' @param alpha_j Juvenile mortality rate.
#' @return A data.table of juvenile agents in stage 1 corresponding to each mother in toLay.
initialise_eggs <- function(toLay, eta_1, eta_2, p_1, alpha_j, mozzie.dt, graveyard){
  #print(paste0(eta_1, " ", eta_2, " ", p_1, " ", alpha_j, " "))
  noMothers <- length(toLay)
  #' Subset mozzie.dt to just include mothers who are laying
  #' this makes it easy to assign variable values on to children
  mothers.dt <- filter(mozzie.dt, ID %in% toLay)
  eggs.dt <- data.table(motherID = numeric(noMothers),
                        fatherID = numeric(noMothers),
                        motherStatus = numeric(noMothers),
                        fatherStatus = numeric(noMothers),
                        age = numeric(noMothers),
                        stage = numeric(noMothers),
                        infProb = numeric(noMothers) ,
                        lat = double(noMothers),
                        long = double(noMothers),
                        clutchSize =  numeric(noMothers),
                        enzyme = double(noMothers),
                        pDeath = double(noMothers),
                        gridID = integer(noMothers))


  # 28/6/20: changing this to use MOTHER ID and NOT THE INDEX YOU GET FROM which()

  #eggs.dt$motherID <- toLay
  #eggs.dt$fatherID <- mozzie.dt$mateID[toLay]

  #'try and use mothers.dt as much as possible as this ensures we
  #'are keeping everything in the right order
  eggs.dt$motherID <- mothers.dt$ID
  eggs.dt$fatherID <- mothers.dt$mateID
  #eggs.dt$fatherID <- mozzie.dt$mateID[mozzie.dt$ID  %in% toLay]

  eggs.dt$age    <- 0 #because they're new!
  eggs.dt$stage  <- 1 #because they're eggs!
 #eggs.dt$lat    <- mozzie.dt$lat[toLay] #mother's position
  #eggs.dt$long   <- mozzie.dt$long[toLay] #mother's position

  # eggs.dt$lat  <- mozzie.dt$lat[mozzie.dt$ID  %in% toLay]
  # eggs.dt$long <- mozzie.dt$long[mozzie.dt$ID  %in% toLay]

  #' Clutch will be initialised at the same position as its mother
  eggs.dt$lat    <- mothers.dt$lat
  eggs.dt$long   <- mothers.dt$long
  eggs.dt$gridID <- mothers.dt$gridID

  # motherStatus is just the mother's infStatus
  eggs.dt$motherStatus <- mothers.dt$infStatus
  eggs.dt$enzyme <- 0 #Since they are new eggs, they haven't had the "chance" to accumulate enzyme yet
  eggs.dt$enzyme[which(eggs.dt$motherStatus == 1)] <- -1.1

  #clutch size determined by mother's infStatus:
  eggs.dt$clutchSize[which(eggs.dt$motherStatus == 1)] <- eta_2
  eggs.dt$clutchSize[which(eggs.dt$motherStatus == 0)] <- eta_1

  #eggs.dt$gridID <- as.integer(mapply(FUN = get_gridID, eggs.dt$lat, eggs.dt$long))

  # This code needs to be removed- fix code so I don't have to do this --------
  # 21-09-20: I actually don't think we need this anymore: I think I fixed it :)
  eggs.dt[eggs.dt=="NULL"] <- 0
  eggs.dt <- na.omit(eggs.dt) #CHANGE THIS it shouldnt be happening
  #eggs.dt <- eggs.dt[complete.cases(eggs.dt$infProb),]

  #eggs.dt$clutchSize[which(mozzie.dt$infStatus[toLay] == 1)] <- eta_2 #finish
  #eggs.dt$clutchSize[which(mozzie.dt$infStatus[toLay] == 0)] <- eta_1 #finish
#  eggs.dt$clutchSize[which(mozzie.dt$infStatus[mozzie.dt$ID  %in% toLay] == 1)] <- eta_2 #finish
#  eggs.dt$clutchSize[which(mozzie.dt$infStatus[mozzie.dt$ID  %in% toLay] == 0)] <- eta_1 #finish


  # CI ----
       #eggs.dt$fatherStatus <- mozzie.dt$infStatus[which(mozzie.dt$ID %in% mozzie.dt$mateID[toLay])] #should give a list of 1s and 0s corresponding to their father's Wolbachia status

       #eggs.dt$fatherStatus <- mozzie.dt$infStatus[mozzie.dt$mateID[toLay]]
       #eggs.dt$motherStatus <- mozzie.dt$infStatus[toLay] #same but for mother

  # 5/7/20 THIS WORKS.... for one number
  #eggs.dt$fatherStatus[1] <- mozzie.dt$infStatus[which(mozzie.dt$ID %in% eggs.dt$fatherID[1])]
  # SO honestly.... just do a loop for now

  # 21-09-20 replace this:
  # eggs.dt$motherStatus <- mozzie.dt$infStatus[mozzie.dt$ID %in% toLay] # mother's Wolbachia status
  #
  # f.IDs       <- eggs.dt$fatherID # List of father IDs
  # alive.f.IDs <- which(f.IDs %in% mozzie.dt$ID) # Which eggs have alive fathers
  # grave.f.IDs <- which(f.IDs %in% graveyard$ID) # Which eggs have fathers in the graveyard
  # end replacement for 21-09-20

  #'below is me adapting the above for optimisation
  #'The idea here is to split eggs.dt into those who have alive dads and thos who
  #'have dead dads. I'll then get the wolbachia status (carry/not carry) for each agent
  #'and then rowbind them together at the end.
  alive.dads <- filter(eggs.dt, fatherID %in% mozzie.dt$ID)
  dead.dads <- filter(eggs.dt, fatherID %in% graveyard$ID)
  #'----- end adaptation

  # We need to do two loops: one for those who have alive fathers and one for those whose fathers are in the
#
#   system.time(
#   for(i in alive.f.IDs){
#     eggs.dt$fatherStatus[i] <- mozzie.dt$infStatus[which(mozzie.dt$ID %in% eggs.dt$fatherID[i])]
#   })
#
#   # Now for those with dads in the graveyard
#   system.time(
#   for(j in grave.f.IDs){
#     eggs.dt$fatherStatus[j] <- graveyard$infStatus[which(graveyard$ID %in% eggs.dt$fatherID[j])]
#   })


  #' EDIT 21-9-20: VECTORISE THE ABOVE CODE
  #' We need to do two loops: one for those who have alive fathers and one for those whose fathers are dead
  #' r is a vector of the father's infStatus for every juvenile clutch in alive.dads.
 system.time(
  r <- foreach (i = 1:length(alive.dads$fatherID), .combine = c) %dopar% {
    mozzie.dt$infStatus[which(mozzie.dt$ID == alive.dads$fatherID[i])]
  })
 #print(paste0(" unique alive dads: ", length(unique(alive.dads$fatherID))))
 #print(paste0("r ", length(r), " ", head(r)))
  alive.dads$fatherStatus <- r
  #' Now for those with dads in the graveyard

 system.time(
  s <- foreach (j = 1:length(dead.dads$fatherID), .combine = c) %dopar% {
    graveyard$infStatus[which(graveyard$ID == dead.dads$fatherID[j])]
  })
 #print(paste0("nrow graveyard: ", nrow(graveyard), "no unique grav: ", length(unique(graveyard$ID))))
 #print(paste0(" unique dead dads: ", length(unique(dead.dads$fatherID))))
 #print(paste0("r ", length(s), " ", head(s)))
 #print(dead.dads)

  dead.dads$fatherStatus <- s

  # smoosh these back together
  l <- list(alive.dads, dead.dads)
  eggs.dt <- rbindlist(l, use.names = TRUE)
  #eggs.dt <- rbind(alive.dads, dead.dads)


  # -----

              # testing out code 5/7/20: remove later ----
              #test <- eggs.dt$fatherID
              #test.l <- which(test %in% mozzie.dt$ID) # which eggs have alive fathers
              #test.g <- which(test %in% graveyard$ID) # which eggs have fathers in the graveyard

              #test.lf <- eggs.dt$fatherID[test.l] #fatherIDs of alive fathers
              #test.gf <- eggs.dt$fatherID[test.g] #fatherIDs of dead fathers

              #which mozzies have the ID of those alive fathers
              #fatherstat.lf <- mozzie.dt$infStatus[test.lf %in% mozzie.dt$fatherID]
              # ok so we need to handle dads that are in the graveyard....
              #fatherstat.l <- mozzie.dt$infStatus[test.l]
              #fatherstat.g <- graveyard$infStatus[test.g]

              #eggs.dt$fatherStatus[test.l] <- mozzie.dt$infStatus[test %in% mozzie.dt$ID]


              #eggs.dt$fatherStatus <- mozzie.dt$infStatus[which(mozzie.dt$ID %in% eggs.dt$fatherID)] #THIS AINT WORKIN
              # End 5/7/20 playing around ----

  # Case I: If mother is a Wolbachia carrier then offspring should be (with probability p_1)
  #eggs.dt$infProb[which(eggs.dt$motherStatus == 1)] <- (1-p_1)
  #CHANGED 4/1 ADDED IF STATEMENT

  #print(paste0("nrow eggs.dt: ", nrow(eggs.dt)))
  #print(paste0("length inf mothers: ", length(which(eggs.dt$motherStatus == 1))))
  #print(p_1)
  if(length(which(eggs.dt$motherStatus == 1)) > 0){
    eggs.dt$infProb[which(eggs.dt$motherStatus == 1)] <- p_1
    eggs.dt$pDeath[which(eggs.dt$motherStatus == 1)]  <- 1.1*alpha_j #natural death rate
  }

  #FIX : HACKY FOR LOOP for ANZIAM 2/2/19

  motheruninf <- which(eggs.dt$motherStatus == 0)


  # This is to remove the confusing if/then statements below
  which.wt <- which(eggs.dt$motherStatus == 0 & eggs.dt$fatherStatus == 0) #Which are Wild Type/no Wolbachia
  eggs.dt$infProb[which.wt] <- 0
  eggs.dt$pDeath[which.wt] <- alpha_j


  which.ci <- which(eggs.dt$motherStatus == 0 & eggs.dt$fatherStatus == 1)

  eggs.dt$infProb[which.ci] <- -1

  eggs.dt$pDeath[which.ci] <- -1

 # print("get up to here")
  #below is deprecated- remove
#  for(i in 1:length(motheruninf)){
#    if(eggs.dt$fatherStatus[motheruninf[i]] == 0){
#      eggs.dt$infProb[i] <- 0
#      eggs.dt$pDeath <- alpha_j
#    }
#    else if(eggs.dt$fatherStatus[motheruninf[i]] == 1){
#      eggs.dt$infProb[i] <- -1
#      eggs.dt$pDeath[i] <- -1
#    }
#    else{
      #error handling case; offspring are just wild CHECK
#      eggs.dt$infProb[i] <- 0
#      eggs.dt$pDeath <- alpha_j
#    }
#  }
  #Get rid of "extra information" in data table that we no longer need
  to.drop <- c("motherStatus","fatherStatus")
  #eggs.dt <- eggs.dt[,c("motherStatus","fatherStatus") := NULL] #this only works if new.eggs.dt is purely a data.table.
  #it seems like some kind of data.frame and data.table chimera, I'm terrified but will ignore it
  eggs.dt <- eggs.dt[ , !(names(eggs.dt) %in% to.drop)]

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
#' \item{enzyme:}{Enzyme Kinetic Score. See init for explanation. Initialised
#'         uniform randomly.}
#' \item{age:}{Age in days. Initialised uniform randomly around what we would expect
#'         young adults to be.}
#' \item{gonoCycle:}{Gonotrophic cycle. Means something different for males and females.
#'            males: number of times they've mated in a day, to be reset daily
#'            females: how many times they've laid a clutch of eggs
#'            we roughly estimate the gonoCycle of females based on age.}
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
#' @param idStart Index starting for new mosquito ID numbers.
#' @param GID Unique identifier for a release site. From data.
#' @return A data.table of N adult agents.
initialise_release <- function(noReleased, noMale, noFemale, lat, long, idStart, GID, propInfRelease, grid.df){
  propInfRelease <- 1 #proportion of released mosquitoes that carry Wolbachia (CHECK: get actual value from WMP)
  noVariables <- 16 #number of variables used to track state of mozzie minus 2 (lat & long). this is so I make sure to remember to change it as necessary
  #18/7: change release.dt so we just have a numeric age for mosquitoes rather than timeBirth/timeAdult/timeDeath
  release.dt <- setNames(data.frame(matrix(ncol = noVariables, nrow = noReleased)), c("ID", "gender", "lat","long","mateID", "enzyme","age","gonoCycle","timeDeath","typeDeath","whereTrapped","motherID","fatherID","infStatus","releaseLoc","gridID"))
#  release.dt <- data.table(ID = integer(noReleased), gender = numeric(noReleased), mateID = numeric(noReleased), enzyme = numeric(noReleased), age = numeric(noReleased), gonoCycle = numeric(noReleased), timeDeath = numeric(noReleased), typeDeath = numeric(noReleased), whereTrapped = numeric(noReleased), motherID = numeric(noReleased), fatherID = numeric(noReleased), infStatus = numeric(noReleased), releaseLoc = numeric(noReleased), long = numeric(noReleased), lat = numeric(noReleased), gridID = integer(noReleased) )


  #print(release.dt)
  release.dt$gender <- 0 #Initialise all as male
  release.dt$gender[1:noFemale] <- 1 #Make #noFemale entries female
  #release.dt$gender <- lapply(release.dt$gender, function(x) x <- rbinom(1,1,1 - pmale)) #probability of male is calculated above. since female mozzies are represented by 1 (a success) we have 1-pmale
  #release.dt$ID     <- (max(mozzie.dt$ID) + 1):((max(mozzie.dt$ID)) + noReleased) #unique ID for each mosquito #CHECK
  release.dt$ID <- idStart:(idStart+noReleased-1) #Unique ID for each mosquito

  #release.dt$age <- lapply(release.dt$age, function (x) x<-round(rtruncnorm(1,mean=20,sd=2,a=14,b=30)))
  release.dt$age        <- as.integer(lapply(release.dt$age, function(x) x <- round(runif(1, min=15, max=25),0)))
  release.dt$motherID   <- -1
  release.dt$fatherID   <- -1
  release.dt$releaseLoc <- GID
  #release.dt$long <- 145.760
  #release.dt$lat <- -16.918
  release.dt$mateID     <- -1
  release.dt$long       <- long
  release.dt$lat        <- lat
  #releaseGrid           <- get_gridID(testlat = lat, testlong = long) #Every mozzie released will be at the same grid
  releaseGrid           <- get_gridID(testlat = lat, testlong = long, gridlong = grid.df$V2, gridlat = grid.df$V1) #Every mozzie released will be at the same grid
  release.dt$gridID     <- releaseGrid
  release.dt$timeDeath  <- -1
  release.dt$typeDeath  <- -1
  release.dt$infStatus    <- lapply(release.dt$infStatus, function(x) x <-rbinom(1,1,propInfRelease))
  #release.dt$infStatus  <- 1
  #infStatus: 1 for wolbachia, 0 for no wolbachia, -1 for CI
  release.dt$whereTrapped <- -1

  release.dt$enzyme <- as.double(lapply(release.dt$enzyme, function(x) x <- runif(1, min = 0, max = 0.5)))
  release.dt$enzyme <- release.dt$enzyme - 2.5
  #CHANGE: want to use initial age distribution to get spread of initial enzymes

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


