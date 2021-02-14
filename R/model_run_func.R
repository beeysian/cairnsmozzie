#' Function that runs the ABM, for use in ABC etc.
#' This version does NOT run microclimates.
#' Input is the ABC parameters
#' Output are vectors of interest over time.
#' @param lambda Shape param for dispersal
#' @param k Parameter for mating
#' @param phi Parameter for trapping
#' @param alpha_a Mortality rate of adult mosquitoes
#' @param alpha_j Mortality rate of juvenile mosquitoes
#' @param eta Clutch size of each egg-laying event.
model_run <- function(lambda, k, phi, alpha_a, alpha_j, eta){

  #' for debugging:
  #' lambda <- param$lambda
  #'  k <- param$k
  # phi <- param$phi
  # alpha_a <- param$alpha_a
  # alpha_j <- param$alpha_j
  # eta <- param$eta_1

  param <- data.frame(lambda = lambda,
                         k = k,
                         phi = phi,
                         alpha_a = alpha_a,
                         alpha_j = alpha_j,
                         eta_1 = eta,
                         eta_2 = eta,
                         pmale = 0.45,
                         p_1 = 0.93)

  # ----------------- Read in data files --------------------------------------------------
  # Boundary data for the boundary of the area used in trial (given as a set of lat/long points)
  boundaryDat <- read.table("inst/exdata/New_pp_Bound.txt", header=TRUE)
  # Bounding box of boundary, used for microclimates
  boundary.bb <- read.table("inst/exdata/testboundary.txt", header=TRUE)
  # Min daily temperatures from Cairns Aero weather station in 2013 provided by BOM
  minTemp     <- read.table("inst/exdata/cairnsAero13min.txt", header=TRUE)
  # Max daily temperatures from Cairns Aero weather station in 2013 provided by BOM
  maxTemp     <- read.table("inst/exdata/cairnsAero13max.txt", header=TRUE)
  # constants used for EKM calculations
  const       <- read.table("inst/exdata/constants.txt", header=TRUE)
  #Data regarding mosquito releases
  releaseDat  <- read.table("inst/exdata/PP_release.txt", header=TRUE)
  # Including all releases within the 110 day period
  releaseDatExtended <- read.table("inst/exdata/releases_combined_final.txt", header = TRUE)
  # Trapping data for mosquitoes
  trappingDat <- read.table("inst/exdata/PP_trapping_data.txt", sep=",", header=TRUE)
  # --------------- End reading in files --------------------------------------------------

  # --------------- Set global constants --------------------------------------------------
  testing         <- 1 # Set to 0 if we aren't doing any testing/optimisation
  max_age         <- 35 #mozzies who reach this age automatically die naturally
  max_juv_age     <- 14 #juveniles who reach this age automatically become adults
  max_pop_size    <- 200000 #if adult population exceeds this size, terminate simulation
  init_prop_inf   <- 0 #initial proportion of wild mozzies carrying Wolbachia- we expect this to be 0
  max_daily_mates <- 2 # maximum number of mating events a MALE mosquito can have in a day (from literature)
  include_microclimates <- 0 #Set to 1 to include microclimates effect, 0 to not include
  gridSize        <- 0.00025 # Size of grid spacing for microclimates
  #noLandTypes     <- 7 # Number of land types considered in microclimate analysis
  noLandTypes     <- 6 # Number of land types considered in microclimate analysis
  #typeTempDevs    <- c(-3, -7, 17, -4, 2, -2)
  trapProb        <- 0.4 # Probability of being trapped when within trap radius
  # --------------- End global constants ---------------------------------------------------

  # --------------- Variable initialisation ------------------------------------------------
  N           <- 15000 #number of adult mosquitoes at time t=0. To be initialised between 4000 to 20000
  Njuv        <- 6000 #number of juvenile CLUSTERS (not individuals) at time t=0
  idStart     <- N+1 #ID numbers of initial adults will be 1:N, then we start at N+1 for any new adult mosquito within simulation
  noTimeSteps <- 110 # Our simulation runs for 110 days
  t           <- 0 #timestep; starts from 0.
  # --------------- End variable initialisation --------------------------------------------

  # --------------- Calculate temperature and growth charts --------------------------------
  #Chart of average daily temperatures for each day of the simulation
  dailyTemps <- temperature_chart(minTemp, maxTemp, noTimeSteps)
  EKMChart   <- initialize_enzyme(dailyTemps, const) # Chart of daily updates for EKS
  # --------------- End chart calculation --------------------------------------------------

  # --------------- Set up traps -----------------------------------------------------------
  # we want to make sure they're active during the sim period
  trapClean       <- trap_clean(trappingDat) # Clean trapping data
  currentTraps    <- trap_empty(trapClean) # Active traps during the simulation period and day they're cleaned
  trapLoc         <- trap_setup(currentTraps) # Dataframe of trap locations in lat/long
  number_of_traps <- length(unique(currentTraps$GID))
  emptyDays       <- unique(currentTraps$DayCollected) # Days when we calculated Wolbachia proportion in graveyard
  daily_trapped   <- list() # Initialising list of trapped mosquitoes
  buf_0_traps     <- currentTraps %>% filter(ToBufDist == 0) # Traps in Buffer Zone 0 only
  # --------------- End trapping data setup ------------------------------------------------

  # --------------- Set up releases --------------------------------------------------------
  releaseDat   <- release_clean(releaseDatExtended) #Cleaning release data
  release.days <- as.data.frame(table(day = releaseDat$Day)) #Grabbing days that releases were made

  # --------------- End of release data wrangling ------------------------------------------

  # --------------- Set up grid ------------------------------------------------------------
  # We do this even if no microclimates used since it's useful for plotting
  grid.df <- make_grid(boundary.bb, gridSize)
  #print(head(grid.df))
  # --------------- End grid setup ---------------------------------------------------------

  # --------------- Set up model output ----------------------------------------------------
  #timesRepeat <- 50 #number of times to repeat simulation
  inf.prop <- seq(from = 1, to = noTimeSteps, by = 1) # Prop. infected over time (full region)
  inf.prop.buf0 <- seq(from = 1, to = noTimeSteps, by = 1) # Prop. infected over time (Buffer 0)
  grav.inf.prop <- seq(from = 1, to = noTimeSteps, by = 1) # Prop. infected in the graveyard
  no.juv.CI <- seq(from = 1, to = noTimeSteps, by = 1) # No. juvenile clutches with CI (cum)
  N.matrix <- seq(from = 1, to = noTimeSteps, by = 1) # N adults over time
  N.vec <- rep(0, times = noTimeSteps)
  juv.N.vec <- rep(0, times = noTimeSteps)
  #graveyard.plot.list <- list()
  # --------------- End model output setup -------------------------------------------------

  # --------------- Model initialisation/construction --------------------------------------
  mozzie.dt     <- as.data.table(initialise_adults(N, param$pmale, boundaryDat, grid.df)) #initial adult population
  juv.dt        <- as.data.table(initialise_juveniles(Njuv, param, grid.df, boundaryDat)) #initial juvenile population
  graveyard     <- initialise_graveyard() # Where entries of dead mosquitoes go.
  juv.graveyard <- initialise_juv_graveyard() # Where entries of eggs with CI go
  # --------------- End initial model construction -----------------------------------------

  # --------------- Optimisation/data checking stuff is below. -----------------------------
  #' For my personal use, feel free to comment out for your use.
  #' Each row will be a timestep, each column will be the time each step takes
  #' PLUS the number of agents.
  optimisation <-  data.frame(matrix(ncol = 10, nrow = noTimeSteps))
  max_gono <- rep(0, noTimeSteps)
  # ---------------  End optimisation framework --------------------------------------------

  t  <- 1 # Since model initialisation is complete
  # Progress bar
  pb <- progress_bar$new(format = " elapsed [:bar] :percent elapsed: :elapsed eta: :eta ", total = noTimeSteps)

  # --------------- Begin model simulation -------------------------------------------------

  for(t in 1:noTimeSteps){
    # -------------- Step 1: Juvenile eclosion ----------------------------------------------
    start <- Sys.time() # debugging

    to.adult <- which(juv.dt$stage == 3 & juv.dt$enzyme > 0.95 & juv.dt$infProb == 0)
    to.adult.inf <- which(juv.dt$stage == 3 & juv.dt$enzyme > 0.95 & juv.dt$infProb == param$p_1)

    #' ------------------------
    #'
    #' ADDED 12/11/20: EGG HATCH RATE OF 0.8 AS PER PERRAN ROSS
    #'
    #' This won't error if length(to.adult) == 0 hopefully
    to.adult <-  sample(to.adult, round(0.8*length(to.adult)))
    to.adult.inf <- sample(to.adult.inf, round(0.7*length(to.adult.inf)))
    to.adult <- c(to.adult, to.adult.inf)
    #' ------------------------
    #'
    if((length(to.adult)) != 0){
      IDvec <- juv.dt$clutchSize[to.adult] # Vector of agent ID indices
      IDvec <- c(0, IDvec)
      IDvec <- IDvec[1:length(to.adult)]
      IDvec <- cumsum(IDvec) + idStart

      #IDvec <- c(0, IDvec)
      #IDvec <- IDvec[1:length(to.adult)]
      # new.adults.dt <- mapply(FUN=juv_to_adult, to.adult, idStart, param$pmale, param$lambda, SIMPLIFY = FALSE)

      #new.adults.dt <- mapply(FUN=juv_to_adult, to.adult, IDvec, param$pmale, param$lambda, SIMPLIFY = FALSE)
      new.adults.dt <- mapply(FUN = juv_to_adult,
                              to.adult, IDvec, param$pmale, param$lambda,
                              MoreArgs = list(juv.dt, grid.df),
                              SIMPLIFY = FALSE)

      new.adults.dt %<>% rbindlist() %>% as.data.table()
      mozzie.dt     <<- rbind(mozzie.dt, new.adults.dt) #CHECK do we need <<- operator or can we just use <- ?
      # Delete row(s) of data table corresponding to the clutches that aged into adult mozzies
      juv.dt        <<- juv.dt[-to.adult]
      idStart       <- idStart + (nrow(new.adults.dt) + 1) # So we ID mosquitoes correctly
      rm(new.adults.dt)
    }

    # ---- Step 1.5: any juveniles moving to the next development stage? ----
    # If enzyme > 0.95, move on to the next development stage and reset enzyme score
    if(length(which(juv.dt$enzyme > 0.95)) > 0){
      # Increment development stage
      juv.dt$stage[which(juv.dt$enzyme > 0.95)]  <- juv.dt$stage[which(juv.dt$enzyme > 0.95)] + 1
      # Reset EkS to 0
      juv.dt$enzyme[which(juv.dt$enzyme > 0.95)] <- 0
    }
    # ---- End Step 1.5 ----
    rm(to.adult)
    rm(to.adult.inf)

    # Exception handling: if adult population size now exceeds max_pop_size, kill the simulation
    if(nrow(mozzie.dt) > max_pop_size){
      print("Termination condition 1: Maximum population size exceeded")
      break
    }
    #  s1 <- s1 +1 #debugging
    end <- Sys.time()
    optimisation[t,1] <- end - start
    rm(start, end)

    # --------------  End Step 1 ------------------------------------------------------------

    # --------------  Step 2: Egg laying ----------------------------------------------------
    start <- Sys.time() # debugging
    #Determine which mothers will lay their eggs today
    #Conditions to be satisfied for a mother to lay a clutch of eggs:
    ## 1) Has a mate
    ## 2) Enzyme Kinetic Score > 1 (should mean that mosquito has a mate)
    ## 3) First egg laying event is at EKS > 1, second at 1.58,
    ##    third at 2.16, etc. for multiples of 0.58 (from Focks et al '93)
    if(length(which(mozzie.dt$gonoCycle > 3))){
      print("gonoCycle too big")
      break
    }
    firstlay  <- mozzie.dt$ID[(which(mozzie.dt$enzyme >= 1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$gender == 1 & mozzie.dt$mateID != -1))]
    secondlay <- mozzie.dt$ID[(which(mozzie.dt$enzyme >= 1.58 & mozzie.dt$gonoCycle == 1 & mozzie.dt$gender == 1 & mozzie.dt$mateID != -1))]
    thirdlay  <- mozzie.dt$ID[(which(mozzie.dt$enzyme >= 2.16 & mozzie.dt$gonoCycle == 2 && mozzie.dt$gender == 1 && mozzie.dt$mateID != -1))]
    to.lay <- c(firstlay, secondlay, thirdlay)

    #' 22-09-20: reduce proportion of females laying eggs per day: CHECK
    to.lay <- sample(to.lay, size = length(to.lay)/5, replace = FALSE)

    # to.lay <- which(mozzie.dt$enzyme[ready.to.lay] >= 1 & mozzie.dt$gonoCycle[ready.to.lay]  == 0)
    #  to.lay <- which(mozzie.dt$gender == 1 & mozzie.dt$mateID != -1 & ((mozzie.dt$enzyme >= 1 & mozzie.dt$gonoCycle == 0) | (mozzie.dt$enzyme >= 1.58 & mozzie.dt$gonoCycle == 1) | (mozzie.dt$enzyme >= 2.16 & mozzie.dt$gonoCycle == 2)))

    if(length(to.lay) != 0){
      # Create new egg agents
      #print(paste0("test p_1: ", param$p_1))
      system.time(new.eggs.dt <- initialise_eggs(to.lay, param$eta_1, param$eta_2, param$p_1, param$alpha_j, mozzie.dt, graveyard))
      # If there are egg clutches with CI, remove them and add to juv.graveyard
      to.CI       <- which(new.eggs.dt$infProb == -1)
      #print(paste0("to CI: ", to.CI))
      if(length(to.CI) > 0){
        # Adding eggs with CI to the juv.graveyard
        l             <- list(juv.graveyard, new.eggs.dt[to.CI, ]) # CHANGED 4/1 ADDED COMMA
        juv.graveyard <- rbindlist(l, use.names = TRUE)
        # Remove CI eggs from new.eggs.dt
        new.eggs.dt   <- new.eggs.dt[!(to.CI),] #CHANGED 4/1 ADDED COMMA
        rm(l)
      }
      # Add new eggs to the juvenile data.table
      juv.dt      <- rbind(juv.dt,new.eggs.dt)
      rm(new.eggs.dt)
      rm(to.CI)

      # You need to increment gonoCycle by 1
      #I think the below works....
      mozzie.dt$gonoCycle[mozzie.dt$ID %in% to.lay] <- mozzie.dt$gonoCycle[mozzie.dt$ID %in% to.lay] + 1
      # mozzie.dt$gonoCycle[which(mozzie.dt$ID[to.lay] %in% mozzie.dt$ID)] <- mozzie.dt$gonoCycle[which(mozzie.dt$ID[to.lay] %in% mozzie.dt$ID)] + 1
    }
    rm(to.lay)
    rm(firstlay)
    rm(secondlay)
    rm(thirdlay)#Remove to.lay
    #  s2 <- s2 +1 #debugging

    end <- Sys.time()
    optimisation[t,2] <- end - start
    rm(start, end)
    # -------------- End Step 2 ------------------------------------------------------------

    # -------------- Step 3: Do we have mosquitoes being trapped (and killed?) -------------
    start <- Sys.time()
    for(tr in 1:number_of_traps){
      to.trap  <- list()
      #temptrap <- find_trapped(trapLoc$V1[tr], trapLoc$V2[tr], param$phi, mozzie.dt)
      temptrap <- find_trapped(trapLoc$Lat[tr], trapLoc$Long[tr], param$phi, mozzie.dt, trapProb)

      if(length(temptrap)> 1 && temptrap[1] != -1){
        # This is slow, work out how to use rbindlist for this in the future.
        # to.trap <- rbind(to.trap, temptrap)
        if(tr == 1){
          to.trap <- temptrap
          mozzie.dt$whereTrapped[mozzie.dt$ID %in% to.trap] <- trapLoc$GID[tr]
          #mozzie.dt[mozzie.dt$ID %in% to.trap,]$whereTrapped <- trapLoc$V3[tr]
        }else{
          to.trap <- append(to.trap, temptrap)
          mozzie.dt$whereTrapped[mozzie.dt$ID %in% to.trap] <- trapLoc$GID[tr]
          #mozzie.dt[mozzie.dt$ID %in% to.trap,]$whereTrapped <- trapLoc$V3[tr]
        }
        # Updating the entries of mozzies that have been trapped
        to.trap <- unlist(to.trap)
        mozzie.dt$typeDeath[mozzie.dt$ID %in% to.trap] <- 1 # Type of death is trapped death
        mozzie.dt$timeDeath[mozzie.dt$ID %in% to.trap] <- t # They died today
        #to.trap   <- unique(to.trap) #in case there are duplicates #CHECK this might be our mystery duplicate IDs
        l         <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.trap,])
        graveyard <- rbindlist(l) # Trapped mozzies go in the graveyard!
        mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.trap),]# Delete entries of mozzies who have been trapped and killed
        rm(l)
      }
      rm(temptrap)
    }
    # The following line of code is for testing: this removes all functional duplicates:
    ## test <- graveyard[!duplicated(graveyard[,!c('whereTrapped')]),]

    # For fun: plotting traps and agent locations
    #bdary <- read.table("inst/exdata/testboundary.txt", sep=" ", header=TRUE) #read in data
    #trapp <- ggplot(trapLoc, aes(x=V1, y=V2)) + geom_point()
    #trapp + geom_point(data=mozzie.dt, aes(x=mozzie.dt$lat, y=mozzie.dt$long, color="red")) +
    #geom_path(data = boundary.bb, aes(x=boundary.bb$lat, y=boundary.bb$long, color="green")) +
    #geom_path(data = boundaryDat, aes(x = boundaryDat$Lat, y=boundaryDat$Long, color="blue"))

    #  s3 <- s3 + 1 #debugging

    # Exception handling: if there are no more adults because they've all been trapped,
    ## kill the simulation
    if(nrow(mozzie.dt) == 0){
      print("Termination condition 3: All adults have died (traps)")
      break
    }
    end <- Sys.time()
    optimisation[t,3] <- end - start
    rm(start, end)
    # -------------- End Step 3 -------------------------------------------------------------

    # -------------- Step 4: Natural mosquito death (juvenile and adult) --------------------
    start <- Sys.time()
    ## Juvenile natural death: remove (alpha_j)% of the population
    ## The subsetting of juvdf$clutchSize is because we don't reduce clutch size of CI-affected eggs:
    ## they're already dead, but we just want to track how many there are so we leave clutchSize alone
    ## ^ The above statement is deprecated because CI-affected eggs now go to juv.graveyard,
    ## but it's still a useful sanity check to have in the code. So it stays
    juv.dt$clutchSize[which(juv.dt$pDeath != -1)] <- mapply(FUN = resize_clutch, juv.dt$clutchSize, param$alpha_j)
    #

    #' ADDED 14/11/20: TESTING MAKING JUV DEATH RATE HIGHER FOR THOSE WITH WB
    #juv.dt$clutchSize[which(juv.dt$pDeath != -1 & juv.dt$infProb == param$p_1)] <- mapply(FUN = resize_clutch, juv.dt$clutchSize[which(juv.dt$pDeath != -1 & juv.dt$infProb == param$p_1)], 1.01*param$alpha_j)

    # juv.dt$clutchSize[which(juv.dt$pDeath != -1)] <- mapply(FUN = resize_clutch,
    #                                                          juv.dt$clutchSize[which(juv.dt$pDeath != -1)],
    #                                                          juv.dt$pDeath[which(juv.dt$pDeath != -1)])

    # Remove juvenile agents where clutchSize == 0.
    # We don't need the information but let's put them in the graveyard anyway.
    # EDIT 12/3/20: don't add empty clutches to the graveyard and see if it improves computation.
    # Uncomment the first two lines under the if statement if we want to bring this back in.
    dead.clutch <- which(juv.dt$clutchSize == 0)
    #print(paste0("no dead clutch: ", length(dead.clutch)))
    #print(dead.clutch)
    if(length(dead.clutch) > 0){
      #  l      <- list(juv.graveyard, juv.dt[dead.clutch,]) #List of dead mozzies to remove
      #  juv.graveyard <- rbindlist(l, use.names = TRUE)
      juv.dt <- juv.dt[-(dead.clutch),]
      # rm(l)
    }
    rm(dead.clutch)

    # Adult (old age) death:
    ## Adults die automatically when they reach max_age
    ## Vector of mosquitoes that will die of old age:
    to.old.die <- which(mozzie.dt$age == max_age & mozzie.dt$timeDeath == -1)
    #print(paste0("no to.old.die: ", length(to.old.die)))
    if((length(to.old.die)) != 0){
     # print(paste0("to.old.die: ", to.old.die))
      #Track day that they died (timeDeath) and set typeDeath = 2 to denote natural (not trapped) death
      mozzie.dt$typeDeath[to.old.die] <- 2
      mozzie.dt$timeDeath[to.old.die] <- t

      #' add to graveyard
      old.to.graveyard <- filter(mozzie.dt, typeDeath == 2)
      l <- list(graveyard, old.to.graveyard)
      graveyard <- rbindlist(l, use.names = TRUE)
      rm(old.to.graveyard, l)

      #'remove from mozzie dataframe
     # mozzie.dt <- mozzie.dt[!(to.old.die),]
      mozzie.dt <- mozzie.dt[-(to.old.die),]
    }
    rm(to.old.die)
    #print(paste0("nrow mozzie after to.old.die:", nrow(mozzie.dt)))
    # Adults can also die as per natural death rate
    # CHECK: do we want the same number killed off every timestep?
    #mozzie.dt <- mozzie.dt[!sample(.N, round(param$alpha_a*nrow(mozzie.dt)))]


    #' THE BELOW LINE IS THE ORIGINAL to.nat.die
    #to.nat.die <- mozzie.dt[sample(.N, round(param$alpha_a*nrow(mozzie.dt)))]

    #' TESTING 05/11/20: RELEASED ADULTS WITH HIGHER DEATH RATE HHHHHHH --------
    to.nat.die.release <- mozzie.dt %>% filter(releaseLoc != -1) %>%
      sample_n(round(param$alpha_a*3*nrow(.)), replace = FALSE)
    to.nat.die <- mozzie.dt %>% filter(releaseLoc == -1) %>%
      sample_n(round(param$alpha_a*nrow(.)), replace = FALSE)

    #' END TESTING 05/11/20 HHHH ----------------------
    if(nrow(to.nat.die) >= 1){
      # mozzie.dt[mozzie.dt$ID %in% to.nat.die$ID,]$typeDeath <- 0
      # mozzie.dt[mozzie.dt$ID %in% to.nat.die$ID,]$timeDeath <- t
      # l <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.nat.die$ID,]) # List of dead mozzies to remove
      # graveyard <- rbindlist(l, use.names = TRUE)
      # mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.nat.die$ID),] # Removing naturally dead mozzies from mozzie.dt
      # rm(to.nat.die)
      # rm(l)
      #
      #print(paste0("no to nat die: ", nrow(to.nat.die)))

      mozzie.dt$typeDeath[which(mozzie.dt$ID %in% to.nat.die$ID)] <- 0
      mozzie.dt$timeDeath[which(mozzie.dt$ID %in% to.nat.die$ID)] <- t
      l <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.nat.die$ID,]) # List of dead mozzies to remove
      graveyard <- rbindlist(l, use.names = TRUE)
      #mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.nat.die$ID),] # Removing naturally dead mozzies from mozzie.dt
      #mozzie.dt <- mozzie.dt[-(mozzie.dt$ID %in% to.nat.die$ID),] # Removing naturally dead mozzies from mozzie.dt
      mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.nat.die$ID),]

      rm(to.nat.die)
      rm(l)
    }

    if(nrow(to.nat.die.release) >= 1){
     # print(paste0("no released to nat die: ", nrow(to.nat.die.release)))
      # mozzie.dt[mozzie.dt$ID %in% to.nat.die.release$ID,]$typeDeath <- 0
      # mozzie.dt[mozzie.dt$ID %in% to.nat.die.release$ID,]$timeDeath <- t
      # l <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.nat.die.release$ID,]) # List of dead mozzies to remove
      # graveyard <- rbindlist(l, use.names = TRUE)
      # mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.nat.die.release$ID),] # Removing naturally dead mozzies from mozzie.dt
      # rm(to.nat.die.release)
      # rm(l)
      #
      mozzie.dt$typeDeath[mozzie.dt$ID %in% to.nat.die.release$ID] <- 0
      mozzie.dt$timeDeath[mozzie.dt$ID %in% to.nat.die.release$ID] <- t
      l <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.nat.die.release$ID,]) # List of dead mozzies to remove
      graveyard <- rbindlist(l, use.names = TRUE)
      mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.nat.die.release$ID),] # Removing naturally dead mozzies from mozzie.dt
      rm(to.nat.die.release)
      rm(l)
    }

    # Exception handling: if all adults have died, kill the simulation
    if(nrow(mozzie.dt) == 0){
      print("Termination condition 2: All adults have died (natural)")
      #print(nrow(mozzie.dt))
      break
    }
    # Exception handling: if all juveniles have died, kill the simulation
    if(nrow(juv.dt) == 0){
      print("Termination condition 4: All juveniles have died (natural)")
      break
    }
    #  s4 <- s4 + 1 #debugging
    end <- Sys.time()
    optimisation[t,4] <- end - start
    rm(start, end)
    # -------------- End Step 4 -------------------------------------------------------------

    # -------------- Step 6 Take 2: a faster way to update microclimates --------------------
    start <- Sys.time()
    if(include_microclimates == 0){
      #' For adults/mozzie.dt: everyone has the same EKS update!
      EKS_today <- calculate_daily_EKS((dailyTemps[t] + 273.15), 4, const)
      mozzie.dt$enzyme <- mozzie.dt$enzyme + as.numeric(EKS_today)

      #'For juveniles, it's a little bit more involved:
      for(s in 1:3){
        #juv.dt$stage[which(juv.dt$stage == s)] <- mozzie.dt$enzyme[which(juv.dt$stage == s)] +
        #                                          calculate_daily_EKS((dailyTemps[t] + 273.15), s)
        juv_EKS_today <- as.numeric(calculate_daily_EKS((dailyTemps[t] + 273.15), s, const))
        #juv.dt$enzyme[which(juv.dt$stage == s)] <- juv.dt$enzyme[which(juv.dt$stage == s)] +
        #                                           as.numeric(juv_EKS_today)

        juv.dt$enzyme[which(juv.dt$stage == s)] <- sapply(juv.dt$enzyme[which(juv.dt$stage == s)],
                                                          function(x){x <- x + as.numeric(juv_EKS_today)},
                                                          simplify = TRUE)
        rm(juv_EKS_today)
      }

    }else if(include_microclimates == 1){
      #' Get the temperature deviate for each mozzie
      mozzie.w.tempdev <- left_join(mozzie.dt, landtype.df, by = "gridID")
      juv.w.tempdev    <- left_join(juv.dt, landtype.df, by = "gridID")

      #'Apply the temperature deviate to today's temperature and then convert to Kelvin
      #'This is the temperature that will then be fed into the EKM
      mozzie.w.tempdev$TemperatureDeviate <- mozzie.w.tempdev$TemperatureDeviate +
        dailyTemps[t] + 273.15
      juv.w.tempdev$TemperatureDeviate    <- juv.w.tempdev$TemperatureDeviate +
        dailyTemps[t] + 273.15
      #'Find EKS score
      #'first: for adults, stage = 4
      #' Find the EKS for each adult for today
      EKS_today <- calculate_daily_EKS(mozzie.w.tempdev$TemperatureDeviate, 4, const)
      #' Add to the current enzyme
      mozzie.dt$enzyme <- mozzie.dt$enzyme + as.numeric(EKS_today)

      #' now: for juveniles
      for(s in 1:3){
        #juv.dt$stage[which(juv.dt$stage == s)] <- mozzie.dt$enzyme[which(juv.dt$stage == s)] +
        #                                          calculate_daily_EKS((dailyTemps[t] + 273.15), s)
        juv_EKS_today <- calculate_daily_EKS(juv.w.tempdev$TemperatureDeviate[which(juv.w.tempdev$stage == s)], s)
        #' DOBULE CHECK the below works....
        juv.dt$enzyme[which(juv.dt$stage == s)] <- sapply(juv.dt$enzyme[which(juv.dt$stage == s)],
                                                          function(x){x <- x + as.numeric(juv_EKS_today)},
                                                          simplify = TRUE)

        #juv.dt$enzyme[which(juv.dt$stage == s)] <- juv.dt$enzyme[which(juv.dt$stage == s)] +
        #                                           as.numeric(juv_EKS_today)
        rm(juv_EKS_today)
      }
    }
    end <- Sys.time()
    optimisation[t,6] <- end - start
    rm(start, end)
    rm(EKS_today)
    #'
    #'
    #'
    # -------------- End Step 6 Take 2 :) ---------------------------------------------------

    # -------------- Step 8: Update numerical 'age' of agent --------------------------------
    start <- Sys.time()
    juv.dt$age    <- mapply('+', juv.dt$age, 1)
    mozzie.dt$age <- mapply('+', mozzie.dt$age, 1)
    #  s8 <- s8 +1 #debugging
    end <- Sys.time()
    optimisation[t,8] <- end - start
    rm(start, end)
    # -------------- End Step 8 -------------------------------------------------------------

    # -------------- Step 9: Mosquito releases ----------------------------------------------
    start <- Sys.time()
    today.releases <- list()
    #Check if it's one of the days where there is a mosquito release
    if(t %in% release.days$day){
      # For each day where a release was held, there are multiple releases in multiple locations
      # So we need to iterate through all of them.
      # Ideally, we should be able to replace this for loop with some vectorised code.
      temp.releases <- releaseDat[which(releaseDat$Day == t),] # Rows of releaseDat corresponding to today's releases
      rownames(temp.releases) <- 1:nrow(temp.releases)
      # We will run initialise_releases for each line of today.releases
      for(i in 1:nrow(temp.releases)){
        if(as.integer(temp.releases$NoMozzie[i] > 0)){
         # today.releases <- rbind(today.releases, initialise_release(as.integer(temp.releases$NoMozzie[i]),
         #                                                            temp.releases$NoMale[i],
        #                                                             temp.releases$NoFemale[i],
        #                                                             temp.releases$Lat[i],
        #                                                             temp.releases$Long[i],
        #                                                             idStart,
        #                                                             temp.releases$GID[i]), param$p_1, grid.df)
          today.releases <- rbind(today.releases, initialise_release(as.integer(temp.releases$NoMozzie[i]),
                                                                     temp.releases$NoMale[i],
                                                                     temp.releases$NoFemale[i],
                                                                     temp.releases$Lat[i],
                                                                     temp.releases$Long[i],
                                                                     idStart,
                                                                     temp.releases$GID[i], param$p_1, grid.df))
          idStart        <- idStart + (as.integer(temp.releases$NoMozzie[i]) + 1) #So we keep track of the right mosquito ID number
        }
      }
      # Mozzies all need to have a dispersal event at time of release
      # TO DO: update gridID as well!!!! * CHECK
      # Idea: perhaps try to do this in constructor? Not sure which one is faster
      today.releases.disp <- mapply(FUN = random_dispersal, today.releases$lat, today.releases$long, param$lambda) #DOES THIS WORK? LOL
      today.releases.disp <- do.call(rbind, today.releases.disp)
      today.releases.lat  <- today.releases.disp[c(TRUE, FALSE)]
      today.releases.long <- today.releases.disp[c(FALSE, TRUE)]
      today.releases$lat  <- today.releases.lat
      today.releases$long <- today.releases.long

      l         <- list(mozzie.dt, today.releases)
      mozzie.dt <- rbindlist(l, use.names = TRUE)

      rm(l)
      rm(today.releases)
      rm(today.releases.disp)
      rm(today.releases.lat)
      rm(today.releases.long)
    }

    # s9 <- s9 + 1 #debugging
    end <- Sys.time()
    optimisation[t,9] <- end - start
    rm(start, end)
    # -------------- End Step 9 -------------------------------------------------------------

    # -------------- Step 5: Mosquito mating ------------------------------------------------
    # Step 5 try 2: find_mate isn't working anymore since passing mozzie.dt doesn't actually pass the data.table
    # Which makes everything really frustrating. So I'll just do a loop for now :/
    # This is more or less all copied and pasted from the find_mate function.
    # This is a really obvious area to improve if we need to speed up simulation.
    # We only need to find mates for:
    ## Females that don't have mates and haven't laid eggs yet (first gono cycle)
    ## AND has EKS >= 1 (From Focks)
    start <- Sys.time()
    #' dataframe for optimisation:
    mate.opto <- data.frame(matrix(ncol = 10, nrow = noTimeSteps))

    startm <- Sys.time()
    to.mate <- which((mozzie.dt$gender == 1 & mozzie.dt$mateID == -1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$enzyme >= 1))
    endm <- Sys.time()
    mate.opto[t,1] <- endm - startm

    #' MATING ATTEMPT NUMBER 3...... I S2G THIS HAS TO BE THE LAST TIME
    #' we do a vectorised for loop. the find_a_mate function works on one agent at a time.

    #' Get a list of females who can mate.
    to.mate <- which((mozzie.dt$gender == 1 & mozzie.dt$mateID == -1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$enzyme >= 1))
    #' There's a ~30% chance that those females who can mate, will mate.
    to.mate <- sample(to.mate, size = length(to.mate)/5, replace = FALSE)

    system.time(
      mates_vector <- foreach (m = 1:length(to.mate), .combine=c) %dopar% {
        #' apply function to find mates
        find_a_mate(to.mate[m], param$k, max_daily_mates, mozzie.dt)
        #' Update female's mate ID
        #mozzie.dt$mateID[to.mate[m]]
        #' Update male's gonoCycle
      })
    # Update female's mate ID
    mozzie.dt$mateID[to.mate] <- mates_vector

    #PERHAPS JUST GET RID OF MALE'S GONOCYCLE? count the max number of times a male mates
    # below is for testing
    max_gono[t] <- max(tabulate(mates_vector))


    rm(to.mate, mates_vector)
    #  s5 <- s5 +1 #debugging
    end <- Sys.time()
    optimisation[t,7] <- end - start
    rm(start, end)
    # -------------- End Step 5 -------------------------------------------------------------

    # -------------- End of timestep procedures ---------------------------------------------
    optimisation[t, 10] <- nrow(mozzie.dt)
    pb$tick() # This is for the progress bar
    # -------------- End of optimisation stuff ----------------------------------------------

    # -------------- Post-simulation analytics and data collection --------------------------
    # The below is some quick and dirty calculations to get Wolbachia proportion over time
    # To do: make the below into a function
    if(length(which(mozzie.dt$infStatus == 0)) == 0){
      inf.prop[t] <- 1
    }else if(length(which(mozzie.dt$infStatus == 1)) == 0){
      inf.prop[t] <- 0
    }else{
      inf.prop[t] <- length(which(mozzie.dt$infStatus == 1))/nrow(mozzie.dt)
    }

    # proportion of trapped mosquitoes with wb over time
    if(nrow(graveyard) > 0){
      if(length(which(graveyard$infStatus == 0 & graveyard$typeDeath == 1 & graveyard$timeDeath == t)) == 0){
        grav.inf.prop[t] <- 1
      }else if(length(which(graveyard$infStatus == 1 & graveyard$typeDeath == 1 & graveyard$timeDeath == t)) == 0){
        grav.inf.prop[t] <- 0
      }else{
        grav.inf.prop[t] <- length(which(graveyard$infStatus == 1 & graveyard$typeDeath == 1 & graveyard$timeDeath == t))/length(which(graveyard$infStatus == 1 & graveyard$timeDeath == t))
      }
    }else{
      grav.inf.prop[t] <- 0
    }


    # number of juvenile clutches with CI over time
    if(nrow(juv.graveyard) > 0){
      no.juv.CI[t] <- length(which(juv.graveyard$infProb == -1))
    }

    juv.N.vec[t] <- nrow(juv.dt)
    # number of adult agents over time
    N.vec[t] <- nrow(mozzie.dt)

    #print(t)
    #print(nrow(mozzie.dt))
    # print(paste0("nrow mozzie: ", nrow(mozzie.dt)))
    # print(paste0("unique mozzie: ", length(unique(mozzie.dt$ID))))
    # print(nrow(juv.dt))
    # print(paste0("nrow grav: ", nrow(graveyard)))
    # print(paste0("unique grav ID: ", length(unique(graveyard$ID))))
    # print(paste0("unique grav entry: ", length(unique(graveyard))))
   # write.csv(graveyard, paste0("graveyard_",t,".csv"))
    #' Set up stuff to get graveyard data

    # ------------------------------------------------
    #
    # Set up graveyard/data collection stuff
    #
    # ------------------------------------------------

    trap_dates <- as.Date(unique(currentTraps$DateCollected))
    trap_dates <- trap_dates[which(trap_dates > as.Date("2013-01-03"))]

    daysvec <- seq(from = 1, to = 110, by = 1)
    datesvec <- seq(from = as.Date("2013-01-03"), to = as.Date("2013-04-21"), "days")

    uniquedates   <- unique(trap_dates) # The dates when traps are cleaned
    trap_days <- which(datesvec %in% uniquedates) # The days in the simulation these correspond to

    N_trapped_total   <- rep(0, length(trap_days))
    NWb_trapped_total <- rep(0, length(trap_days))
    NWT_trapped_total <- rep(0, length(trap_days))
    prop_wb_total     <- rep(0, length(trap_days))
    #' below is prop wb from buffer zone 0
    prop_0_wb         <- rep(0, length(trap_days))

    for(te in 1:length(trap_days)){

      if(te == 1){
        #filtered_trap   <- graveyard %>% filter(timeDeath <= trap_days[te],
        #                                        typeDeath == 1)
        #filtered_trap_0 <- graveyard %>% filter(timeDeath <= trap_days[te],
        #                                        typeDeath == 1,
        #                                        whereTrapped > -1)
        filtered_trap   <- graveyard %>% filter(timeDeath <= trap_days[te],
                                                whereTrapped > -1)
        #filtered_trap_0 <- graveyard %>% filter(timeDeath <= trap_days[te],
        #                                       whereTrapped > -1)
        filtered_trap_0 <- filtered_trap[which(filtered_trap$whereTrapped %in% buf_0_traps$GID),]

        #print(head(filtered_trap_0))
      }else{

        filtered_trap   <- graveyard %>% filter(timeDeath <= trap_days[te],
                                                timeDeath >= trap_days[te - 1],
                                                whereTrapped > -1)
        filtered_trap_0 <- filtered_trap[which(filtered_trap$whereTrapped %in% buf_0_traps$GID),]
        # filtered_trap_0 <- graveyard %>% filter(timeDeath <= trap_days[te],
        #                                         timeDeath >= trap_days[te - 1],
        #                                         typeDeath == 1,
        #                                         whereTrapped > -1)
      }

      N_trapped_total[te]   <- nrow(filtered_trap)
      NWb_trapped_total[te] <- length(which(filtered_trap$infStatus == 1))
      NWT_trapped_total[te] <- length(which(filtered_trap$infStatus == 0))
      prop_wb_total[te] <- NWb_trapped_total[te]/ N_trapped_total[te]
      prop_0_wb[te] <- length(which(filtered_trap_0$infStatus == 1))/nrow(filtered_trap_0)
      #print(prop_0_wb[te])
      #print(length(which(filtered_trap_0$infStatus == 1)))
      #print(nrow(filtered_trap_0))
    }

    graveyard_plot_data <- as.data.frame(cbind(trap_days,
                                               N_trapped_total,
                                               NWb_trapped_total,
                                               NWT_trapped_total,
                                               prop_wb_total,
                                               prop_0_wb))
    #print(graveyard_plot_data)

    #' end graveyard data
    #'
    # ------------------------------------------------
    #
    # NB the above wasn't working so
    #
    # ------------------------------------------------


  }

# ----------------- End of simulation ------------------------------------------------------
  #return_list <- list(inf.prop, grav.inf.prop, no.juv.CI, N.vec, juv.N.vec)
  return_list <- list(trap_days,
                      N_trapped_total,
                      NWb_trapped_total,
                      NWT_trapped_total,
                      prop_wb_total,
                      prop_0_wb)

  # save graveyard for testing
  graveyard$infStatus<- sapply(graveyard$infStatus, function(x) paste0(unlist(x), collapse = "\n"))
  save_graveyard <- as.data.frame(graveyard)
  #colnames(grav_inf_prop_matrix) <-  c("Day", (paste0("S",seq(1,timesRepeat, by = 1))))
 # write_csv(save_graveyard, paste0(lubridate::hour(Sys.time()), "h_graveyard.csv"))

  #write_csv(graveyard, "graveyard.csv")

  #return(return_list)

  #try returning the graveyard instead:
  return(save_graveyard)
}






#' God i don't know hopefully this will work.
#' This version does NOT run microclimates.
#' Input is the ABC parameters
#' Output are vectors of interest over time.
#' @param lambda Shape param for dispersal
#' @param k Parameter for mating
#' @param phi Parameter for trapping
#' @param alpha_a Mortality rate of adult mosquitoes
#' @param alpha_j Mortality rate of juvenile mosquitoes
#' @param eta Clutch size of each egg-laying event.
model_run_new <- function(lambda, k, phi, alpha_a, alpha_j, eta){

  boundaryDat <- read.table("inst/exdata/New_pp_Bound.txt", header=TRUE) #Boundary data for the boundary of the area used in trial (given as a set of lat/long points)
  boundary.bb <- read.table("inst/exdata/testboundary.txt", header=TRUE) # Bounding box of boundary, used for microclimates
  #minTemp     <- read.table("inst/exdata/cairnsAero_minTemp.txt", header=TRUE) # Min daily temperatures from Cairns Aero weather station in 2013 provided by BOM
  #maxTemp     <- read.table("inst/exdata/cairnsAero_maxTemp.txt", header=TRUE) # Max daily temperatures from Cairns Aero weather station in 2013 provided by BOM
  minTemp     <- read.table("inst/exdata/cairnsAero13min.txt", header=TRUE) # Min daily temperatures from Cairns Aero weather station in 2013 provided by BOM
  maxTemp     <- read.table("inst/exdata/cairnsAero13max.txt", header=TRUE) # Max daily temperatures from Cairns Aero weather station in 2013 provided by BOM
  const       <- read.table("inst/exdata/constants.txt", header=TRUE) #constants used for EKM calculations
  releaseDat  <- read.table("inst/exdata/PP_release.txt", header=TRUE) #Data regarding mosquito releases
  releaseDatExtended <- read.table("inst/exdata/releases_combined_final.txt", header = TRUE) # Including all releases within the 110 day period
  trappingDat <- read.table("inst/exdata/PP_trapping_data.txt", sep=",", header=TRUE) #Trapping data for mosquitoes
  # End reading in files

  # Set global constants ----
  testing         <- 1 # Set to 0 if we aren't doing any testing/optimisation
  max_age         <- 35 #mozzies who reach this age automatically die naturally
  max_juv_age     <- 14 #juveniles who reach this age automatically become adults
  max_pop_size    <- 200000 #if adult population exceeds this size, terminate simulation
  init_prop_inf   <- 0 #initial proportion of wild mozzies carrying Wolbachia- we expect this to be 0
  max_daily_mates <- 2 #maximum number of mating events a MALE mosquito can have in a day (from literature)
  include_microclimates <- 0 #Set to 1 to include microclimates effect, 0 to not include
  gridSize        <- 0.00025 # Size of grid spacing for microclimates
  #noLandTypes     <- 7 # Number of land types considered in microclimate analysis
  noLandTypes     <- 6 # Number of land types considered in microclimate analysis
  #typeTempDevs    <- c(-3, -7, 17, -4, 10, -3, 2) # Temperature deviates for each land type- FIX (should be read in)
  typeTempDevs    <- c(-3, -7, 17, -4, 2, -2)
  trapProb        <- 0.4 # Probability of being trapped when within trap radius
  # End global constants

  # For testing: alternative values of parameters ----
  #lambda <- rgamma(1,shape=0.044,scale=3) #CHECK/CHANGE THIS omg i don't know how the Gamma distribution works
  #lambda <- rgamma(1,3,22) #this is the alternative for lambda given in priors.pdf. Seems to give reasonable results with lat/long conversion so tempted to go with
  #pmale <- runif(1,0.3,0.7) #can also set to be 0.45
  #OR pmale <- rnorm(1,mean=0.5,sd=0.02) ??? check with Carla which one
  # k <-  rgamma(1,shape=0.2,scale=0.06) #CHANGE THIS #rate at which probability of mating drops wrt distance of potential mate
  #k <- rgamma(1,1,1/35) #rate at which probability of mating drops wrt distance of potential mate

  #WMP say that the following two values of eta are wrong/too high: so will put lower ones below
  #eta_1 <- round(rnorm(1,mean=120,sd=7),0) #no. of offspring produced by Wb non-carrier mothers (values from literature)
  #eta_2 <- round(rnorm(1,mean=70,sd=5),0) #no. of offspring produced by Wb carrier mothers (values from literature)

  #To re-iterate, these are MADE UP and need validifying
  #eta_1 <- round(rnorm(1,mean=40,sd=7),0) #no. of offspring produced by Wb non-carrier mothers
  #eta_2 <- round(rnorm(1,mean=30,sd=5),0) #no. of offspring produced by Wb carrier mothers
  #p_1 = runif(1,0,0.05)

  # Parameter generation for RABC ----
  param <- data.frame(lambda = rgamma(1,3,22),
                      pmale = 0.45, #proportion of male mosquitoes in wild population: FIXED
                      k = rgamma(1,shape=1,scale=1/35),
                      phi = runif(1,0.00001,0.001), #hyperparameter for probability that mozzies are trapped by any given trap
                      a = rtruncnorm(1, mean=0.02, sd=0.01, a=0, b=0.5), #Scale parameter regarding probability of natural death (Gompertz model)
                      b = rtruncnorm(1, mean=0.22, sd=0.05, a=0, b=0.5), #Scale parameter regarding probability of natural death (Gompertz model)
                      alpha_a = 0.1, #adult mortality rate, approx. from literature (0.2 is made up)
                      alpha_j = 0.3, #juvenile mortality rate: this is made up
                      p_1 = 1, #Probability of complete maternal transmission (virtually 100% according to Ross '16 & '17)
                      eta_1 = round(rnorm(1, mean=40, sd=7), 0), #no. of offspring produced by Wb non-carrier mothers
                      eta_2 = 0) #this was changed 8/1/20: Dutra15 and Ant18 imply eta_1 = eta_2
  #eta_2 = round(rnorm(1,mean=30,sd=5),0)) #no. of offspring produced by Wb carrier mothers
  param$eta_2 <- param$eta_1 # Dutra15 and Ant18 imply eta_1 = eta_2
  # End of parameter generation

  # Parameter generation for non- RABC ----
  param <- data.frame(lambda = rgamma(1,3,22),
                      pmale = 0.45,
                      k = rgamma(1, shape = 1, scale = 1/35),
                      a = rtruncnorm(1,mean=0.02,sd=0.01,a=0,b=0.5), #Scale parameter regarding probability of natural death (Gompertz model)
                      b = rtruncnorm(1,mean=0.22,sd=0.05,a=0,b=0.5),
                      phi = 0.00007,
                      alpha_a = 0.11, #from the book that peter ryan sent me
                      alpha_j = 0.3,
                      p_1 = 1,
                      eta_1 = 20,
                      eta_2 = 20)
  # End of parameter generation

  #' Option to select either high, low, or medium fixed values for each parameter
  #' comment this all out after the plots are done.
  #' low = -1, medium = 0, high = 1.
  highlowparams <- 0

  if(highlowparams == 0){
    #' Mid parameters:
    #param <- data.frame(lambda = 3/22, pmale = 0.45,
    #                    k = 1/35, a = 0.02, b = 0.22, phi = 0.00009,
    #                    alpha_a = 0.11, alpha_j = 0.3,
    #                    p_1 = 1, eta_1 = 20, eta_2 = 20)
    param <- data.frame(lambda = 3/22, pmale = 0.45,
                        k = 1/1135, a = 0.02, b = 0.22, phi = 0.00010,
                        alpha_a = 0.1, alpha_j = 0.3,
                        p_1 = 0.93, eta_1 = 20, eta_2 = 18)

    #' End mid parameters
  }else if(highlowparams == -1){
    #' Low parameters:
    param <- data.frame(lambda = (3/22 - 6/(22^2)), pmale = 0.45,
                        k = (1/35 - 2/(35^2)), a = 0.02, b = 0.22, phi = 0.00002,
                        alpha_a = 0.05, alpha_j = 0.1,
                        p_1 = 1, eta_1 = 10, eta_2 = 10)
    #' End low parameters
  }else if(highlowparams == 1){
    #' High parameters:
    param <- data.frame(lambda = (3/22 + 6/(22^2)), pmale = 0.45,
                        k = (1/35 + 2/(35^2)), a = 0.02, b = 0.22, phi = 0.0001,
                        alpha_a = 0.16, alpha_j = 0.5,
                        p_1 = 1, eta_1 = 30, eta_2 = 30)
    #' End high parameters
  }

  # Variable initialisation ----
  N           <- 15000 #number of adult mosquitoes at time t=0. To be initialised between 4000 to 20000
  Njuv        <- 6000 #number of juvenile CLUSTERS (not individuals) at time t=0
  idStart     <- N+1 #ID numbers of initial adults will be 1:N, then we start at N+1 for any new adult mosquito within simulation
  noTimeSteps <- 110 # Our simulation runs for 110 days
  t           <- 0 #timestep; starts from 0.
  # End variable initialisation

  # Calculate temperature and growth charts ----
  #Chart of average daily temperatures for each day of the simulation
  dailyTemps <- temperature_chart(minTemp, maxTemp, noTimeSteps)
  # If we're running microclimates, set them up ----
  if(include_microclimates == 1){
    print("Simulation includes microclimates")
    mc.setup <- initialise_enzyme_wmicroclim(boundary.bb, dailyTemps, typeTempDevs, noTimeSteps, noLandTypes, gridSize, const)
    # Extract EKM chart
    EKMChart <- list(mc.setup[[1]], mc.setup[[2]], mc.setup[[3]], mc.setup[[4]])
    # Retain other useful info
    landtype.df           <- data.frame(cbind(mc.setup$Latitude, mc.setup$Longitude, mc.setup$LandType, mc.setup$TemperatureDeviate))
    temp.grid.IDs         <- seq(1, nrow(landtype.df))
    landtype.df           <- data.frame(cbind(landtype.df, temp.grid.IDs))
    colnames(landtype.df) <- c("Latitude","Longitude","LandType", "TemperatureDeviate", "gridID")
    rm(temp.grid.IDs)
  } else if(include_microclimates == 0){
    print("No microclimates in simulation")
    EKMChart <- initialize_enzyme(dailyTemps,const) # Chart of daily updates for EKS
  }
  # End microclimate setup
  # End chart calculation


  #' This section sets up alternative temperature deviate setup methods ----
  #' There are three alternative TD setups for testing purposes:
  #' 1) TDs assigned uniformly at random
  #' 2) TDs assigned proportionally randomly to how frequently they appear
  #'    in the "realistic" microclimates setup
  #' 3) The grid is partitioned into the different land types, where the size
  #'    of each partition is proportional to how frequently they appear in the
  #'    "realistic" microclimates setup. (The order that the partitions appear
  #'    in is randomly generated).
  #' There are three function that do the above. The "tempdev_type" variable can
  #' be changed to 0 (keep realistic), 1 (option 1), 2 (option 2) and 3 (option 3).
  #' The functions assign land type, then we match it to the associated temperature
  #' deviate below.
  tempdev_type <- 0
  if(include_microclimates == 1 && tempdev_type == 1){
    #type 1
    new_landtype <- uniform_random_tempdev(landtype.df, noLandTypes, typeTempDevs)
  }else if(include_microclimates == 1 && tempdev_type == 2){
    #type 2
    new_landtype <- proportional_random_tempdev(landtype.df, noLandTypes, typeTempDevs)
  }else if(include_microclimates == 1 && tempdev_type == 3){
    #type 3
    new_landtype <- partition_tempdev(landtype.df, noLandTypes, typeTempDevs)
  }

  if(include_microclimates == 1 && tempdev_type > 0){
    landtypenumbers <- seq(1, noLandTypes, by = 1)
    tempLookup <- as.data.frame(cbind(landtypenumbers, typeTempDevs))
    colnames(tempLookup) <- c("LandType", "TemperatureDeviate")
    # do a MF left join
    #test <- landtype.df
    landtype.df$LandType <- new_landtype

    #I can't do this properly....
    landtype.df <- left_join(x = landtype.df, y = tempLookup, by = "LandType")
    colnames(landtype.df) <- c("Latitude", "Longitude", "LandType", "oldtempdev", "gridID", "TemperatureDeviate")
    landtype.df <- subset(landtype.df, select = -c(oldtempdev))
  }
  #' End alternative temperature deviate methods ----


  # Set up trap locations and clean trapping data ----
  #trapClean       <- trap_clean(trappingDat) # Clean trapping data
  #trapLoc         <- trap_setup(trapClean) # Dataframe of trap locations in lat/long
  #number_of_traps <- length(trapLoc$V1) # Number of traps in simulation
  #currentTraps    <- trap_empty(trapClean) # Active traps during the simulation period and day they're cleaned
  #emptyDays       <- unique(currentTraps$DayCollected) # Days when we calculated Wolbachia proportion in graveyard
  #daily_trapped   <- list() # Initialising list of trapped mosquitoes
  # End trapping data setup ----

  # 11/3/20: Set up traps a little differently: ----
  # we want to make sure they're active during the sim period
  trapClean       <- trap_clean(trappingDat) # Clean trapping data
  currentTraps    <- trap_empty(trapClean) # Active traps during the simulation period and day they're cleaned
  trapLoc         <- trap_setup(currentTraps) # Dataframe of trap locations in lat/long
  number_of_traps <- length(unique(currentTraps$GID))
  emptyDays       <- unique(currentTraps$DayCollected) # Days when we calculated Wolbachia proportion in graveyard
  daily_trapped   <- list() # Initialising list of trapped mosquitoes
  # End trapping data setup ----


  # Wrangling of release data used in simulation ----
  #releaseDat   <- release_clean(releaseDat) #Cleaning release data
  #release.days <- as.data.frame(table(day = releaseDat$Day)) #Grabbing days that releases were made

  releaseDat   <- release_clean(releaseDatExtended) #Cleaning release data
  release.days <- as.data.frame(table(day = releaseDat$Day)) #Grabbing days that releases were made

  # End of release data wrangling

  # Set up grid ----
  # We do this even if no microclimates used since it's useful for plotting
  grid.df <- make_grid(boundary.bb, gridSize)
  # End grid setup

  # ---- SET UP DOING MULTIPLE SIMULATIONS ----
  timesRepeat <- 1 #number of times to repeat simulation
  inf.prop.matrix <- seq(from = 1, to = noTimeSteps, by = 1)
  grav.inf.prop.matrix <- seq(from = 1, to = noTimeSteps, by = 1)
  no.juv.CI.matrix <- seq(from = 1, to = noTimeSteps, by = 1)
  N.matrix <- seq(from = 1, to = noTimeSteps, by = 1)
  graveyard.plot.list <- list()

  #' ---- Set WD for saving graveyard ----
  setwd("~/Documents/masters tomfoolery/modeloutput/nomicroclim/")


  for(timesrep in 1:timesRepeat){

    # Model initialisation/construction ----
    mozzie.dt     <- as.data.table(initialise_adults(N, param$pmale, boundaryDat, grid.df)) #initial adult population
    juv.dt        <- as.data.table(initialise_juveniles(Njuv, grid.df, boundaryDat)) #initial juvenile population
    graveyard     <- initialise_graveyard() # Where entries of dead mosquitoes go.
    juv.graveyard <- initialise_juv_graveyard() # Where entries of eggs with CI go
    # End initial model construction

    # Initialising data output ----
    # Vector to keep track of Wolbachia infection proportion over time
    inf.prop      <- rep(0, noTimeSteps) # Proportion of Wb infection over time: alive adults
    grav.inf.prop <- rep(0, noTimeSteps) # Proportion of Wb infection from trapped mozzies.
    no.juv.CI     <- rep(0, noTimeSteps) # Number of juvenile clutches with CI over time
    N.vec         <- rep(0, noTimeSteps) # Number of adult agents over time
    juv.N.vec     <- rep(0, noTimeSteps)
    #graveyard.plot.list <- list()
    # End data output initialisation ----

    # Below is for debugging/identifying where things mucked up ----
    #s1 <- 0
    #s2 <- 0
    #s3 <- 0
    #s4 <- 0
    #s5 <- 0
    #s6 <- 0
    #s7 <- 0
    #s8 <- 0
    #s9 <- 0
    # end debugging variables ----

    # Framework for optimisation: how much time does each step take? ----
    #' Each row will be a timestep, each column will be the time each step takes
    #' PLUS the number of agents.
    optimisation <-  data.frame(matrix(ncol = 10, nrow = noTimeSteps))
    max_gono <- rep(0, noTimeSteps)
    # End optimisation framework

    t  <- 1 # Since model initialisation is complete
    # Progress bar
    pb <- progress_bar$new(format = " elapsed [:bar] :percent elapsed: :elapsed eta: :eta ", total = noTimeSteps)

    for(t in 1:noTimeSteps){
      # Step 1: Juvenile emergence -----
      start <- Sys.time()

      to.adult <- which(juv.dt$stage >= 3 & juv.dt$enzyme > 0.95 & juv.dt$infProb == 0)

      to.adult.inf <- which(juv.dt$stage >= 3 & juv.dt$enzyme > 0.95 & juv.dt$infProb == param$p_1)

      #to.nat.die <- mozzie.dt[sample(.N, round(param$alpha_a*nrow(mozzie.dt)))]

      #' ------------------------
      #'
      #' ADDED 12/11/20: EGG HATCH RATE OF 0.8 AS PER PERRAN ROSS
      #'
      #' This won't error if length(to.adult) == 0 hopefully
      to.adult <-  sample(to.adult, round(0.8*length(to.adult)))
      to.adult.inf <- sample(to.adult.inf, round(0.7*length(to.adult.inf)))
      to.adult <- c(to.adult, to.adult.inf)
      #' ------------------------
      #'

      if((length(to.adult)) != 0){
        IDvec <- juv.dt$clutchSize[to.adult] # Vector of agent ID indices
        IDvec <- c(0, IDvec)
        IDvec <- IDvec[1:length(to.adult)]
        IDvec <- cumsum(IDvec) + idStart

        #IDvec <- c(0, IDvec)
        #IDvec <- IDvec[1:length(to.adult)]
        # new.adults.dt <- mapply(FUN=juv_to_adult, to.adult, idStart, param$pmale, param$lambda, SIMPLIFY = FALSE)
        #new.adults.dt <- mapply(FUN=juv_to_adult, to.adult, IDvec, param$pmale,
        #                        param$lambda, juv.dt, grid.df, SIMPLIFY = FALSE)

        new.adults.dt <- mapply(FUN = juv_to_adult,
                                to.adult, IDvec, param$pmale, param$lambda,
                                MoreArgs = list(juv.dt, grid.df),
                                SIMPLIFY = FALSE)

        new.adults.dt %<>% rbindlist() %>% as.data.table()
        mozzie.dt     <<- rbind(mozzie.dt, new.adults.dt) #CHECK do we need <<- operator or can we just use <- ?
        # Delete row(s) of data table corresponding to the clutches that aged into adult mozzies
        juv.dt        <<- juv.dt[-to.adult]
        idStart       <- idStart + (nrow(new.adults.dt) + 1) # So we ID mosquitoes correctly
        rm(new.adults.dt)
      }

      # Step 1.5: any juveniles moving to the next development stage? ----
      # If enzyme > 0.95, move on to the next development stage and reset enzyme score
      if(length(which(juv.dt$enzyme > 0.95)) > 0){
        # Increment development stage
        juv.dt[which(juv.dt$enzyme > 0.95)]$stage  <- juv.dt[which(juv.dt$enzyme > 0.95)]$stage + 1
        # Reset EkS to 0
        juv.dt[which(juv.dt$enzyme > 0.95)]$enzyme <- 0
      }
      # End Step 1.5

      rm(to.adult)
      rm(to.adult.inf)

      # Exception handling: if adult population size now exceeds max_pop_size, kill the simulation
      if(nrow(mozzie.dt) > max_pop_size){
        print("Termination condition 1: Maximum population size exceeded")
        break
      }
      #  s1 <- s1 +1 #debugging
      end <- Sys.time()
      optimisation[t,1] <- end - start
      rm(start, end)

      # End Step 1



      # Step 2: Egg laying ----
      start <- Sys.time()
      #Determine which mothers will lay their eggs today
      #Conditions to be satisfied for a mother to lay a clutch of eggs:
      ## 1) Has a mate
      ## 2) Enzyme Kinetic Score > 1 (should mean that mosquito has a mate)
      ## 3) First egg laying event is at EKS > 1, second at 1.58,
      ##    third at 2.16, etc. for multiples of 0.58 (from Focks et al '93)
      if(length(which(mozzie.dt$gonoCycle > 3))){
        print("gonoCycle too big")
        break
      }

      # below is the "real" to.lay line
      #to.lay <- mozzie.dt$ID[which(mozzie.dt$gender == 1 & mozzie.dt$mateID != -1 & ((mozzie.dt$enzyme >= 1 & mozzie.dt$gonoCycle == 0) || (mozzie.dt$enzyme >= 1.58 & mozzie.dt$gonoCycle == 1) || (mozzie.dt$enzyme >= 2.16 & mozzie.dt$gonoCycle == 2)))]

      # to.lay <- mozzie.dt$ID[(which(mozzie.dt$enzyme >= 1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$gender == 1 & mozzie.dt$mateID != -1)) |
      #                       (which(mozzie.dt$enzyme >= 1.58 & mozzie.dt$gonoCycle == 1 & mozzie.dt$gender == 1 & mozzie.dt$mateID != -1)) |
      #                      (which(mozzie.dt$enzyme >= 2.16 & mozzie.dt$gonoCycle == 2 & mozzie.dt$gender == 1 & mozzie.dt$mateID != -1)) ]

      firstlay  <- mozzie.dt$ID[(which(mozzie.dt$enzyme >= 1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$gender == 1 & mozzie.dt$mateID != -1))]
      secondlay <- mozzie.dt$ID[(which(mozzie.dt$enzyme >= 1.58 & mozzie.dt$gonoCycle == 1 & mozzie.dt$gender == 1 & mozzie.dt$mateID != -1))]
      thirdlay  <- mozzie.dt$ID[(which(mozzie.dt$enzyme >= 2.16 & mozzie.dt$gonoCycle == 2 && mozzie.dt$gender == 1 && mozzie.dt$mateID != -1))]
      to.lay <- c(firstlay, secondlay, thirdlay)

      #' 22-09-20: reduce proportion of females laying eggs per day: CHECK
      to.lay <- sample(to.lay, size = length(to.lay)/5, replace = FALSE)

      # to.lay <- which(mozzie.dt$enzyme[ready.to.lay] >= 1 & mozzie.dt$gonoCycle[ready.to.lay]  == 0)
      #  to.lay <- which(mozzie.dt$gender == 1 & mozzie.dt$mateID != -1 & ((mozzie.dt$enzyme >= 1 & mozzie.dt$gonoCycle == 0) | (mozzie.dt$enzyme >= 1.58 & mozzie.dt$gonoCycle == 1) | (mozzie.dt$enzyme >= 2.16 & mozzie.dt$gonoCycle == 2)))

      if(length(to.lay) != 0){
        # Create new egg agents
        #system.time(new.eggs.dt <- initialise_eggs(to.lay, param$eta_1, param$eta_2, param$p_1, param$alpha_j))
        system.time(new.eggs.dt <- initialise_eggs(to.lay, param$eta_1, param$eta_2, param$p_1, param$alpha_j, mozzie.dt, graveyard))
        # If there are egg clutches with CI, remove them and add to juv.graveyard
        to.CI       <- which(new.eggs.dt$infProb == -1)
        if(length(to.CI) > 0){
          # Adding eggs with CI to the juv.graveyard
          l             <- list(juv.graveyard, new.eggs.dt[to.CI])
          juv.graveyard <- rbindlist(l, use.names = TRUE)
          # Remove CI eggs from new.eggs.dt
          new.eggs.dt   <- new.eggs.dt[!(to.CI)]
          rm(l)
        }
        # Add new eggs to the juvenile data.table
        juv.dt      <- rbind(juv.dt,new.eggs.dt)
        rm(new.eggs.dt)
        rm(to.CI)

        # You need to increment gonoCycle by 1
        #I think the below works....
        mozzie.dt$gonoCycle[mozzie.dt$ID %in% to.lay] <- mozzie.dt$gonoCycle[mozzie.dt$ID %in% to.lay] + 1
        # mozzie.dt$gonoCycle[which(mozzie.dt$ID[to.lay] %in% mozzie.dt$ID)] <- mozzie.dt$gonoCycle[which(mozzie.dt$ID[to.lay] %in% mozzie.dt$ID)] + 1
      }
      rm(to.lay)
      rm(firstlay)
      rm(secondlay)
      rm(thirdlay)#Remove to.lay
      #  s2 <- s2 +1 #debugging

      end <- Sys.time()
      optimisation[t,2] <- end - start
      rm(start, end)
      # End Step 2

      # Step 3: Do we have mosquitoes being trapped (and killed?) ----
      start <- Sys.time()
      for(tr in 1:number_of_traps){
        to.trap  <- list()
        #temptrap <- find_trapped(trapLoc$V1[tr], trapLoc$V2[tr], param$phi, mozzie.dt)
        temptrap <- find_trapped(trapLoc$Lat[tr], trapLoc$Long[tr], param$phi, mozzie.dt, trapProb)

        if(length(temptrap)> 1 && temptrap[1] != -1){
          # This is slow, work out how to use rbindlist for this in the future.
          # to.trap <- rbind(to.trap, temptrap)
          if(tr == 1){
            to.trap <- temptrap
            mozzie.dt$whereTrapped[mozzie.dt$ID %in% to.trap] <- trapLoc$GID[tr]
            #mozzie.dt[mozzie.dt$ID %in% to.trap,]$whereTrapped <- trapLoc$V3[tr]
          }else{
            to.trap <- append(to.trap, temptrap)
            mozzie.dt$whereTrapped[mozzie.dt$ID %in% to.trap] <- trapLoc$GID[tr]
            #mozzie.dt[mozzie.dt$ID %in% to.trap,]$whereTrapped <- trapLoc$V3[tr]
          }
          # Updating the entries of mozzies that have been trapped
          break;
          to.trap <- unlist(to.trap)
          mozzie.dt$typeDeath[mozzie.dt$ID %in% to.trap] <- 1 # Type of death is trapped death
          mozzie.dt$timeDeath[mozzie.dt$ID %in% to.trap] <- t # They died today
          #mozzie.dt[mozzie.dt$ID %in% to.trap,]$typeDeath <- 1 # Type of death is trapped death
          #mozzie.dt[mozzie.dt$ID %in% to.trap,]$timeDeath <- t # They died today
          #to.trap   <- unique(to.trap) #in case there are duplicates #CHECK this might be our mystery duplicate IDs
          l         <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.trap,])
          graveyard <- rbindlist(l) # Trapped mozzies go in the graveyard!
          mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.trap),]# Delete entries of mozzies who have been trapped and killed
          rm(l)
        }
        rm(temptrap)
      }
      # The following line of code is for testing: this removes all functional duplicates:
      ## test <- graveyard[!duplicated(graveyard[,!c('whereTrapped')]),]

      # For fun: plotting traps and agent locations
      #bdary <- read.table("inst/exdata/testboundary.txt", sep=" ", header=TRUE) #read in data
      #trapp <- ggplot(trapLoc, aes(x=V1, y=V2)) + geom_point()
      #trapp + geom_point(data=mozzie.dt, aes(x=mozzie.dt$lat, y=mozzie.dt$long, color="red")) +
      #geom_path(data = boundary.bb, aes(x=boundary.bb$lat, y=boundary.bb$long, color="green")) +
      #geom_path(data = boundaryDat, aes(x = boundaryDat$Lat, y=boundaryDat$Long, color="blue"))

      #  s3 <- s3 + 1 #debugging

      # Exception handling: if there are no more adults because they've all been trapped,
      ## kill the simulation
      if(nrow(mozzie.dt) == 0){
        print("Termination condition 3: All adults have died (traps)")
        break
      }
      end <- Sys.time()
      optimisation[t,3] <- end - start
      rm(start, end)
      # End Step 3

      # Step 4: Natural mosquito death (juvenile and adult) ----
      start <- Sys.time()
      ## Juvenile natural death: remove (alpha_j)% of the population
      ## The subsetting of juvdf$clutchSize is because we don't reduce clutch size of CI-affected eggs:
      ## they're already dead, but we just want to track how many there are so we leave clutchSize alone
      ## ^ The above statement is deprecated because CI-affected eggs now go to juv.graveyard,
      ## but it's still a useful sanity check to have in the code. So it stays


      juv.dt$clutchSize[which(juv.dt$pDeath != -1)] <- mapply(FUN = resize_clutch, juv.dt$clutchSize, param$alpha_j)
      #

      #' ADDED 14/11/20: MAKING JUV DEATH RATE HIGHER FOR THOSE WITH WB
      #juv.dt$clutchSize[which(juv.dt$pDeath != -1 & juv.dt$infProb == param$p_1)] <- mapply(FUN = resize_clutch, juv.dt$clutchSize[which(juv.dt$pDeath != -1 & juv.dt$infProb == param$p_1)], 1.01*param$alpha_j)

      # juv.dt$clutchSize[which(juv.dt$pDeath != -1)] <- mapply(FUN = resize_clutch,
      #                                                          juv.dt$clutchSize[which(juv.dt$pDeath != -1)],
      #                                                          juv.dt$pDeath[which(juv.dt$pDeath != -1)])

      # Remove juvenile agents where clutchSize == 0.
      # We don't need the information but let's put them in the graveyard anyway.
      # EDIT 12/3/20: don't add empty clutches to the graveyard and see if it improves computation.
      # Uncomment the first two lines under the if statement if we want to bring this back in.
      dead.clutch <- which(juv.dt$clutchSize == 0)
      if(length(dead.clutch) > 0){
        #  l      <- list(juv.graveyard, juv.dt[dead.clutch,]) #List of dead mozzies to remove
        #  juv.graveyard <- rbindlist(l, use.names = TRUE)
        juv.dt <- juv.dt[!(dead.clutch)]
        # rm(l)
      }
      rm(dead.clutch)

      # Adult (old age) death:
      ## Adults die automatically when they reach max_age
      ## Vector of mosquitoes that will die of old age:
      to.old.die <- which(mozzie.dt$age == max_age & mozzie.dt$timeDeath == -1)
      if((length(to.old.die)) != 0){
        #Track day that they died (timeDeath) and set typeDeath = 2 to denote natural (not trapped) death
        mozzie.dt$typeDeath[to.old.die] <- 2
        mozzie.dt$timeDeath[to.old.die] <- t

        #' add to graveyard
        old.to.graveyard <- filter(mozzie.dt, typeDeath == 2)
        l <- list(graveyard, old.to.graveyard)
        graveyard <- rbindlist(l, use.names = TRUE)
        rm(old.to.graveyard, l)

        #'remove from mozzie dataframe
        mozzie.dt <- mozzie.dt[!to.old.die,]
      }
      rm(to.old.die)

      # Adults can also die as per natural death rate
      # CHECK: do we want the same number killed off every timestep?
      #mozzie.dt <- mozzie.dt[!sample(.N, round(param$alpha_a*nrow(mozzie.dt)))]
      #to.nat.die <- mozzie.dt[sample(.N, round(param$alpha_a*nrow(mozzie.dt)))]

      #' TESTING 05/11/20: RELEASED ADULTS WITH HIGHER DEATH RATE HHHHHHH --------
      to.nat.die.release <- mozzie.dt %>% filter(releaseLoc != -1) %>%
        sample_n(round(param$alpha_a*3*nrow(.)), replace = FALSE)
      to.nat.die <- mozzie.dt %>% filter(releaseLoc == -1) %>%
        sample_n(round(param$alpha_a*nrow(.)), replace = FALSE)

      #' END TESTING 05/11/20 HHHH ----------------------
      if(length(to.nat.die) >= 1){
        mozzie.dt[mozzie.dt$ID %in% to.nat.die$ID,]$typeDeath <- 0
        mozzie.dt[mozzie.dt$ID %in% to.nat.die$ID,]$timeDeath <- t
        l <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.nat.die$ID,]) # List of dead mozzies to remove
        graveyard <- rbindlist(l, use.names = TRUE)
        mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.nat.die$ID),] # Removing naturally dead mozzies from mozzie.dt
        rm(to.nat.die)
        rm(l)
      }

      if(length(to.nat.die.release) >= 1){
        mozzie.dt[mozzie.dt$ID %in% to.nat.die.release$ID,]$typeDeath <- 0
        mozzie.dt[mozzie.dt$ID %in% to.nat.die.release$ID,]$timeDeath <- t
        l <- list(graveyard, mozzie.dt[mozzie.dt$ID %in% to.nat.die.release$ID,]) # List of dead mozzies to remove
        graveyard <- rbindlist(l, use.names = TRUE)
        mozzie.dt <- mozzie.dt[!(mozzie.dt$ID %in% to.nat.die.release$ID),] # Removing naturally dead mozzies from mozzie.dt
        rm(to.nat.die.release)
        rm(l)
      }

      # Exception handling: if all adults have died, kill the simulation
      if(nrow(mozzie.dt) == 0){
        print("Termination condition 2: All adults have died (natural)")
        break
      }
      # Exception handling: if all juveniles have died, kill the simulation
      if(nrow(juv.dt) == 0){
        print("Termination condition 4: All juveniles have died (natural)")
        break
      }
      #  s4 <- s4 + 1 #debugging
      end <- Sys.time()
      optimisation[t,4] <- end - start
      rm(start, end)
      # End Step 4 ----

      # Step 5: Mosquito mating ----
      #SMB 8/9/20: commented this out because I have two mating steps? uncomment up to 506 if we add back
      #   start <- Sys.time()
      #   # Step 5 try 2: find_mate isn't working anymore since passing mozzie.dt doesn't actually pass the data.table
      #   # Which makes everything really frustrating. So I'll just do a loop for now :/
      #   # This is more or less all copied and pasted from the find_mate function.
      #   # This is a really obvious area to improve if we need to speed up simulation.
      #   # We only need to find mates for:
      #   ## Females that don't have mates and haven't laid eggs yet (first gono cycle)
      #   ## AND has EKS >= 1 (From Focks)
      #   to.mate <- which((mozzie.dt$gender == 1 & mozzie.dt$mateID == -1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$enzyme >= 1))
      #   if(length(to.mate) != 0){
      #     for(m in 1:length(to.mate)){
      #       # Grab the position of the female that's mating
      #       to.mate.pos    <- as.numeric(c("lat"  = mozzie.dt$lat[to.mate[m]],
      #                                      "long" = mozzie.dt$long[to.mate[m]]))
      #       # Draw a square around her with width 2*k: this is the search area
      #       to.mate.bdary  <- c("latmin"  = to.mate.pos[1]-param$k, "latmax"  = to.mate.pos[1]+param$k,
      #                           "longmin" = to.mate.pos[2]-param$k, "longmax" = to.mate.pos[2]+param$k)
      #       # Find possible mates within this search area
      #       possible.mates <- which(mozzie.dt$lat >= to.mate.bdary["latmin"] & mozzie.dt$lat <= to.mate.bdary["latmax"] & mozzie.dt$long >= to.mate.bdary["longmin"] & mozzie.dt$long <= to.mate.bdary["longmax"] & mozzie.dt$gender == 0)
      #       no.bachelors   <- length(possible.mates)
      #       if(no.bachelors == 0){
      #         # There were no suitable mates within the search area
      #         new.mate <- -1
      #         #print(paste0(to.mate[m], " no mate found"))
      #       }else if(no.bachelors == 1){
      #         if(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates){
      #           # Even though there is only one male, he has mated too much today
      #           new.mate <- -1
      #          # print(paste0(femID, " 1 mate found, unsuitable"))
      #         }else{
      #           #New mate is just the only male they found in the search area
      #           new.mate <- as.integer(possible.mates)
      #           mozzie.dt$gonoCycle[new.mate] <- mozzie.dt$gonoCycle[new.mate] + 1
      #           #print(paste0(femID, " 1 mate found"))
      #         }
      #       }else{
      #         if(length(which(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates)) == 0){
      #           # More than one suitable mate was found in the search area.
      #           # Randomly permutes the list of possible mates and then picks the one at the top of the pile
      #           new.mate <- sample(possible.mates)[1]
      #           # Increment number of mates of male for this timestep by 1
      #           mozzie.dt$gonoCycle[new.mate] <- mozzie.dt$gonoCycle[new.mate] + 1
      #           #print(paste0(femID, " multiple mates, all ok"))
      #         }else{
      #           # Remove males who have mated more than "max_daily_mates" number of times in a day
      #           possible.mates <- possible.mates[-(which(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates))]
      #           # Now we have one more condition to check for:
      #           ## If we drop some males from "possible.mates" we might end up dropping them all,
      #           ## in which case, we return -1
      #           if(length(possible.mates) == 0){
      #             #print(paste0(femID, " multiple mates, had to drop some AND then ended up with none"))
      #             new.mate <- 1
      #           }else{
      #             # Randomly permutes the list of possible mates and then picks the one at the top of the pile
      #             new.mate <- sample(possible.mates)[1]
      #             mozzie.dt$gonoCycle[new.mate] <- mozzie.dt$gonoCycle[new.mate] + 1 #increment number of mates by 1
      #            # print(paste0(femID, " multiple mates, had to drop some"))
      #           }
      #         }
      #       }
      #       mozzie.dt$mateID[to.mate[m]] <- new.mate
      #       rm(to.mate.pos)
      #       rm(to.mate.bdary)
      #       rm(possible.mates)
      #       rm(no.bachelors)
      #       new.mate <- -1
      #     }
      #   }
      #   rm(to.mate)
      #   # For male agents, gonoCycle has a different purpose:
      #   # it tracks the number of times per day that the male has mated.
      #   # So we now reset gonoCycle for the next day for all male agents:
      #   mozzie.dt$gonoCycle[which(mozzie.dt$gender == 0)] <- 0
      # #  s5 <- s5 +1 #debugging
      #   end <- Sys.time()
      #   optimisation[t,5] <- end - start
      #   rm(start, end)
      # End Step 5

      # Step 6: Updating Enzyme Kinetic Score for each agent ----
      #   start <- Sys.time()
      #   # Step 6 Option 1: No microclimate effects
      #   if(include_microclimates == 0){
      #     ##Update EKS for juvenile agents
      #     juv.enzyme.update <- mapply(FUN = update_enzyme, juv.dt$stage, t)
      #     #This is a super hacky method....
      #     juv.dt$enzyme <- as.numeric(juv.dt$enzyme) + as.numeric(juv.enzyme.update)
      #
      #     ##Update EKS for adult agents
      #     ##Adult mosquitoes always update with the same EKM formula, hence are always in "stage 4"
      #     mozzie.enzyme.update <- mapply(FUN = update_enzyme, 4, t)
      #     mozzie.dt$enzyme <- as.numeric(mozzie.dt$enzyme) + as.numeric(mozzie.enzyme.update)
      #
      #     rm(juv.enzyme.update)
      #     rm(mozzie.enzyme.update)
      #   }else if(include_microclimates == 1){
      #     # Step 6 Option 2: Microclimate effects included
      #     # Remember: access is EKMChart[[stage]][land type, ][timestep]
      #     # Also remember: the above type is a LIST, we need to convert to a double
      #     juv.enzyme.update   <- unlist(lapply(juv.dt$gridID, grid_to_land_type))
      #     mozzie.enzyme.update <- unlist(lapply(mozzie.dt$gridID, grid_to_land_type))
      #
      #     # NB system.time() of the following for 10000 mozzies is 4.491s! Not sure how to make this quicker.
      #     mozzie.enzyme.update <- mapply(FUN = update_enzyme_microclim, stage = 4, timestep = t, landType = mozzie.enzyme.update)
      #     juv.enzyme.update    <- mapply(FUN = update_enzyme_microclim, stage = juv.dt$stage, timestep = t, landType = juv.enzyme.update)
      #     # test <- lapply(adult.enzyme.update, update_enzyme_microclim, stage = 4, timestep = t) #Slightly faster but needs rejigging from a list
      #
      #     # Updating the EKS of each agent
      #     mozzie.dt$enzyme <- as.numeric(mozzie.dt$enzyme) + as.numeric(mozzie.enzyme.update)
      #     juv.dt$enzyme    <- as.numeric(juv.dt$enzyme) + as.numeric(juv.enzyme.update)
      #
      #     rm(juv.enzyme.update)
      #     rm(mozzie.enzyme.update)
      #   }
      # #  s6 <- s6 + 1 #debugging
      #   end <- Sys.time()
      #   optimisation[t,6] <- end - start
      #   rm(start, end)
      # # End Step 6

      #' Begin Step 6 Take 2: a faster way to update microclimates ----
      start <- Sys.time()
      if(include_microclimates == 0){
        #' For adults/mozzie.dt: everyone has the same EKS update!
        EKS_today <- calculate_daily_EKS((dailyTemps[t] + 273.15), 4, const)
        mozzie.dt$enzyme <- mozzie.dt$enzyme + as.numeric(EKS_today)

        #'For juveniles, it's a little bit more involved:
        for(s in 1:3){
          #juv.dt$stage[which(juv.dt$stage == s)] <- mozzie.dt$enzyme[which(juv.dt$stage == s)] +
          #                                          calculate_daily_EKS((dailyTemps[t] + 273.15), s)
          juv_EKS_today <- as.numeric(calculate_daily_EKS((dailyTemps[t] + 273.15), s, const))
          #juv.dt$enzyme[which(juv.dt$stage == s)] <- juv.dt$enzyme[which(juv.dt$stage == s)] +
          #                                           as.numeric(juv_EKS_today)


          juv.dt$enzyme[which(juv.dt$stage == s)] <- sapply(juv.dt$enzyme[which(juv.dt$stage == s)],
                                                            function(x){x <- x + as.numeric(juv_EKS_today)},
                                                            simplify = TRUE)
          rm(juv_EKS_today)
        }

      }else if(include_microclimates == 1){
        #' Get the temperature deviate for each mozzie
        mozzie.w.tempdev <- left_join(mozzie.dt, landtype.df, by = "gridID")
        juv.w.tempdev    <- left_join(juv.dt, landtype.df, by = "gridID")

        #'Apply the temperature deviate to today's temperature and then convert to Kelvin
        #'This is the temperature that will then be fed into the EKM
        mozzie.w.tempdev$TemperatureDeviate <- mozzie.w.tempdev$TemperatureDeviate +
          dailyTemps[t] + 273.15
        juv.w.tempdev$TemperatureDeviate    <- juv.w.tempdev$TemperatureDeviate +
          dailyTemps[t] + 273.15
        #'Find EKS score
        #'first: for adults, stage = 4
        #' Find the EKS for each adult for today
        EKS_today <- calculate_daily_EKS(mozzie.w.tempdev$TemperatureDeviate, 4, const)
        #' Add to the current enzyme
        mozzie.dt$enzyme <- mozzie.dt$enzyme + as.numeric(EKS_today)

        #' now: for juveniles
        for(s in 1:3){
          #juv.dt$stage[which(juv.dt$stage == s)] <- mozzie.dt$enzyme[which(juv.dt$stage == s)] +
          #                                          calculate_daily_EKS((dailyTemps[t] + 273.15), s)
          if(length(which(juv.w.tempdev$stage == s)) > 0){
            juv_EKS_today <- calculate_daily_EKS(juv.w.tempdev$TemperatureDeviate[which(juv.w.tempdev$stage == s)], s, const)
            #' DOBULE CHECK the below works....
            #' commented the below out 7/1/20
            #juv.dt$enzyme[which(juv.dt$stage == s)] <- sapply(juv.dt$enzyme[which(juv.dt$stage == s)],
            #                                                  function(x){x <- x + as.numeric(juv_EKS_today)},
            #                                                  simplify = TRUE)

            juv.dt$enzyme[which(juv.dt$stage == s)] <- juv.dt$enzyme[which(juv.dt$stage == s)] +
              as.numeric(juv_EKS_today)
            #juv.dt$enzyme[which(juv.dt$stage == s)] <- juv.dt$enzyme[which(juv.dt$stage == s)] +
            #                                           as.numeric(juv_EKS_today)
            rm(juv_EKS_today)
          }
        }
      }
      end <- Sys.time()
      optimisation[t,6] <- end - start
      rm(start, end)
      rm(EKS_today)
      #'
      #'
      #'
      #' End Step 6 Take 2 :)


      # Step 7: Mosquito dispersal ----
      # A decision was made to have 1 dispersal event: upon juvenile -> adult emergence.
      # This is built into the juv_to_adult function called in Step 1.
      # The below code is legacy code, and can be uncommented if we wish to add another dispersal event
      # You just need to have an appropriate condition to make a to.migrate list.

      # ---- the below was commented out on 27/2/20 ----
      #  if(length(to.migrate) != 0){
      #CHECK this works with length(toMigrate) > 1
      #    update.disp <- mapply(FUN = random_dispersal, mozzie.dt$lat[to.migrate], mozzie.dt$long[to.migrate], param$lambda)
      #    update.disp <- do.call(rbind, update.disp)

      # Suppressing warnings for this assignment at the moment:
      ## This operation coerces a list to a double, going from a precision of 47 dp to 45 dp.
      ## However for lat/long any value to more than 10 dp is fairly nonsensical (sub milimetre) so we don't care
      #   suppressWarnings(mozzie.dt[to.migrate, "lat" := update.disp[1,]])
      #    suppressWarnings(mozzie.dt[to.migrate, "long" := update.disp[2,]])
      #   rm(update.disp)
      #  }
      #rm(to.migrate)
      # ---- uncomment the above if we want to add a second migration event
      #  s7 <- s7 +1 #debugging
      # End Step 7


      # Step 8: Update numerical 'age' of agent ----
      start <- Sys.time()
      juv.dt$age    <- mapply('+', juv.dt$age, 1)
      mozzie.dt$age <- mapply('+', mozzie.dt$age, 1)
      #  s8 <- s8 +1 #debugging
      end <- Sys.time()
      optimisation[t,8] <- end - start
      rm(start, end)
      # End Step 8

      # Step 9: Mosquito releases ----
      start <- Sys.time()
      today.releases <- list()
      #Check if it's one of the days where there is a mosquito release
      if(t %in% release.days$day){
        # For each day where a release was held, there are multiple releases in multiple locations
        # So we need to iterate through all of them.
        # Ideally, we should be able to replace this for loop with some vectorised code.
        temp.releases <- releaseDat[which(releaseDat$Day == t),] # Rows of releaseDat corresponding to today's releases
        rownames(temp.releases) <- 1:nrow(temp.releases)
        # We will run initialise_releases for each line of today.releases
        for(i in 1:nrow(temp.releases)){
          if(as.integer(temp.releases$NoMozzie[i] > 0)){
            #today.releases <- rbind(today.releases, initialise_release(as.integer(temp.releases$NoMozzie[i]), temp.releases$NoMale[i], temp.releases$NoFemale[i], temp.releases$Lat[i], temp.releases$Long[i], idStart, temp.releases$GID[i]), param$p_1)
            today.releases <- rbind(today.releases, initialise_release(as.integer(temp.releases$NoMozzie[i]),
                                                                       temp.releases$NoMale[i],
                                                                       temp.releases$NoFemale[i],
                                                                       temp.releases$Lat[i],
                                                                       temp.releases$Long[i],
                                                                       idStart,
                                                                       temp.releases$GID[i], param$p_1, grid.df))

            idStart        <- idStart + (as.integer(temp.releases$NoMozzie[i]) + 1) #So we keep track of the right mosquito ID number
          }
        }
        # Mozzies all need to have a dispersal event at time of release
        # TO DO: update gridID as well!!!! * CHECK
        # Idea: perhaps try to do this in constructor? Not sure which one is faster
        today.releases.disp <- mapply(FUN = random_dispersal, today.releases$lat, today.releases$long, param$lambda) #DOES THIS WORK? LOL
        today.releases.disp <- do.call(rbind, today.releases.disp)
        today.releases.lat  <- today.releases.disp[c(TRUE, FALSE)]
        today.releases.long <- today.releases.disp[c(FALSE, TRUE)]
        today.releases$lat  <- today.releases.lat
        today.releases$long <- today.releases.long

        l         <- list(mozzie.dt, today.releases)
        mozzie.dt <- rbindlist(l, use.names = TRUE)

        rm(l)
        rm(today.releases)
        rm(today.releases.disp)
        rm(today.releases.lat)
        rm(today.releases.long)
      }

      # s9 <- s9 + 1 #debugging
      end <- Sys.time()
      optimisation[t,9] <- end - start
      rm(start, end)
      # End Step 9 ----
      start <- Sys.time()

      # Step 5: Mosquito mating ----
      # Step 5 try 2: find_mate isn't working anymore since passing mozzie.dt doesn't actually pass the data.table
      # Which makes everything really frustrating. So I'll just do a loop for now :/
      # This is more or less all copied and pasted from the find_mate function.
      # This is a really obvious area to improve if we need to speed up simulation.
      # We only need to find mates for:
      ## Females that don't have mates and haven't laid eggs yet (first gono cycle)
      ## AND has EKS >= 1 (From Focks)

      #' dataframe for optimisation:
      mate.opto <- data.frame(matrix(ncol = 10, nrow = noTimeSteps))

      startm <- Sys.time()
      to.mate <- which((mozzie.dt$gender == 1 & mozzie.dt$mateID == -1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$enzyme >= 1))
      endm <- Sys.time()
      mate.opto[t,1] <- endm - startm

      # if(length(to.mate) != 0){
      #   for(m in 1:length(to.mate)){
      #     # Grab the position of the female that's mating
      #     startm <- Sys.time()
      #     to.mate.pos    <- as.numeric(c("lat"  = mozzie.dt$lat[to.mate[m]],
      #                                    "long" = mozzie.dt$long[to.mate[m]]))
      #     # Draw a square around her with width 2*k: this is the search area
      #     to.mate.bdary  <- c("latmin"  = to.mate.pos[1]-param$k, "latmax"  = to.mate.pos[1] + param$k,
      #                         "longmin" = to.mate.pos[2]-param$k, "longmax" = to.mate.pos[2] + param$k)
      #     # Find possible mates within this search area
      #     possible.mates <- which(mozzie.dt$lat >= to.mate.bdary["latmin"] & mozzie.dt$lat <= to.mate.bdary["latmax"] & mozzie.dt$long >= to.mate.bdary["longmin"] & mozzie.dt$long <= to.mate.bdary["longmax"] & mozzie.dt$gender == 0)
      #     no.bachelors   <- length(possible.mates)
      #     endm <- Sys.time()
      #     mate.opto[t,2] <- endm - startm
      #
      #     startm <- Sys.time()
      #     if(no.bachelors == 0){
      #       # There were no suitable mates within the search area
      #       new.mate <- -1
      #       #print(paste0(to.mate[m], " no mate found"))
      #     }else if(no.bachelors == 1){
      #       if(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates){
      #         # Even though there is only one male, he has mated too much today
      #         new.mate <- -1
      #         # print(paste0(femID, " 1 mate found, unsuitable"))
      #       }else{
      #         #New mate is just the only male they found in the search area
      #         new.mate <- as.integer(possible.mates)
      #         mozzie.dt$gonoCycle[new.mate] <- mozzie.dt$gonoCycle[new.mate] + 1
      #         #print(paste0(femID, " 1 mate found"))
      #       }
      #     }else{
      #       if(length(which(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates)) == 0){
      #         # More than one suitable mate was found in the search area.
      #         # Randomly permutes the list of possible mates and then picks the one at the top of the pile
      #         new.mate <- sample(possible.mates)[1]
      #         # Increment number of mates of male for this timestep by 1
      #         mozzie.dt$gonoCycle[new.mate] <- mozzie.dt$gonoCycle[new.mate] + 1
      #         #print(paste0(femID, " multiple mates, all ok"))
      #       }else{
      #         # Remove males who have mated more than "max_daily_mates" number of times in a day
      #         possible.mates <- possible.mates[-(which(mozzie.dt[possible.mates]$gonoCycle >= max_daily_mates))]
      #         # Now we have one more condition to check for:
      #         ## If we drop some males from "possible.mates" we might end up dropping them all,
      #         ## in which case, we return -1
      #         if(length(possible.mates) == 0){
      #           #print(paste0(femID, " multiple mates, had to drop some AND then ended up with none"))
      #           new.mate <- 1
      #         }else{
      #           # Randomly permutes the list of possible mates and then picks the one at the top of the pile
      #           new.mate <- sample(possible.mates)[1]
      #           mozzie.dt$gonoCycle[new.mate] <- mozzie.dt$gonoCycle[new.mate] + 1 #increment number of mates by 1
      #           # print(paste0(femID, " multiple mates, had to drop some"))
      #         }
      #       }
      #       endm <- Sys.time()
      #       mate.opto[t,3] <- endm - startm
      #     }
      #     mozzie.dt$mateID[to.mate[m]] <- new.mate
      #     rm(to.mate.pos)
      #     rm(to.mate.bdary)
      #     rm(possible.mates)
      #     rm(no.bachelors)
      #     new.mate <- -1
      #   }
      # }
      # rm(to.mate)
      # # For male agents, gonoCycle has a different purpose:
      # # it tracks the number of times per day that the male has mated.
      # # So we now reset gonoCycle for the next day for all male agents:
      # mozzie.dt$gonoCycle[which(mozzie.dt$gender == 0)] <- 0

      #' MATING ATTEMPT NUMBER 3...... I S2G THIS HAS TO BE THE LAST TIME
      #' we do a vectorised for loop. the find_a_mate function works on one agent at a time.

      #' Get a list of females who can mate.
      to.mate <- which((mozzie.dt$gender == 1 & mozzie.dt$mateID == -1 & mozzie.dt$gonoCycle == 0 & mozzie.dt$enzyme >= 1))
      #' There's a ~30% chance that those females who can mate, will mate.
      to.mate <- sample(to.mate, size = length(to.mate)/5, replace = FALSE)

      system.time(
        mates_vector <- foreach (m = 1:length(to.mate), .combine=c) %dopar% {
          #' apply function to find mates
          find_a_mate(to.mate[m], param$k, max_daily_mates, mozzie.dt)
          #' Update female's mate ID
          #mozzie.dt$mateID[to.mate[m]]
          #' Update male's gonoCycle
        })
      # Update female's mate ID
      mozzie.dt$mateID[to.mate] <- mates_vector

      #PERHAPS JUST GET RID OF MALE'S GONOCYCLE? count the max number of times a male mates
      # below is for testing
      max_gono[t] <- max(tabulate(mates_vector))


      rm(to.mate, mates_vector)
      #  s5 <- s5 +1 #debugging
      end <- Sys.time()
      optimisation[t,7] <- end - start
      rm(start, end)
      # End Step 5

      # End of timestep procedures ----
      optimisation[t, 10] <- nrow(mozzie.dt)
      pb$tick() # This is for the progress bar

      # Post-simulation analytics and data collection ----
      # The below is some quick and dirty calculations to get Wolbachia proportion over time
      # To do: make the below into a function
      if(length(which(mozzie.dt$infStatus == 0)) == 0){
        inf.prop[t] <- 1
      }else if(length(which(mozzie.dt$infStatus == 1)) == 0){
        inf.prop[t] <- 0
      }else{
        inf.prop[t] <- length(which(mozzie.dt$infStatus == 1))/nrow(mozzie.dt)
      }

      # proportion of trapped mosquitoes with wb over time
      if(nrow(graveyard) > 0){
        if(length(which(graveyard$infStatus == 0 & graveyard$typeDeath == 1 & graveyard$timeDeath == t)) == 0){
          grav.inf.prop[t] <- 1
        }else if(length(which(graveyard$infStatus == 1 & graveyard$typeDeath == 1 & graveyard$timeDeath == t)) == 0){
          grav.inf.prop[t] <- 0
        }else{
          grav.inf.prop[t] <- length(which(graveyard$infStatus == 1 & graveyard$typeDeath == 1 & graveyard$timeDeath == t))/length(which(graveyard$infStatus == 1 & graveyard$timeDeath == t))
        }
      }else{
        grav.inf.prop[t] <- 0
      }


      # number of juvenile clutches with CI over time
      if(nrow(juv.graveyard) > 0){
        no.juv.CI[t] <- length(which(juv.graveyard$infProb == -1))
      }

      juv.N.vec[t] <- nrow(juv.dt)
      # number of adult agents over time
      N.vec[t] <- nrow(mozzie.dt)
    }

    #---- End of time loop ----
    print(timesrep)


    #' Set up stuff to get graveyard data


    trap_dates <- as.Date(unique(currentTraps$DateCollected))
    trap_dates <- trap_dates[which(trap_dates > as.Date("2013-01-03"))]

    daysvec <- seq(from = 1, to = 110, by = 1)
    datesvec <- seq(from = as.Date("2013-01-03"), to = as.Date("2013-04-21"), "days")

    uniquedates   <- unique(trap_dates) # The dates when traps are cleaned
    trap_days <- which(datesvec %in% uniquedates) # The days in the simulation these correspond to

    N_trapped_total   <- rep(0, length(trap_days))
    NWb_trapped_total <- rep(0, length(trap_days))
    NWT_trapped_total <- rep(0, length(trap_days))
    prop_wb_total     <- rep(0, length(trap_days))
    for(te in 1:length(trap_days)){

      if(te == 1){
        filtered_trap <- graveyard %>% filter(timeDeath <= trap_days[te],
                                              typeDeath == 1)
      }else{

        filtered_trap <- graveyard %>% filter(timeDeath <= trap_days[te],
                                              timeDeath >= trap_days[te - 1],
                                              typeDeath == 1)
      }

      N_trapped_total[te]   <- nrow(filtered_trap)
      NWb_trapped_total[te] <- length(which(filtered_trap$infStatus == 1))
      NWT_trapped_total[te] <- length(which(filtered_trap$infStatus == 0))
      prop_wb_total[te] <- NWb_trapped_total[te]/ N_trapped_total[te]
    }

    graveyard_plot_data <- as.data.frame(cbind(trap_days,
                                               N_trapped_total,
                                               NWb_trapped_total,
                                               NWT_trapped_total,
                                               prop_wb_total))

    #' end graveyard data

    #' capture data
    grav.inf.prop.matrix <- cbind(grav.inf.prop.matrix, grav.inf.prop)
    inf.prop.matrix <- cbind(inf.prop.matrix, inf.prop)
    no.juv.CI.matrix <- cbind(no.juv.CI.matrix, no.juv.CI)
    N.matrix <- cbind(N.matrix, N.vec)
    graveyard.plot.list[[timesrep]] <- graveyard_plot_data

    #' save graveyard as a csv
    graveyard$infStatus<- sapply(graveyard$infStatus, function(x) paste0(unlist(x), collapse = "\n"))
    save_graveyard <- as.data.frame(graveyard)
    #colnames(grav_inf_prop_matrix) <-  c("Day", (paste0("S",seq(1,timesRepeat, by = 1))))
    write_csv(save_graveyard, paste0(lubridate::hour(Sys.time()), "h_", timesrep, "_iter_graveyard.csv"))


  } # End of multiple simulations

  return(graveyard_plot_data)
}
