#' Runs the 'juvenile emergence' step.
#' Determines which juvenile agents are to age into adult agents,
#'   and calls \code{juv_to_adult} on them.
#' IMPORTANT: this should update \code{mozzie.dt} and \code{juv.dt} on a global scale.
#'
#' @param juv.dt The juvenile data.table.
juv_emergence <- function(juv.dt){
  toAdult <- which(juv.dt$stage == 3 & juv.dt$enzyme > 0.95)
  if((length(toAdult)) != 0){
    newAdults.dt <- mapply(FUN=juv_to_adult, toAdult, idStart, pmale, lambda, SIMPLIFY = FALSE)
    newAdults.dt %<>% rbindlist() %>% as.data.table()
    mozzie.dt    <<- rbind(mozzie.dt, newAdults.dt)

    #delete row(s) of data table corresponding to the clutches that aged into adult mozzies
    juv.dt       <<- juv.dt[-toAdult]
    rm(newAdults.dt)
  }
  rm(toAdult)
}


#---Step 1.5: any juveniles moving to the next development stage?

#if enzyme > 0.95, move on to the next development stage and reset enzyme score

#' Determines the juvenile agents that develop into the next stage, and
#'   updates the development stage of these juveniles.
#' Once \code{enzyme} exceeds 0.95, a juvenile agent develops to the next stage.
#' Recall that Stage 1: egg, 2: larvae, 3: pupae.
#' Any pupae that are to develop into the next stage, i.e. adulthood, should
#'   have been handled by \code{juv_emergence}.
#' IMPORTANT: this should update \code{juv.dt} on a global scale.
#' @param juv.dt The juvenile data.table.
juv_stage_update <- function(juv.dt){
  if(length(which(juv.dt$enzyme > 0.95))){
    #which(juv.dt$enzyme > 0.95)
    juv.dt[which(juv.dt$enzyme > 0.95)]$stage  <<- juv.dt[which(juv.dt$enzyme > 0.95)]$stage + 1

    #this gives HELLA warnings- try mapply or lapply or something?
    juv.dt[which(juv.dt$enzyme > 0.95)]$enzyme <<- 0
    #juv.dt[which(juv.dt$enzyme > 0.95)]$stage <- lapply(juv.dt[which(juv.dt$enzyme > 0.95)]$stage, function(x) x <- juv.dt$enzyme[x] <- 0)
  }
}

#' Determines which mothers are ready to lay eggs and calls \code{initialise_eggs}.
#' Conditions to be satisfied for a mother to lay a clutch of eggs:
#' \enumerate{
#'  \item {Has a mate}
#'  \item {Enzyme Kinetic Score > 1 and gonoCycle = 0 OR}
#'  \item {EKS > 1.58 and gonoCycle = 1 OR}
#'  \item {EKS > 2.16 and gonoCycle = 2.}
#' }
#'
#' IMPORTANT: This should update \code{juv.dt} on a global scale.
#' @param mozzie.dt data.table of adult agents.
produce_eggs <- function(mozzie.dt){
  toLay <- which(mozzie.dt$gender == 1 & mozzie.dt$mateID != -1 & ((mozzie.dt$enzyme >= 1 & mozzie.dt$gonoCycle == 0) | (mozzie.dt$enzyme >= 1.58 & mozzie.dt$gonoCycle == 1) | (mozzie.dt$enzyme >= 2.16 & mozzie.dt$gonoCycle == 2)))
  if(length(toLay) != 0){
    newEggsdt <- initialise_eggs(toLay)
    juv.dt    <- rbind(juv.dt, newEggsdt)
    rm(newEggsdt)
  }
  rm(toLay)
}

#---Step 3: Do we have mosquitoes being trapped (and killed?)---#
##--- End Step 3


#' Runs the 'juvenile emergence' step.
#' Determines which juvenile agents are to age into adult agents,
#'   and calls \code{juv_to_adult} on them.
#' IMPORTANT: this should update \code{mozzie.dt} and \code{juv.dt} on a global scale.
#'
#' @param juv.dt The juvenile data.table.
