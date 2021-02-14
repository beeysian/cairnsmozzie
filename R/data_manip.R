#' Calculates the average daily temperature over the simulation period.
#' Temperature calculations are in Celcius.
#'
#' Function splits the data into separate months, lays it end to end,
#' and then calculates the mean temperature for each day of the
#' simulation period.
#'
#' Note that due to months not having the same number of days, we remove
#' any -1 entries which correspond to days like Feb 30.
#'
#' This function is deterministic and should yield the same result for
#' the same input.
#'
#' @param mindf Raw data of all minimum daily temperatures, read into a dataframe.
#' @param maxdf Raw data of all maximum daily temperatures, read into a dataframe.
#' @param noTimeSteps Number of timesteps in the simulation.
#' @return The mean temperature for each day of the simulation period, expressed as
#' a dataframe of length \code{noTimeSteps}.
temperature_chart <- function(mindf, maxdf, noTimeSteps){

  # Old stuff- leaving it in in case I break something
  #oct <- rowMeans(cbind(mindf$Oct,maxdf$Oct))
  #nov <- rowMeans(cbind(mindf$Nov,maxdf$Nov))
  #dec <- rowMeans(cbind(mindf$Dec,maxdf$Dec))
  #dailyTemps <- c(oct,nov,dec)
  #dailyTemps <- dailyTemps[which(dailyTemps!=-1)]
  #dailyTemps <- dailyTemps[1:(noTimeSteps+1)]

  jan <- rowMeans(cbind(mindf$Jan, maxdf$Jan))
  feb <- rowMeans(cbind(mindf$Feb, maxdf$Feb))
  mar <- rowMeans(cbind(mindf$Mar, maxdf$Mar))
  apr <- rowMeans(cbind(mindf$Apr, maxdf$Apr))
  dailyTemps <- c(jan, feb, mar, apr)
  dailyTemps <- dailyTemps[which(dailyTemps!=-1)]
  dailyTemps <- dailyTemps[1:(noTimeSteps)]

  return(dailyTemps)
}

#' Initialises the Enzyme Kinetic model lookup table for the model run.
#' The EKM is how we dictate the development rate of agents.
#' The EKM model is the same one used in Focks 1993.
#' Input data is in CELSIUS, EKM works with KELVIN.
#'
#' The model uses different sets of experimentally-derived parameters
#' for each four development stage: Egg, Larvae, Pupae, Gonotrophic.
#' Hence there are four equations.
#'
#' @section Notation and constants:
#' KELV_CONV: Conversion constant between Kelvin <-> Celcius.
#' R: Universal gas constant.
#' All other constants are experimentally-determined constants
#' which can be found in Focks 1993.
#' Suffixes E, L, P, G for the constants refer to each of the
#' development stages.
#'
#' @section "Why do you multiply by 24?":
#' We multiply by 24 since timesteps are in days, and EKM model
#' is in timestep unit of hours.
#'
#' @param dailyTemps Dataframe of daily average temperatures as output from
#' function temperature_chart
#' @return Chart of EKM updates for each development stage for each day.
initialize_enzyme <- function(dailyTemps, const){
  kelvinTemps <- dailyTemps + const$KELV_CONV

  EKM_egg    <- 24*((const$RHO_E*(kelvinTemps/298)*exp((const$HA_E/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_E/const$R)*(((1/const$THALF_E)-(1/kelvinTemps))))))
  EKM_larval <- 24*((const$RHO_L*(kelvinTemps/298)*exp((const$HA_L/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_L/const$R)*(((1/const$THALF_L)-(1/kelvinTemps))))))
  EKM_pupal  <- 24*((const$RHO_P*(kelvinTemps/298)*exp((const$HA_P/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_P/const$R)*(((1/const$THALF_P)-(1/kelvinTemps))))))
  EKM_gono   <- 24*((const$RHO_G*(kelvinTemps/298)*exp((const$HA_G/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_G/const$R)*(((1/const$THALF_G)-(1/kelvinTemps))))))

  EKM_chart  <- as.data.table(cbind(EKM_egg,EKM_larval,EKM_pupal,EKM_gono))

  return(EKM_chart)
}

#' New method of calculating Enzyme Kinetic Score for daily update (12-09-20)
#' This function takes in a vector of temperatures in Kelvin,
#' after a temperature deviate has been applied (if applicable).
#' Equation is from Focks 93.
#' @param kelvinTemps Vector of agent temperatures, in Kelvin
#' @param stage Development stage of agent, where:
#' \enumerate{
#' \item Egg
#' \item Larvae
#' \item Pupae
#' \item Adult
#' }
#' @return A vector of daily enzyme kinetic scores for each agent.
calculate_daily_EKS <- function(kelvinTemps, stage, const){
  if(stage == 1){
    EKS <- 24*((const$RHO_E*(kelvinTemps/298)*exp((const$HA_E/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_E/const$R)*(((1/const$THALF_E) - (1/kelvinTemps))))))
  }else if(stage == 2){
    EKS <- 24*((const$RHO_L*(kelvinTemps/298)*exp((const$HA_L/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_L/const$R)*(((1/const$THALF_L) - (1/kelvinTemps))))))
  }else if(stage == 3){
    EKS <- 24*((const$RHO_P*(kelvinTemps/298)*exp((const$HA_P/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_P/const$R)*(((1/const$THALF_P) - (1/kelvinTemps))))))
  }else{
    EKS <- 24*((const$RHO_G*(kelvinTemps/298)*exp((const$HA_G/const$R)*((1/298) - (1/kelvinTemps))))/(1 + exp((const$HH_G/const$R)*(((1/const$THALF_G)-(1/kelvinTemps))))))
  }

return(EKS)
}



#' Makes plot data for one graveyard dataframe.
#' Because I'm tired of copying and pasting the same code over and over...
#'
#' @param grav.data A graveyard data.frame.
make_plot_data <- function(grav.data){

  plot.data <- grav.data %>% mutate(trapEmptyDay = case_when(timeDeath <= 8 ~ 6,
                                                                timeDeath <= 13 ~ 13,
                                                                timeDeath <= 20 ~ 20,
                                                                timeDeath <= 27 ~ 27,
                                                                timeDeath <= 34 ~ 34,
                                                                timeDeath <= 41 ~ 41,
                                                                timeDeath <= 48 ~ 48,
                                                                timeDeath <= 55 ~ 55,
                                                                timeDeath <= 62 ~ 62,
                                                                timeDeath <= 69 ~ 69,
                                                                timeDeath <= 76 ~ 76,
                                                                timeDeath <= 83 ~ 83,
                                                                timeDeath <= 91 ~ 91,
                                                                timeDeath <= 97 ~ 97,
                                                                timeDeath <= 104 ~ 104,
                                                                timeDeath > 104 ~ 110)) %>%
    group_by(trapEmptyDay, infStatus) %>%
    summarise(trapped = n()) %>%
    mutate(prop = trapped/sum(trapped))

  plot.data  <- plot.data %>% filter(infStatus == 1)


return(plot.data)

}





