#' Plots positions of mosquitoes in simulation region.
#'
#' Due to the number of agents in the simulation, this is not
#' a good summary plot. Try a density plot instead.
#' @param mozzie.dt data.table of adult agents
#' @return A plot of mozzie positions.
plot_positions <- function(mozzie.dt){
  g <- ggplot(data = mozzie.dt, aes(x = lat, y = long)) + geom_point() +
       ggtitle("Position of adult agents")
  return(g)
}


#' Plots proportion of mosquitoes trapped that have Wolbachia over time for a buffer zone.
#'
#' The main buffers we should care about is 0-100m, 100-200m, then 200m and beyond.
#'
#' @param graveyardTrapped Subset of graveyard for trapped mosquitoes only. This is modified graveyard
#' data with only columns whereTrapped, infStatus, timeDeath, Buf and count
#' @param buf Buffer zone to plot.
#' @return A plot of trapped Wolbachia proportions over time for a given buffer zone.
plot_proportions_by_buf <- function(graveyardTrapped, buf){
  by.buffer        <- graveyardTrapped[which(graveyardTrapped$Buf == buf), ]
  buf.totals       <- aggregate(by.buffer$count, by = list(Date = by.buffer$timeDeath), FUN = sum)
  buf.by.infStatus <- aggregate(by.buffer$count, by = list(Date = by.buffer$timeDeath,  infStatus = by.buffer$infStatus), FUN = sum)
  # The issue with the above line is that since Wb-carrying mosquitoes aren't trapped until Day 9,
  ## so there are no entries with infStatus = 1 for Days 1-9.
  ## obviously this needs to be changed if release data is different.
  nulldays <- as.data.frame(cbind(seq(1:8), rep(1,8), rep(0,8)))
  colnames(nulldays) <- c("Date", "infStatus", "x")
  wb.cases <- buf.by.infStatus[which(buf.by.infStatus$infStatus == 1), ]
  wb.cases <- rbind(nulldays, wb.cases)
  props <- wb.cases$x/buf.totals$x
  wb.cases <- cbind(wb.cases, props)


  # The above throws an error if there's another day where no mozzies are caught
  day.setup <- as.data.frame(cbind(seq(1:90), rep(1,90), rep(0,90)))
  colnames(day.setup) <- c("Date", "infStatus", "x")
  day.setup$x[day.setup$Date %in% wb.cases$Date] <- wb.cases$x + day.setup$x[day.setup$Date %in% wb.cases$Date]
  props <- day.setup$x/buf.totals$x



  for(i in 1:90){

  }


  g <- ggplot(data = wb.cases, aes(x = Date, y = props)) + geom_line() +
       ggtitle("Proportion of trapped mosquitoes that carry Wolbacchia \n for buffer =",buf)

  #g <- ggplot(data = mozzie.dt, aes(x = lat, y = long)) + geom_point()
  return(g)
}
