#' library(mvtnorm)
#' Function to run SMC ABC
#' Created by Carla Chen and Sarah Belet
#' @param tol tolerance level
#' @param Np Number of particles
#' @param alpha proportion of particles dropped each step
#' @param x thingo vector
#' @param obs observations vector
#' @param c_v Probability that a particle is not moved at least once
ABC_rep <- function(tol = 101, Np, alpha, x, obs, c_v){
  #ini.tol: initial tolerance level
  #Np: number of particles
  #tol: tolerance level
  #alpha: proportion of particles dropped
  #c_v: is the probability that particle is not moved at least once

  #' ---- Initial setup ---- #
  Nd   <- floor(Np*alpha)
  para <- data.frame(matrix(-99, nrow = Np, ncol = 7))
  # colnames(para) <- c("lambda",
  #                        "k",
  #                        "phi", #hyperparameter for probability that mozzies are trapped by any given trap
  #                        "alpha", #CHANGE #this is for natural deathac
  #                        "p_1", #parameter for CI (prop. of uninfected larvae given both parents are infected)
  #                        "eta_1", #no. of offspring produced by Wb non-carrier mothers
  #                        "eta_2", #no. of offspring produced by Wb carrier mothers
  #                        "rho")
  colnames(para) <- c("lambda", #dispersal
                      "k", #trapping
                      "phi", #hyperparameter for probability that mozzies are trapped by any given trap
                      "alpha_a", #adult mortality
                      "alpha_j", # juvenile mortality
                      "eta",
                      "rho")

  #' Set initial tolerance
  ini.tol <- 100000

  #' ---- Rejection ABC part --------------
  for (i in 1:Np){
    flag  <- 1
    p     <- para[i,]
    p$rho <- 100001
    while(p$rho > ini.tol || flag == 1){
      #p$lambda = rgamma(1,3,22)
      #p$k = rgamma(1,shape=1,scale=1/35)
      #p$phi = runif(1,0,0.1) #hyperparameter for probability that mozzies are trapped by any given trap
      #p$alpha = runif(1,0.001,0.3)  #CHANGE #this is for natural death
      #p$p_1 = runif(1,0,0.05) #parameter for CI (prop. of uninfected larvae given both parents are infected)
      #p$eta_1 = round(rnorm(1,mean=40,sd=7),0) #no. of offspring produced by Wb non-carrier mothers
      #p$eta_2 = round(rnorm(1,mean=30,sd=5),0) #no. of offspring produced by Wb carrier mothers
      #start test
      #mu<-0.12+1.75*x1+2.5*x2-1.4*x3+6.7*x4-0.5*x5-1.23*x6+2.54*x7

      # p$lambda = rgamma(1,3,22)
      # p$k      = rgamma(1,shape = 1, scale = 1/35)
      # p$phi    = runif(1, 0.00001,0.001) #hyperparameter for probability that mozzies are trapped by any given trap
      # p$alpha  = runif(1, 0.05, 0.16)  #CHANGE #this is for natural death
      # p$p_1    = runif(1,-1,1) #parameter for CI (prop. of uninfected larvae given both parents are infected)
      # p$eta_1  = rnorm(1,mean=1,sd=5) #no. of offspring produced by Wb non-carrier mothers
      # p$eta_2  = rnorm(1,mean=5,sd=5) #no. of offspring produced by Wb carrier mothers


      p$lambda = rgamma(1,3,22)
      p$k      = rgamma(1,shape = 1, scale = 1/35)
      p$phi    = runif(1, 0.00001,0.001) #hyperparameter for probability that mozzies are trapped by any given trap
      p$alpha_a  = runif(1, 0.05, 0.16)
      p$alpha_j  = runif(1, 0.1, 0.4)
      p$eta  = rnorm(1, mean = 20, sd = 5) #no. of offspring produced by Wb non-carrier mothers


      # Sarah'a model here run with parameters- spits out the proportion of infected individuals at a specified time
      #flag -

      # as of 14/2, this runs the model with the output being the graveyard.
      model_grav <- mozzie::model_run(p$lambda, p$k, p$phi, p$alpha_a, p$alpha_j, p$eta_1)

      pred <- 0.12 + p$lambda*x[,1] + p$k*x[,2] + p$phi*x[,3] + p$alpha*x[,4] +
                     p$p_1*x[,5] + p$eta_1*x[,6] + p$eta_2*x[,7]
      flag <- 0
      # Summary statistics? proportion of animals infected | trapped

      # end test
      # read in the trap dat

      p$rho <- sum((sim - pred)^2)
    }
    para[i,] <- p
  }
  #set initial R value

  R <- 100
  e_max <- max(para$rho)

  #replenishment part
  while(e_max > tol){
    # sort particles- from large to small discrepancy value
    para <- para[order(para$rho, decreasing = TRUE),]

    # replenish particles with the largest Nd values
    e_next   <- para$rho[Nd + 1]
    keep     <- para[(Nd + 1):nrow(para), 1:(ncol(para) - 1)]
    prop_sig <- var(keep)
    i_acc    <- 0
    for (j in 1:Nd){
      #proposal mean
      p_mean <- keep[sample(nrow(keep), 1), -ncol(para)]

      #Theta^{N_\alpha+1}
      curr <- para[j, -ncol(para)]

      #initial condition making sure while-loop runs
      k <- 1
      repeat{
        #Theta**
        prop <- rmvnorm(1, as.matrix(p_mean), prop_sig)

        #Metropolis-Hasting
        #numerator
        #log-sum
        #nu_p1<- log(dgamma(curr$lambda,3,22))+log(dgamma(curr$k,shape=1,scale=1/35))+log(1/0.1)+log(1/(0.3-0.001))+log(0.05)+log(dnorm(curr$eta_1,mean=40,sd=7))+log(dnorm(curr$eta_2,mean=30,sd=5))
        #test----------
        nu_p1 <- log(dnorm(curr$lambda, 3, 5)) + log(dnorm(curr$k, 3, 5)) +
                 log(1/20) + log(1/(5)) + log(1/2) + log(dnorm(curr$eta_1, mean = -1, sd = 5)) +
                 log(dnorm(curr$eta_2, mean = 5, sd = 5))
        #test------------------------
        nu_p2 <- log(dmvnorm(curr, mean = prop, prop_sig))

        #denominator
        #de_p1<-log(dgamma(prop$lambda,3,22))+log(dgamma(prop$k,shape=1,scale=1/35))+log(1/0.1)+log(1/(0.3-0.001))+log(0.05)+log(dnorm(prop$eta_1,mean=40,sd=7),0)+log(dnorm(prop$eta_2,mean=30,sd=5),0)
        #test
        de_p1 <- log(dgamma(prop[ , 1], 3, 22)) +
                 log(dnorm(prop[ , 2], 3, 5)) +
                 log(1/20) + log(1/5) + log(1/2) + log(dnorm(prop[ , 6], mean = -1, sd = 5)) +
                 log(dnorm(prop[ ,7], mean = 5, sd = 5))
        #--------------test

        de_p2 <- log(dmvnorm(prop, mean = as.matrix(curr), prop_sig))

        #proposal
        MH <- min(0, nu_p1 + nu_p2 - de_p1 - de_p2)
        u  <- log(runif(1, 0, 1))
        if(u < MH){
          #run Sarah's code with proposal - prop
          #test------------------------------------
          #test
          #mu<-0.12+1.75*x1+2.5*x2-1.4*x3+6.7*x4-0.5*x5-1.23*x6+2.54*x7
          #flag -
          pred <- 0.12 + prop[,1]*x[,1] +
                  prop[,2]*x[,2] +
                  prop[,3]*x[,3] +
                  prop[,4]*x[,4] +
                  prop[,5]*x[,5] +
                  prop[,6]*x[,6] +
                  prop[,7]*x[,7]
          #test---------------------------------

          rho <- sum((obs - pred)^2)

          if(prop[ , ncol(prop)] <= e_next){
            para[j,] <- cbind(prop, rho)
            i_acc    <- i_acc + 1
            break;
          }#if
        }#if
        k <- k + 1
        if(k > R){break}
      }#k
    }#j
    p_acc <- i_acc/(R*Nd)
    R     <- ceiling(log(c_v)/log(1 - p_acc))
    e_max <- max(para$rho)
  }
  return(para)
}



