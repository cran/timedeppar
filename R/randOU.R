#' Draw from an Ornstein-Uhlenbeck process
#' 
#' This function draws a realization of an Ornstein-Uhlenbeck process 
#' with a random start value (drawn from the marginal distribution),
#' conditional of the start value, or conditional on both start and end values.
#' The function includes the option of obtaining a lognormal marginal by 
#' exponential tranformation.
#' 
#' @param  mean   asymptotic mean of the process.
#' @param  sd     asymptotic standard deviation of the process
#' @param  gamma  rate coefficient for return to the mean
#' @param  t      vector of time points at which the process should be sampled
#'                (note: the value at t[1] will be the starting value yini,
#'                the value at t[length(t)] the end value yend if these are specified)
#' @param  yini   start value of the process 
#'                (NA indicates random with asymptotic mean and sd)
#' @param  yend   end value of the process 
#'                (NA indicates no conditioning at the end)
#' @param  log    indicator whether the log of the variable should be an Ornstein-Uhlenbeck 
#'                process (log=TRUE) rather than the variable itself
#'                (mean and sd are interpreted in original units also for log=TRUE)
#' 
#' @return a data frame with t and y columns for time and for the realization of the Ornstein-Uhlenbeck process
#' 
#' @examples 
#' plot(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000),type="l",ylim=2.5*c(-1,1))
#' abline(h=0)
#' lines(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000),col="red")
#' lines(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000),col="blue")
#' lines(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000),col="green")
#' 
#' plot(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000,yini=0,yend=0),type="l",ylim=2.5*c(-1,1))
#' abline(h=0)
#' lines(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000,yini=0,yend=0),col="red")
#' lines(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000,yini=0,yend=0),col="blue")
#' lines(randOU(mean=0,sd=1,gamma=1,t=0:1000/1000,yini=0,yend=0),col="green")


randOU <- function(mean=0,sd=1,gamma=1,t=0:1000/1000,yini=NA,yend=NA,log=FALSE)
{
  # consistency checks:
  
  n <- length(t)
  if ( min(diff(t)) <= 0 )
  {
     warning("t must be strictly increasing")
     return(NA)
  }
  if ( sd <= 0 )
  {
    warning("sd must be positive")
    return(NA)
  }
  if ( gamma <= 0 )
  {
    warning("gamma must be positive")
    return(NA)
  }
  
  res <- NA
  if ( log ) 
  {
    # if log=TRUE calculate meanlog and sdlog (mean and sd are in original units),
    # call function with these parameters, and exponentiate the result:
    
    if ( mean <= 0 )
    {
      warning("if log=TRUE, mean must be positive")
      return(NA)
    }
    if ( !is.na(yini) )
    {
      if ( yini <= 0 )
      {
        warning("if log=TRUE, yini must be NA or positive")
        return(NA)
      }
    }
    if ( !is.na(yend) )
    {
      if ( yend <= 0 )
      {
        warning("if log=TRUE, yend must be NA or positive")
        return(NA)
      }
    }
    meanlog <- log(mean/(sqrt(1+sd*sd/(mean*mean))))
    sdlog   <- sqrt(log(1+sd*sd/(mean*mean)))
    res <- randOU(mean  = meanlog,
                  sd    = sdlog,
                  gamma = gamma,
                  t     = t,
                  yini  = log(yini),
                  yend  = log(yend),
                  log   = FALSE)
    res$y <- exp(res$y)
  }
  else
  {  
    # draw and return sample:
  
    y <- rep(NA,n)
    if ( is.na(yini) ) 
    {
      y[1] <- rnorm(n=1,mean=mean,sd=sd) 
    }
    else 
    {
      y[1] <- yini  # initial point
    }
    if ( n == 1 )
    {
      if ( !is.na(yend) )
      {
        if ( is.na(yini) )
        {
          y[1] <- yend
        }
        else
        {
          if ( yini != yend )
          {
            warning("if conditioned on both ends and of length 1, the two values must be the same")
            return(NA)
          }
        }
      }
    }
    else
    {
      if ( is.na(yend) )
      {
        for ( i in 2:n )  # iterative sampling conditional on previous point
        {
          y[i] <- rnorm(n    = 1,
                        mean = mean+(y[i-1]-mean)*exp(-gamma*(t[i]-t[i-1])),
                        sd   = sd*sqrt(1-exp(-2*gamma*(t[i]-t[i-1]))))
        }
      }
      else
      {
        y[n] <- yend  # end point
        if ( n > 2 )
        {
          for ( i in 2:(n-1) )  # iterative sampling conditional on previous point and end point
          {
            y[i] <- rnorm(n    = 1,
                          mean = mean+(y[i-1]-mean)*exp(-gamma*(t[i]-t[i-1]))*(1-exp(-2*gamma*(t[n]-t[i])))/
                                                    (1-exp(-2*gamma*(t[n]-t[i-1]))) +
                                      (y[n]-mean)  *exp(-gamma*(t[n]-t[i]))*(1-exp(-2*gamma*(t[i]-t[i-1])))/
                                                    (1-exp(-2*gamma*(t[n]-t[i-1]))),
                          sd   = sd*sqrt((1-exp(-2*gamma*(t[i]-t[i-1])))*(1-exp(-2*gamma*(t[n]-t[i])))/
                                         (1-exp(-2*gamma*(t[n]-t[i-1])))))
          }
        }
      }
    }
    res <- data.frame(t=t,y=y)
  }
  return(res)
}

