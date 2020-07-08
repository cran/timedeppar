#' Draw indices for a random split of a vector into intervals of the same mean length
#' 
#' This function draws indices for a random split of a vector into sub-vectors.
#' 
#' @param  n.grid       number of grid points to divide into intervals
#' @param  n.interval   number of intervals
#' @param  method       method for random splitting:\cr
#'                      \code{modunif}  modification of uniform intervals\cr
#'                      \code{random}   random split (higher variability in inverval lengths)\cr
#'                      \code{weighted} random split with weights; 
#'                          non-normalized weights must be specified by the argument \code{weights}\cr
#' @param  weights      weights for choosing interval boundaries for method \code{weighted}; 
#'                      vector of length i2-i1+1 (does not need to be normalized and will be
#'                      ignored for all methods except for method \code{weighted})
#' @param  offset       offset to shift subset of potential interval boundaries to draw from.
#'                      To guarantee different intervals on subsequent calls, offset should be
#'                      increased by one between subsequent calls for the same variable.
#' @param  min.internal minimum number of internal points between interval boundary points
#' 
#' @return the function returns an index vector of length n+1 with the endpoint indices of
#'         the random intervals. 
#'
#' @examples 
#' randsplit(100,10)
#' randsplit(100,10)
#' randsplit(100,10,method="random")
#' randsplit(100,10,method="weighted",weights=1:100)
#' for ( i in 1:10 ) print(randsplit(100,10,method="weighted",weights=1:100,offset=i))

randsplit <- function(n.grid,n.interval,method=c("modunif","random","weighted"),
                      weights=numeric(0),offset=0,min.internal=2)
{
  n.interval <- round(n.interval,0); if ( n.interval <= 0 ) n.interval <- 1
  min.internal <- max(1,round(min.internal,0))
  if ( n.interval == 1 ) return(c(1,n.grid))
  
  if ( method[1] == "modunif")
  {
    n.max <- max(floor((n.grid-1)/4),1)  # at least three points in between for equal spacing
    if ( n.interval > n.max ) stop(paste("*** randsplit: number of splits",n.interval-1,"too large for",n.grid,"indices; maximum is",n.max))
    l.mean <- (n.grid-1)/n.interval
    endpoints <- c(1,round(1+1:(n.interval-1)*l.mean,0),n.grid)
    for ( i in 2:n.interval ) endpoints[i] <- round(runif(1,endpoints[i]-l.mean/3,endpoints[i]+l.mean/3),0)
  }
  else
  {
    if ( n.grid < (n.interval+1)*(min.internal+1) ) stop(paste("*** randsplit: n.grid[",n.grid,"] < (n.interval[",n.interval,"]+1)",
                                                               "*(min.internal[",min.internal,"]+1)",sep=""))
    if ( method[1] == "random" )
    {
      offset <- offset %% (min.internal+1)
      subset <- seq(from=1+(min.internal+1)+offset,to=n.grid-(min.internal+1),by=min.internal+1)
      if ( length(subset) == 1 )
      {
        if ( n.interval == 2 ) samp <- subset
        else                   stop("*** randsplit: internal error")
      }
      else
      {
        samp <- sort(sample(subset,size=n.interval-1))
      }
      endpoints <- c(1,samp,n.grid)
    }
    else
    {
      if ( method[1] == "weighted" )
      {
        if ( length(weights)>0 & length(weights) != n.grid )
        {
          weights <- rep(1,n.grid)
          warning("randsplit: weight vector with incorrect length, weights ignored")
        }
        if ( length(weights) == 0 ) weights <- rep(1,n.grid)
        offset <- offset %% (min.internal+1)
        subset <- seq(from=1+(min.internal+1)+offset,to=n.grid-(min.internal+1),by=min.internal+1)
        if ( length(subset) == 1 )
        {
          if ( n.interval == 2 ) samp <- subset
          else                   stop("*** randsplit: internal error")
        }
        else
        {
          subweights <- weights[subset]
          samp <- sort(sample(subset,size=n.interval-1,prob=subweights))
        }
        endpoints <- c(1,samp,n.grid)
      }
      else
      {
        stop(paste("unknown splitting method: \"",method[1],"\"",sep=""))
      }
    }
  }
  return(endpoints)
}

