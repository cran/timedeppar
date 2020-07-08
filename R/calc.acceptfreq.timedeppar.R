#' Calculate apparent acceptance frequency of time-dependent parameter
#' 
#' This function calculates the apparent acceptance frequency from a potentially thinned Markov chain sample.
#' 
#' @param  x                  results from the function \code{infer.timedeppar} of class \code{timedeppar}.
#' @param  n.burnin           number of (unthinned) burnin points of the Markov chain to omit from analysis.
#' 
#' @return list two-column matrices with time and apparent acceptance frequencies

calc.acceptfreq <- function(x,n.burnin) UseMethod("calc.acceptfreq")

calc.acceptfreq.timedeppar <- function(x,n.burnin=0)
{
  res <- x
  
  # check input:
  
  if ( class(res) != "timedeppar" )
  {
    warning("argument not of class timedeppar")
    return()
  }
  
  # get thin and n.adapt:
  
  thin    <- 1; if ( !is.null(res$control) ) { if(!is.null(res$control$thin))    thin    <- res$control$thin } 
  n.adapt <- 0; if ( !is.null(res$control) ) { if(!is.null(res$control$n.adapt)) n.adapt <- min(res$control$n.adapt,res$control$n.iter-1) }

  # calculate apparent acceptance frequencies:

  if ( is.null(res$sample.param.timedep) )    stop("*** calc.acceptfreq: sample of time dependent parameters not found")
  if ( length(res$sample.param.timedep) < 1 ) stop("*** calc.acceptfreq: sample of time dependent parameters not found")
  acceptfreq <- list()
  for ( i in 1:length(res$sample.param.timedep) )
  {
    acceptfreq[[names(res$sample.param.timedep)[i]]] <- matrix(NA,ncol=2,nrow=ncol(res$sample.param.timedep[[i]]))
    acceptfreq[[names(res$sample.param.timedep)[i]]][,1] <- res$sample.param.timedep[[i]][1,]
    sample <- signif(res$sample.param.timedep[[i]][-1,],digits=6)
    ind.range <- 1:nrow(sample)
    if ( floor(max(n.adapt,n.burnin)/thin)+1 < nrow(sample) ) ind.range <- (floor(max(n.adapt,n.burnin)/thin)+1):nrow(sample)
    for ( j in 1:ncol(sample) ) acceptfreq[[names(res$sample.param.timedep)[i]]][j,2] <- (length(unique(sample[ind.range,j]))-1)/length(ind.range)
    colnames(acceptfreq[[names(res$sample.param.timedep)[i]]]) <- c("time","acceptfreq")
  }

  return(acceptfreq)
}
