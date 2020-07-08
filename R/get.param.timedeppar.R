#' Extract parameter list and process parameter vectors from an object of type \code{timedeppar}
#' 
#' This function extracts an element of the stored Markov chain from an object of type \code{timedeppar} 
#' produced by the function \code{\link{infer.timedeppar}} and converts it to a format that facilitates
#' re-evaluation of the posterior, evaluation of the underlying model, and sampling from the Ornstein-Uhlenbeck
#' process.
#' In particular, the constant and time-dependent parameters are provided in the same list format as supplied
#' as \code{param.ini} to the function \code{\link{infer.timedeppar}}.
#' In addition, the parameters of the Ornstein-Uhlenbeck process(es) of the time-dependent parameters are
#' provided in different formats (see details below under \code{Value}.
#' 
#' @param  x                  results from the function \code{\link{infer.timedeppar}} of class \code{timedeppar}.
#' @param  ind.sample         index of the stored (potentially thinned) Markov chain defining which parameters 
#'                            to reconstruct.
#'                            Default is to extract the parameters corresponding to the maximum posterior 
#'                            solution.
#' 
#' @return list with the following elements:\cr
#'         \code{param}: list of constant and time-dependent model parameters,\cr
#'         \code{param.ou.estim}: vector of estimated process parameters of all time-dependent parameters,\cr
#'         \code{param.ou.fixed}: vector of fixed process parameters of all time-dependent parameters,\cr
#'         \code{param.ou}: matrix of Ornstein-Uhlenbeck parameters for all time-dependent parameters;
#'         columns are mean, sd and gamma of the processes, rows are the time-dependent parameter(s).\cr
#'         \code{logpdf}: corresponding lopdf values (posterior, intermediate densities, and priors),\cr
#'         \code{ind.timedeppar}: indices of \code{param} at which parameters are time-dependent.\cr
#'         \code{ind.sample}: sample index.\cr
#'         \code{ind.chain}: corresponding index of the Markov chain.\cr


get.param <- function(x,ind.sample) UseMethod("get.param")

get.param.timedeppar <- function(x,ind.sample=NA)
{
  res <- x
  if ( class(res) != "timedeppar" )
  {
    warning("the first argument or get.param.timedeppar must be of class timedeppar")
    return(NA)
  }
  if ( is.null(res$param.ini) | 
       is.null(res$sample.param.const) | 
       is.null(res$sample.param.timedep) | 
       is.null(res$sample.param.ou) |
       is.null(res$sample.logpdf) |
       is.null(res$control) )
  {
    warning("get.param: the first argument must contain objects param.ini, sample.param.const, sample.param.timedep, sample.param.ou, sample.logpdf and control")
    return(NA)
  }
  if ( is.na(ind.sample) ) ind.sample <- which.max(res$sample.logpdf[,"logposterior"])
  if ( ind.sample < 1 | ind.sample > max(nrow(res$sample.param.const),nrow(res$sample.paramou)) )
  {
    warning("get.param: the argument ind.sample must either be NA or between 1 and the sample size of the stored (and potentially thinned) Markov chain")
    return(NA)
  }
  
  param <- res$param.ini
  ind.timedeppar <- numeric(0)
  for ( i in 1:length(param) )
  {
    parname <- names(param)[i]
    if ( is.matrix(param[[i]]) )
    {
      ind.timedeppar <- c(ind.timedeppar,i)
      if ( ! parname %in% names(res$sample.param.timedep) )
      {
        warning(paste("get.param: sample for time-dependent parameter",parname,"not found"))
        return(NA)
      }
      param[[i]][,2] <- as.numeric(res$sample.param.timedep[[parname]][1+ind.sample,])
    }
    else
    {
      if ( ! parname %in% colnames(res$sample.param.const) )
      {
        warning(paste("get.param: sample for constant parameter",parname,"not found"))
        return(NA)
      }
      param[[i]] <- as.numeric(res$sample.param.const[ind.sample,names(param)[i]])
    }
  }
  
  param.ou.estim <- res$sample.param.ou[ind.sample,]

  param.ou.fixed <- res$param.ou.fixed
  
  logpdf <- res$sample.logpdf[ind.sample,]
  
  param.ou <- matrix(nrow=length(ind.timedeppar),ncol=3)
  colnames(param.ou) <- c("mean","sd","gamma")
  if ( length(ind.timedeppar) > 0 )
  {
    rownames(param.ou) <- names(res$sample.param.timedep)
    for ( i in 1:length(ind.timedeppar) )
    {
      parname <-names(res$sample.param.timedep)[i]
      if ( paste(parname,"mean",sep="_") %in% names(param.ou.estim) )  param.ou[parname,"mean"]  <- param.ou.estim[paste(parname,"mean",sep="_")]
      else                                                             param.ou[parname,"mean"]  <- param.ou.fixed[paste(parname,"mean",sep="_")]
      if ( paste(parname,"sd",sep="_") %in% names(param.ou.estim) )    param.ou[parname,"sd"]    <- param.ou.estim[paste(parname,"sd",sep="_")]
      else                                                             param.ou[parname,"sd"]    <- param.ou.fixed[paste(parname,"sd",sep="_")]
      if ( paste(parname,"gamma",sep="_") %in% names(param.ou.estim) ) param.ou[parname,"gamma"] <- param.ou.estim[paste(parname,"gamma",sep="_")]
      else                                                             param.ou[parname,"gamma"] <- param.ou.fixed[paste(parname,"gamma",sep="_")]
    }
  }
    
  return(list(param          = param,
              param.ou.estim = param.ou.estim,
              param.ou.fixed = param.ou.fixed,
              param.ou       = param.ou,
              logpdf         = logpdf,
              ind.timedeppar = ind.timedeppar,
              ind.sample     = ind.sample,
              ind.chain      = (ind.sample-1)*res$control$thin))
}
