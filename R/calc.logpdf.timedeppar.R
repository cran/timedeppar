#' Calculate log pdf values (prior, internal, posterior) from an object of type \code{timedeppar}
#' 
#' This function calculated log priors, log pdf of Ornstein-Uhlenbeck time dependent parameters, and
#' log posterior from an object of type \code{timedeppar} produced by the function \code{\link{infer.timedeppar}}.
#' 
#' @param  x                  results from the function \code{\link{infer.timedeppar}} of class \code{timedeppar}.
#' @param  param              list of parameter lists and vectors extracted from an object of class
#'                            \code{timedeppar} using the function \code{\link{get.param}}.
#' @param  verbose            boolean indicator for writing the results to the console
#'                            (default is not to do it).
#' 
#' @return numeric vecctor of values of the different log pdfs.


calc.logpdf <- function(x,param,verbose) UseMethod("calclogpdf")

calc.logpdf.timedeppar <- function(x,param,verbose=FALSE)
{
  res <- x
  
  if ( is.null(res$dot.args) ) 
  {
    warning("results from newer version of infer.timedeppar needed")
    return(NA)
  }
  
  param.names    <- names(param$param)
  param.ou.names <- names(param$param.ou.estim)
  ind.timedeppar <- numeric(0)
  ind.constpar   <- numeric(0)
  for ( i in 1:length(param$param) )
  {
    if ( is.matrix(param$param[[i]]) | is.data.frame(param$param[[i]]) ) ind.timedeppar <- c(ind.timedeppar,i)
    else                                                                 ind.constpar   <- c(ind.constpar,i)
  }
  
  logposterior <- NA
  
  loglikeliobs <- do.call(res$loglikeli,c(list(param$param),res$dot.args))
  
  logprior.const <- 0
  if ( length(ind.constpar) > 0 ) logprior.const <- as.numeric(res$param.logprior(unlist(param$param[ind.constpar])))
  
  logpdf <- c(logposterior      = logposterior,
              loglikeliobs      = loglikeliobs,
              logprior_constpar = logprior.const)
  
  if ( length(ind.timedeppar) > 0 )
  {
    logpriorou <- rep(0,length(ind.timedeppar))
    logpdfou   <- rep(0,length(ind.timedeppar))
    for ( i in 1:length(ind.timedeppar) )
    {
      name <- param.names[ind.timedeppar[i]]
      par.ou  <- numeric(0)
      if ( paste(name,"mean",sep="_") %in% param.ou.names )
      {
        par.mean  <- param$param.ou.estim[[paste(name,"mean",sep="_")]]
        names(par.mean) <- paste(name,"mean",sep="_")
        par.ou <- c(par.ou,par.mean)
      }
      else
      {
        par.mean  <- param$param.ou.fixed[paste(name,"mean",sep="_")]
      }
      if ( paste(name,"sd",sep="_") %in% param.ou.names )
      {
        par.sd    <- param$param.ou.estim[[paste(name,"sd",sep="_")]]
        names(par.sd) <- paste(name,"sd",sep="_")
        par.ou <- c(par.ou,par.sd)
      }
      else
      {
        par.sd    <- param$param.ou.fixed[paste(name,"sd",sep="_")]
      }
      if ( paste(name,"gamma",sep="_") %in% param.ou.names )
      {
        par.gamma <- param$param.ou.estim[[paste(name,"gamma",sep="_")]]
        names(par.gamma) <- paste(name,"gamma",sep="_")
        par.ou <- c(par.ou,par.gamma)
      }
      else
      {
        par.gamma <- param$param.ou.fixed[paste(name,"gamma",sep="_")]
      }
      par.log <- FALSE; if ( name %in% names(res$param.log) ) { if(res$param.log[name]) par.log <- TRUE }
      
      if ( length(par.ou) > 0 ) logpriorou[i] <- res$param.ou.logprior(par.ou)
      logpdfou[i] <- logpdfOU(param$param[[ind.timedeppar[i]]][,1],param$param[[ind.timedeppar[i]]][,2],
                              mean=par.mean,sd=par.sd,gamma=par.gamma,cond=0,log=par.log)
      names(logpriorou)[i] <- paste("logprior_oupar",name,sep="_")
      names(logpdfou)[i]   <- paste("logpdfou_timedeppar",name,sep="_")
    }
    logpdf <- c(logpdf,logpriorou,logpdfou)
  }
  
  logpdf["logposterior"] <- sum(logpdf[-1])
  
  if ( verbose )
  {
    par.all <- numeric(0)
    if ( length(ind.constpar) > 0 ) par.all <- unlist(param$param[ind.constpar])
    par.all <- c(par.all,param$param.ou.estim,param$param.ou.fixed)
    print(par.all)
    print(logpdf)
  }
  
  return(logpdf)
}


