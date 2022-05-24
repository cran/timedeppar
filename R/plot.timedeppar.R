#' Plot results of time-dependent parameter estimation
#' 
#' This function plot Markov chains and marginal densities of constant parameters, distributions of
#' time dependent parameters, and Markov chains and marginal densities of time-dependent parameters at
#' selected points in time.
#' 
#' @param  x                  results from the function \code{\link{infer.timedeppar}} of class \code{timedeppar} or
#'                            list of such results for comparing multiple chains (in the latter case you have to call
#'                            explicitly \code{plot.timedeppar} rather than being able to do the generic call \code{plot}
#'                            as the list of results is not an object of class \code{timedeppar}).
#' @param  type               vector of plot types:\cr
#'                            \code{"traces"} or \code{"marginals"}: traces and 1d marginals of Markov chains of 
#'                                              constant parameters, of time-dependent parameters at certain time points
#'                                              (see argument \code{chains.at}), chains of log posterior and
#'                                              log observational likelihood values. 
#'                                              For selected outputs only, specify
#'                                              \code{traces.constpar}, \code{traces.timedeppar}, 
#'                                              \code{traces.logposterior}.\cr
#'                            \code{"summary"}: print summary of acceptance rates and maximum log posterior and 
#'                                              log likelihood values.\cr
#'                            \code{"pairs"}: scatterplot matrix of posterior sample of constant parameters.\cr
#'                            \code{"time-series"}: uncertainty range and median time series of
#'                                                  time-dependent parameters.\cr
#'                            \code{"accept"}: time series of apparent acceptance frequencies (at the level of thinning).\cr
#'                            \code{"realizations"}: realizations of time-dependent parameters to check for burnin..\cr
#'                            \code{"diagnostics"}: plot diagnostics for inference of time-dependent parameters
#'                                                  (note that the plot file could become very large 
#'                                                  to follow the inference steps).
#' @param  chains.at          vector of time points at which chains and marginals of time-dependent parameters 
#'                            should be plotted if \code{"traces"} or \code{"marginals"} is contained in the vector argument
#'                            \code{type} (default: none [numeric(0)]).
#' @param  labels             optional named vector of expressions to label variables in the plots (names of the expression
#'                            have to correspond to the variable names as used by the program, expressions can have special
#'                            symbols, e.g. \code{expression(a=alpha,b=beta,c1=gamma[1])}).
#' @param  units              optional named vector of expressions to add units to variables in the plots (names of the expression
#'                            have to correspond to the variable names as used by the program, expressions can have special
#'                            symbols, e.g. \code{expression(a=m^3/s,b=h^-1,c1=m)}).
#' @param  prob.band          probability defining the width of the uncertainty bands plotted for output variables
#'                            (default value: 0.9)
#' @param  max.diag.plots     maximum number of diagnostic plots of inference steps
#' @param  xlim.ts            optional range of time values for time-series plot
#' @param  n.burnin           number of Markov chain points to omit for density and pairs plots
#'                            (number of omitted points is max(control$n.adapt,n.burnin)).
#' @param  nrow               number of plot rows per page (except for pairs plot).
#' @param  nrow.constpar      number of plot rows per page for traces and marginals (default is nrow).
#' @param  nrow.timedeppar    number of plot rows per page for time-dependent parameters (default is nrow).
#' @param  nrow.diagnostics   number of plot rows per page for diagnostics plots (default is nrow).
#' @param  ...                additional arguments passed to the plotting function.


plot.timedeppar <- function(x,type=c("traces","marginals","summary","pairs","time-series","accept"),
                            chains.at=numeric(0),labels=NA,units=NA,
                            prob.band=0.9,max.diag.plots=100,xlim.ts=numeric(0),
                            n.burnin=0,
                            nrow=4,nrow.constpar=NA,nrow.timedeppar=NA,nrow.diagnostics=NA,...)
{
  message(paste("plot.timedeppar (timedeppar ",
                as.character(packageVersion("timedeppar"))," ",
                as.character(packageDate("timedeppar")),"): plotting inference results: ",
                paste(type,collapse=","),sep=""))
  
  # save and restore graphics parameters:
  # -------------------------------------
  
  par.old <- par(no.readonly=TRUE)
  on.exit(par(par.old))
  
  # preliminary checks:
  # -------------------
  
  res <- x
  if ( ! inherits(res[[1]],"timedeppar") ) res <- list(res)

  # check input:
  
  for ( i in 1:length(res) )
  {
    if ( ! inherits(res[[i]],"timedeppar") )
    {
      warning("first argument must be of class timedeppar or a list of objects of class timedeppar")
      return()
    }
    if ( !is.null(res[[i]]$errmsg) )
    {
      warning("there were errors when producing the first argument")
      for ( j in length(res[[i]]$errmsg) ) warning(res[[i]]$errmsg[j])
      return()
    }
  }
  
  if ( prob.band <= 0 | prob.band >= 1  )
  {
    warning("illegal value of argument prob.band (needs to be between 0 and 1)")
    return()
  }
  
  # plot parameters for traces and marginals:
  # -----------------------------------------
  
  mgp <- c(1.8,0.6,0)          # distance of axes labels and numbers
  mar <- c(2.9,2.9,0.4,1.0)    # c(bottom, left, top, right)
  
  # variables of general use:
  # -------------------------
  
  n.res          <- length(res)    # number if results from infer.timedeppar contained in list
  thin           <- rep(1,n.res)   # thin parameter for all results
  n.adapt        <- rep(0,n.res)   # length of adaptation phase for all results
  ind.sample.all <- list()         # sample indices available for all results (for plotting chains with adaptationn and burnin)
  ind.sample.use <- list()         # sample indices to be used for analyses (after adaptation and burnin)
  iter.max       <- 0              # maximum number of iterations available across all results
  for ( i in 1:n.res )
  {
    if ( !is.null(res[[i]]$control) ) { if(!is.null(res[[i]]$control$thin))    thin[i]    <- res[[i]]$control$thin } 
    if ( !is.null(res[[i]]$control) ) { if(!is.null(res[[i]]$control$n.adapt)) n.adapt[i] <- res[[i]]$control$n.adapt }
    if ( "sample.logpdf" %in% names(res[[i]]) )
    {
      iter.max <- max(iter.max,(nrow(res[[i]]$sample.logpdf)-1)*thin[i])
      ind.sample.all[[i]] <- 1:nrow(res[[i]]$sample.logpdf)
      ind.sample.use[[i]] <- ind.sample.all[[i]]
      if ( floor(max(n.adapt[i],n.burnin)/thin[i])+1 < nrow(res[[i]]$sample.logpdf) )
      {
        ind.sample.use[[i]] <- (floor(max(n.adapt[i],n.burnin)/thin[i])+1):nrow(res[[i]]$sample.logpdf)
      }
      else
      {
        warning("analysis is done using data from the adaptation or burnin phase")
      }
    }
    else
    {
      warning(paste("no sample present in input",i))
      return()
    }
  }
  
  # plot traces and marginals of constant parameters:
  # -------------------------------------------------
  
  if ( "traces" %in% type | "marginals" %in% type | "traces.constpar" %in% type | "marginals.constpar" %in% type )
  {
    # check presence of constant parameters:
    
    if ( is.null(res[[1]]$sample.param.const) & is.null(res[[1]]$sample.param.ou) )
    {
      warning("no sample of constant parameters found to plot traces and marginals")
    }
    else
    {
      if ( ncol(res[[1]]$sample.param.const)+ncol(res[[1]]$sample.param.ou) == 0 )
      {
        warning("no sample of constant parameters found to plot traces and marginals")
      }
      else
      {
        # collect data:
    
        sample.constpar.all <- list()
        var.names <- character(0)
        for ( i in 1:n.res )
        {
          sample.constpar.all[[i]] <- cbind(res[[i]]$sample.param.const,res[[i]]$sample.param.ou)
          var.names <- c(var.names,colnames(sample.constpar.all[[i]]))
        }
        var.names <- unique(var.names)
    
        # plot traces and marginals of constant parameters:

        par.def <- par(no.readonly=TRUE)
        par(mfrow=c(ifelse(is.na(nrow.constpar),nrow,nrow.constpar),2),mgp=mgp,mar=mar)
        for ( var.name in var.names )
        {
          # get var.lim:
          var.lim <- c(NA,NA)
          for ( i in 1:n.res )
          {
            if ( var.name %in% colnames(sample.constpar.all[[i]]) )
            {
              var.lim[1] <- min(var.lim[1],min(sample.constpar.all[[i]][ind.sample.use[[i]],var.name],na.rm=T),na.rm=T)
              var.lim[2] <- max(var.lim[2],max(sample.constpar.all[[i]][ind.sample.use[[i]],var.name],na.rm=T),na.rm=T)
            }
          }
              
          # plot traces:
          plot(numeric(0),numeric(0),type="n",xlim=c(0,iter.max),ylim=var.lim,xaxs="i",yaxs="i",
               xlab="iterations",ylab=get.label(var.name,labels,units))
          for ( i in 1:n.res )
          {
            if ( var.name %in% colnames(sample.constpar.all[[i]]))
            {
              lines((ind.sample.all[[i]]-1)*thin[i],sample.constpar.all[[i]][,var.name],col=i)
            }
            abline(v=n.adapt[i],col="red",lwd=2)
          }
          abline(v=n.burnin,col="blue",lty="dashed",lwd=2)
          
          # plot densities:
          densities <- list()
          max.dens  <- NA
          for ( i in 1:n.res )
          {
            if ( var.name %in% colnames(sample.constpar.all[[i]]) )
            {
              densities[[i]] <- density(sample.constpar.all[[i]][ind.sample.use[[i]],var.name])
              max.dens <- max(max.dens,densities[[i]]$y,na.rm=TRUE)
            }
            else
            {
              densities[[i]] <- NA
            }
          }
          plot(numeric(0),numeric(0),type="n",xlim=var.lim,ylim=c(0,1.1*max.dens),yaxs="i",yaxs="i",
               xlab=get.label(var.name,labels,units),ylab="density")
          for ( i in 1:n.res )
          {
            if ( var.name %in% colnames(sample.constpar.all[[i]]))
            {
              lines(densities[[i]],col=i)
            }
          }
        }
        par(par.def)
      }
    }
  }
    
  # plot traces and marginals of time-dependent parameters at fixed time points:
  # ----------------------------------------------------------------------------
    
  if ( "traces" %in% type | "marginals" %in% type | "traces.timedeppar" %in% type | "marginals.timedeppar" %in% type )
  {
    # check presence of time-dependent parameters:
    
    if ( length(chains.at) > 0 )
    {
      if ( is.null(res[[1]]$sample.param.timedep) )
      {
        warning("no sample of time-dependent parameters found to plot traces and marginals")
      }
      else
      {
        if ( length(res[[1]]$sample.param.timedep) == 0 )
        {
          warning("no sample of time-dependent parameters found to plot traces and marginals")
        }
        else
        {
          par.def <- par(no.readonly=TRUE)
          par(mfrow=c(ifelse(is.na(nrow.constpar),nrow,nrow.constpar),2),mgp=mgp,mar=mar)

          # collect iterations, indices and variable names:
          
          var.names <- character(0)
          for ( i in 1:n.res )
          {
            var.names <- c(var.names,names(res[[i]]$sample.param.timedep))
          }
          var.names <- unique(var.names)

          # plot traces and marginals of time-dependent parameters:
          
          for ( var.name in var.names )
          {
            sample.param.timedep.interpol <- list()
            for ( i in 1:n.res )
            {
              if ( var.name %in% names(res[[i]]$sample.param.timedep) )
              {
                sample.param.timedep.interpol[[i]] <-
                  matrix(NA,nrow=nrow(res[[i]]$sample.param.timedep[[var.name]])-1,ncol=length(chains.at))
                for ( k in 1:(nrow(res[[i]]$sample.param.timedep[[var.name]])-1) )
                  sample.param.timedep.interpol[[i]][k,] <- approx(x    = res[[i]]$sample.param.timedep[[var.name]][1,],
                                                                   y    = res[[i]]$sample.param.timedep[[var.name]][1+k,],
                                                                   xout = chains.at)$y
              }
              else
              {
                sample.param.timedep.interpol[[i]] <- NA
              }
            }
            for ( j in 1:length(chains.at) )
            {
              # get bounds:
              var.lim <- c(NA,NA)
              for ( i in 1:n.res )
              {
                if ( var.name %in% names(res[[i]]$sample.param.timedep) )
                {
                  var.lim[1] <- min(var.lim[1],min(sample.param.timedep.interpol[[i]][ind.sample.use[[i]],j],na.rm=T),na.rm=T)
                  var.lim[2] <- max(var.lim[2],max(sample.param.timedep.interpol[[i]][ind.sample.use[[i]],j],na.rm=T),na.rm=T)
                }
              }

              # plot traces:
              plot(numeric(0),numeric(0),type="n",xaxs="i",yaxs="i",xlim=c(0,iter.max),ylim=var.lim,
                   xlab="iterations",ylab=get.label(var.name,labels,units,t2=paste(" at t=",signif(chains.at[j],3)," ",sep="")))
              for ( i in 1:n.res )
              {
                if ( var.name %in% names(res[[i]]$sample.param.timedep) )
                {
                  lines((ind.sample.all[[i]]-1)*thin[i],
                        sample.param.timedep.interpol[[i]][ind.sample.all[[i]],j],
                        col=i)
                }
                abline(v=n.adapt[i],col="red",lwd=2)
              }
              abline(v=n.burnin,col="blue",lty="dashed",lwd=2)

              # plot densities:
              densities <- list()
              max.dens <- NA
              for ( i in 1:n.res )
              {
                if ( var.name %in% names(res[[i]]$sample.param.timedep) )
                {
                  densities[[i]] <- density(sample.param.timedep.interpol[[i]][ind.sample.use[[i]],j])
                  max.dens <- max(max.dens,densities[[i]]$y,na.rm=TRUE)
                }
                else
                {
                  densities[[i]] <- NA
                }
              }
              plot(numeric(0),numeric(0),type="n",xlim=var.lim,ylim=c(0,1.1*max.dens),yaxs="i",yaxs="i",
                   xlab=get.label(var.name,labels,units,t2=paste(" at t=",signif(chains.at[j],3)," ",sep="")),ylab="density")
              for ( i in 1:n.res )
              {
                if ( var.name %in% names(res[[i]]$sample.param.timedep) )
                {
                  lines(densities[[i]],col=i)
                }
              }
            }
          }
          par(par.def)
        }
      }
    }
  }

  # plot traces of log posterior and log observational likelihood:
  # --------------------------------------------------------------
  
  if ( "traces" %in% type | "marginals" %in% type | "traces.logposterior" %in% type | "marginals.logposterior" %in% type )
  {
    # check presence of logpdf:
    
    if ( !"sample.logpdf" %in% names(res[[1]]) )
    {
      warning("no sample for log pdf values found to plot traces")
    }
    else
    {
      par.def <- par(no.readonly=TRUE)
      par(mfrow=c(ifelse(is.na(nrow.constpar),nrow,nrow.constpar),2),mgp=mgp,mar=mar)
      
      # plot traces of log posterior:
      
      ylim <- c(NA,NA)
      for ( i in 1:n.res )
      {
        if ( "sample.logpdf" %in% names(res[[i]]) )
        {
          ylim[1] <- min(ylim[1],min(res[[i]]$sample.logpdf[ind.sample.use[[i]],"logposterior"],na.rm=T),na.rm=T)
          ylim[2] <- max(ylim[2],max(res[[i]]$sample.logpdf[ind.sample.use[[i]],"logposterior"],na.rm=T),na.rm=T)
        }
      }
      plot(numeric(0),numeric(0),type="n",xaxs="i",yaxs="i",
           xlab="iterations",ylab="log posterior",xlim=c(0,iter.max),ylim=ylim)
      for ( i in 1:n.res )
      {
        if ( "sample.logpdf" %in% names(res[[i]]) )
        {
          lines((ind.sample.all[[i]]-1)*thin[i],res[[i]]$sample.logpdf[,"logposterior"],col=i)
        }
        abline(v=n.adapt[i],col="red",lwd=2)
      }
      abline(v=n.burnin,col="blue",lty="dashed",lwd=2)
      
      # plot traces of log likeliobs:
      
      ylim <- c(NA,NA)
      for ( i in 1:n.res )
      {
        if ( "sample.logpdf" %in% names(res[[i]]) )
        {
          ylim[1] <- min(ylim[1],min(res[[i]]$sample.logpdf[ind.sample.use[[i]],"loglikeliobs"],na.rm=T),na.rm=T)
          ylim[2] <- max(ylim[2],max(res[[i]]$sample.logpdf[ind.sample.use[[i]],"loglikeliobs"],na.rm=T),na.rm=T)
        }
      }
      plot(numeric(0),numeric(0),type="n",xaxs="i",yaxs="i",
           xlab="iterations",ylab="log likeliobs",xlim=c(0,iter.max),ylim=ylim)
      for ( i in 1:n.res )
      {
        if ( "sample.logpdf" %in% names(res[[i]]) )
        {
          lines((ind.sample.all[[i]]-1)*thin[i],res[[i]]$sample.logpdf[,"loglikeliobs"],col=i)
        }
        abline(v=n.adapt[i],col=i,lwd=2)
      }
      abline(v=n.burnin,col="blue",lty="dashed",lwd=2)
      
      par(par.def)
    }
  }
  
  # plot summary of acceptance rates and max. log posterior and log likelihood values:
  # ----------------------------------------------------------------------------------

  if ( "summary" %in% type )
  {
    par.def <- par(no.readonly=TRUE)
    par(mfrow=c(ifelse(is.na(nrow.constpar),nrow,nrow.constpar),2),mgp=mgp,mar=mar)
    for ( i in 1:n.res )
    {
      if ( !is.null(res[[i]]$acceptfreq.constpar) & !is.null(res[[i]]$acceptfreq.timedeppar) & !is.null(res[[i]]$acceptfreq.oupar) )
      {
        cex.text <- 0.8
        
        x1 <- 0.1
        x2 <- 0.5
        x3 <- 0.8
        y1 <- 0.8
        y2 <- 0.6
        y3 <- 0.4
        y4 <- 0.2
        
        plot(numeric(0),numeric(0),type="n",xlim=c(0,1),ylim=c(0,1),
             xlab="",ylab="",xaxt="n",yaxt="n")
        if ( n.res > 1 ) text(x=x1,y=y1,labels=paste("chain",i),pos=4,cex=cex.text,col=i)
        text(x=x2,y=y1,labels="num. par.:",cex=cex.text)
        text(x=x3,y=y1,labels="acc. freq. (%):",cex=cex.text)
        text(x=x1,y=y2,labels="const. par:",pos=4,cex=cex.text)
        text(x=x2,y=y2,labels=ncol(res[[i]]$sample.param.const),cex=cex.text)
        if ( !is.na(res[[i]]$acceptfreq.constpar) ) text(x=x3,y=y2,labels=round(100*res[[i]]$acceptfreq.constpar,1),cex=cex.text)
        text(x=x1,y=y3,labels="timedep. par:",pos=4,cex=cex.text)
        text(x=x2,y=y3,labels=length(res[[i]]$sample.param.timedep),cex=cex.text)
        if ( length(res[[i]]$acceptfreq.timedeppar) > 0 )
        {
          if ( !is.na(res[[i]]$acceptfreq.timedeppar[1]) )
          {
            text(x=x3,y=y3,labels=paste(round(100*res[[i]]$acceptfreq.timedeppar,1),collapse=", "),cex=cex.text)
          }
        }
        text(x=x1,y=y4,labels="proc. par:",pos=4,cex=cex.text)
        text(x=x2,y=y4,labels=ncol(res[[i]]$sample.param.ou),cex=cex.text)
        if ( length(res[[i]]$acceptfreq.oupar) > 0 )
        {
          if ( !is.na(res[[i]]$acceptfreq.oupar[1]) )
          {
            text(x=x3,y=y4,labels=paste(round(100*res[[i]]$acceptfreq.oupar,1),collapse=", "),cex=cex.text)
          }
        }
        
        ind.maxpost <- which.max(res[[i]]$sample.logpdf[,"logposterior"])
        max.logpost <- res[[i]]$sample.logpdf[ind.maxpost,"logposterior"]
        ind.maxlikeli <- which.max(res[[i]]$sample.logpdf[,"loglikeliobs"])
        max.loglikeli <- res[[i]]$sample.logpdf[ind.maxlikeli,"loglikeliobs"]
        plot(numeric(0),numeric(0),type="n",xlim=c(0,1),ylim=c(0,1),
             xlab="",ylab="",xaxt="n",yaxt="n")
        text(x=x1,y=y1,labels=paste("max log posterior =",signif(max.logpost,6)),pos=4,cex=cex.text)
        text(x=x1,y=y2,labels=paste("   at iteration",(ind.maxpost-1)*thin),pos=4,cex=cex.text)
        text(x=x1,y=y3,labels=paste("max log likelihood =",signif(max.loglikeli,6)),pos=4,cex=cex.text)
        text(x=x1,y=y4,labels=paste("   at iteration",(ind.maxlikeli-1)*thin),pos=4,cex=cex.text)
        
      }
    }
    par(par.def)
  }
      
  # plot scatterplot matrix of sample of constant parameters:
  # ---------------------------------------------------------
  
  if ( "pairs" %in% type )
  {
    # check presence of constant parameters:
    
    if ( is.null(res[[1]]$sample.param.const) & is.null(res[[1]]$sample.param.ou) )
    {
      warning("no sample of constant parameters found to plot 2d marginals")
    }
    else
    {
      if ( ncol(res[[1]]$sample.param.const)+ncol(res[[1]]$sample.param.ou) == 0 )
      {
        warning("no sample of constant parameters found to plot 2d marginals")
      }
      else
      {
        # collect data:

        sample.constpar.all <- cbind(res[[1]]$sample.param.const,res[[1]]$sample.param.ou)[ind.sample.use[[1]],,drop=FALSE]
        if ( ncol(sample.constpar.all) < 2 )
        {
          warning("plot of 2d marginals not possible for a single parameter")
        }
        else
        {
          col <- rep(1,nrow(sample.constpar.all))
          if ( n.res > 1 )
          {
            for ( i in 2:n.res )
            {
              sample.constpar.all <- 
                rbind(sample.constpar.all,
                      cbind(res[[i]]$sample.param.const,res[[i]]$sample.param.ou)[ind.sample.use[[i]],])
              col <- c(col,rep(i,length(ind.sample.use[[i]])))
            }
          }
          
          panel.cor <- function(x,y,...)
          {
            ind <- which(!is.na(x) & !is.na(y))
            corr <- cor(x[ind],y[ind],method="spearman")
            text(x=0.5*(min(x,na.rm=T)+max(x,na.rm=T)),y=0.5*(min(y,na.rm=T)+max(y,na.rm=T)),label=round(corr,2),cex=2)
          }
          panel.dat <- function(x,y,...)
          {
            points(x,y,...)
          }
          plot.timedeppar.labels.global <- labels
          plot.timedeppar.units.global  <- units
          panel.diag <- function(x,y,labels,cex,font,...) 
          {
            text(x=mean(range(x)),y=mean(range(y)),
                 labels=get.label(var.name=labels,labels=plot.timedeppar.labels.global,units=plot.timedeppar.units.global),
                 cex=2)
          }
          
          par.def <- par(no.readonly=TRUE)
          par(mfrow=c(1,1))
          pairs(sample.constpar.all,lower.panel=panel.cor,upper.panel=panel.dat,text.panel=panel.diag,
                pch=19,cex=0.3,col=col)
          par(par.def)
        }
      }
    }
  }
        
  # plot results for time-dependent parameters:
  # -------------------------------------------
  
  if ( "time-series" %in% type )
  {
    if ( is.null(res[[1]]$sample.param.timedep) )
    {
      warning("no sample of time-dependent parameters found to plot traces and marginals")
    }
    else
    {
      if ( length(res[[1]]$sample.param.timedep) == 0 )
      {
        warning("no sample of time-dependent parameters found to plot traces and marginals")
      }
      else
      {
        par.def <- par(no.readonly=TRUE)
        par(mfrow=c(ifelse(is.na(nrow.timedeppar),nrow,nrow.timedeppar),1),mgp=mgp,mar=mar)
        
        # collect iterations, indices and variable names:
        
        var.names <- character(0)
        for ( i in 1:n.res )
        {
          var.names <- c(var.names,names(res[[i]]$sample.param.timedep))
        }
        var.names <- unique(var.names)
        
        # plot traces and marginals of time-dependent parameters:
        
        for ( var.name in var.names )
        {
          quantiles <- list()
          xlim      <- c(NA,NA)
          ylim      <- c(NA,NA)
          for ( i in 1:n.res )
          {
            if ( var.name %in% names(res[[i]]$sample.param.timedep) )
            {
              quantiles[[i]]       <- matrix(NA,nrow=4,ncol=ncol(res[[i]]$sample.param.timedep[[var.name]]))
              quantiles[[i]][1,]   <- res[[i]]$sample.param.timedep[[var.name]][1,]
              # for ( j in 1:ncol(res[[i]]$sample.param.timedep[[var.name]]) )
              # {
              #   quantiles[[i]][2:4,j] <- quantile(res[[i]]$sample.param.timedep[[var.name]][ind.sample.use[[i]]+1,j],
              #                                     probs=c(0.5*(1-prob.band),0.5,1-0.5*(1-prob.band)))
              # }
              quantiles[[i]][2:4,] <- apply(res[[i]]$sample.param.timedep[[var.name]][ind.sample.use[[i]]+1,],2,quantile,
                                            probs=c(0.5*(1-prob.band),0.5,1-0.5*(1-prob.band)))
              xlim[1] <- min(xlim[1],min(quantiles[[i]][1,],na.rm=T),na.rm=T)
              xlim[2] <- max(xlim[2],max(quantiles[[i]][1,],na.rm=T),na.rm=T)
              ylim[1] <- min(ylim[1],min(quantiles[[i]][2,],na.rm=T),na.rm=T)
              ylim[2] <- max(ylim[2],max(quantiles[[i]][4,],na.rm=T),na.rm=T)
            }
            else
            {
              quantiles[[i]] <- NA
            }
          }
          plot(numeric(0),numeric(0),type="n",xaxs="i",yaxs="i",xlim=xlim,ylim=ylim,
               xlab=get.label("time",labels,units),ylab=get.label(var.name,labels,units))
          for ( i in 1:n.res )
          {
            if ( var.name %in% names(res[[i]]$sample.param.timedep) )
            {
              polygon(c(quantiles[[i]][1,],rev(quantiles[[i]][1,]),quantiles[[i]][1,1]),
                      c(quantiles[[i]][2,],rev(quantiles[[i]][4,]),quantiles[[i]][2,1]),
                      col=adjustcolor(i,alpha.f=0.3),border=NA)
            }
          }
          for ( i in 1:n.res )
          {
            if ( var.name %in% names(res[[i]]$sample.param.timedep) )
            {
              lines(quantiles[[i]][1,],quantiles[[i]][3,],col=i,lwd=1)
            }
          }
        }
        par(par.def)
      }
    }
  }

  # plot apparent acceptance frequencies:
  # -------------------------------------
  
  if ( "accept" %in% type | "diagnostics" %in% type )
  {
    par.def <- par(no.readonly=TRUE)
    par(mfrow=c(ifelse(is.na(nrow.diagnostics),nrow,nrow.constpar),1),mgp=mgp,mar=mar)
    for ( i in 1:n.res )
    {
      if ( "sample.param.timedep" %in% names(res[[i]]) )
      {
        if ( length(res[[i]]$sample.param.timedep) > 0 )
        {
          acceptfreq <- calc.acceptfreq(res[[i]],n.burnin=n.burnin)
          for ( j in 1:length(res[[i]]$sample.param.timedep) )
          {
            t <- res[[i]]$sample.param.timedep[[j]][1,]
            xlim <- range(acceptfreq[[j]][,1]); if ( length(xlim.ts)==2 ) xlim <- xlim.ts
            plot(numeric(0),numeric(0),xlim=xlim,ylim=c(0,100),xaxs="i",yaxs="i",
                 #main=paste("Apparent accept. freq. of sample (thin=",thin[i],") of ",names(res[[i]]$sample.param.timedep)[j],sep=""),
                 xlab=get.label("time",labels,units),ylab=get.label(names(res[[i]]$sample.param.timedep)[j],labels,units=NA,t1="acc. freq. ",t2=" [%]"))
            lines(acceptfreq[[j]][,1],100*acceptfreq[[j]][,2])    
            if ( n.res > 1 ) text(x=xlim[2]-0.9*diff(xlim),y=80,label=paste("chain",i),col=i)
            
            if ( !is.null(res[[i]]$control$splitmethod) )
            {
              if ( res[[i]]$control$splitmethod == "weighted" | res[[i]]$control$splitmethod == "autoweights" )
              {
                weights <- numeric(0)
                if ( !is.null(res[[i]]$control$interval.weights) )
                {
                  if ( is.list(res[[i]]$control$interval.weights) )
                  {
                    weights <- res[[i]]$control$interval.weights[[names(res[[i]]$sample.param.timedep)[j]]]
                  }
                  else
                  {
                    weights <- res[[i]]$control$interval.weights
                  }
                }
                if ( length(weights) > 0 )
                {
                  lines(t,weights/max(weights)*100,lty="dashed",col="red")
                }
              }
            }
          }
        }
      } 
    }
    par(par.def)
  }
  
  # plot realizations of time series:
  # -------------------------------------
  
  if ( "realizations" %in% type | "diagnostics" %in% type )
  {
    par.def <- par(no.readonly=TRUE)
    par(mfrow=c(ifelse(is.na(nrow.diagnostics),nrow,nrow.constpar),1),mgp=mgp,mar=mar)
    for ( i in 1:n.res )
    {
      if ( "sample.param.timedep" %in% names(res[[i]]) )
      {
        if ( length(res[[i]]$sample.param.timedep) > 0 )
        {
          for ( j in 1:length(res[[i]]$sample.param.timedep) )
          {
            t <- res[[i]]$sample.param.timedep[[j]][1,]
            q <- matrix(NA,nrow=length(t),ncol=3)
            for ( k in 1:length(t) ) q[k,] <- quantile(res[[i]]$sample.param.timedep[[j]][ind.sample.use[[i]]+1,k],probs=c(0.5*(1-prob.band),0.5,1-0.5*(1-prob.band)))
            xlim <- range(t); if ( length(xlim.ts)==2 ) xlim <- xlim.ts
            ylim <- range(q)
            plot(numeric(0),numeric(0),xlim=xlim,ylim=ylim,xaxs="i",yaxs="i",
                 #main=paste("Realizations of parameter",names(res[[i]]$sample.param.timedep)[i]),
                 xlab=get.label("time",labels,units),ylab=get.label(names(res[[i]]$sample.param.timedep)[j],labels,units))
            if ( n.res > 1 ) text(x=xlim[2]-0.9*diff(xlim),y=ylim[2]-0.2*diff(ylim),label=paste("chain",i),col=i)
            n.curves <- 11
            #cols <- rainbow(n.curves)
            #cols <- heat.colors(n.curves,rev=TRUE)
            cols <- topo.colors(n.curves,rev=TRUE)
            #cols <- terrain.colors(n.curves)
            indices <- ind.sample.use[[i]][round((length(ind.sample.use[[i]])-1)*(0:(n.curves-1))/(n.curves-1)) + 1]
            for ( k in 1:n.curves )
            {
              lines(t,res[[i]]$sample.param.timedep[[j]][indices[k]+1,],lwd=1,col=cols[k])
            }
            legend("topright",paste("iter =",thin[i]*(indices-1)),col=cols,lwd=1,cex=0.7)
            lines(t,q[,2],lwd=1)
            lines(t,q[,1],lwd=1,lty="dashed")
            lines(t,q[,3],lwd=1,lty="dashed")
          }
        }
      }
    }
    par(par.def)
  }

  # plot results for diagnostics:
  # -----------------------------
  
  if ( "diagnostics" %in% type )
  {
    par.def <- par(no.readonly=TRUE)
    par(mfrow=c(ifelse(is.na(nrow.diagnostics),nrow,nrow.constpar),1))
    for ( i in 1:n.res )
    {
      if ( "sample.diag" %in% names(res[[i]]) )
      {
        if ( length(res[[i]]$sample.diag) > 0 )
        {
          ind.range.diag <- 1:(nrow(res[[i]]$sample.param.timedep[[1]])-2)
          if ( floor(max(n.adapt[i],n.burnin)/thin[i])+1 < nrow(res[[i]]$sample.param.timedep[[1]])-1 ) ind.range.diag <- (floor(max(n.adapt[i],n.burnin)/thin[i])+1):(nrow(res[[i]]$sample.param.timedep[[1]])-2)
          else warning("diagnostics plot is using data from the adaptation or burnin phase")
          par.def <- par(no.readonly=TRUE)
          par(mfrow=c(ifelse(is.na(nrow.timedeppar),nrow,nrow.timedeppar),1))
          
          # plot interval lengths:
          
          for ( j in 1:length(res[[i]]$sample.diag) )
          {
            t <- res[[i]]$sample.param.timedep[[j]][1,]
            q <- matrix(NA,nrow=length(t),ncol=3)
            for ( k in 1:length(t) ) q[k,] <- quantile(res[[i]]$sample.diag[[j]]$internalpts[ind.range.diag,k],
                                                       probs=c(0.5*(1-prob.band),0.5,1-0.5*(1-prob.band)),na.rm=TRUE)
            xlim <- range(t); if ( length(xlim.ts)==2 ) xlim <- xlim.ts
            ylim <- c(0,max(q))
            plot(numeric(0),numeric(0),xlim=xlim,ylim=ylim,xaxs="i",yaxs="i",
                 main=get.label(names(res[[i]]$sample.param.timedep)[j],labels,
                                t1=paste("Median and ",100*prob.band,"% band of internal points for ",sep="")),
                 xlab=get.label("time",labels,units),ylab="internal interval points")
            if ( n.res > 1 ) text(x=xlim[2]-0.9*diff(xlim),y=ylim[2]-0.2*diff(ylim),label=paste("chain",i),col=i)
            polygon(c(t,rev(t),t[1]),c(q[,1],rev(q[,3]),q[1,1]),col=grey(0.8),border=NA)
            lines(t,q[,2],lwd=1)
          }
          
          # plot inference steps:
          
          ind.diag <- ind.range.diag; if ( length(ind.diag) > max.diag.plots ) ind.diag <- ind.range.diag[round((length(ind.range.diag)-1)/(max.diag.plots-1)*0:(max.diag.plots-1))+1]
          for ( k in ind.diag )
          {
            for ( j in 1:length(res[[i]]$sample.diag) )
            {
              t <- res[[i]]$sample.param.timedep[[j]][1,]
              xlim <- range(t); if ( length(xlim.ts)==2 ) xlim <- xlim.ts
              ylim <- range(res[[i]]$sample.param.timedep[[j]][ind.diag+2,])
              ylogacceptratio <- ylim[2] + 0.05*(ylim[2]-ylim[1])
              dylogacceptratio <-  0.04*(ylim[2]-ylim[1])
              ylim[2] <- ylim[2] + 0.15*(ylim[2]-ylim[1])
              plot(numeric(0),numeric(0),xlim=xlim,ylim=ylim,xaxs="i",yaxs="i",
                   main=get.label(names(res[[i]]$sample.param.timedep)[j],labels,t1="Parameter ",t3=paste(", iteration ",thin[i]*k,sep="")),
                   xlab=get.label("time",labels,units),ylab=get.label(names(res[[i]]$sample.param.timedep)[j],labels,units))
              if ( n.res > 1 ) text(x=xlim[2]-0.9*diff(xlim),y=ylim[2]-0.2*diff(ylim),label=paste("chain",i),col=i)
              abline(h=ylogacceptratio+dylogacceptratio*(-2:2),lty=c("dotted","dashed","solid","dashed","dotted"))
              ind.nodes <- which(is.na(res[[i]]$sample.diag[[j]]$accepted[k,]))
              for ( kk in ind.nodes ) abline(v=t[kk])
              kk0 <- 0
              for ( kk in c(ind.nodes,length(t)+1) )
              {
                lines(0.5*(t[kk-1]+t[kk0+1])*c(1,1),
                      c(ylogacceptratio,ylogacceptratio-res[[i]]$sample.diag[[j]]$logacceptratio[k,kk-1]/log(10)*dylogacceptratio),
                      col=ifelse(res[[i]]$sample.diag[[j]]$accepted[k,kk-1],"green","red"),lwd=2)
                kk0 <- kk
              }
              lines(t,res[[i]]$sample.param.timedep[[j]][k+1,])
              col.accept <- rep("red",length(t))
              for ( kk in 1:(length(ind.nodes)+1) )
              {
                ind <- c(1,ind.nodes,length(t))[kk]:c(1,ind.nodes,length(t))[kk+1]
                lines(t[ind],res[[i]]$sample.diag[[j]]$proposal[k,ind],
                      col=ifelse(any(res[[i]]$sample.diag[[j]]$accepted[k,ind],na.rm=TRUE),"green","red"))
              }
              cex.r <- 0.75
              mtext("  r" ,side=4,line=0.2,las=1,at=ylogacceptratio+3*dylogacceptratio,cex=cex.r)
              mtext("0.01",side=4,line=0.2,las=1,at=ylogacceptratio+2*dylogacceptratio,cex=cex.r)
              mtext(" 0.1",side=4,line=0.2,las=1,at=ylogacceptratio+1*dylogacceptratio,cex=cex.r)
              mtext("  1" ,side=4,line=0.2,las=1,at=ylogacceptratio                   ,cex=cex.r)
              mtext(" 10" ,side=4,line=0.2,las=1,at=ylogacceptratio-1*dylogacceptratio,cex=cex.r)
              mtext("100" ,side=4,line=0.2,las=1,at=ylogacceptratio-2*dylogacceptratio,cex=cex.r)
              # legend("bottomright",
              #        legend = c("accepted","rejected","r = 0","r = 0.1, 10","r = 0.01, 100"),
              #        col    = c("green","red","black","black","black"),
              #        lty    = c("solid","solid","solid","dashed","dotted"))
            }
          }
        }
      }
    }
    par(par.def)
  }
}


