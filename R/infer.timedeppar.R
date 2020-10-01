#' Jointly infer constant and time-dependent parameters of a dynamic model given time-series data
#' 
#' This function draws a Markov Chain from the posterior of constant and time-dependent parameters
#' (following Ornstein-Uhlenbeck processes) of a dynamic model.
#' The dynamic model is specified by a function that calculates the log likelihood for given time-series 
#' data.
#' The Ornstein-Uhlenbeck processes of time-dependent processes are characterized by there mean (\code{mean}),
#' standard deviation (\code{sd}), and a rate parameter (\code{gamma}) that quantifies temporal correlation.
#' 
#' @param  loglikeli            function that calculates the log likelihood of the model for given constant or
#'                              time-dependent parameters and given observational data.\cr
#'                              The parameters are passed as a named list in the first argument of the function.
#'                              The list elements are either scalar values representing constant parameters
#'                              or two-column matrices with columns for time points and values for time-dependent 
#'                              parameters.
#'                              If the argument \code{loglikeli.keepstate} is \code{FALSE} no further arguments
#'                              are needed (but can be provided, see below).
#'                              In this case, the function should return the log likelihood as a single double
#'                              value.\cr
#'                              If the argument \code{loglikeli.keepstate} is \code{TRUE}, the second argument
#'                              provides the time range over which a time-dependent parameter was changed or NA
#'                              if the full simulation time has to be evaluated, and the third argument
#'                              provides the state of the functon at the last successful call. 
#'                              This allows the function to only calculate and return modifications to that
#'                              previous state.
#'                              In this case, the function has to return a list with the log likelihood value
#'                              as its first element and the current state of the function as the second argument.
#'                              This state can be an R variable of an arbitrary data type.
#'                              The version from the last accepted MCMC step will be returnde at the next call.\cr
#'                              Further arguments provided to \code{infer.timedeppar} will be passed to this function.
#' @param  loglikeli.keepstate  boolean to indicate which kind of interface to the likelihood function is used.
#'                              See argument \code{loglikeli} for details.
#' @param  param.ini            named list of initial vallues of parameters to be estimated.
#'                              scalar initial values for constant parameters,
#'                              two-column matrices for time and parameter values for time-dependent parameters
#'                              (values of time-dependent parameters may be NA and are then drawn from the 
#'                              Ornstein-Uhlenbeck process).
#'                              The list \code{param.ini} needs to be a legal and complete first element of the 
#'                              function passed by the argument \code{loglikeli}.
#'                              For each time-dependent parameter with name \code{<name>} 
#'                              initial values or fixed value of the parameters
#'                              \code{<name>_mean}, \code{<name>_sd} and \code{<name>_gamma} must be provided
#'                              in the arguments \code{param.ou.ini} or \code{param.ou.fixed}, respectively.
#'                              These parameters represent the mean, the asymptotic standard deviation, and 
#'                              the rate parameter of the Ornstein-Uhlenbeck process.
#'                              If these parameters are given in the argument \code{param.ou.ini}, 
#'                              they are used as initial condition of the inference process and the parameters 
#'                              are estimated, 
#'                              if they are given in the argument \code{param.ou.fixed},
#'                              they are assumed to be given and are kept fixed.
#' @param  param.range          named list of ranges (2 element vectors with minimum and maximum) of parameters 
#'                              that are constrained (non-logarithmic for all parameters)
#' @param  param.log            named vector of logicals indicating if inference should be done on the log scale 
#'                              (parameters are still given and returned on non-log scales).
#'                              For time-dependent parameters, selecting this option implies the use of a lognormal
#'                              marginal for the Ornstein-Uhlenbeck process.
#'                              This means that the parameter is modelled as exp(Ornstein-Uhlenbeck), but mean
#'                              and standard devistion of the process are still on non-log scales.
#' @param  param.logprior       function to calculate the (joint) log prior of all estimaged constant parameters 
#'                              of the model.
#'                              The function gets as its argument a named vector of the values of the estimaged  
#'                              constant parameters to allow the function to identify for which parameters a
#'                              joint prior is required in the current setting).
#' @param  param.ou.ini         named vector of initial values of parameters of the Ornstein-Uhlenbeck processes of
#'                              time-dependent parameters; see description of argument \code{param.ini}.
#' @param  param.ou.fixed       named vector of values of parameters of the Ornstein-Uhlenbeck processes of 
#'                              time-dependent parameters that are kept fixed rather than being extimated.
#'                              If all process parameters are kept fixed, these names are 
#'                              \code{<name>_mean}, \code{<name>_sd} and \code{<name>_gamma} for each
#'                              time-dependent parameter with name \code{<name>}; see description of the
#'                              argument \code{param.ini}.
#'                              The values specified in \code{param.fixed} are ignored if the parameter is also
#'                              given in the argument \code{param.ini}; in this case it is estimated.
#' @param  param.ou.logprior    function to calculate the (joint) log prior of all estimated parameters of
#'                              the Ornstein-Uhlenbeck processes of a single time-dependent parameter.
#'                              The function gets as its argument a named vector of the values of the process 
#'                              parameters to estimated.
#'                              These names are a subset of 
#'                              \code{<name>_mean}, \code{<name>_sd} and \code{<name>_gamma}; 
#'                              see description of the argument \code{param.ini}.
#'                              The function has to work for each time-dependent parameter by being sensitive
#'                              to the parameter names.
#' @param  task                 Which task to perform (default value: "start"):\cr
#'                              \code{"start"}: start an inference process from scratch based on the arguments of
#'                              the function.
#'                              The argument \code{res.infer.timedeppar} is ignored.\cr
#'                              \code{"continue"}: continue a Markov chain from a previous call to 
#'                              \code{\link{infer.timedeppar}}.
#'                              The results of a previous call have to be provided as the argument 
#'                              \code{res.infer.timedeppar} in the form of an object of type \code{timedeppar}.
#'                              To guarantee convergence of the chain, all numerical specifications including the
#'                              final state of the chain are taken from the object provided by 
#'                              \code{res.infer.timedeppar} and the actual arguments of the function are ignored 
#'                              except \code{n.iter} which specifies the number of iterations to be added 
#'                              to the chain.\cr
#'                              \code{"restart"}: A new chain is started from the last point of a previous chain
#'                              except if the argument \code{param.ini} is provided.
#'                              The results of a previous call have to be provided as the argument 
#'                              \code{res.infer.timedeppar} in the form of an object of type \code{timedeppar}.
#'                              Likelihood, prior pdf functions, initial proposal covariance matrices, 
#'                              scales and control parameters are taken from 
#'                              the previous chain unless explicitly provided.
#' @param  n.iter               number of iterations of the Markov chain to be performed (default value: 10000).
#' @param  cov.prop.const.ini   scaled covariance matrix of the proposal distribution for the Metropolis step of
#'                              constant parameters.
#'                              The proposal distribution of the Metropolis step is a normal distribution centered
#'                              at the last point of the chain with a covariance matrix equal to 
#'                              \code{scale.prop.const^2 * cov.prop.const}.
#'                              Note that if \code{param.log} is \code{TRUE} for a parameter, then the proposal
#'                              is evaluated at the log scale of the parameter.
#'                              During the adaptation phase, the covariance matrix is periodically adapted to the
#'                              covariance matrix of the current sample and the scale to get a reasonable
#'                              acceptance rate.
#'                              After the adaptation phase, both variables are kept constant to guarantee 
#'                              convergence.
#' @param  cov.prop.ou.ini      list of scaled covariance matrices of the proposal distributions for the Metropolis 
#'                              step of the parameters of the Ornstein-Uhlenbeck processes of time-dependent
#'                              parameters.
#'                              The proposal distribution of the Metropolis step for the process parameters of the
#'                              time-dependent parameter i is a normal distribution centered
#'                              at the last point of the chain with a covariance matrix equal to 
#'                              \code{cov.prop.ou[[i]] * scale.prop.ou[i]^2}.
#'                              Note that if \code{param.log} is \code{TRUE} for a parameter, then the proposal
#'                              for the mean is evaluated at the log scale of the parameter.
#'                              This is anyway the case for the standard deviation and the rate parameter of the
#'                              Ornstein-Uhlenbeck process.
#'                              During the adaptation phase, the covariance matrices are periodically adapted to the
#'                              covariance matrix of the current sample and the scale to get a reasonable
#'                              acceptance rate.
#'                              After the adaptation phase, both variables are kept constant to guarantee 
#'                              convergence.      
#' @param  scale.prop.const.ini scale factor for the covariance matrix of the Metropolis step for constant parameters
#'                              with a proposal distribution equal to a normal distribution centered at the previous 
#'                              point of the chain and a covariance matrix equal to 
#'                              \code{scale.prop.const^2 * cov.prop.const}.
#'                              During the adaptation phase, the covariance matrix is periodically adapted to the
#'                              covariance matrix of the current sample and the scale to get a reasonable
#'                              acceptance rate.
#'                              After the adaptation phase, both variables are kept constant to guarantee 
#'                              convergence.      
#' @param  scale.prop.ou.ini    vector of scale factors for the covariance matrices of the Metropolis step for parameters
#'                              of Ornstein-Uhlenbeck processes with a proposal distribution equal to a normal 
#'                              distribution centered at the previous point of the chain and a covariance matrix 
#'                              for the time dependent parameter i equal to \code{cov.prop.ou[[i]] * scale.prop.ou[i]^2}.
#'                              During the adaptation phase, the covariance matrix is periodically adapted to the
#'                              covariance matrix of the current sample and the scale to get a reasonable
#'                              acceptance rate.
#'                              After the adaptation phase, both variables are kept constant to guarantee 
#'                              convergence.      
#' @param  control              list of control parameters of the algorithm:
#'                              \itemize{
#'                              \item
#'                              \code{n.interval}: number of sub-intervals into which the time domain is splitted
#'                                                 to infer the time-dependent parameters; either scalar for universal
#'                                                 choice for all parameters or named vector for parameter-specific
#'                                                 choices
#'                                                 (default value: 50; this number must be increased if the acceptance
#'                                                 rates of the time-dependent parameters are very low, it can be decreased
#'                                                 if they are high);
#'                              \item
#'                              \code{min.internal}: minimum number of internal points in an interval 
#'                                                  (default value: 1; may be increased if time resolution is high).
#'                              \item
#'                              \code{splitmethod}: method used for random splitting of time domain into sub-intervals.
#'                                                  Possible values:
#'                                                  \code{"modunif"}:     modification of uniform intervals;
#'                                                  \code{"random"}:      random split (higher variability in inverval lengths);
#'                                                  \code{"weighted"}:    weighted random split leading to shorter intervals
#'                                                                        where the acceptance frequency is low;
#'                                                  \code{"autoweights"}: use weighted random split but adjusts weights 
#'                                                                        adaptively.
#'                                                  (default value: \code{"modunif"}).
#'                              \item
#'                              \code{interval.weights}: numerical vector or named list of numerical vectors 
#'                                                       (by time-dependent parameter) of weights for sampling interval 
#'                                                       boundaries (the length(s) of the vector(s) must be equal to the
#'                                                       time series length in the parameter specification).
#'                                                       The weight vectors do not have to be normalized.
#'                                                       The weights are used if the parameter \code{splitmethod} is
#'                                                       equal to \code{"weighted"} or as optional initial weights
#'                                                       if \code{splitmethod} is equal to \code{"autoweights"}.
#'                              \item
#'                              \code{n.autoweighting}: number of past iterations to consider for weight calculation
#'                                                      for \code{splitmethod} \code{"autoweights"}
#'                                                      (default value: 1000). 
#'                                                      Note that the calculation of weights only starts after 
#'                                                      \code{n.autoweights} iterations and that only stored points 
#'                                                      are considerd so that the number of points considered is
#'                                                      equal to \code{n.autoweighting}/\code{thin}.
#'                              \item
#'                              \code{offset.weighting}: offset used to caluclate weights from apparent acceptance
#'                                                       frequencies for \code{splitmethod} \code{"autoweights"}
#'                                                       (default value: 0.05).
#'                              \item
#'                              \code{n.widening}: number grid points used to widen areas of high weight
#'                                                 for \code{splitmethod} \code{"autoweights"}
#'                                                 (default value: 10).
#'                              \item
#'                              \code{n.timedep.perstep}: number of updates of the time-dependent parameter(s) before
#'                                                        updating the constant parameters (default value: 1).
#'                              \item
#'                              \code{n.const.perstep}: number of Markov chain steps for the constant parameters
#'                                                      to be performed between updating the time-dependent parameters
#'                                                      (default value: 1).
#'                              \item
#'                              \code{n.adapt}: number of iterations of the Markov chain during which adaptation is made 
#'                                              (default value: 2000; only during this phase, the covariance matrix 
#'                                              and the scaling factors are adapted).
#'                              \item
#'                              \code{n.adapt.scale}: number of iterations after which the acceptance rate is checked for
#'                                                   potentially adapting the scaling factor (default value: 30).
#'                              \item
#'                              \code{n.adapt.cov}: number of iterations of the Markov chain, after which the covariance matrix 
#'                                                  of the proposal distribution is adapted 
#'                                                  (default value: 900; 0 means no adaptation of the covariance matrix;
#'                                                  note that after \code{control$n.adapt} iterations 
#'                                                  adaptation is turned off; for this reason, after the last multiple of
#'                                                  \code{n.adapt.cov} below \code{n.adapt} there should be sufficient 
#'                                                  iterations left to adapt the scaling factors).
#'                              \item
#'                              \code{f.reduce.cor}: factor by which sample correlations are reduced when constructing
#'                                                   the covariance matrix of the proposal distribution (default value: 0.90).
#'                              \item
#'                              \code{f.accept.decscale}: acceptance rate below which the proposal scaling factor is decreased
#'                                                        during the adaptation phase (default value: 0.05).
#'                              \item
#'                              \code{f.accept.incscale}: acceptance rate above which the proposal scaling factor is increased
#'                                                        during the adaptation phase (default value: 0.30).
#'                              \item
#'                              \code{f.max.scalechange}: max. factor for changing proposal distribution scale from reference
#'                                                        (default value: 10; reference is either initial value or modified 
#'                                                        value when the covariance matrix was adapted).
#'                              \item
#'                              \code{f.sample.cov.restart}: fraction of previous samples to be used to calculate the 
#'                                                           covariance matrix of proposal distribution when restarting 
#'                                                           inference (default value: 0.3;
#'                                                           the last part of the samples is used).
#'                              \item
#'                              \code{thin}: thinning for storing Markov chain results (default value: 1).
#'                              \item
#'                              \code{n.save}: number of iterations after which the results are (periodically) saved 
#'                                             (default value: 1000).
#'                              \item
#'                              \code{save.diag}: save diagnostic information about acceptance ratio, acceptance, and 
#'                                                interval lengths for inference of the time-dependent parameters.
#'                              }
#' @param  res.infer.timedeppar results of a previous call to this function.
#'                              These results are ignored if the argument \code{task} is equal to \code{"start"}, but it is
#'                              needed for the tasks \code{"continue"} and \code{"restart"}.     
#' @param  verbose              integer parameter indicating the level of progress reporting:\cr
#'                              0: no reporting;\cr
#'                              1: reporting of thinned and accepted Markov Chain steps and of adapted proposal covariance matrices;\cr
#'                              2: reporting of proposals and accepted steps before thinning.
#' @param  file.save            if non-empty string, the intermediate results are saved to this file as variable \code{res}
#'                              in a workspace after every \code{control$n.save} iterations 
#'                              (the extension \code{.RData} will be appended to the file name).
#' @param  ...                  additional parameters passed to the function \code{loglikeli}.
#' 
#' @return class of type \code{timedeppar} with the following elements:
#'         \itemize{
#'         \item
#'         \code{package}: package timedeppar: version and date,
#'         \item
#'         \code{func}: function called (infer.timedeppar),
#'         \item
#'         \code{date}: date of call,
#'         \item
#'         \code{dot.args}: arguments passed to the likelihood function (included for reproducibility of results),
#'         \item
#'         \code{task}: task that was performed (\code{start}, \code{restart} or \code{continue}),
#'         \item
#'         \code{file}: name of file to which output was written,
#'         \item
#'         \code{param.ini}: initial values of likelihood parameters (constant and time-dependent),
#'         \item
#'         \code{param.ou.ini}: initial values of Ornstein-Uhlenbeck process parameters that are estimated,
#'         \item
#'         \code{param.ou.fixed}: values of Ornstein-Uhlenbeck process parameters that are not estimated,
#'         \item
#'         \code{loglikeli}: function that was passed to calculate the log likelihood of the observations,
#'         \item
#'         \code{loglikeli.keepstate}: boolean indicating whether or not the state from the previous run should be kept 
#'                                     (this allows only partial time evaluation when only part of the input was replaced),
#'         \item
#'         \code{param.logprior}: function that was passed to calculate the joint log prior of the constant likelihood parameters,
#'         \item
#'         \code{param.ou.logprior}: function that was passed to calculate the joint log prior of the estimated Ornstein-Uhlenbeck process parameters
#'               (in case of multiple Ornstein-Uhlenbeck processes the function has to return the prior for the correct process;
#'               this can be identified by the names of the argument),
#'         \item
#'         \code{param.range}: parameter ranges,
#'         \item
#'         \code{param.log}: named logical vector of indicators for log inference,
#'         \item
#'         \code{control}: named list of control parameters as passed to the call (or read from a previous call),
#'         \item
#'         \code{n.iter}: number of iterations peformed (note that the size of the sample will be n.iter/control$thin),
#'         \item
#'         \code{sample.diag}: list of samples of proposals, log acceptance ratios, and interval lengths of time-dependent
#'               parameters (only available if the control variable \code{save.diag} is set to \code{TRUE}),
#'         \item
#'         \code{sample.param.timedep}: list of samples of time dependent parameters (first row contains time points),
#'         \item
#'         \code{sample.param.ou}: sample of Ornstein-Uhlenbeck process parameters,
#'         \item
#'         \code{sample.param.const}: sample of constant parameters,
#'         \item
#'         \code{sample.logpdf}: sample of prior, Ornstein-Uhlenbeck and posterior pdf,
#'         \item
#'         \code{acceptfreq.constpar}: acceptance frequency of constant parameters after adaptation phase,
#'         \item
#'         \code{acceptfreq.oupar}: acceptance frequencies of Ornstein-Uhlenbeck process parameters after adaptation phase,
#'         \item
#'         \code{acceptfreq.timedeppar}: acceptance frequencies of time-depenent parameters,
#'         \item
#'         \code{param.maxpost}: parameters at the maximum posterior (constant and time-dependent parameters),
#'         \item
#'         \code{param.ou.maxpost}: Ornstein-Uhlenbeck process parameters at the maximum posterior,
#'         \item
#'         \code{cov.prop.const}: final covariance matrix used for proposal distribution of constant parameters,
#'         \item
#'         \code{cov.prop.ou}: list of final covariance matrices used for proposal distribution of Ornstein-Uhlenbeck process paramemters,
#'         \item
#'         \code{scale.prop.const}: final scale of proposal distribution of constant parameters,
#'         \item
#'         \code{scale.prop.ou}: final scale of proposal distribution of Ornstein-Uhlenbeck process parameters,
#'         \item
#'         \code{sys.time}: run time used for the previous inference job.
#'         }
#'         
#' @seealso
#' 
#' \code{\link{plot.timedeppar}} for visualizing results.\cr
#' \code{\link{calc.acceptfreq}} for calculating (apparent) acceptance frequencies.\cr
#' \code{\link{calc.logpdf}} for calculating log pdf values (prior, internal, posterior) from the results.\cr
#' \code{\link{get.param}} for extracting individual parameters from the Markov chain.\cr
#' \code{\link{get.parsamp}} for extracting subsamples of the Markov chain.\cr
#' \code{\link{readres.timedeppar}} for reading saved results from a previous run.\cr
#' \code{\link{randOU}} for sampling from an Ornstein-Uhlenbeck process.\cr
#' \code{\link{logpdfOU}} for calculating the probability density of a sample from an Ornstein-Uhlenbeck process.
#' 
#' @references 
#' Reichert, P.
#' timedeppar: An R package for inferring stochastic, time-dependent model parameters
#' in preparation, 2020.\cr\cr
#' Reichert, P., Ammann, L. and Fenicia, F.
#' Potential and challenges of investigating intrinsic uncertainty of hydrological models with stochastic, time-dependent parameters.
#' in preparation, 2020.\cr\cr
#' Reichert, P. and Mieleitner, J. 
#' Analyzing input and structural uncertainty of nonlinear dynamic models with stochastic, time-dependent parameters. 
#' \emph{Water Resources Research}, 45, W10402, 2009. 
#' \url{https://doi.org/10.1029/2009WR007814}\cr\cr
#' Tomassini, L., Reichert, P., Kuensch, H.-R. Buser, C., Knutti, R. and Borsuk, M.E. 
#' A smoothing algorithm for estimating stochastic, continuous-time model parameters 
#' and its application to a simple climate model. 
#' \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} 58, 679-704, 2009. 
#' \url{https://doi.org/10.1111/j.1467-9876.2009.00678.x}
#' 
#' @examples 
#' # re-infer an artificially produced Ornstein-Uhlenbeck process with noise
#' # (didactical example for package use; increase length of chain for convergence)
#' 
#' # generate data:
#' 
#' set.seed(123)
#' data  <- randOU(mean=0,sd=1,gamma=10,t=seq(from=0,to=2,length.out=101)) 
#' data$yobs <- data$y + rnorm(nrow(data),mean=0,sd=0.2)
#' 
#' # define observational likelihood:
#' 
#' loglikeli <- function(param,data)
#' {
#'   # get parameter y at time points of observations::
#'   
#'   y <- param$y
#'   if ( is.matrix(y) | is.data.frame(y) ) y <- approx(x=y[,1],y=y[,2],xout=data[,1])$y
#'   
#'   # calculate likelihood:
#'   loglikeli <- sum(dnorm(data[,"yobs"],mean=y,sd=0.2))
#'   
#'   # return result:
#'   return(loglikeli)
#' }
#' 
#' # do a few MCMC steps for inference of y and epsilon for didactical purpose
#' # (much longer chains would be needed for convergence, try n.iter=10000, n.adapt=2000):
#' 
#' res <- infer.timedeppar(loglikeli      = loglikeli,
#'                         param.ini      = list(y=data[,1:2]),
#'                         param.ou.ini   = c(y_mean=0,y_sd=1),
#'                         param.ou.fixed = c(y_gamma=10),
#'                         n.iter         = 150,
#'                         control        = list(n.interval = 10,
#'                                               n.adapt    = 75),
#'                         data           = data)
#' 
#' plot(res,
#'      labels=expression(y       = y,
#'                        y_mean  = mu[y],
#'                        y_sd    = sigma[y],
#'                        y_gamma = gamma[y]))


infer.timedeppar <- function(loglikeli            = NULL,
                             loglikeli.keepstate  = FALSE,
                             param.ini            = list(),
                             param.range          = list(),
                             param.log            = logical(0),
                             param.logprior       = NULL,
                             param.ou.ini         = numeric(0),
                             param.ou.fixed       = numeric(0),
                             param.ou.logprior    = NULL,
                             task                 = c("start","continue","restart"),
                             n.iter               = NA,
                             cov.prop.const.ini   = NA,
                             cov.prop.ou.ini      = NA,
                             scale.prop.const.ini = NA,
                             scale.prop.ou.ini    = NA,
                             control              = list(),
                             res.infer.timedeppar = list(),
                             verbose              = 0,
                             file.save            = "",
                             ...)
{
  start.time <- proc.time()
  res <- list()
  class(res)   <- "timedeppar"
  res$package  <- paste("timedeppar",
                        as.character(packageVersion("timedeppar")),
                        as.character(packageDate("timedeppar"))) 
  res$func     <- "infer.timedeppar"
  res$date     <- Sys.Date()
  res$dot.args <- as.list(match.call(expand.dots=FALSE))[["..."]]
  if ( length(res$dot.args) > 0 )
  {
    for ( i in 1:length(res$dot.args) ) { if ( !is.null(res$dot.args[[i]]) ) res$dot.args[[i]] <- eval(res$dot.args[[i]]) }
  }  
  res$task     <- task[1]
  res$file     <- file.save
    
  # check input:
  # ------------
  
  errmsg <- character(0)
  
  # check argment types:
  
  if ( is.null(loglikeli) )
  {
    if ( task[1] == "start" )
    {
      errmsg <- c(errmsg,"*** infer.timedeppar: the argument loglikeli must be a function")
    }
  }
  else
  {
    if ( !is.function(loglikeli) )
    {
      errmsg <- c(errmsg,"*** infer.timedeppar: the argument loglikeli must be a function")
    }
  }
  if ( !is.list(param.ini) )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: the argument param.ini must be a list")
  }
  if ( !is.list(param.range) )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: the argument param.range must be a list")
  }
  if ( !is.logical(param.log) )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: the argument param.log must be a logical vector")
  }
  if ( !is.null(param.logprior) )
  {
    if ( !is.function(param.logprior) )
    {
      errmsg <- c(errmsg,"*** infer.timedeppar: the argument param.logprior must be a function")
    }
  }
  if ( !is.numeric(param.ou.ini) )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: the argument param.ou.ini must be a numerical vector")
  }
  if ( !is.numeric(param.ou.fixed) )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: the argument param.ou.fixed must be a numerical vector")
  }
  if ( !is.null(param.ou.logprior) )
  {
    if ( !is.function(param.ou.logprior) )
    {
      errmsg <- c(errmsg,"*** infer.timedeppar: the argument param.ou.logprior must be a function")
    }
  }
  if ( !is.character(task) )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: the argument task must be of type character")
  }
  if ( !is.list(control) )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: the argument control must be a list")
  }
  
  # check availability and plausibility of arguments:
  
  if ( ! task[1] %in% c("start","continue","restart") )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: task must be \"start\", \"continue\", or \"restart\"")
  }
  
  if ( task[1] == "start" )
  {
    if ( is.null(loglikeli) )
    {
      errmsg <- c(errmsg,
                  paste("*** infer.timedeppar: function loglikeli must be provided for task \"start\""))
    }
    if ( length(param.ini) == 0 | !is.list(param.ini) )
    {
      errmsg <- c(errmsg,
                  paste("*** infer.timedeppar: param.ini must be a named list of parameters;",
                        "constant parameters must be represented by an initival value,",
                        "time-dependent parameters by a two column matrix for time and values"))
    }
  }
  else
  {
    if ( length(res.infer.timedeppar) == 0 | !is.list(res.infer.timedeppar) | class(res.infer.timedeppar)!="timedeppar" )
    {
      errmsg <- c(errmsg,
                  paste("*** infer.timedeppar: for tasks \"continue\" or task \"restart\" old results must be provided",
                        "as argument \"res.infer.timedeppar\""))
    }
    else
    {
      if ( is.null(res.infer.timedeppar$param.ini) |
           is.null(res.infer.timedeppar$control) |
           is.null(res.infer.timedeppar$cov.prop.const) |
           is.null(res.infer.timedeppar$cov.prop.ou) |
           is.null(res.infer.timedeppar$scale.prop.const) |
           is.null(res.infer.timedeppar$scale.prop.ou) |
           is.null(res.infer.timedeppar$sample.param.timedep) | 
           is.null(res.infer.timedeppar$sample.param.ou) |
           is.null(res.infer.timedeppar$sample.param.const) |
           is.null(res.infer.timedeppar$sample.param.const) )
      {
        errmsg <- c(errmsg,
                    paste("*** infer.timedeppar: illegal results passed as argument \"res.infer.timedeppar\":",
                          "no element(s)",
                          "\"param.ini\",",
                          "\"control\",",
                          "\"cov.prop.const\",",
                          "\"cov.prop.ou\",",
                          "\"scale.prop.const\",",
                          "\"scale.prop.ou\",",
                          "\"sample.param.timedep\",",
                          "\"sample.param.ou\",",
                          "\"sample.param.const\", or",
                          "\"sample.param.const\""))
      }
    }
  }
  
  
  if ( length(errmsg) > 0 )
  {
    res$errmsg <- errmsg
    res$sys.time <- proc.time() - start.time
    for ( i in 1:length(errmsg) ) warning(errmsg[i])
    return(res)
  }
  
  # analyze structure of likelihood parameter list (which constant, which time-dependent):
  # --------------------------------------------------------------------------------------
  
  errmsg <- character(0)
  param <- param.ini; if ( task[1] != "start" ) param <- res.infer.timedeppar$param.ini
  param.names <- names(param)
  param.ou <- param.ou.ini; if ( task[1] != "start" ) param.ou <- res.infer.timedeppar$param.ou.ini
  param.ou.names <- names(param.ou)
  if ( task[1] == "continue" | (task[1]=="restart" & length(param.log)==0) )
  {
    param.log <- res.infer.timedeppar$param.log
  }
  if ( task[1] == "continue" | (task[1]=="restart" & length(param.range)==0) ) 
  {
    param.range <- res.infer.timedeppar$param.range
  }
  ind.timedeppar      <- numeric(0)
  ind.constpar        <- numeric(0)
  param.log.loc       <- rep(FALSE,length(param))
  names(param.log.loc) <- param.names
  for ( i in 1:length(param) )
  {
    if ( param.names[i] %in% names(param.log) )
    {
      if ( param.log[param.names[i]] ) param.log.loc[i] <- TRUE
    }
    if ( is.matrix(param[[i]]) | is.data.frame(param[[i]]))
    {
      ind.timedeppar <- c(ind.timedeppar,i)
      if ( anyNA(param[[i]][,1]) )
      {
        errmsg <- c(errmsg,paste("*** infer.timedeppar: time points of time-dependent parameter ",param.names[i],
                                 " contain NAs",sep=""))
      }
      else
      {
        if ( any(diff(param[[i]][,1])<=0) )
        {
          errmsg <- c(errmsg,paste("*** infer.timedeppar: time points of time-dependent parameter ",param.names[i],
                                   " are not strictly increasing",sep=""))
        }
      }
      if ( param.names[i] %in% names(param.range) )
      {
        range <- param.range[[param.names[i]]]
        if ( ! anyNA(param[[i]][,2]) )  # otherwise the parameter value will be drawn below from an OU process
        {
          if ( min(param[[i]][,2]) < range[1] | max(param[[i]][,2]) > range[2] )
          {
            errmsg <- c(errmsg,paste("*** infer.timedeppar: initial value of parameter \"",param.names[i],"\" is out of range.",sep=""))
          }
        }
      }
    }
    else
    {
      ind.constpar   <- c(ind.constpar,i)
      if ( is.na(param[[i]]) )
      {
        errmsg <- c(errmsg,paste("*** infer.timedeppar: initial value of parameter \"",param.names[i],"\" is not available.",sep=""))
      }
      else
      {
        if ( param.names[i] %in% names(param.range) )
        {
          range <- param.range[[param.names[i]]]
          if ( param[[i]] < range[1] | param[[i]] > range[2] )
          {
            errmsg <- c(errmsg,paste("*** infer.timedeppar: initial value of parameter \"",param.names[i],"\" is out of range.",sep=""))
          }
        }
      }
    }
  }
  if ( length(param.ou) > 0 )
  {
    for ( i in 1:length(param.ou) )
    {
      if ( is.na(param.ou[i]) )
      {
        errmsg <- c(errmsg,paste("*** infer.timedeppar: initial value of parameter \"",param.ou.names[i],"\" is not available.",sep=""))
      }
      else
      {
        if ( param.ou.names[i] %in% names(param.range) )
        {
          range <- param.range[[param.ou.names[i]]]
          if ( param.ou[i] < range[1] | param.ou[i] > range[2] )
          {
            errmsg <- c(errmsg,paste("*** infer.timedeppar: initial value of parameter \"",param.ou.names[i],"\" is out of range.",sep=""))
          }
        }
      }
    }
  }
  if ( length(errmsg) > 0 )
  {
    res$errmsg <- errmsg
    res$sys.time <- proc.time() - start.time
    for ( i in 1:length(errmsg) ) warning(errmsg[i])
    return(res)
  }
  
  func.version <- paste("infer.timedeppar (timedeppar ",
                        as.character(packageVersion("timedeppar"))," ",
                        as.character(packageDate("timedeppar")),"):",sep="")
  if ( task[1] == "start" )
  {
    message(paste(func.version,"starting new Markov Chain"))
  }
  else
  {
    if ( task[1] == "continue" )
    {
      message(paste(func.version,"continuing existing Markov chain"))
    }
    else
    {
      message(paste(func.version,"restarting Markov Chain from existing chain"))
    }
  }
  message(paste("  number of constant parameters:         ",length(ind.constpar)))
  message(paste("  number of time-dependent parameters:   ",length(ind.timedeppar)))
  if ( length(ind.timedeppar) > 0 )
  {
    message(paste("  number of estimated process parameters:",length(param.ou)))
  }
  
  # initialize variables:
  # ---------------------
  
  thin                 <-     1
  n.interval           <-    50
  min.internal         <-     1
  splitmethod          <- "modunif"
  interval.weights     <- numeric(0)
  n.autoweighting      <-  1000
  offset.weighting     <-     0.05
  n.widening           <-    10
  n.timedep.perstep    <-     1
  n.const.perstep      <-     1
  n.adapt              <-  2000
  n.adapt.cov          <-   900
  f.reduce.cor         <-     0.90
  f.accept.incscale    <-     0.30
  f.accept.decscale    <-     0.05
  f.max.scalechange    <-    10
  f.sample.cov.restart <-     0.3
  n.adapt.scale        <-    30
  n.save               <-  1000
  save.diag            <- FALSE
  
  if ( !is.null(control$thin) )                 thin                 <- control$thin
  if ( !is.null(control$n.interval) )           n.interval           <- control$n.interval
  if ( !is.null(control$min.internal) )         min.internal         <- control$min.internal
  if ( !is.null(control$splitmethod) )          splitmethod          <- control$splitmethod
  if ( !is.null(control$interval.weights) )     interval.weights     <- control$interval.weights
  if ( !is.null(control$n.autoweighting) )      n.autoweighting      <- control$n.autoweighting
  if ( !is.null(control$offset.weighting) )     offset.weighting     <- control$offset.weighting
  if ( !is.null(control$n.widening) )           n.widening           <- control$n.widening
  if ( !is.null(control$n.timedep.perstep) )    n.timedep.perstep    <- control$n.timedep.perstep
  if ( !is.null(control$n.const.perstep) )      n.const.perstep      <- control$n.const.perstep
  if ( !is.null(control$n.adapt) )              n.adapt              <- control$n.adapt
  if ( !is.null(control$n.adapt.scale) )        n.adapt.scale        <- control$n.adapt.scale
  if ( !is.null(control$n.adapt.cov) )          n.adapt.cov          <- control$n.adapt.cov
  if ( !is.null(control$f.reduce.cor) )         f.reduce.cor         <- control$f.reduce.cor
  if ( !is.null(control$f.accept.incscale) )    f.accept.incscale    <- control$f.accept.incscale
  if ( !is.null(control$f.accept.decscale) )    f.accept.decscale    <- control$f.accept.decscale
  if ( !is.null(control$f.max.scalechange) )    f.max.scalechange    <- control$f.max.scalechange
  if ( !is.null(control$f.sample.cov.restart) ) f.sample.cov.restart <- control$f.sample.cov.restart
  if ( !is.null(control$n.save) )               n.save               <- control$n.save
  if ( !is.null(control$save.diag) )            save.diag            <- control$save.diag
  
  cov.prop.const    <- cov.prop.const.ini
  scale.prop.const  <- scale.prop.const.ini
  cov.prop.ou       <- cov.prop.ou.ini
  scale.prop.ou     <- scale.prop.ou.ini
  
  if ( task[1] == "start" )
  {
    res$param.ini           <- param.ini
    res$param.ou.ini        <- param.ou.ini
    res$param.ou.fixed      <- param.ou.fixed
    res$loglikeli           <- loglikeli
    res$loglikeli.keepstate <- loglikeli.keepstate
    if ( is.null(param.logprior) )    res$param.logprior    <- function(x) return(0)
    else                              res$param.logprior    <- param.logprior
    if ( is.null(param.ou.logprior) ) res$param.ou.logprior <- function(x) return(0)
    else                              res$param.ou.logprior <- param.ou.logprior
    
    if ( is.na(n.iter) ) n.iter <- 10000
    
    if ( !is.matrix(cov.prop.const) )
    {
      cov.prop.const   <- diag(1,length(ind.constpar))
      if ( length(ind.constpar) > 0 )
      {
        colnames(cov.prop.const) <- param.names[ind.constpar]
        rownames(cov.prop.const) <- param.names[ind.constpar]
        for ( i in 1:length(ind.constpar) )
        {
          if ( param[[param.names[ind.constpar[i]]]] != 0 & !param.log.loc[ind.constpar[i]] ) 
          {
            cov.prop.const[i,i] <- (0.1*param[[param.names[ind.constpar[i]]]])^2
          }
        }
      }
    }
    if ( is.na(scale.prop.const) ) scale.prop.const <- 0.05
  }
  if ( task[1] == "restart" )
  {
    if ( !is.null(res.infer.timedeppar$control$thin) )                 thin                 <- res.infer.timedeppar$control$thin
    if ( !is.null(res.infer.timedeppar$control$n.interval) )           n.interval           <- res.infer.timedeppar$control$n.interval
    if ( !is.null(res.infer.timedeppar$control$min.internal) )         min.internal         <- res.infer.timedeppar$control$min.internal
    if ( !is.null(res.infer.timedeppar$control$splitmethod) )          splitmethod          <- res.infer.timedeppar$control$splitmethod
    if ( !is.null(res.infer.timedeppar$control$interval.weights) )     interval.weights     <- res.infer.timedeppar$control$interval.weights
    if ( !is.null(res.infer.timedeppar$control$n.autoweighting) )      n.autoweighting      <- res.infer.timedeppar$control$n.autoweighting
    if ( !is.null(res.infer.timedeppar$control$offset.weighting) )     offset.weighting     <- res.infer.timedeppar$control$offset.weighting
    if ( !is.null(res.infer.timedeppar$control$n.widening) )           n.widening           <- res.infer.timedeppar$control$n.widening
    if ( !is.null(res.infer.timedeppar$control$n.timedep.perstep) )    n.timedep.perstep    <- res.infer.timedeppar$control$n.timedep.perstep
    if ( !is.null(res.infer.timedeppar$control$n.const.perstep) )      n.const.perstep      <- res.infer.timedeppar$control$n.const.perstep
    if ( !is.null(res.infer.timedeppar$control$n.adapt) )              n.adapt              <- res.infer.timedeppar$control$n.adapt
    if ( !is.null(res.infer.timedeppar$control$n.adapt.scale) )        n.adapt.scale        <- res.infer.timedeppar$control$n.adapt.scale
    if ( !is.null(res.infer.timedeppar$control$n.adapt.cov) )          n.adapt.cov          <- res.infer.timedeppar$control$n.adapt.cov
    if ( !is.null(res.infer.timedeppar$control$f.reduce.cor) )         f.reduce.cor         <- res.infer.timedeppar$control$f.reduce.cor
    if ( !is.null(res.infer.timedeppar$control$f.accept.incscale) )    f.accept.incscale    <- res.infer.timedeppar$control$f.accept.incscale
    if ( !is.null(res.infer.timedeppar$control$f.accept.decscale) )    f.accept.decscale    <- res.infer.timedeppar$control$f.accept.decscale
    if ( !is.null(res.infer.timedeppar$control$f.max.scalechange) )    f.max.scalechange    <- res.infer.timedeppar$control$f.max.scalechange
    if ( !is.null(res.infer.timedeppar$control$f.sample.cov.restart) ) f.sample.cov.restart <- res.infer.timedeppar$control$f.sample.cov.restart
    if ( !is.null(res.infer.timedeppar$control$n.save) )               n.save               <- res.infer.timedeppar$control$n.save
    if ( !is.null(res.infer.timedeppar$control$save.diag) )            save.diag            <- res.infer.timedeppar$control$save.diag
    
    if ( !is.null(control$thin) )                 thin                 <- control$thin
    if ( !is.null(control$n.interval) )           n.interval           <- control$n.interval
    if ( !is.null(control$min.internal) )         min.internal         <- control$min.internal
    if ( !is.null(control$splitmethod) )          splitmethod          <- control$splitmethod
    if ( !is.null(control$interval.weights) )     interval.weights     <- control$interval.weights
    if ( !is.null(control$n.autoweighting) )      n.autoweighting      <- control$n.autoweighting
    if ( !is.null(control$offset.weighting) )     offset.weighting     <- control$offset.weighting
    if ( !is.null(control$n.widening) )           n.widening           <- control$n.widening
    if ( !is.null(control$n.timedep.perstep) )    n.timedep.perstep    <- control$n.timedep.perstep
    if ( !is.null(control$n.const.perstep) )      n.const.perstep      <- control$n.const.perstep
    if ( !is.null(control$n.adapt) )              n.adapt              <- control$n.adapt
    if ( !is.null(control$n.adapt.scale) )        n.adapt.scale        <- control$n.adapt.scale
    if ( !is.null(control$n.adapt.cov) )          n.adapt.cov          <- control$n.adapt.cov
    if ( !is.null(control$f.reduce.cor) )         f.reduce.cor         <- control$f.reduce.cor
    if ( !is.null(control$f.accept.incscale) )    f.accept.incscale    <- control$f.accept.incscale
    if ( !is.null(control$f.accept.decscale) )    f.accept.decscale    <- control$f.accept.decscale
    if ( !is.null(control$f.max.scalechange) )    f.max.scalechange    <- control$f.max.scalechange
    if ( !is.null(control$f.sample.cov.restart) ) f.sample.cov.restart <- control$f.sample.cov.restart
    if ( !is.null(control$n.save) )               n.save               <- control$n.save
    if ( !is.null(control$save.diag) )            save.diag            <- control$save.diag
    
    if ( is.na(n.iter) ) n.iter <- res.infer.timedeppar$n.iter
    
    if ( length(res$param.ini) == 0 )
    {
      # start with final parameter values, not with maximum posterior ones:
      
      if ( length(ind.constpar) > 0 ) 
      {
        for ( i in ind.constpar ) param[[param.names[i]]] <- as.numeric(res.infer.timedeppar$sample.param.const[nrow(res.infer.timedeppar$sample.param.const),param.names[i]])
      }
      if ( length(ind.timedeppar) > 0 )
      {
        param.ou <- res.infer.timedeppar$sample.param.ou[nrow(res.infer.timedeppar$sample.param.ou),]
        for ( i in ind.timedeppar ) 
        {
          param[[param.names[i]]] <- cbind(res.infer.timedeppar$sample.param.timedep[[param.names[i]]][1,],
                                           res.infer.timedeppar$sample.param.timedep[[param.names[i]]][nrow(res.infer.timedeppar$sample.param.timedep[[param.names[i]]]),])
        }
      }
    }
    res$param.ini         <- param
    res$param.ou.ini      <- param.ou
    if ( length(param.ou.fixed) == 0 ) param.ou.fixed <- res.infer.timedeppar$param.ou.fixed
    res$param.ou.fixed    <- param.ou.fixed
    if ( !is.null(loglikeli) )
    {
      res$loglikeli           <- loglikeli
      res$loglikeli.keepstate <- loglikeli.keepstate
    }
    else
    {
      res$loglikeli           <- res.infer.timedeppar$loglikeli
      if ( is.null(res.infer.timedeppar$loglikeli.keepstate) ) res$loglikeli.keepstate <- loglikeli.keepstate  # for compatibility with old files
      else                                                     res$loglikeli.keepstate <- res.infer.timedeppar$loglikeli.keepstate
    }
    if ( !is.null(param.logprior) ) res$param.logprior <- param.logprior
    else                            res$param.logprior <- res.infer.timedeppar$param.logprior
    if ( !is.null(param.ou.logprior) ) res$param.ou.logprior <- param.ou.logprior
    else                               res$param.ou.logprior <- res.infer.timedeppar$param.ou.logprior

    if ( is.na(n.iter) ) n.iter <- res.infer.timedeppar$n.iter
    
    if ( is.na(f.sample.cov.restart) | f.sample.cov.restart <= 0 )
    {
      if ( !is.matrix(cov.prop.const) ) cov.prop.const   <- res.infer.timedeppar$cov.prop.const
      if ( is.na(scale.prop.const) )    scale.prop.const <- res.infer.timedeppar$scale.prop.const
      if ( !is.list(cov.prop.ou) )      cov.prop.ou      <- res.infer.timedeppar$cov.prop.ou
      if ( is.na(scale.prop.ou[1])    ) scale.prop.ou    <- res.infer.timedeppar$scale.prop.ou
    }
    else
    {
      n.sample <- nrow(res.infer.timedeppar$sample.logpdf)
      ind.sample <- max(1,ceiling((1-f.sample.cov.restart)*n.sample)):n.sample
      if ( length(ind.constpar) > 0 )
      {
        sample.log <-  res.infer.timedeppar$sample.param.const[ind.sample,,drop=FALSE]
        for ( i in 1:length(ind.constpar) ) 
        { 
          if ( param.log.loc[ind.constpar[i]] ) sample.log[,i] <- log(sample.log[,i,drop=FALSE]) 
          if ( var(sample.log[,i]) == 0 ) 
          {
            warning("infer.timedeppar: unable to define covariance matrix of proposal of constant parameters from old sample")
            return(res)
          }
        }
        vol.old <- sqrt(prod(diag(res.infer.timedeppar$cov.prop.const)))  # expect to be more robust without considering correlation
        cov.prop.const <- cov(sample.log)
        if ( ncol(cov.prop.const) > 1 ) cov.prop.const <- f.reduce.cor*cov.prop.const + (1-f.reduce.cor)*diag(diag(cov.prop.const))
        colnames(cov.prop.const) <- param.names[ind.constpar]
        rownames(cov.prop.const) <- param.names[ind.constpar]
        vol.new <- sqrt(prod(diag(cov.prop.const)))  # expect to be more robust without considering correlation
        scale.prop.const <- res.infer.timedeppar$scale.prop.const * (vol.old/vol.new)^(1/nrow(cov.prop.const))
      }
       
      if ( length(ind.timedeppar) > 0 )
      {
        ind.oupar <- list()
        poss.names <- character(0)
        for ( i in ind.timedeppar )
        {
          ind.oupar[[param.names[i]]] <- numeric(0)
          poss.names <- c(poss.names,paste(param.names[i],c("mean","sd","gamma"),sep="_"))
          if ( paste(param.names[i],"mean",sep="_") %in% param.ou.names ) 
          {
            ind.oupar[[param.names[i]]] <- c(ind.oupar[[param.names[i]]],match(paste(param.names[i],"mean",sep="_"),param.ou.names))
          }
          if ( paste(param.names[i],"sd",sep="_") %in% param.ou.names ) 
          {
            ind.oupar[[param.names[i]]] <- c(ind.oupar[[param.names[i]]],match(paste(param.names[i],"sd",sep="_"),param.ou.names))
          }
          if ( paste(param.names[i],"gamma",sep="_") %in% param.ou.names ) 
          {
            ind.oupar[[param.names[i]]] <- c(ind.oupar[[param.names[i]]],match(paste(param.names[i],"gamma",sep="_"),param.ou.names))
          }
        }

        cov.prop.ou <- list()
        for ( i in 1:length(ind.timedeppar) )
        {
          if ( length(ind.oupar[[i]]) > 0 )
          {
            name <- param.names[ind.timedeppar[i]]
            sample.log <- res.infer.timedeppar$sample.param.ou[ind.sample,ind.oupar[[i]],drop=FALSE]
            for ( ii in 1:length(ind.oupar[[i]]) )
            {
              if ( ! ( !param.log.loc[ind.timedeppar[i]] & param.ou.names[ind.oupar[[i]][ii]]==paste(name,"mean",sep="_") ) )
              {
                sample.log[,ii] <- log(sample.log[,ii,drop=FALSE])
              }
              if ( var(sample.log[,ii]) == 0 ) 
              {
                warning(paste("infer.timedeppar: unable to define covariance matrix of proposal of time-dependent parameter",name,"from old sample"))
                return(res)
              }
            }
            vol.old <- sqrt(prod(diag(res.infer.timedeppar$cov.prop.ou[[i]])))  # expect to be more robust without considering correlation
            cov.prop.ou[[i]] <- cov(sample.log)
            if ( ncol(cov.prop.ou[[i]]) > 1 ) cov.prop.ou[[i]] <- f.reduce.cor*cov.prop.ou[[i]] + (1-f.reduce.cor)*diag(diag(cov.prop.ou[[i]]))
            colnames(cov.prop.ou[[i]]) <- param.ou.names[ind.oupar[[i]]]
            rownames(cov.prop.ou[[i]]) <- param.ou.names[ind.oupar[[i]]]
            vol.new <- sqrt(prod(diag(cov.prop.ou[[i]])))  # expect to be more robust without considering correlation
            scale.prop.ou[i] <- res.infer.timedeppar$scale.prop.ou[i] * (vol.old/vol.new)^(1/nrow(cov.prop.ou[[i]])) 
          }
          else
          {
            cov.prop.ou[[i]] <- matrix(nrow=0,ncol=0)
            scale.prop.ou[i] <- res.infer.timedeppar$scale.prop.ou 
          }
        }
      }
    }
  }
  if ( task[1] == "continue" )
  {
    if ( !is.null(res.infer.timedeppar$control$thin) )                 thin                 <- res.infer.timedeppar$control$thin
    if ( !is.null(res.infer.timedeppar$control$n.interval) )           n.interval           <- res.infer.timedeppar$control$n.interval
    if ( !is.null(res.infer.timedeppar$control$min.internal) )         min.internal         <- res.infer.timedeppar$control$min.internal
    if ( !is.null(res.infer.timedeppar$control$splitmethod) )          splitmethod          <- res.infer.timedeppar$control$splitmethod
    if ( !is.null(res.infer.timedeppar$control$interval.weights) )     interval.weights     <- res.infer.timedeppar$control$interval.weights
    if ( !is.null(res.infer.timedeppar$control$n.autoweighting) )      n.autoweighting      <- res.infer.timedeppar$control$n.autoweighting
    if ( !is.null(res.infer.timedeppar$control$offset.weighting) )     offset.weighting     <- res.infer.timedeppar$control$offset.weighting
    if ( !is.null(res.infer.timedeppar$control$n.widening) )           n.widening           <- res.infer.timedeppar$control$n.widening
    if ( !is.null(res.infer.timedeppar$control$n.timedep.perstep) )    n.timedep.perstep    <- res.infer.timedeppar$control$n.timedep.perstep
    if ( !is.null(res.infer.timedeppar$control$n.const.perstep) )      n.const.perstep      <- res.infer.timedeppar$control$n.const.perstep
    if ( !is.null(res.infer.timedeppar$control$n.adapt) )              n.adapt              <- res.infer.timedeppar$control$n.adapt
    if ( !is.null(res.infer.timedeppar$control$n.adapt.scale) )        n.adapt.scale        <- res.infer.timedeppar$control$n.adapt.scale
    if ( !is.null(res.infer.timedeppar$control$n.adapt.cov) )          n.adapt.cov          <- res.infer.timedeppar$control$n.adapt.cov
    if ( !is.null(res.infer.timedeppar$control$f.reduce.cor) )         f.reduce.cor         <- res.infer.timedeppar$control$f.reduce.cor
    if ( !is.null(res.infer.timedeppar$control$f.accept.incscale) )    f.accept.incscale    <- res.infer.timedeppar$control$f.accept.incscale
    if ( !is.null(res.infer.timedeppar$control$f.accept.decscale) )    f.accept.decscale    <- res.infer.timedeppar$control$f.accept.decscale
    if ( !is.null(res.infer.timedeppar$control$f.max.scalechange) )    f.max.scalechange    <- res.infer.timedeppar$control$f.max.scalechange
    if ( !is.null(res.infer.timedeppar$control$f.sample.cov.restart) ) f.sample.cov.restart <- res.infer.timedeppar$control$f.sample.cov.restart
    if ( !is.null(res.infer.timedeppar$control$n.save) )               n.save               <- res.infer.timedeppar$control$n.save
    if ( !is.null(res.infer.timedeppar$control$save.diag) )            save.diag            <- res.infer.timedeppar$control$save.diag
    
    #if ( !is.null(control$thin) )                 thin              <- control$thin
    #if ( !is.null(control$n.interval) )           n.interval        <- control$n.interval
    #if ( !is.null(control$min.internal) )         min.internal      <- control$min.internal
    #if ( !is.null(control$splitmethod) )          splitmethod       <- control$splitmethod
    #if ( !is.null(control$interval.weights) )     interval.weights  <- control$interval.weights
    #if ( !is.null(control$n.autoweighting) )      n.autoweighting   <- control$n.autoweighting
    #if ( !is.null(control$offset.weighting) )     offset.weighting  <- control$offset.weighting
    #if ( !is.null(control$n.widening) )           n.widening        <- control$n.widening
    if ( !is.null(control$n.timedep.perstep) )    n.timedep.perstep <- control$n.timedep.perstep
    #if ( !is.null(control$n.const.perstep) )      n.const.perstep   <- control$n.const.perstep
    if ( !is.null(control$n.adapt) )              n.adapt           <- control$n.adapt
    if ( !is.null(control$n.adapt.scale) )        n.adapt.scale     <- control$n.adapt.scale
    if ( !is.null(control$n.adapt.cov) )          n.adapt.cov       <- control$n.adapt.cov
    if ( !is.null(control$f.reduce.cor) )         f.reduce.cor      <- control$f.reduce.cor
    if ( !is.null(control$f.accept.incscale) )    f.accept.incscale <- control$f.accept.incscale
    if ( !is.null(control$f.accept.decscale) )    f.accept.decscale <- control$f.accept.decscale
    if ( !is.null(control$f.max.scalechange) )    f.max.scalechange <- control$f.max.scalechange
    if ( !is.null(control$f.sample.cov.restart) ) f.sample.cov.restart <- control$f.sample.cov.restart
    if ( !is.null(control$n.save) )               n.save            <- control$n.save
    if ( save.diag ) { if ( !is.null(control$save.diag) ) save.diag <- control$save.diag }

    if ( is.na(n.iter) ) n.iter <- res.infer.timedeppar$n.iter

    if ( length(ind.constpar) > 0 ) 
    {
      for ( i in ind.constpar ) param[[param.names[i]]] <- as.numeric(res.infer.timedeppar$sample.param.const[nrow(res.infer.timedeppar$sample.param.const),param.names[i]])
    }
    if ( length(ind.timedeppar) > 0 )
    {
      param.ou <- res.infer.timedeppar$sample.param.ou[nrow(res.infer.timedeppar$sample.param.ou),]
      for ( i in ind.timedeppar ) 
      {
        param[[param.names[i]]] <- cbind(res.infer.timedeppar$sample.param.timedep[[param.names[i]]][1,],
                                         res.infer.timedeppar$sample.param.timedep[[param.names[i]]][nrow(res.infer.timedeppar$sample.param.timedep[[param.names[i]]]),])
      }
    }
    
    res$param.ini           <- res.infer.timedeppar$param.ini
    res$param.ou.ini        <- res.infer.timedeppar$param.ou.ini
    param.ou.fixed          <- res.infer.timedeppar$param.ou.fixed
    res$param.ou.fixed      <- param.ou.fixed
    res$loglikeli           <- res.infer.timedeppar$loglikeli
    if ( is.null(res.infer.timedeppar$loglikeli.keepstate) ) res$loglikeli.keepstate <- loglikeli.keepstate  # for compatibility with old files
    else                                                     res$loglikeli.keepstate <- res.infer.timedeppar$loglikeli.keepstate
    res$param.logprior      <- res.infer.timedeppar$param.logprior
    res$param.ou.logprior   <- res.infer.timedeppar$param.ou.logprior

    if ( is.na(n.iter) ) n.iter <- res.infer.timedeppar$n.iter
    
    cov.prop.const        <- res.infer.timedeppar$cov.prop.const
    cov.prop.ou           <- res.infer.timedeppar$cov.prop.ou
    scale.prop.const      <- res.infer.timedeppar$scale.prop.const
    scale.prop.ou         <- res.infer.timedeppar$scale.prop.ou
  }
  res$param.range       <- param.range
  res$param.log         <- param.log
  if ( length(ind.constpar) > 0 & verbose > 0 )
  {
    cat("* proposal distribution of constant parameters:\n")
    cat("correlation matrix:\n")
    print(cov2cor(cov.prop.const))
    cat("standard deviations:\n")
    print(sqrt(diag(cov.prop.const)))
    cat("log:\n")
    print(param.log.loc[ind.constpar])
    cat("scaling factor:\n",scale.prop.const,"\n")
  }
  
  res$control <- list(thin                 = thin,
                      n.interval           = n.interval,
                      min.internal         = min.internal,
                      splitmethod          = splitmethod,
                      interval.weights     = interval.weights,
                      n.autoweighting      = n.autoweighting,
                      offset.weighting     = offset.weighting,
                      n.widening           = n.widening,
                      n.timedep.perstep    = n.timedep.perstep,
                      n.const.perstep      = n.const.perstep,
                      n.adapt              = n.adapt,
                      n.adapt.scale        = n.adapt.scale,
                      n.adapt.cov          = n.adapt.cov,
                      f.reduce.cor         = f.reduce.cor,
                      f.accept.incscale    = f.accept.incscale,
                      f.accept.decscale    = f.accept.decscale,
                      f.max.scalechange    = f.max.scalechange,
                      f.sample.cov.restart = f.sample.cov.restart,
                      n.save               = n.save,
                      save.diag            = save.diag)
  
  # check completeness of Ornstein-Uhlenbeck process parameters (either to be estimated or fixed)
  # and initialize time-dependent parameters:
  # ---------------------------------------------------------------------------------------------

  n.interval.all <- n.interval
  ind.oupar <- list()
  if ( length(ind.timedeppar) > 0 )
  {
    poss.names <- character(0)
    errmsg <- character(0)
    for ( i in ind.timedeppar )
    {
      n.interval.i <- n.interval[1]
      if ( length(n.interval) > 1 )
      {
        if ( !is.na(match(param.names[i],names(n.interval))) )
        {
          n.interval.i <- n.interval[param.names[i]]
        }
        else
        {
          errmsg <- c(errmsg,paste("*** infer.timedeppar: control parameter n.interval is not defined for variable",
                                   param.names[i],"(and also not universally)"))
        }
      }
      if ( n.interval.i > max(floor((nrow(param[[i]])-1)/4),1) )
      {
        errmsg <- c(errmsg,paste("*** infer.timedeppar: control parameter n.interval (=",n.interval.i,">",
                                 max(floor((nrow(param[[i]])-1)/4),1),
                                 ") is too large (number of time points = ",nrow(param[[i]]),") for parameter ",
                                 param.names[i],sep=""))
      }
      ind.oupar[[param.names[i]]] <- numeric(0)
      poss.names <- c(poss.names,paste(param.names[i],c("mean","sd","gamma"),sep="_"))
      if ( paste(param.names[i],"mean",sep="_") %in% param.ou.names ) 
      {
        par.mean <- param.ou[[paste(param.names[i],"mean",sep="_")]]
        ind.oupar[[param.names[i]]] <- c(ind.oupar[[param.names[i]]],match(paste(param.names[i],"mean",sep="_"),param.ou.names))
      }
      else
      {
        if ( paste(param.names[i],"mean",sep="_") %in% names(param.ou.fixed) ) 
        {
          par.mean <- param.ou.fixed[paste(param.names[i],"mean",sep="_")]
        }
        else
        {
          errmsg <- c(errmsg,paste("*** infer.timedeppar: initial or fixed value for parameter",paste(param.names[i],"mean",sep="_"),"is missing"))
        }
      }
      if ( paste(param.names[i],"sd",sep="_") %in% param.ou.names ) 
      {
        par.sd <- param.ou[[paste(param.names[i],"sd",sep="_")]]
        ind.oupar[[param.names[i]]] <- c(ind.oupar[[param.names[i]]],match(paste(param.names[i],"sd",sep="_"),param.ou.names))
      }
      else
      {
        if ( paste(param.names[i],"sd",sep="_") %in% names(param.ou.fixed) ) 
        {
          par.sd <- param.ou.fixed[paste(param.names[i],"sd",sep="_")]
        }
        else
        {
          errmsg <- c(errmsg,paste("*** infer.timedeppar: initial or fixed value for parameter",paste(param.names[i],"sd",sep="_"),"is missing"))
        }
      }
      if ( paste(param.names[i],"gamma",sep="_") %in% param.ou.names ) 
      {
        par.gamma <- param.ou[[paste(param.names[i],"gamma",sep="_")]]
        ind.oupar[[param.names[i]]] <- c(ind.oupar[[param.names[i]]],match(paste(param.names[i],"gamma",sep="_"),param.ou.names))
      }
      else
      {
        if ( paste(param.names[i],"gamma",sep="_") %in% names(param.ou.fixed) ) 
        {
          par.gamma <- param.ou.fixed[paste(param.names[i],"gamma",sep="_")]
        }
        else
        {
          errmsg <- c(errmsg,paste("*** infer.timedeppar: initial or fixed value for parameter",paste(param.names[i],"gamma",sep="_"),"is missing"))
        }
      }
      if ( length(errmsg) == 0 )  # initialize timd-dependent parameters:
      {
        if ( anyNA(param[[i]][,2]) )
        {
          param[[i]] <- randOU(mean  = par.mean,
                               sd    = par.sd,
                               gamma = par.gamma,
                               t     = param[[i]][,1],
                               yini  = NA,
                               yend  = NA,
                               log   = param.log.loc[i])
        }
        if ( param.names[i] %in% names(param.range) )
        {
          val   <- param[[i]][,2]
          range <- param.range[[param.names[i]]]
          param[[i]][,2] <- ifelse(val<range[1],range[1],ifelse(val>range[2],range[2],val))
        }
      }
    }
    n.interval.all <- rep(n.interval[1],length(ind.timedeppar))
    if ( length(n.interval) > 1 ) n.interval.all <- n.interval[param.names[ind.timedeppar]]  # existence checked above
    if ( sum(! param.ou.names %in% poss.names ) > 0 )
    {
      errmsg <- c(errmsg,paste("*** infer.timedeppar: undefined name(s)",paste(param.ou.names[! param.ou.names %in% poss.names]),collapse=", ",
                               "in argument param.ou.ini"))
    }
    res$param.ini <- param  # update with initialized time-dependent parameters

    if ( task[1] == "start" )
    {
      if ( !is.list(cov.prop.ou) )
      {
        cov.prop.ou <- list()
        if ( is.na(scale.prop.ou[1]) ) scale.prop.ou <- rep(0.05,length(ind.timedeppar))
        for ( i in 1:length(ind.timedeppar) )
        {
          cov.prop.ou[[i]] <- diag(1,length(ind.oupar[[i]]))
          if ( length(ind.oupar[[i]]) > 0 )
          {
            ii <- match(paste(param.names[ind.timedeppar[i]],"mean",sep="_"),param.ou.names[ind.oupar[[i]]])
            if ( ! is.na(ii) )
            {
              if ( ! param.log.loc[param.names[ind.timedeppar[i]]] & param.ou[ind.oupar[[i]][ii]] != 0 )
              {
                cov.prop.ou[[i]][ii,ii] <- (0.1*param.ou[ind.oupar[[i]][ii]])^2
              }
            }
            colnames(cov.prop.ou[[i]]) <- param.ou.names[ind.oupar[[i]]]
            rownames(cov.prop.ou[[i]]) <- param.ou.names[ind.oupar[[i]]]
          }
        }
      }
    }
    if ( length(cov.prop.ou) != length(ind.timedeppar) )
    {
      errmsg <- c(errmsg,"*** infer.timedeppar: incorrect length of list \"cov.prop.ou\"")
    }
    for ( i in 1:length(ind.timedeppar) )
    {
      if( nrow(cov.prop.ou[[i]])!=length(ind.oupar[[i]]) | ncol(cov.prop.ou[[i]])!=length(ind.oupar[[i]]) )
      {
        errmsg <- c(errmsg,"*** infer.timedeppar: incorrect dimension of covariance matrix in list \"cov.prop.ou\"")
      }
    }
    if ( length(scale.prop.ou) != length(ind.timedeppar) )
    {
      errmsg <- c(errmsg,"*** infer.timedeppar: incorrect length of vector \"scale.prop.ou\"")
    }
    if ( length(interval.weights) != 0 )
    {
      if ( is.list(interval.weights) )  # specific weights per time-dependent parameter
      {
        for ( i in ind.timedeppar )
        {
          if ( length(interval.weights[[param.names[i]]]) != nrow(param[[i]]) )
          {
            errmsg <- c(errmsg,"*** infer.timedeppar: length of weights vector does not match defintion of parameter",param.names[i])
          }
          if ( splitmethod == "autoweights" )
          {
            interval.weights[[param.names[i]]] <- interval.weights[[param.names[i]]]/max(interval.weights[[param.names[i]]])
          }
        }
      }
      else  # universal weights for all parameters
      {
        for ( i in ind.timedeppar )
        {
          if ( length(interval.weights) != nrow(param[[i]]) )
          {
            errmsg <- c(errmsg,"*** infer.timedeppar: length of weights vector does not match defintion of parameter",param.names[i])
          }
        }
        if ( splitmethod == "autoweights" )
        {
          weights <- list()
          for ( i in ind.timedeppar ) weights[[param.names[i]]] <- interval.weights/max(interval.weights)
          interval.weights <- weights
        }
      }
    }
    else
    {
      if ( splitmethod == "weighted" )
      {
        errmsg <- c(errmsg,"*** infer.timedeppar: weights must be provided for control$splitmethod = \"weighted\"")
      }
      if ( splitmethod == "autoweights" )  # generate initial uniform weights for option "autoweights"
      {
        interval.weights <- list()
        for ( i in ind.timedeppar )
        {
          interval.weights[[param.names[i]]] <- rep(1,nrow(param[[i]]))
        }
      }
    }
    if ( length(errmsg) > 0 )
    {
      for ( i in 1:length(errmsg) ) warning(errmsg[i])
      res$errmsg   <- errmsg
      res$sys.time <- proc.time()-start.time
      return(res)
    }
    for ( i in 1:length(ind.timedeppar) )
    {
      if ( length(ind.oupar[[i]]) > 0 )
      {
        if ( verbose > 0 )
        {
          cat("* proposal distribution of process parameters of time-dependent parameter",param.names[ind.timedeppar[i]],":\n")
          cat("correlation matrix:\n")
          print(cov2cor(cov.prop.ou[[i]]))
          cat("standard deviations:\n")
          print(sqrt(diag(cov.prop.ou[[i]])))
          if ( param.log.loc[param.names[ind.timedeppar[i]]]) cat("logarithmic Ornstein-Uhlenbeck process\n")
          else                                                cat("non-logarithmic Ornstein-Uhlenbeck process\n")
          cat("scaling factor:\n",scale.prop.ou[i],"\n")
        }
      }
    }
  }
  scale.prop.const.ref <- scale.prop.const
  scale.prop.ou.ref    <- scale.prop.ou
  
  # calculate initial log prior, logpdfOU and log likelihood:
  # ---------------------------------------------------------
  
  # log prior
  
  errmsg <- character(0)
  logprior.const <- 0
  if ( length(ind.constpar) > 0 ) logprior.const <- res$param.logprior(unlist(param[ind.constpar]))
  if ( !is.finite(logprior.const) ) errmsg <- c(errmsg,"*** infer.timedeppar: initial value of prior of constant parameters is not finite")
  logprior.ou <- numeric(0)
  if ( length(ind.timedeppar) > 0 )
  {
    logprior.ou <- rep(0,length(ind.timedeppar))
    names(logprior.ou) <- param.names[ind.timedeppar]
    for ( i in ind.timedeppar )
    {
      if ( length(ind.oupar[[param.names[i]]]) > 0 )
      {
        logprior.ou[param.names[i]] <- res$param.ou.logprior(param.ou[ind.oupar[[param.names[i]]]])
        if ( !is.finite(logprior.ou[param.names[i]]) ) errmsg <- c(errmsg,paste("*** infer.timedeppar: initial value of prior of parameters of Ornstein-Uhlenbeck parameters for variable",
                                                                                param.names[i],"is not finite"))
      }
    }
  }
  
  # pdf of process:
    
  logpdfou <- numeric(0)
  if ( length(ind.timedeppar) > 0 )
  {
    logpdfou <- rep(0,length(ind.timedeppar))
    for ( i in 1:length(ind.timedeppar) )
    {
      if ( paste(param.names[ind.timedeppar[i]],"mean",sep="_") %in% param.ou.names )  par.mean  <- param.ou[[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]]
      else                                                                             par.mean  <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]
      if ( paste(param.names[ind.timedeppar[i]],"sd",sep="_") %in% param.ou.names )    par.sd    <- param.ou[[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]]
      else                                                                             par.sd    <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]
      if ( paste(param.names[ind.timedeppar[i]],"gamma",sep="_") %in% param.ou.names ) par.gamma <- param.ou[[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]]
      else                                                                             par.gamma <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]
      logpdfou[i] <- logpdfOU(param[[ind.timedeppar[i]]][,1],param[[ind.timedeppar[i]]][,2],mean=par.mean,sd=par.sd,gamma=par.gamma,cond=0,log=param.log.loc[ind.timedeppar[i]])
    }
  }
  
  # observational likelihood:
  
  if ( !res$loglikeli.keepstate )
  {
    loglikeliobs <- res$loglikeli(param,...)
  }
  else
  {
    loglikeli.list <- res$loglikeli(param,NA,NA,...)
    loglikeliobs   <- loglikeli.list[[1]]
    loglikelistate <- loglikeli.list[[2]]
  }
  if ( !is.finite(loglikeliobs) ) errmsg <-c(errmsg,"*** infer.timedeppar: initial value of likelihood function is not finite")

  if ( length(errmsg) > 0 )
  {
    errmsg <- c(errmsg,"*** infer.timedeppar: unable to start Markov chain")
    for ( i in 1:length(errmsg) ) warning(errmsg[i])
    res$errmsg   <- errmsg
    res$sys.time <- proc.time()-start.time
    return(res)
  }

  # allocate and initialize sample matrices:
  # ----------------------------------------

  n.sample <- floor((n.iter+0.1)/thin)   # samples to be stored (in addition to initial state)

  sample.param.timedep <- list()
  if ( length(ind.timedeppar) > 0 )
  {
    for ( i in ind.timedeppar )  # rows: points in time, initial state, n.sample sample points
    {
      sample.param.timedep[[param.names[i]]] <- rbind(t(param[[param.names[i]]][,1]),
                                                      matrix(NA,nrow=n.sample+1,ncol=nrow((param[[param.names[i]]]))))
    }
  }
  sample.param.ou              <- matrix(NA,nrow=n.sample+1,ncol=length(param.ou))  # rows: initial state, n.sample sample points
  colnames(sample.param.ou)    <- param.ou.names
  sample.param.const           <- matrix(NA,nrow=n.sample+1,ncol=length(ind.constpar))  # rows: initial state, n.sample sample points
  colnames(sample.param.const) <- param.names[ind.constpar]
  sample.logpdf                <- matrix(NA,nrow=n.sample+1,ncol=2+ifelse(length(ind.constpar)>0,1,0)+2*length(ind.timedeppar))  # rows: initial state, n.sample sample points
  cnames                       <- c("logposterior","loglikeliobs")
  if ( length(ind.constpar) > 0 ) cnames <- c(cnames,"logprior_constpar")
  if ( length(ind.timedeppar) > 0 ) cnames <- c(cnames,
                                                paste("logprior_oupar",param.names[ind.timedeppar],sep="_"),
                                                paste("logpdfou_timedeppar",param.names[ind.timedeppar],sep="_"))
  colnames(sample.logpdf)      <- cnames
  if ( save.diag & length(ind.timedeppar) > 0 )
  {
    sample.diag <- list()
    for ( i in ind.timedeppar )  # rows: sample points (initial state not needed for diagnostics, time points not stored in diag arrays)
    {
      diag <- list()
      diag[["proposal"]]       <- matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]])))
      diag[["logacceptratio"]] <- matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]])))
      diag[["accepted"]]       <- matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]])))
      diag[["internalpts"]]    <- matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]])))
      sample.diag[[param.names[i]]] <- diag
    }
  }
  i.sample <- 1
  if ( task[1] == "continue" )
  {
    i.sample <- floor((res.infer.timedeppar$n.iter+0.1)/res.infer.timedeppar$control$thin) + 1
    if ( length(ind.timedeppar) >  0 )
    {
      for ( i in ind.timedeppar )
      {
        sample.param.timedep[[param.names[i]]] <- rbind(res.infer.timedeppar$sample.param.timedep[[param.names[i]]],
                                                        matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]]))))
        if ( save.diag )
        {
          sample.diag[[param.names[i]]][["proposal"]]       <- rbind(res.infer.timedeppar$sample.diag[[param.names[i]]][["proposal"]],
                                                                     matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]]))))
          sample.diag[[param.names[i]]][["logacceptratio"]] <- rbind(res.infer.timedeppar$sample.diag[[param.names[i]]][["logacceptratio"]],
                                                                     matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]]))))
          sample.diag[[param.names[i]]][["accepted"]]       <- rbind(res.infer.timedeppar$sample.diag[[param.names[i]]][["accepted"]],
                                                                     matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]]))))
          sample.diag[[param.names[i]]][["internalpts"]]    <- rbind(res.infer.timedeppar$sample.diag[[param.names[i]]][["internalpts"]],
                                                                     matrix(NA,nrow=n.sample,ncol=nrow((param[[param.names[i]]]))))
        }
      }
    }
    sample.param.const   <- rbind(res.infer.timedeppar$sample.param.const,sample.param.const[-1,,drop=FALSE])
    sample.param.ou      <- rbind(res.infer.timedeppar$sample.param.ou   ,sample.param.ou[-1,,drop=FALSE])
    sample.logpdf        <- rbind(res.infer.timedeppar$sample.logpdf     ,sample.logpdf[-1,,drop=FALSE])
  }
  if ( length(ind.timedeppar) > 0 )
  {
    for ( i in ind.timedeppar )
    {
      sample.param.timedep[[param.names[i]]][1+i.sample,] <- param[[param.names[i]]][,2]
    }
    if ( length(param.ou) > 0 )
    {
      sample.param.ou[i.sample,] <- param.ou
    }
  }
  if ( length(ind.constpar) > 0 )
  {
    sample.param.const[i.sample,] <- unlist(param[ind.constpar])
  }
  sample.logpdf[i.sample,"logposterior"] <- logprior.const + sum(logprior.ou)+sum(logpdfou)+loglikeliobs
  sample.logpdf[i.sample,"loglikeliobs"] <- loglikeliobs
  if ( length(ind.constpar) > 0 ) sample.logpdf[i.sample,"logprior_constpar"] <- logprior.const
  if ( length(ind.timedeppar) > 0 )
  {
    sample.logpdf[i.sample,paste("logprior_oupar",param.names[ind.timedeppar],sep="_")] <- logprior.ou 
    sample.logpdf[i.sample,paste("logpdfou_timedeppar",param.names[ind.timedeppar],sep="_")] <- logpdfou 
  }

  # loop over elements of Markov chain:
  # -----------------------------------

  n.accept.timedeppar           <- rep(NA,length(ind.timedeppar))
  n.accept.oupar                <- rep(NA,length(ind.timedeppar))
  n.accept.constpar             <- NA

  n.accept.oupar.sinceupdate    <- rep(0,length(ind.timedeppar))
  n.iter.oupar.sinceupdate      <- 0
  n.accept.constpar.sinceupdate <- 0
  n.iter.constpar.sinceupdate   <- 0
  
  adaptation.completed          <- FALSE

  if ( task[1] == "continue" )
  {
    n.iter.afteradapt   <- res.infer.timedeppar$n.iter - res.infer.timedeppar$control$n.adapt
    n.accept.timedeppar <- res.infer.timedeppar$acceptfreq.timedeppar*n.iter.afteradapt*n.interval.all
    n.accept.oupar      <- res.infer.timedeppar$acceptfreq.oupar*n.iter.afteradapt*n.const.perstep
    n.accept.constpar   <- res.infer.timedeppar$acceptfreq.constpar*n.iter.afteradapt*n.const.perstep
  }
  names(n.accept.timedeppar) <- param.names[ind.timedeppar]
  names(n.accept.oupar)      <- param.names[ind.timedeppar]
  start.iter <- 1
  if ( task[1] == "continue" ) 
  {
    n.iter <- (nrow(res.infer.timedeppar$sample.logpdf)-1)*thin + n.iter
    start.iter <- (nrow(res.infer.timedeppar$sample.logpdf)-1)*thin + 1
  }

  # write initial log posterior values:
  
  if ( verbose > 0 )
  {
    cat("iter = ",ifelse(task[1]=="continue",start.iter-1,0),": ",sep="")
    if ( length(ind.constpar) > 0 )
    {
      cat(paste(paste(param.names[ind.constpar],signif(unlist(param[ind.constpar]),4),sep="="),collapse=", "))
      if ( length(ind.timedeppar) > 0 ) cat(", ")
    }
    if ( length(ind.timedeppar) > 0 )
    {
      cat(paste(paste(param.ou.names,signif(param.ou,4),sep="="),collapse=", "))
    }
    cat("\n  ",
        "logpost=",logprior.const + sum(logprior.ou)+sum(logpdfou)+loglikeliobs,
        ", logprior=",logprior.const + sum(logprior.ou),
        ", logpdfou=",sum(logpdfou),
        ", loglikeli=",loglikeliobs,"\n",sep="")
  }
  
  for ( k in start.iter:n.iter )
  {
    # inference step for time dependent parameters:

    if ( length(ind.timedeppar) > 0 )
    {
      for ( ii in 1:n.timedep.perstep )
      {
        for ( i in ind.timedeppar )
        {
          n.interval.i <- n.interval[1]
          if ( length(n.interval) > 1 )
          {
            n.interval.i <- n.interval[param.names[i]]  # existence has been tested above
          }
          if ( paste(param.names[i],"mean",sep="_") %in% param.ou.names )  par.mean  <- param.ou[[paste(param.names[i],"mean",sep="_")]]
          else                                                             par.mean  <- param.ou.fixed[paste(param.names[i],"mean",sep="_")]
          if ( paste(param.names[i],"sd",sep="_") %in% param.ou.names )    par.sd    <- param.ou[[paste(param.names[i],"sd",sep="_")]]
          else                                                             par.sd    <- param.ou.fixed[paste(param.names[i],"sd",sep="_")]
          if ( paste(param.names[i],"gamma",sep="_") %in% param.ou.names ) par.gamma <- param.ou[[paste(param.names[i],"gamma",sep="_")]]
          else                                                             par.gamma <- param.ou.fixed[paste(param.names[i],"gamma",sep="_")]
          #if ( splitmethod != "weighted" & splitmethod != "autoweights" )  # non-weighted interval generation
          if ( splitmethod == "modunif" | splitmethod == "random" )  # non-weighted interval generation
          {
            split <- randsplit(n.grid=nrow(param[[i]]),n.interval=n.interval.i,method=splitmethod)
          }
          else  # aweighted interval generation
          {
            # update weights if needed:
            
            if ( splitmethod == "autoweights" & k %% n.autoweighting == 0 & i.sample > 10 )
            {
              sample.ind <- max(2,i.sample-floor(n.autoweighting/thin)+1):i.sample
              sample <- round(sample.param.timedep[[param.names[i]]][sample.ind,],6)
              n.t <- ncol(sample)
              acceptfreq <- rep(NA,n.t)
              for ( j in 1:n.t ) acceptfreq[j] <- (length(unique(sample[,j]))-1)/(n.t-1)
              afmin <- min(acceptfreq)
              afmax <- max(acceptfreq)
              w <- rep(1,n.t)
              if ( afmax > afmin )
              {
                af <- (acceptfreq-afmin)/(afmax-afmin)
                w1 <- offset.weighting*(1+offset.weighting*af)/(offset.weighting+af)
                for ( j in 1:n.t ) w[j] <- max(w1[max(1,j-n.widening):min(n.t,j+n.widening)])  # widen ranges of large values
              }
              
              interval.weights[[param.names[i]]] <- ifelse(w>interval.weights[[param.names[i]]],w,0.5*(w+interval.weights[[param.names[i]]]))
              
              # renormalize and bound below:
              interval.weights[[param.names[i]]] <- interval.weights[[param.names[i]]]/max(interval.weights[[param.names[i]]])
              interval.weights[[param.names[i]]] <- ifelse(interval.weights[[param.names[i]]]<offset.weighting,offset.weighting,interval.weights[[param.names[i]]])
            }
            
            # extract weights and generate intervals:
            
            weights <- numeric(0)
            if ( is.list(interval.weights) )  # specific weights per time-dependent parameter
            {
              weights <- interval.weights[[param.names[i]]]
            }
            else  # universal weights for all parameters
            {
              weights <- interval.weights
            }
            split <- randsplit(n.grid=nrow(param[[i]]),n.interval=n.interval.i,method="weighted",
                               weights=weights,offset=k+ii,min.internal=min.internal)
          }
          
          # loop over intervals to generate and accept or reject proposals:
          
          n.accept.timedeppar.singlerun <- 0
          n.accept.with.ratio.lt.1      <- 0
          for ( j in 1:n.interval.i )
          {
            param.prop <- param
            yini <- ifelse(j==1,NA,param.prop[[i]][split[j],2])
            yend <- ifelse(j==n.interval.i,NA,param.prop[[i]][split[j+1],2])
            param.prop[[i]][split[j]:split[j+1],2]  <-
              randOU(mean  = par.mean,
                     sd    = par.sd,
                     gamma = par.gamma,
                     t     = param.prop[[i]][split[j]:split[j+1],1],
                     yini  = yini,
                     yend  = yend,
                     log   = param.log.loc[i])[,2]
            if ( param.names[i] %in% names(param.range) )
            {
              val   <- param.prop[[i]][,2]
              range <- param.range[[param.names[i]]]
              param.prop[[i]][,2] <- ifelse(val<range[1],range[1],ifelse(val>range[2],range[2],val))
            }
            t.mod <- NA; if ( j > 1 ) t.mod <- param.prop[[i]][c(split[j],split[j+1]),1]
            if ( !res$loglikeli.keepstate )
            {
              loglikeliobs.prop <- res$loglikeli(param.prop,...)
            }
            else
            {
              loglikeli.list <- res$loglikeli(param.prop,t.mod,loglikelistate,...)
              loglikeliobs.prop   <- loglikeli.list[[1]]
              loglikelistate.prop <- loglikeli.list[[2]]
            }
            r <- exp(loglikeliobs.prop-loglikeliobs)
            if ( verbose > 1 )
            {
              if ( j == 1 ) cat("  sampling time-dependent parameter \"",param.names[i],"\"\n",sep="")
              cat("    interval ",j,": loglikeliobs = ",loglikeliobs.prop," (",loglikeliobs,") r =",r,sep="")
            }
            r.random <- runif(1)
            if ( save.diag & k %% thin == 0 )
            {
              internalind <- (split[j]+1):(split[j+1]-1)
              if ( j==1 )            internalind <- split[j]:(split[j+1]-1)
              if ( j==n.interval.i ) internalind <- (split[j]+1):split[j+1]
              sample.diag[[param.names[i]]][["proposal"]][i.sample,split[j]:split[j+1]] <- param.prop[[i]][split[j]:split[j+1],2]
              sample.diag[[param.names[i]]][["logacceptratio"]][i.sample,internalind]   <- loglikeliobs.prop-loglikeliobs
              sample.diag[[param.names[i]]][["accepted"]][i.sample,internalind]         <- r > r.random
              sample.diag[[param.names[i]]][["internalpts"]][i.sample,internalind]      <- length(internalind)
            }
            if ( r > r.random )
            {
              param[[i]][,2]                      <- param.prop[[i]][,2]
              loglikeliobs                        <- loglikeliobs.prop
              if ( res$loglikeli.keepstate ) loglikelistate <- loglikelistate.prop
              n.accept.timedeppar[param.names[i]] <- n.accept.timedeppar[param.names[i]] + 1
              n.accept.timedeppar.singlerun <- n.accept.timedeppar.singlerun + 1
              if ( r < 1 ) n.accept.with.ratio.lt.1 <- n.accept.with.ratio.lt.1 + 1
              if ( verbose > 1 ) cat(" - proposal accepted\n")
            }
            else
            {
              if ( verbose > 1 ) cat(" - proposal rejected\n")
            }
          }
          if ( verbose > 1 )
          {
            cat("  proposal for time-dependent parameter \"",param.names[i],
                "\" was accepted in ",n.accept.timedeppar.singlerun," (",n.accept.with.ratio.lt.1," with r.acc<1) out of ",n.interval.i," intervals\n",sep="")
          }
        }
      }

      # inference step for Ornstein-Uhlenbeck parameters:
    
      if ( length(param.ou) > 0 )
      {
        for ( j in 1:n.const.perstep )
        {
          for ( i in 1:length(ind.timedeppar) )
          {
            if ( length(ind.oupar[[i]]) > 0 )
            {
              if ( paste(param.names[ind.timedeppar[i]],"mean",sep="_") %in% param.ou.names )  par.mean  <- param.ou[[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]]
              else                                                                             par.mean  <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]
              if ( paste(param.names[ind.timedeppar[i]],"sd",sep="_") %in% param.ou.names )    par.sd    <- param.ou[[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]]
              else                                                                             par.sd    <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]
              if ( paste(param.names[ind.timedeppar[i]],"gamma",sep="_") %in% param.ou.names ) par.gamma <- param.ou[[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]]
              else                                                                             par.gamma <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]
              logpdfou[i] <- logpdfOU(param[[ind.timedeppar[i]]][,1],param[[ind.timedeppar[i]]][,2],mean=par.mean,sd=par.sd,gamma=par.gamma,cond=0,log=param.log.loc[ind.timedeppar[i]])
              prior.log.corr      <- 0
              prior.log.corr.prop <- 0
              param.withinrange <- TRUE
              param.ou.prop <- param.ou[ind.oupar[[i]]]
              step <- rmvnorm(1,mean=rep(0,length(ind.oupar[[i]])),sigma=cov.prop.ou[[i]]*scale.prop.ou[i]^2)
              for ( ii in 1:length(ind.oupar[[i]]) )
              {
                name <- param.names[ind.timedeppar[i]]
                if ( !param.log.loc[ind.timedeppar[i]] & param.ou.names[ind.oupar[[i]][ii]]==paste(name,"mean",sep="_") )
                {
                  param.ou.prop[ii]     <- param.ou[ind.oupar[[i]][ii]] + step[ii]
                }
                else
                {
                  param.ou.prop[ii]   <- exp( log(param.ou[ind.oupar[[i]][ii]]) + step[ii] )
                  prior.log.corr      <- prior.log.corr      + log(param.ou[ind.oupar[[i]][ii]])
                  prior.log.corr.prop <- prior.log.corr.prop + log(param.ou.prop[ii])
                }
                if ( param.ou.names[ind.oupar[[i]][ii]] %in% names(param.range) )
                {
                  range <- param.range[[param.ou.names[ind.oupar[[i]][ii]]]]
                  if ( param.ou.prop[ii] <= range[1] | param.ou.prop[ii] >= range[2] )
                  {
                    param.withinrange <- FALSE
                    break
                  }
                }
              }
              if ( verbose > 1 )
              {
                cat("  proposal for process parameters for time-dependent parameter \"",param.names[ind.timedeppar[i]],"\" (scale = ",scale.prop.ou[i],"):\n",sep="")
                cat("   ",paste(paste(param.ou.names[ind.oupar[[i]]],param.ou.prop,sep="="),collapse=", "),"\n")
              }
              if ( param.withinrange )
              {
                logprior.ou.prop <- res$param.ou.logprior(param.ou.prop)
                if ( verbose > 1 )
                {
                  cat("    ","logprior = ",logprior.ou.prop," (",logprior.ou[i],")",sep="")
                  cat("  prior.logcorr = ",prior.log.corr.prop," (",prior.log.corr,")",sep="")
                }
                if ( is.finite(logprior.ou.prop) )
                {
                  if ( paste(param.names[ind.timedeppar[i]],"mean",sep="_") %in% param.ou.names )  par.mean  <- param.ou.prop[[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]]
                  else                                                                             par.mean  <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]
                  if ( paste(param.names[ind.timedeppar[i]],"sd",sep="_") %in% param.ou.names )    par.sd    <- param.ou.prop[[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]]
                  else                                                                             par.sd    <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]
                  if ( paste(param.names[ind.timedeppar[i]],"gamma",sep="_") %in% param.ou.names ) par.gamma <- param.ou.prop[[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]]
                  else                                                                             par.gamma <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]
                  logpdfou.prop <- logpdfOU(param[[ind.timedeppar[i]]][,1],param[[ind.timedeppar[i]]][,2],mean=par.mean,sd=par.sd,gamma=par.gamma,cond=0,log=param.log.loc[ind.timedeppar[i]])
                  if ( verbose > 1 )
                  {
                    cat(" logpdfOU = ",logpdfou.prop," (",logpdfou[i],")",sep="")
                  }
                  r <- exp( (logpdfou.prop + logprior.ou.prop + prior.log.corr.prop) -
                              (logpdfou[i]   + logprior.ou[i]   + prior.log.corr) )
                  if ( r > runif(1) )
                  {
                    param.ou[ind.oupar[[i]]]      <- param.ou.prop
                    logprior.ou[i]                <- logprior.ou.prop
                    logpdfou[i]                   <- logpdfou.prop
                    n.accept.oupar[i]             <- n.accept.oupar[i] + 1
                    n.accept.oupar.sinceupdate[i] <- n.accept.oupar.sinceupdate[i] + 1
                    if ( verbose > 1 )
                    {
                      cat("\n    proposal accepted (r = ",r,")\n",sep="")
                    }
                  }
                  else
                  {
                    if ( verbose > 1 )
                    {
                      cat("\n    proposal rejected (r = ",r,")\n",sep="")
                    }
                  }
                }
                else
                {
                  warning(paste("log prior of process parameters =",logprior.ou.prop,
                                "at",paste(paste(param.ou.names[ind.oupar[[i]]],signif(param.ou.prop),sep="="),collapse=", "),"\n"))
                  if ( verbose > 1 )
                  {
                    cat("    proposal rejected as log prior is not finite\n")
                  }
                }
              }
              else
              {
                if ( verbose > 1 )
                {
                  cat("    proposal rejected as it is out of range\n")
                }
              }
            }
          }
          n.iter.oupar.sinceupdate <- n.iter.oupar.sinceupdate + 1
        }
      }
      else   # calculate logpdfou in case of given Ornstein-Uhlenbeck parameters:
      {
        for ( i in 1:length(ind.timedeppar) )
        {
          if ( paste(param.names[ind.timedeppar[i]],"mean",sep="_") %in% param.ou.names )  par.mean  <- param.ou[[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]]
          else                                                                             par.mean  <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"mean",sep="_")]
          if ( paste(param.names[ind.timedeppar[i]],"sd",sep="_") %in% param.ou.names )    par.sd    <- param.ou[[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]]
          else                                                                             par.sd    <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"sd",sep="_")]
          if ( paste(param.names[ind.timedeppar[i]],"gamma",sep="_") %in% param.ou.names ) par.gamma <- param.ou[[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]]
          else                                                                             par.gamma <- param.ou.fixed[paste(param.names[ind.timedeppar[i]],"gamma",sep="_")]
          logpdfou[i] <- logpdfOU(param[[ind.timedeppar[i]]][,1],param[[ind.timedeppar[i]]][,2],mean=par.mean,sd=par.sd,gamma=par.gamma,cond=0,log=param.log.loc[ind.timedeppar[i]])
        }
      }
    }

    # inference step for constant parameters:
    
    if ( length(ind.constpar) > 0 )
    {
      for ( j in 1:n.const.perstep )
      {
        prior.log.corr      <- 0
        prior.log.corr.prop <- 0
        param.withinrange <- TRUE
        param.prop <- param
        step <- rmvnorm(1,mean=rep(0,length(ind.constpar)),sigma=cov.prop.const*scale.prop.const^2)
        for ( i in 1:length(ind.constpar) )
        {
          if ( !param.log.loc[ind.constpar[i]] ) 
          {
            param.prop[[ind.constpar[i]]] <- param[[ind.constpar[i]]] + step[i]
          }
          else
          {
            param.prop[[ind.constpar[i]]] <- exp( log(param[[ind.constpar[i]]]) + step[i] )
            prior.log.corr                <- prior.log.corr      + log(param[[ind.constpar[i]]])
            prior.log.corr.prop           <- prior.log.corr.prop + log(param.prop[[ind.constpar[i]]])
          }
          if ( param.names[ind.constpar[i]] %in% names(param.range) )
          {
            range <- param.range[[param.names[ind.constpar[i]]]]
            if ( param.prop[[ind.constpar[i]]] <= range[1] | param.prop[[ind.constpar[i]]] >= range[2] )
            {
              param.withinrange <- FALSE
              break
            }
          }
        }
        if ( verbose > 1 )
        {
          cat("  proposal for constant parameters (scale = ",scale.prop.const,"):\n",sep="")
          cat("   ",paste(paste(param.names[ind.constpar],param.prop[ind.constpar],sep="="),collapse=", "),"\n")
        }
        if ( param.withinrange )
        {
          logprior.const.prop <- res$param.logprior(unlist(param.prop[ind.constpar]))
          if ( verbose > 1 )
          {
            cat("    ","logprior = ",logprior.const.prop," (",logprior.const,")",sep="")
            cat("  prior.logcorr = ",prior.log.corr.prop," (",prior.log.corr,")",sep="")
          }
          if ( is.finite(logprior.const.prop) )
          {
            if ( !res$loglikeli.keepstate )
            {
              loglikeliobs.prop <- res$loglikeli(param.prop,...)
            }
            else
            {
              loglikeli.list <- res$loglikeli(param.prop,NA,NA,...)
              loglikeliobs.prop   <- loglikeli.list[[1]]
              loglikelistate.prop <- loglikeli.list[[2]]
            }
            if ( verbose > 1 )
            {
              cat(" loglikeliobs = ",loglikeliobs.prop, " (",loglikeliobs,")",sep="")
            }
            if ( is.finite(loglikeliobs.prop) )
            {
              r <- exp( (loglikeliobs.prop + logprior.const.prop + prior.log.corr.prop) -
                          (loglikeliobs      + logprior.const      + prior.log.corr) )
              if ( r > runif(1) )
              {
                param                         <- param.prop
                logprior.const                <- logprior.const.prop
                loglikeliobs                  <- loglikeliobs.prop
                if ( res$loglikeli.keepstate ) loglikelistate <- loglikelistate.prop
                n.accept.constpar             <- n.accept.constpar + 1
                n.accept.constpar.sinceupdate <- n.accept.constpar.sinceupdate + 1
                if ( verbose > 1 )
                {
                  cat("\n    proposal accepted (r = ",r,")\n",sep="")
                }
              }
              else
              {
                if ( verbose > 1 )
                {
                  cat("\n    proposal rejected (r = ",r,")\n",sep="")
                }
              }
            }
            else
            {
              warning(paste("log likelihood =",loglikeliobs.prop,
                            "at",paste(paste(param.names[ind.constpar],signif(unlist(param.prop[ind.constpar]),3),sep="="),collapse=", ")))
              if ( length(ind.timedeppar) > 0 ) warning(paste(" and timedep parameter(s) ",paste(param.names[ind.timedeppar],collapse=", "),sep=""))
              if ( verbose > 1 )
              {
                cat("    proposal rejected as log likelihood is not finite\n")
              }
            }
          }
          else
          {
            warning(paste("log prior of model parameters =",logprior.const.prop,
                          "at",paste(paste(param.names[ind.constpar],signif(unlist(param.prop[ind.constpar]),3),sep="="),collapse=", ")))
            if ( verbose > 1 )
            {
              cat("    proposal rejected as log prior is not finite\n")
            }
          }
        }
        else
        {
          if ( verbose > 1 )
          {
            cat("    proposal rejected as it is out of range\n")
          }
        }
        n.iter.constpar.sinceupdate <- n.iter.constpar.sinceupdate + 1
      }
    }

    # save new element of chain (consider thinning):

    if ( k %% thin == 0 )
    {
      i.sample <- i.sample + 1
      if ( length(ind.constpar) > 0 )
      {
        for ( i in 1:length(ind.constpar) ) sample.param.const[i.sample,i] <- param[[ind.constpar[i]]]
      }
      if ( length(ind.timedeppar) > 0 )
      {
        sample.param.ou[i.sample,] <- param.ou
        for ( i in ind.timedeppar )
        {
          sample.param.timedep[[param.names[i]]][1+i.sample,] <- param[[param.names[i]]][,2]
        }
      }
      sample.logpdf[i.sample,"logposterior"] <- logprior.const + sum(logprior.ou)+sum(logpdfou)+loglikeliobs
      sample.logpdf[i.sample,"loglikeliobs"] <- loglikeliobs
      if ( length(ind.constpar) > 0 ) sample.logpdf[i.sample,"logprior_constpar"] <- logprior.const
      if ( length(ind.timedeppar) > 0 )
      {
        sample.logpdf[i.sample,paste("logprior_oupar",param.names[ind.timedeppar],sep="_")] <- logprior.ou 
        sample.logpdf[i.sample,paste("logpdfou_timedeppar",param.names[ind.timedeppar],sep="_")] <- logpdfou 
      }

      if ( verbose > 0 )
      {
        cat("iter = ",k,": ",sep="")
        if ( length(ind.constpar) > 0 )
        {
          cat(paste(paste(param.names[ind.constpar],signif(unlist(param[ind.constpar]),4),sep="="),collapse=", "))
          if ( length(ind.timedeppar) > 0 ) cat(", ")
        }
        if ( length(ind.timedeppar) > 0 )
        {
          cat(paste(paste(param.ou.names,signif(param.ou,4),sep="="),collapse=", "))
        }
        cat("\n  ",
            "logpost=",logprior.const + sum(logprior.ou)+sum(logpdfou)+loglikeliobs,
            ", logprior=",logprior.const + sum(logprior.ou),
            ", logpdfou=",sum(logpdfou),
            ", loglikeli=",loglikeliobs,"\n",sep="")
      } 
    }

    # adapt covariance matrix and step size:

    if ( k < n.adapt )
    {
      if ( n.adapt.cov > 0 )   # n.adapt.cov == 0 means no covariance adaptation
      {
        if ( k >= n.adapt.cov & k %% n.adapt.cov == 0 )
        {
          # adapt covariance matrices:
          
          errmsg <- character(0)
          if ( length(ind.constpar) > 0 )
          {
            sample.log <- sample.param.const[ceiling(0.25*i.sample):i.sample,,drop=FALSE]
            for ( i in 1:length(ind.constpar) ) 
            { 
              if ( param.log.loc[ind.constpar[i]] ) sample.log[,i] <- log(sample.log[,i,drop=FALSE]) 
              if ( var(sample.log[,i]) == 0 ) 
              {
                errmsg <- c(errmsg,paste("*** infer.timedeppar: unable to advance parameter",param.names[ind.constpar[i]]))
              }
            }
            if ( length(errmsg) == 0 )  # update only if no errors, otherwise inference will be stopped below
            {
              vol.old <- sqrt(prod(diag(cov.prop.const)))  # expect to be more robust without considering correlation
              cov.prop.const <- cov(sample.log)
              if ( ncol(cov.prop.const) > 1 ) cov.prop.const <- f.reduce.cor*cov.prop.const + (1-f.reduce.cor)*diag(diag(cov.prop.const))
              colnames(cov.prop.const) <- param.names[ind.constpar]
              rownames(cov.prop.const) <- param.names[ind.constpar]
              vol.new <- sqrt(prod(diag(cov.prop.const)))  # expect to be more robust without considering correlation
              scale.prop.const <- scale.prop.const * (vol.old/vol.new)^(1/nrow(cov.prop.const))
              scale.prop.const.ref <- scale.prop.const
              if ( verbose > 0 )
              {
                cat("* new proposal distribution of constant parameters:\n")
                cat("correlation matrix:\n")
                print(cov2cor(cov.prop.const))
                cat("standard deviations:\n")
                print(sqrt(diag(cov.prop.const)))
                cat("log:\n")
                print(param.log.loc[ind.constpar])
                cat("new (reference) scaling factor:\n",scale.prop.const,"\n")
              }
            }
          }
          if ( length(ind.timedeppar) > 0 )
          {
            for ( i in 1:length(ind.timedeppar) )
            {
              if ( length(ind.oupar[[i]]) > 0 )
              {
                name <- param.names[ind.timedeppar[i]]
                sample.log <- sample.param.ou[ceiling(0.25*i.sample):i.sample,ind.oupar[[i]],drop=FALSE]
                for ( ii in 1:length(ind.oupar[[i]]) )
                {
                  if ( ! ( !param.log.loc[ind.timedeppar[i]] & param.ou.names[ind.oupar[[i]][ii]]==paste(name,"mean",sep="_") ) )
                  {
                    sample.log[,ii] <- log(sample.log[,ii,drop=FALSE])
                  }
                  if ( var(sample.log[,ii]) == 0 ) 
                  {
                    errmsg <- c(errmsg,paste("*** infer.timedeppar: unable to advance parameter",param.ou.names[ind.oupar[[i]][ii]]))
                  }
                }
                if ( length(errmsg) == 0 )  # update only if no errors, otherwise inference will be stopped below
                {
                  vol.old <- sqrt(prod(diag(cov.prop.ou[[i]])))  # expect to be more robust without considering correlation
                  cov.prop.ou[[i]] <- cov(sample.log)
                  if ( ncol(cov.prop.ou[[i]]) > 1 ) cov.prop.ou[[i]] <- f.reduce.cor*cov.prop.ou[[i]] + (1-f.reduce.cor)*diag(diag(cov.prop.ou[[i]]))
                  colnames(cov.prop.ou[[i]]) <- param.ou.names[ind.oupar[[i]]]
                  rownames(cov.prop.ou[[i]]) <- param.ou.names[ind.oupar[[i]]]
                  vol.new <- sqrt(prod(diag(cov.prop.ou[[i]])))  # expect to be more robust without considering correlation
                  scale.prop.ou[i] <- scale.prop.ou[i] * (vol.old/vol.new)^(1/nrow(cov.prop.ou[[i]])) 
                  scale.prop.ou.ref[i] <- scale.prop.ou[i]
                  if ( verbose > 0 )
                  {
                    cat("* new proposal distribution of process parameters of time-dependent parameter",param.names[ind.timedeppar[i]],":\n")
                    cat("correlation matrix:\n")
                    print(cov2cor(cov.prop.ou[[i]]))
                    cat("standard deviations:\n")
                    print(sqrt(diag(cov.prop.ou[[i]])))
                    if ( param.log.loc[param.names[ind.timedeppar[i]]]) cat("logarithmic Ornstein-Uhlenbeck process\n")
                    else                                                cat("non-logarithmic Ornstein-Uhlenbeck process\n")
                    cat("new (reference) scaling factor:\n",scale.prop.ou[i],"\n")
                  }
                }
              }
            }
          }
          if ( length(errmsg) > 0 )
          {
            errmsg <- c(errmsg,"*** infer.timedeppar: inference stopped ***")
            for ( i in 1:length(errmsg) ) warning(errmsg[i])
            res$errmsg   <- errmsg
            res$sys.time <- proc.time()-start.time
            return(res)
          }
        }
      }

      # adapt step size:

      if ( length(ind.constpar) > 0 )
      {
        if ( n.iter.constpar.sinceupdate >= n.adapt.scale )
        {
          if ( n.accept.constpar.sinceupdate < f.accept.decscale*n.adapt.scale )
          {
            # reduce step size:
  
            scale.prop.const <- 0.8*scale.prop.const
            if ( verbose > 0 )
            {
              cat("* acc. rate of constant parameters: ",n.accept.constpar.sinceupdate,"/",n.adapt.scale,
                  "=",signif(n.accept.constpar.sinceupdate/n.adapt.scale,3)," < ",f.accept.decscale,": ",sep="")
            }
            if ( scale.prop.const >= scale.prop.const.ref/f.max.scalechange )
            {
              if ( verbose > 0 )
              {
                cat("scaling factor decreased to ",signif(scale.prop.const,3),"\n",sep="")
              }
            }
            else
            {
              scale.prop.const <- scale.prop.const.ref/f.max.scalechange
              if ( verbose > 0 )
              {
                cat("scaling factor at minimum: ",signif(scale.prop.const,3),"\n",sep="")
              }
            }
          }
          if ( n.accept.constpar.sinceupdate > f.accept.incscale*n.adapt.scale )
          {
            # increase step size:

            scale.prop.const <- 1.2*scale.prop.const
            if ( verbose > 0 )
            {
              cat("* acc. rate of constant parameters: ",n.accept.constpar.sinceupdate,"/",n.adapt.scale,
                  "=",signif(n.accept.constpar.sinceupdate/n.adapt.scale,3)," > ",f.accept.incscale,": ",sep="")
            }
            if ( scale.prop.const <= scale.prop.const.ref*f.max.scalechange )
            {
              if ( verbose > 0 )
              {
                cat("scaling factor increased to ",signif(scale.prop.const,3),"\n",sep="")
              }
            }
            else
            {
              scale.prop.const <- scale.prop.const.ref*f.max.scalechange
              if ( verbose > 0 )
              {
                cat("scaling factor at maximum: ",signif(scale.prop.const,3),"\n",sep="")
              }
            }
          }
          n.iter.constpar.sinceupdate   <- 0
          n.accept.constpar.sinceupdate <- 0
        }
      }
      if ( length(ind.timedeppar) > 0 )
      {
        if ( n.iter.oupar.sinceupdate >= n.adapt.scale )
        {
          for ( i in 1:length(ind.timedeppar) )
          {
            name <- param.names[ind.timedeppar[i]]
            if ( n.accept.oupar.sinceupdate[i] < f.accept.decscale*n.adapt.scale )
            {
              # reduce step size:
              
              scale.prop.ou[i] <- 0.8*scale.prop.ou[i]
              if ( verbose > 0 )
              {
                cat("* acc. rate of process parameters for ",name,": ",n.accept.oupar.sinceupdate[i],"/",n.adapt.scale,
                    "=",signif(n.accept.oupar.sinceupdate[i]/n.adapt.scale,3)," < ",f.accept.decscale,": ",sep="")
              }
              if ( scale.prop.ou[i] >= scale.prop.ou.ref[i]/f.max.scalechange )
              {
                if ( verbose > 0 )
                {
                  cat("scaling factor decreased to ",signif(scale.prop.ou[i],3),"\n",sep="")
                }
              }
              else
              {
                scale.prop.ou[i] <- scale.prop.ou.ref[i]/f.max.scalechange
                if ( verbose > 0 )
                {
                  cat("scaling factor at minimum: ",signif(scale.prop.ou[i],3),"\n",sep="")
                }
              }
            }
            if ( n.accept.oupar.sinceupdate[i] > f.accept.incscale*n.adapt.scale )
            {
              # increase step size:
              
              scale.prop.ou[i] <- 1.2*scale.prop.ou[i]
              if ( verbose > 0 )
              {
                cat("* acc. rate of process parameters for ",name,": ",n.accept.oupar.sinceupdate[i],"/",n.adapt.scale,
                    "=",signif(n.accept.oupar.sinceupdate[i]/n.adapt.scale,3)," > ",f.accept.incscale,": ",sep="")
              }
              if ( scale.prop.ou[i] <= scale.prop.ou.ref[i]*f.max.scalechange )
              {
                if ( verbose > 0 )
                {
                  cat("scaling factor increased to ",signif(scale.prop.ou[i],3),"\n",sep="")
                }
              }
              else
              {
                scale.prop.ou[i] <- scale.prop.ou.ref[i]*f.max.scalechange
                if ( verbose > 0 )
                {
                  cat("scaling factor at maximum: ",signif(scale.prop.ou[i],3),"\n",sep="")
                }
              }
            }
            n.accept.oupar.sinceupdate[i] <- 0
          }
          n.iter.oupar.sinceupdate        <- 0
        }
      }
    }
    else
    {
      if ( adaptation.completed == FALSE )
      {
        n.accept.timedeppar        <- rep(0,length(ind.timedeppar))
        names(n.accept.timedeppar) <- param.names[ind.timedeppar]
        n.accept.oupar             <- rep(0,length(ind.timedeppar))
        names(n.accept.oupar)      <- param.names[ind.timedeppar]
        if ( length(ind.constpar) > 0 ) n.accept.constpar <- 0
        adaptation.completed       <- TRUE
      }
    }

    # save most important variables for recovery of aborted run:

    if ( nchar(file.save) > 0 )
    {
      if ( k %% n.save == 0 )
      {
        n.iter.afteradapt     <- k-n.adapt
        acceptfreq.timedeppar <- n.accept.timedeppar/(n.iter.afteradapt*n.interval.all)
        acceptfreq.oupar      <- n.accept.oupar/(n.iter.afteradapt*n.const.perstep)
        acceptfreq.constpar   <- n.accept.constpar/(n.iter.afteradapt*n.const.perstep)
        ind.maxpost           <- which.max(sample.logpdf[,"logposterior"])
        param.maxpost         <- NA
        param.ou.maxpost      <- NA
        if ( length(ind.maxpost) > 0 )
        {
          param.maxpost         <- param
          param.ou.maxpost      <- param.ou
          if ( length(ind.constpar) > 0 )
          {
            for ( i in 1:length(ind.constpar) )
            {
              param.maxpost[[ind.constpar[i]]] <- as.numeric(sample.param.const[ind.maxpost,i])
            }
          }
          if ( length(ind.timedeppar) > 0 )
          {
            for ( i in ind.timedeppar )
            {
              param.maxpost[[param.names[i]]][,2] <- sample.param.timedep[[param.names[i]]][1+ind.maxpost,]  # first row is time
            }
            param.ou.maxpost <- sample.param.ou[ind.maxpost,]
          }
        }

        res$n.iter <- k
        if ( save.diag  & length(ind.timedeppar) > 0 ) res$sample.diag <- sample.diag
        res$sample.param.timedep <- sample.param.timedep
        if (  length(ind.timedeppar) > 0 )
        {
          for ( i in ind.timedeppar )
          {
            res$sample.param.timedep[[param.names[i]]] <- sample.param.timedep[[param.names[i]]][1:(1+i.sample),,drop=FALSE]
            if ( save.diag )
            {
              res$sample.diag[[param.names[i]]][["proposal"]]       <- sample.diag[[param.names[i]]][["proposal"]][1:(i.sample-1),,drop=FALSE]
              res$sample.diag[[param.names[i]]][["logacceptratio"]] <- sample.diag[[param.names[i]]][["logacceptratio"]][1:(i.sample-1),,drop=FALSE]
              res$sample.diag[[param.names[i]]][["accepted"]]       <- sample.diag[[param.names[i]]][["accepted"]][1:(i.sample-1),,drop=FALSE]
              res$sample.diag[[param.names[i]]][["internalpts"]]    <- sample.diag[[param.names[i]]][["internalpts"]][1:(i.sample-1),,drop=FALSE]
            }
          }
        }
        res$sample.param.ou       <- sample.param.ou[1:i.sample,,drop=FALSE]
        res$sample.param.const    <- sample.param.const[1:i.sample,,drop=FALSE]
        res$sample.logpdf         <- sample.logpdf[1:i.sample,,drop=FALSE]
        res$acceptfreq.constpar   <- acceptfreq.constpar
        res$acceptfreq.oupar      <- acceptfreq.oupar
        res$acceptfreq.timedeppar <- acceptfreq.timedeppar
        res$param.maxpost         <- param.maxpost
        res$param.ou.maxpost      <- param.ou.maxpost
        res$cov.prop.const        <- cov.prop.const
        res$cov.prop.ou           <- cov.prop.ou
        res$scale.prop.const      <- scale.prop.const
        res$scale.prop.ou         <- scale.prop.ou
        res$sys.time              <- proc.time()-start.time
        if ( splitmethod == "autoweights" ) res$control$interval.weights <- interval.weights  # update weights!
        save(res,file=paste(file.save,"RData",sep="."))
      }
    }
  }

  # compile and return results:
  # ---------------------------

  n.iter.afteradapt <- n.iter-n.adapt
  acceptfreq.timedeppar <- n.accept.timedeppar/(n.iter.afteradapt*n.interval.all)
  acceptfreq.oupar      <- n.accept.oupar/(n.iter.afteradapt*n.const.perstep)
  acceptfreq.constpar   <- n.accept.constpar/(n.iter.afteradapt*n.const.perstep)
  ind.maxpost <- which.max(sample.logpdf[,"logposterior"])
  param.maxpost <- NA
  param.ou.maxpost <- NA
  if ( is.finite(ind.maxpost) )
  {
    param.maxpost <- param
    param.ou.maxpost <- param.ou
    if ( length(ind.constpar) > 0 )
    {
      for ( i in 1:length(ind.constpar) )
      {
        param.maxpost[[ind.constpar[i]]] <- as.numeric(sample.param.const[ind.maxpost,i])
      }
    }
    if ( length(ind.timedeppar) > 0 )
    {
      for ( i in ind.timedeppar )
      {
        param.maxpost[[param.names[i]]][,2] <- sample.param.timedep[[param.names[i]]][1+ind.maxpost,]  # first row is time
      }
      param.ou.maxpost <- sample.param.ou[ind.maxpost,]
    }
  }

  if ( task[1] != "continue" ) message(paste(n.iter,"iterations completed"))
  else                         message(paste(n.iter-res.infer.timedeppar$n.iter,"iterations added"))
  if ( length(ind.constpar) > 0 )
  {
    message(paste("  acceptance frequency of constant parameters:            ",signif(acceptfreq.constpar,3)))
  }
  if ( length(ind.timedeppar) > 0 )
  {
    message(paste("  acceptance frequencies of time-dependent parameters:    ",paste(signif(acceptfreq.timedeppar,3),collapse=", ")))
    if ( length(param.ou) > 0 )
    {
      message(paste("  acceptance frequencies of Ornstein-Uhlenbeck parameters:",paste(signif(acceptfreq.oupar,3),collapse=", ")))
    }
  }
  else
  {
    message("  no time-dependent parameters")
  }

  res$n.iter                <- n.iter
  if ( save.diag  & length(ind.timedeppar) > 0 ) res$sample.diag <- sample.diag
  res$sample.param.timedep  <- sample.param.timedep
  res$sample.param.ou       <- sample.param.ou
  res$sample.param.const    <- sample.param.const
  res$sample.logpdf         <- sample.logpdf
  res$acceptfreq.constpar   <- acceptfreq.constpar
  res$acceptfreq.oupar      <- acceptfreq.oupar
  res$acceptfreq.timedeppar <- acceptfreq.timedeppar
  res$param.maxpost         <- param.maxpost
  res$param.ou.maxpost      <- param.ou.maxpost
  res$cov.prop.const        <- cov.prop.const
  res$cov.prop.ou           <- cov.prop.ou
  res$scale.prop.const      <- scale.prop.const
  res$scale.prop.ou         <- scale.prop.ou
  res$sys.time              <- proc.time()-start.time
  if ( splitmethod == "autoweights" ) res$control$interval.weights <- interval.weights  # update weights!
  
  if ( nchar(file.save) > 0 ) save(res,file=paste(file.save,"RData",sep="."))
  
  return(res)
}


