#' Reads an object of type \code{timedeppar} saved to a file by \code{\link{infer.timedeppar}}
#' 
#' This function read a workspace stored by the function \code{\link{infer.timedeppar}} and returns
#' the object of type \code{timedeppar} that contains the intermediate or final results of the 
#' inference process
#' 
#' @param  file         file name of the workspace to be read.
#' 
#' @return object of type \code{timedeppar} containing the intermediate or final results of the inference process 
#'         or a list of length zero if the file was not found of no object called \code{res} of type \code{timedeppar}
#'         was found in the workspace.

readres.timedeppar <- function(file)
{
  env.res = new.env()
  if ( ! file.exists(file) ) { warning(paste("readres.timedeppar: file",file,"not found")); return(list()) }
  load(file,envir=env.res)
  if ( ! "res" %in% names(env.res) ) { warning("readres.timedeppar: no results found"); return(list()) }
  if ( class(env.res$res) != "timedeppar" ) { warning("readres.timedeppar: no results found"); return(list()) }
  return(env.res$res)
}
