#' Construct a plot label from expressions for variables and units
#' 
#' This function produces an expression to label plots from expressions for variables and units and from character strings. 
#' 
#' @param  var.name           name of the variable.
#' @param  labels             named vector of expressions encoding variable names with greek letters, sub- and superscripts.
#' @param  units              named vector of expressions encoding units with sub- and superscripts.
#' @param  t1                 optional text (see below).
#' @param  t2                 optional text (see below).
#' @param  t3                 optional text (see below).
#' 
#' @return expression of a label of the form: t1 label t2 [unit] t3\cr
#'         label and unit are extracted from the vectors \code{labels} and \code{units}
#'         by using the compontents named \code{var.name}


get.label <- function(var.name,labels=NA,units=NA,t1="",t2="",t3="")
{
  if ( var.name %in% names(labels) )
  {
    if ( var.name %in% names(units) ) label <- substitute(t1*label*t2*" ["*unit*"]"*t3,
                                                          lapply(list(label=labels[var.name],unit=units[var.name],
                                                                      t1=as.expression(t1),t2=as.expression(t2),t3=as.expression(t3)),
                                                                 "[[",1))
    else                              label <- substitute(t1*label*t2*t3,
                                                          lapply(list(label=labels[var.name],
                                                                      t1=as.expression(t1),t2=as.expression(t2),t3=as.expression(t3)),
                                                                 "[[",1))
  }
  else
  {
    if ( var.name %in% names(units) ) label <- substitute(t1*label*t2*" ["*unit*"]"*t3,
                                                          lapply(list(label=as.expression(var.name),unit=units[var.name],
                                                                      t1=as.expression(t1),t2=as.expression(t2),t3=as.expression(t3)),
                                                                 "[[",1))
    else                              label <- substitute(t1*label*t2*t3,
                                                          lapply(list(label=as.expression(var.name),
                                                                      t1=as.expression(t1),t2=as.expression(t2),t3=as.expression(t3)),
                                                                 "[[",1))
  }
  return(label)
}

