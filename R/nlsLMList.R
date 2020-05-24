###                  Create a list of nls 'LM' objects
###
### Copyright 1997-2003  Jose C. Pinheiro,
###                      Douglas M. Bates <bates@stat.wisc.edu>
### Copyright 2006-2016 The R Core team
###
### Modified by Fernando Miguez to use nlsLM for fitting in addition to (optionally) nls (2020-01-08)
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#

#' @title Create a list of nls objects with the option of using nlsLM in addition to nls
#' @name nlsLMList
#' @description This function is a copy of 'nlsList' from the 'nlme' package modified
#' to use the 'nlsLM' function in addition to (optinally) 'nls'. By changing the algorithm argument it is possible
#' to use 'nls' as well
#' @param model either a nonlinear model formula, with the response on the left of a ~ operator and an expression involving parameters, covariates, and a grouping factor separated by the | operator on the right, or a selfStart function. 
#' @param data a data frame
#' @param start list with starting values
#' @param control control list, see \code{\link{nls}}
#' @param level an optional integer specifying the level of grouping to be used when multiple nested levels of grouping are present.
#' @param subset subset of rows to use
#' @param na.action a function that indicates what should happen when the data contain NAs. The default action (na.fail) causes nlsList to print an error message and terminate if there are any incomplete observations.
#' @param algorithm choice of algorithm. Default is 'LM' which uses 'nlsLM' from the \CRANpkg{minpack.lm} package. Other options are: \dQuote{default}, \dQuote{port} and \dQuote{plinear} (nls).
#' @param pool an optional logical value that is preserved as an attribute of the returned value. This will be used as the default for pool in calculations of standard deviations or standard errors for summaries.
#' @param warn.nls logical indicating if nls errors (all of which are caught by tryCatch) should be signalled as a “summarizing” warning.
#' @details See function \code{\link[nlme]{nlsList}} and \code{\link[minpack.lm]{nlsLM}}. This function is a copy of nlsList but with minor changes to use LM instead as the default algorithm. The authors of the original function are Pinheiro and Bates.
#' @author Jose C. Pinheiro and Douglas M. Bates \email{bates@@stat.wisc.edu} wrote the original \code{\link[nlme]{nlsList}}. Fernando E. Miguez made minor changes to use \code{\link[minpack.lm]{nlsLM}} in addition to (optionally) \code{\link{nls}}. R-Core maintains copyright after 2006.
#' @export
#' 

nlsLMList <-
  ## A list of nls objects
  function(model, data, start, control, level, subset, na.action = na.fail,
           algorithm = c("LM","default","port","plinear"),
           pool = TRUE, warn.nls = NA) # Deprecation: will be 'TRUE'
      UseMethod("nlsLMList")

nlsLMList.selfStart <-
  function (model, data, start, control, level, subset, na.action = na.fail,
            algorithm = c("LM","default","port","plinear"),
            pool = TRUE, warn.nls = NA) # Deprecation: will be 'TRUE'
{
    if(algorithm == "LM"){
      if(!requireNamespace("minpack.lm", quietly = TRUE)){
        warning("minpack.lm package is required for this algorithm in this function")
        return(NULL)
      }
    }
    
  mCall <- as.list(match.call())[-1]
  if (!inherits(data, "groupedData")) {
    stop("second argument must be a groupedData object")
  }
  marg <- substitute(model)
  if (mode(marg) != "name") {
    stop("cannot use an anonymous function for the model")
  }
	## Build up a call to the model function
  m <- call(as.character(marg))
  args <- lapply(names(formals(eval(marg))), as.name)
  args[[1]] <- getCovariateFormula(data)[[2]]
  m[1 + seq_along(args)] <- args
  form <- formula(data)
  form[[3]][[2]] <- m
  mCall$model <- form
  do.call("nlsLMList.formula", mCall)
}

#' @title Formula method for nls 'LM' list method
#' @name nlsLMList.formula
#' @description formula method for nlsLMList
#' @param model see \code{\link{nlsList}}
#' @param data see \code{\link{nlsList}}
#' @param start see \code{\link{nlsList}}
#' @param control see \code{\link{nls}}
#' @param level see \code{\link{nlsList}}
#' @param subset see \code{\link{nlsList}}
#' @param na.action see \code{\link{nlsList}}
#' @param algorithm choice of algorithm default is 'LM' which uses 'nlsLM' from the minpack.lm package.
#' @param pool see \code{\link{nlsList}}
#' @param warn.nls see \code{\link{nlsList}}
#' @export
nlsLMList.formula <-
  function(model, data, start = NULL, control, level, subset,
           na.action = na.fail, 
           algorithm = c("LM","default","port","plinear"),
           pool = TRUE,
           warn.nls = NA) # Deprecation: will be 'TRUE'
{
  if (!missing(level) && length(level) > 1)
    stop("multiple levels not allowed")
  ## Deprecation: options(show.error.messages = FALSE) should continue to work for now
  algorithm <- match.arg(algorithm)
  if(is.na(warn.nls <- as.logical(warn.nls)))
    warn.nls <- !identical(FALSE, getOption("show.error.messages"))
  Call <- match.call()
  if (!missing(subset)) {
    data <-
      data[eval(asOneSidedFormula(Call[["subset"]])[[2]], data),, drop = FALSE]
  }
  if (!is.data.frame(data)) data <- as.data.frame(data)
  data <- na.action(data)
  if (is.null(grpForm <- getGroupsFormula(model))) {
    if (inherits(data, "groupedData")) {
      if (missing(level))
        level <- length(getGroupsFormula(data, asList = TRUE))
      groups <- getGroups(data, level = level)[drop = TRUE]
      grpForm <- getGroupsFormula(data)
    } else {
      stop("'data' must be a \"groupedData\" object if 'formula' does not include groups")
    }
  } else {
    if (missing(level))
      level <- length(getGroupsFormula(model, asList = TRUE))
    model <- eval(substitute(Y ~ RHS,
			     list(Y  = model[[2]],
				  RHS= getCovariateFormula(model)[[2]])))
    groups <- getGroups(data, form = grpForm, level = level)[drop = TRUE]
  }
  if (is.null(start) && is.null(attr(data, "parameters"))) {
    ## no starting values
    ## checking for old-style selfStart functions
    FUN <- eval(model[[3]][[1]])
    if (is.function(FUN) && class(FUN) != "selfStart" &&
        !is.null(attr(FUN, "initial"))) {
      stop("old-style self-starting model functions\nare no longer supported.\nNew selfStart functions are available.\nUse\n  SSfpl instead of fpl,\n  SSfol instead of first.order.log,\n  SSbiexp instead of biexp,\n  SSlogis instead of logistic.\nIf writing your own selfStart model, see\n  \"help(selfStart)\"\nfor the new form of the \"initial\" attribute.")
    }
  }

  controlvals <- nls.control()
  if(!missing(control)) controlvals[names(control)] <- control
  val <- lapply(split(data, groups),
		function(dat)
                  tryCatch({
                    data <- as.data.frame(dat)
                    if (is.null(start)) {
                      if(algorithm == "LM"){
                        minpack.lm::nlsLM(model, data = data, control = controlvals)  
                      }else{
                        nls(model, data = data, control = controlvals, algorithm = algorithm)  
                      }
                    } else {
                      if(algorithm == "LM"){
                        minpack.lm::nlsLM(model, data = data, control = controlvals, start = start)  
                      }else{
                        nls(model, data = data, control = controlvals, start = start, algorithm = algorithm)  
                      }
                    }
                  }, error = function(e) e))
  val <- utils::warnErrList(val, warn = warn.nls)
  if (inherits(data, "groupedData")) {
    ## saving labels and units for plots
    attr(val, "units") <- attr(data, "units")
    attr(val, "labels") <- attr(data, "labels")
    attr(val, "outer") <- attr(data, "outer")
  }

  structure(val, class = c("nlsList", "lmList"),
            call = Call,
            dims = list(N = nrow(data), M = length(val)),
            groups = ordered(groups, levels = names(val)),
            origOrder = match(unique(as.character(groups)), names(val)),
            pool = pool,
            groupsForm = grpForm)
}


