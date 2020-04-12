### FEM: I include this because it was in the original file gnls.R
###
### Copyright 1997-2003  Jose C. Pinheiro,
###                      Douglas M. Bates <bates@stat.wisc.edu>
### Copyright 2006-2015  The R Core team
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

#' Simulate gnls based on \sQuote{nlme::predict.gnls} function 
#' 
#' @title Simulate fitted values from an object of class \sQuote{gnls}
#' @name simulate_gnls
#' @param object object of class \sQuote{gnls}
#' @param psim parameter simulation level, 0: for fitted values, 1: for simulation from 
#' fixed parameters (assuming a fixed vcov matrix), 2: for simulation considering the 
#' uncertainty in vcov (not implemented yet)
#' @param na.action default \sQuote{na.fail}
#' @param naPattern missing value pattern
#' @details It uses function \sQuote{MASS::mvrnorm} to generate new values for the coefficients
#' of the model using the Variance-Covariance matrix \sQuote{vcov}
#' @export
#' 

simulate_gnls <- function(object, psim = 1, na.action = na.fail, naPattern = NULL, ...){
    ##
    ## method for predict() designed for objects inheriting from class gnls
    ##
    if(!inherits(object, "gnls")) stop("This function is only for 'gnls' objects")
  
    mCall <- object$call
    
    ndata <- eval(object$call$data)
    
    mfArgs <- list(formula =
                     nlme::asOneFormula(formula(object),
                                  mCall$params, naPattern,
                                  omit = c(names(object$plist), "pi",
                                           deparse(nlme::getResponseFormula(object)[[2]]))),
                   data = ndata, na.action = na.action,
                   drop.unused.levels = TRUE)
    
    dataMod <- do.call("model.frame", mfArgs)
    
    ## making sure factor levels are the same as in contrasts
    contr <- object$contrasts
    
    for(i in names(dataMod)) {
      if (inherits(dataMod[,i], "factor") &&
          !is.null(contr[[i]]) && is.matrix(contr[[i]]) ) {
        levs <- levels(dataMod[,i])
        levsC <- dimnames(contr[[i]])[[1]]
        if (any(wch <- is.na(match(levs, levsC)))) {
          stop(sprintf(ngettext(sum(wch),
                                "level %s not allowed for %s",
                                "levels %s not allowed for %s"),
                       paste(levs[wch], collapse = ",")), domain = NA)
        }
        attr(dataMod[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
       }
    }
    N <- nrow(dataMod)
    ##
    ## evaluating the naPattern expression, if any
    ##
    naPat <- if(is.null(naPattern)){
      rep(TRUE, N)
    }else{
      as.logical(eval(stats::asOneSidedFormula(naPattern)[[2]], dataMod))
    } 
    ##
    ## Getting  the plist for the new data frame
    ##
    plist <- object$plist
    pnames <- names(plist)
    if(is.null(params <- eval(object$call$params))){
      params <- list(formula(paste0(paste(pnames, collapse = "+"), "~ 1")))
    }else{ 
      if (!is.list(params)) params <- list(params)
    }
      
    params <- unlist(lapply(params, function(pp){
        if(is.name(pp[[2]])){
          list(pp)
        }else{
        ## multiple parameters on left hand side
        eval(parse(text = paste("list(",
                                paste(paste(all.vars(pp[[2]]), deparse(pp[[3]]), sep = "~"),
                                      collapse = ","),
                                ")")))
      }
    }), recursive=FALSE)
      
    names(params) <- pnames
    
    ## This section is by FEM 2020-04-11
    if(psim == 0){
      prs <- coef(object)  
    }
    
    if(psim == 1){
      prs <- MASS::mvrnorm(n = 1, mu = coef(object), Sigma = vcov(object))
    }
    
    if(psim == 2){
      ## This option will, in the future, sample Sigma first
      stop("not implemented yet")
    }
    
    for(nm in pnames) {
      if (!is.logical(plist[[nm]])) {
        form1s <- stats::asOneSidedFormula(params[[nm]][[3]])
        plist[[nm]] <- model.matrix(form1s, model.frame(form1s, dataMod))
      }
    }
    
    modForm <- nlme::getCovariateFormula(object)[[2]]
    val <- eval(modForm, data.frame(dataMod,
                                    nlme:::getParsGnls(plist, object$pmap, prs, N)))[naPat]
    names(val) <- row.names(ndata)
    lab <- "Simulated values"
    if (!is.null(aux <- attr(object, "units")$y)) {
      lab <- paste(lab, aux)
    }
    attr(val, "label") <- lab
    
    return(val)
}
