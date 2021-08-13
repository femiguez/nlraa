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

#' This function is based on \code{\link[nlme]{predict.gnls}} function 
#' 
#' @title Simulate fitted values from an object of class \code{\link[nlme]{gnls}}
#' @description Simulate values from an object of class gnls. Unequal variances, 
#' as modeled using the \sQuote{weights} option are supported, and there is experimental
#' code for dealing with the \sQuote{correlation} structure. 
#' @name simulate_gnls
#' @param object object of class \code{\link[nlme]{gnls}}
#' @param psim parameter simulation level, 0: for fitted values, 1: for simulation from 
#' fixed parameters (assuming a fixed vcov matrix), 2: for simulation considering the 
#' uncertainty in the residual standard error (sigma), this returns data which
#' will appear similar to the observed values
#' @param na.action default \sQuote{na.fail}. See \code{\link[nlme]{predict.gnls}}
#' @param naPattern missing value pattern. See \code{\link[nlme]{predict.gnls}}
#' @param data the data argument is needed when using this function inside user defined functions.
#' It should be identical to the data used to fit the model.
#' @param ... additional arguments (it is possible to supply a newdata this way)
#' @return It returns a vector with simulated values with length equal to the number of rows 
#' in the original data
#' @details It uses function \code{\link[MASS]{mvrnorm}} to generate new values for the coefficients
#' of the model using the Variance-Covariance matrix \code{\link{vcov}}. This variance-covariance matrix 
#' refers to the one for the parameters \sQuote{beta}, not the one for the residuals.
#' @seealso \code{\link[nlme]{predict.gnls}}
#' @export
#' @examples 
#' \donttest{
#' require(nlme)
#' data(barley, package = "nlraa")
#' 
#' fit.gnls <- gnls(yield ~ SSlinp(NF, a, b, xs), data = barley)
#' 
#' sim <- simulate_gnls(fit.gnls)
#' }

simulate_gnls <- function(object, psim = 1, na.action = na.fail, naPattern = NULL, 
                          data = NULL, ...){
    ##
    ## method for predict() designed for objects inheriting from class gnls
    ##
    if(!inherits(object, "gnls")) stop("This function is only for 'gnls' objects")
  
    mCall <- object$call
    
    ## ndata <- eval(object$call$data)
    ## Is this more robust?
    args <- list(...)
    if(!is.null(args$newdata)){
      ndata <- args$newdata
      if(!is.null(object$modelStruct$corStruct) && psim == 2)
        stop("At this point 'newdata' is not compatible with psim = 2 and correlated residuals",
             call. = FALSE)
    }else{
      if(is.null(data)){
        ndata <- try(nlme::getData(object), silent = TRUE)
        if(inherits(ndata, "try-error") || is.null(ndata)) 
          stop("'data' argument is required. It is likely you are using simulate_gnls inside another function")
      }else{
        if(object$dims$N != nrow(data)){
          stop("Number of rows in data argument does not match the original data \n
              The data argument should only be used to pass the same data.frame \n 
              used to fit the model",
               call. = FALSE)
        }
        ndata <- data
      } 
    } 
    
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
    
    ##------ This section is by FEM 2020-04-11 ------##
    if(psim == 0){
      prs <- coef(object)  
    }
    
    if(psim == 1){
      ## Only sample from the vector of coefficiencts
      prs <- MASS::mvrnorm(n = 1, mu = coef(object), Sigma = vcov(object))
    }
    
    if(psim == 2){
      prs <- MASS::mvrnorm(n = 1, mu = coef(object), Sigma = vcov(object))
      ## Sample standardized residuals and scale by the 
      ## residual standard error
      ## N is the number of rows in the data
      ## this works for uncorrelated errors
      if(is.null(object$modelStruct$corStruct)){
        if(is.null(args$newdata) || is.null(object$modelStruct$varStruct)){
          rsds.std <- stats::rnorm(N, 0, 1)
          rsds <- rsds.std * attr(residuals(object), "std") ## This last term is 'sigma'        
        }else{
          ## This chunk of code was added 2021-08-13. It allows for newdata, varStruct and psim = 2
          rsds.std <- stats::rnorm(nrow(ndata), 0, 1)
          rsds <- rsds.std * predict_varFunc(object, newdata = ndata)
        }
      }else{
        ## This was added 2020-05-04
        ## This extracts the variance covariance of the error matrix
        ## This is potentially a huge matrix
        var.cov.err <- var_cov(object, sparse = TRUE, data = ndata)
        ## This generates residuals considering the variance covariance of the error matrix
        ## It is computationally demanding
        ## I'm going to try using the mean of the errors for the mean (2020-05-26)
        ## I'm thinking now (2020-06-03) that the Cholesky factorization is better
        ## rsds <- MASS::mvrnorm(mu = residuals(object), Sigma = var.cov.err) - old, pre 2020-06-03
        ## Update: The Cholesky method is 10x faster than the MASS::mvrnorm method
        ## I have no idea if they produce similar simulations
        chol.var.cov.err <- Matrix::chol(var.cov.err)
        rsds <- Matrix::as.matrix(chol.var.cov.err %*% rnorm(nrow(chol.var.cov.err)))
      }
    }
    ##------ End FEM section ----------------------##
    
    for(nm in pnames) {
      if (!is.logical(plist[[nm]])) {
        form1s <- stats::asOneSidedFormula(params[[nm]][[3]])
        plist[[nm]] <- model.matrix(form1s, model.frame(form1s, dataMod))
      }
    }
    
    modForm <- nlme::getCovariateFormula(object)[[2]]
    val <- eval(modForm, data.frame(dataMod,
                                    nlraa_getParsGnls(plist, object$pmap, prs, N)))[naPat]
    
    ## ------------ FEM 04-14-2020 ---------------- ##
    ## If the simulation level is 2 we also add residuals
    if(psim == 2) val <- as.vector(val + rsds)
    ## --------------------------------------------##
      
    names(val) <- row.names(ndata)
    lab <- "Simulated values"
    if (!is.null(aux <- attr(object, "units")$y)) {
      lab <- paste(lab, aux)
    }
    attr(val, "label") <- lab
    
    return(val)
}

##' This function is not exported by nlme so I need to include it here
##' As with predict.gnls it was originally written by Pinheiro and Bates
##' and copied verbatim from the nlme package.
##' @name nlraa_getParsGnls
##' @param plist parameter list
##' @param pmap parameter map
##' @param beta beta coefficients
##' @param N I presume the number of observations
##' @noRd
nlraa_getParsGnls <- function(plist, pmap, beta, N)
{
  pars <- array(0, c(N, length(plist)), list(NULL, names(plist)))
  for (nm in names(plist)) {
    pars[, nm] <-
      if (is.logical(p <- plist[[nm]]))
        beta[pmap[[nm]]]
    else
      p %*% beta[pmap[[nm]]]
  }
  pars
}