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

#' This function is based on \code{\link[nlme]{predict.gls}} function 
#' 
#' @title Simulate fitted values from an object of class \code{\link[nlme]{gls}}
#' @description Simulate values from an object of class gls. Unequal variances, 
#' as modeled using the \sQuote{weights} option are supported, and there is experimental
#' code for dealing with the \sQuote{correlation} structure. This generates just one simulation
#' from these type of models. To generate multiple simulations use \code{\link{simulate_lme}}
#' @name simulate_gls
#' @param object object of class \code{\link[nlme]{gls}}
#' @param psim parameter simulation level, 0: for fitted values, 1: for simulation from 
#' fixed parameters (assuming a fixed vcov matrix), 2: for simulation considering the 
#' uncertainty in the residual standard error (sigma), this returns data which
#' will appear similar to the observed values
#' @param na.action default \sQuote{na.fail}. See \code{\link[nlme]{predict.gls}}
#' @param naPattern missing value pattern. See \code{\link[nlme]{predict.gls}}
#' @param data the data argument is needed when using this function inside user defined functions.
#' It should be identical to the data used to fit the model.
#' @param ... additional arguments (it is possible to supply a newdata this way)
#' @return It returns a vector with simulated values with length equal to the number of rows 
#' in the original data
#' @details It uses function \code{\link[MASS]{mvrnorm}} to generate new values for the coefficients
#' of the model using the Variance-Covariance matrix \code{\link{vcov}}. This variance-covariance matrix 
#' refers to the one for the parameters \sQuote{beta}, not the one for the residuals.
#' @seealso \code{\link[nlme]{predict.gls}} \code{\link{simulate_lme}}
#' @export
#' @examples 
#' \donttest{
#' require(nlme)
#' data(Orange)
#' 
#' fit.gls <- gls(circumference ~ age, data = Orange, 
#'                weights = varPower())
#' 
#' ## Visualize covariance matrix
#' fit.gls.vc <- var_cov(fit.gls)
#' image(log(fit.gls.vc[,ncol(fit.gls.vc):1]))
#' 
#' sim <- simulate_gls(fit.gls)
#' }

simulate_gls <- function(object, psim = 1, na.action = na.fail, naPattern = NULL, data = NULL, ...){

  if(!inherits(object, "gls")) stop("This function is only for 'gls' objects")
  
  ##
  ## method for predict() designed for objects inheriting from class gls
  ##
  args <- list(...)
  if(!is.null(args$newdata)){
    newdata <- args$newdata
    if(!is.null(object$modelStruct$corStruct) && psim == 2)
      stop("At this point 'newdata' is not compatible with psim = 2 and correlated residuals",
           call. = FALSE)
  }else{
    if(is.null(data)){
      newdata <- try(nlme::getData(object), silent = TRUE)
      if(inherits(newdata, "try-error") || is.null(newdata)) 
        stop("'data' argument is required. It is likely you are using simulate_gls inside another function")
    }else{
      if(object$dims$N != nrow(data)){
        ## I could fit the model here to make sure
        stop("Number of rows in data argument does not match the original data \n
              The data argument should only be used to pass the same data.frame \n 
              used to fit the model",
             call. = FALSE)
      }
      newdata <- data
    } 
  } 

  form <- getCovariateFormula(object)
  mfArgs <- list(formula = form, data = newdata, na.action = na.action)
  mfArgs$drop.unused.levels <- TRUE
  dataMod <- do.call(model.frame, mfArgs)
  ## making sure factor levels are the same as in contrasts
  contr <- object$contrasts
  for(i in names(dataMod)) {
    if (inherits(dataMod[,i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMod[,i])
      levsC <- dimnames(contr[[i]])[[1]]
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(sprintf(ngettext(sum(wch),
                              "level %s not allowed for %s",
                              "levels %s not allowed for %s"),
                     paste(levs[wch], collapse = ",")),
             domain = NA)
      }
      attr(dataMod[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
      ##      if (length(levs) < length(levsC)) {
      ##        if (inherits(dataMod[,i], "ordered")) {
      ##          dataMod[,i] <- ordered(as.character(dataMod[,i]), levels = levsC)
      ##        } else {
      ##          dataMod[,i] <- factor(as.character(dataMod[,i]), levels = levsC)
      ##        }
      ##      }
    }
  }
  N <- nrow(dataMod)
  if (length(all.vars(form)) > 0) {
    ##    X <- model.matrix(form, dataMod, contr)
    X <- model.matrix(form, dataMod)
  } else {
    X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
  }
  
  ##----------- FEM 2020-06-01 ---------------##
  if(psim == 0){
    cf <- coef(object)  
  }
  
  if(psim == 1){
    cf <- MASS::mvrnorm(n = 1, mu = coef(object), Sigma = vcov(object))
  }
  
  if(psim == 2){
    ## For more details see the comments in simulate_gnls
    cf <- MASS::mvrnorm(n = 1, mu = coef(object), Sigma = vcov(object))
    if(is.null(object$modelStruct$corStruct)){
      if(is.null(args$newdata) || is.null(object$modelStruct$varStruct)){
        rsds.std <- stats::rnorm(N, 0, 1)
        rsds <- rsds.std * attr(residuals(object), "std") ## This last term is 'sigma'        
      }else{
        rsds.std <- stats::rnorm(nrow(newdata), 0, 1)
        rsds <- rsds.std * predict_varFunc(object, newdata = newdata)
      }
    }else{
      ## This is one way of doing this, but might not be the best
      var.cov.err <- var_cov(object, sparse = TRUE, data = newdata)
      ## rsds <- MASS::mvrnorm(mu = residuals(object), Sigma = var.cov.err)
      ## An alternative is to do a Cholesky decomposition first
      ## Since var.cov.err is sparse using the Matrix package might be
      ## Much better
      chol.var.cov.err <- Matrix::chol(var.cov.err)
      rsds <- chol.var.cov.err %*% rnorm(nrow(chol.var.cov.err))
    }
  }
  
  val <- c(X[, names(cf), drop = FALSE] %*% cf)
  
  ## If psim == 2, we add residuals
  if(psim == 2) val <- as.vector(val + rsds)
  
  lab <- "Predicted values"
  if (!is.null(aux <- attr(object, "units")$y)) {
    lab <- paste(lab, aux)
  }
  structure(val, label = lab)
}

#' Function to predict the residual variance for newdata
#' @name predict_varFunc
#' @param object object of class \sQuote{gls}, \sQuote{gnls}, \sQuote{lme} or \sQuote{nlme}
#' @param newdata new data frame for which to predict the variance based on the structure in \sQuote{object}
#' @noRd
predict_varFunc <- function(object, newdata){
  
  ## If this fails, we know that newdata is not appropriate for this object
  fttd <- predict(object, newdata = newdata)
  
  if(is.null(object$modelStruct$varStruct))
    stop("varStruct should not be null for this function", call. = TRUE)
  
  ## These are the standard deviations, for the original model
  stds <- sigma(object)/nlme::varWeights(object$modelStruct$varStruct)
  vrSt <- object$modelStruct$varStruct
  
  if(!inherits(vrSt, c("varIdent", "varFixed", "varExp", "varPower"))){
    stop("Only varIdent, varFixed, varExp, varPower classes are supported at this time", call. = FALSE)
  }
      
  ## Here I would need a specific method for each variance function
  if(inherits(vrSt, "varFixed")){
    ## I don't think groups are relevant in this case
    ## nlme actually throws a warning saying that grouping factors are ignored
    ans <- sigma(object) * sqrt(newdata[[as.character(formula(vrSt)[2])]])
  }
  
  if(inherits(vrSt, "varIdent")){
    ## Groups are relevant in this case
    if(is.null(getGroups(vrSt)))
      stop("Groups should be present for varIdent", call. = FALSE)
    ## The covariate is ignored for this variance function
    ans <- numeric(nrow(newdata))
    grp.nm <- as.character(getGroupsFormula(vrSt)[[2]])
    if(!grp.nm %in% names(newdata))
      stop("Grouping factor should be present in 'newdata' object", call. = FALSE)
    if(grepl("*", grp.nm, fixed = TRUE))
      stop("This is not supported yet. Please submit this as an issue to github if you need it.")
    for(i in 1:nrow(newdata)){
      crr.grp <- newdata[i, grp.nm] ## This is the current group
      wch.grp.nm <- which(names(varWeights(vrSt)) == crr.grp)[1] ## This is the index for the weight of the current group
      ans[i] <- sigma(object) * (1/varWeights(vrSt))[wch.grp.nm]  ## This computes the std for the current group row
    }
  }
    
  if(inherits(vrSt, "varExp")){
    var_exp_fun <- function(x, t) exp(2 * t * x)
    if(is.null(getGroups(vrSt))){
      if(any(grepl("fitted", as.character(formula(vrSt))))){
        ## Need to fit the model to get the covariate
        cvrt <- fttd
      }else{
        cvrt.nm <- as.character(getCovariateFormula(vrSt))[2]
        if(!grepl(cvrt.nm, names(newdata)))
          stop("Variance covariate should be present in 'newdata' object", call. = FALSE)
        cvrt <- newdata[[cvrt.nm]] ## Is this guranteed to always work?
      }
      ans <- sigma(object) * sqrt(var_exp_fun(cvrt, coef(vrSt)))
    }else{
      ans <- numeric(nrow(newdata))
      grp.nm <- as.character(getGroupsFormula(vrSt)[[2]])
      if(!grp.nm %in% names(newdata))
        stop("Grouping factor should be present in 'newdata' object", call. = FALSE)
      if(grepl("*", grp.nm, fixed = TRUE))
        stop("This is not supported yet. Please submit this as an issue to github if you need it.")
      ## Looping over group levels
      for(i in unique(newdata[[grp.nm]])){
        ## For each group compute the answer
        wch.crr.grp <- which(newdata[[grp.nm]] == i)
        grp.coef <- coef(vrSt)[which(attr(vrSt, "groupNames") == i)] ## This is the coef for the current group
        if(any(grepl("fitted", as.character(formula(vrSt))))){
          cvrt <- fttd[wch.crr.grp]
        }else{
          cvrt <- newdata[wch.crr.grp, as.character(getCovariateFormula(vrSt))[[2]]] ## Is this guranteed to always work?
        }
        ans[wch.crr.grp] <- sigma(object) * sqrt(var_exp_fun(cvrt, grp.coef)) ## This computes the std for the current group row
      }
    }
  }
  
  if(inherits(vrSt, "varPower")){
    var_power_fun <- function(x, delta) abs(x)^(2 * delta)
    if(is.null(getGroups(vrSt))){
      if(any(grepl("fitted", as.character(formula(vrSt))))){
        ## Need to fit the model to get the covariate
        cvrt <- fttd
      }else{
        cvrt.nm <- as.character(getCovariateFormula(vrSt))[2]
        if(!grepl(cvrt.nm, names(newdata)))
          stop("Variance covariate should be present in 'newdata' object", call. = FALSE)
        cvrt <- newdata[[cvrt.nm]]
      }
      ans <- sigma(object) * sqrt(var_power_fun(cvrt, coef(vrSt)))
    }else{
      ans <- numeric(nrow(newdata))
      grp.nm <- as.character(getGroupsFormula(vrSt)[[2]])
      if(!grp.nm %in% names(newdata))
        stop("Grouping factor should be present in 'newdata' object", call. = FALSE)
      if(grepl("*", grp.nm, fixed = TRUE))
        stop("This is not supported yet. Please submit this as an issue to github if you need it.", call. = FALSE)
      ## I don't think this loop is that inefficienct
      for(i in unique(newdata[[grp.nm]])){
        ## For each group compute the answer
        wch.crr.grp <- which(newdata[[grp.nm]] == i)
        grp.coef <- coef(vrSt)[which(attr(vrSt, "groupNames") == i)] ## This is the coef for the current group
        if(any(grepl("fitted", as.character(formula(vrSt))))){
          cvrt <- fttd[wch.crr.grp]
        }else{
          cvrt <- newdata[wch.crr.grp, as.character(getCovariateFormula(vrSt))[[2]]] ## Is this guranteed to always work?
        }
        ans[wch.crr.grp] <- sigma(object) * sqrt(var_power_fun(cvrt, grp.coef)) ## This computes the std for the current group row
      }
    }
  }
  
  ans <- c(as.vector(ans))
  return(ans)
}

