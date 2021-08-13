### FEM: I include this because it was in the original file lme.R
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

#' This function is based on \code{\link[nlme]{predict.lme}} function 
#' 
#' @title Simulate values from an object of class \code{\link[nlme]{lme}}
#' @description Simulate values from an object of class lme. Unequal variances, 
#' as modeled using the \sQuote{weights} option are supported, and there is
#' experimental code for considering the \sQuote{correlation} structure.
#' @name simulate_lme
#' @param object object of class \code{\link[nlme]{lme}} or \code{\link[nlme]{gls}}
#' @param nsim number of samples, default 1
#' @param psim parameter simulation level, 0: for fitted values, 1: for simulation from 
#' fixed parameters (assuming a fixed vcov matrix), 2: for simulation considering the 
#' uncertainty in the residual standard error (sigma), this returns data which
#' will appear similar to the observed values. 3: in addition samples a new set of random effects.
#' @param value whether to return a matrix (default) or an augmented data frame
#' @param data the data argument is needed when using this function inside user defined functions.
#' @param ... additional arguments (it is possible to supply a newdata this way)
#' @return It returns a vector with simulated values with length equal to the number of rows 
#' in the original data
#' @details It uses function \code{\link[MASS]{mvrnorm}} to generate new values for the coefficients
#' of the model using the Variance-Covariance matrix \code{\link{vcov}}. This variance-covariance matrix 
#' refers to the one for the parameters 'beta', not the one for the residuals.
#' @note I find the simulate.merMod in the lme4 pacakge confusing. There is use.u and several versions of re.form.
#' From the documentation it seems that if use.u = TRUE, then the current values of the random effects are used.
#' This would mean that it is equivalent to psim = 2 in this function. Then use.u = FALSE, would be equivalent 
#' to psim = 3. re.form allows for specifying the formula of the random effects.
#' @seealso \code{\link[nlme]{predict.lme}} and \sQuote{simulate.merMod} in the \sQuote{lme4} package.
#' @export
#' @examples 
#' \donttest{
#' require(nlme)
#' data(Orange)
#' 
#' fm1 <- lme(circumference ~ age, random = ~ 1 | Tree, data = Orange)
#' 
#' sims <- simulate_lme(fm1, nsim = 10)
#' 
#' }

simulate_lme <- function(object, nsim = 1, psim = 1,
                         value = c("matrix", "data.frame"), 
                         data = NULL, ...){
  
  ## Error checking
  if(!inherits(object, c("gls","lme"))) stop("object should be of class 'gls' or 'lme'")
  
  value <- match.arg(value)
  
  if(is.null(list(...)$newdata)){
    sim.mat <- matrix(ncol = nsim, nrow = length(fitted(object)))
  }else{
    sim.mat <- matrix(ncol = nsim, nrow = nrow(list(...)$newdata))  
  } 
    
  ## First example for the gls case
  for(i in seq_len(nsim)){
    if(inherits(object, "gls")){
      sim.mat[,i] <- as.vector(simulate_gls(object, psim = psim, data = data, ...))
    }
    if(inherits(object, "lme")){
      sim.mat[,i] <- as.vector(simulate_lme_one(object, psim = psim, data = data, ...))
    }
  }
  
  if(value == "matrix"){
    colnames(sim.mat) <- paste0("sim_", 1:nsim)
    return(sim.mat)  
  }else{
    if(is.null(data)){
      dat <- try(nlme::getData(object), silent = TRUE)
      if(inherits(dat, "try-error") || is.null(dat)) 
        stop("'data' argument is required. It is likely you are using simulate_lme inside another function")
    }else{
      if(is.null(list(...)$newdata)){
        dat <- data  
      }else{
        dat <- list(...)$newdata
      }
    } 
    ## I don't understand why I need to do this
    ## This code is not needed for simulate_nlme
    ndat <- NULL 
    for(i in seq_len(nsim)) ndat <- rbind(ndat, dat)
    ## End ugly code
    adat <- data.frame(ii = as.factor(rep(1:nsim, each = nrow(dat))),
                       ndat,
                       sim.y = c(sim.mat),
                       row.names = 1:c(nsim * nrow(dat)))   
    return(adat)
  }
}

## May be should document this
simulate_lme_one <- function(object, psim = 1, level = Q, asList = FALSE, na.action = na.fail, data = NULL, ...){
    ##
    ## method for predict() designed for objects inheriting from class lme
    ##
  
  if(!missing(level) && length(level) > 1) stop("level must be of length = 1")
  
  Q <- object$dims$Q
  # if (missing(newdata)) {		# will return fitted values
  #   val <- fitted(object, level, asList)
  #   if (length(level) == 1) return(val)
  #   return(data.frame(object[["groups"]][,level[level != 0], drop = FALSE],
  #                     predict = val))
  # }
  if(!missing(level) && psim == 2 && level < Q)
    stop("psim = 2 should only be used for the deepest level of the hierarchy")
  
  ## Data
  args <- list(...)
  if(!is.null(args$newdata)){
    newdata <- args$newdata
    if(!is.null(object$modelStruct$corStruct) && psim > 1)
      stop("At this point 'newdata' is not compatible with psim > 1 and correlated residuals",
           call. = FALSE)
  }else{
    if(is.null(data)){
      newdata <- try(nlme::getData(object), silent = TRUE)
      if(inherits(newdata, "try-error") || is.null(newdata)) 
        stop("'data' argument is required. It is likely you are using simulate_lme_one inside another function")
    }else{
      if(object$dims$N != nrow(data)){
        stop("Number of rows in data argument does not match the original data \n
              The data argument should only be used to pass the same data.frame \n 
              used to fit the model",
             call. = FALSE)
      }
      newdata <- data
    }  
  } 
    
  maxQ <- max(level)			# maximum level for predictions
  nlev <- length(level)
  mCall <- object$call
  fixed <- eval(eval(mCall$fixed)[-2])  # RHS
  Terms <- object$terms
  newdata <- as.data.frame(newdata)
  if (maxQ > 0) {			# predictions with random effects
    whichQ <- Q - (maxQ-1):0
    reSt <- object$modelStruct$reStruct[whichQ]
    lmeSt <- lmeStruct(reStruct = reSt)
    groups <- getGroupsFormula(reSt)
    if (any(is.na(match(all.vars(groups), names(newdata))))) {
      ## groups cannot be evaluated in newdata
      stop("cannot evaluate groups for desired levels on 'newdata'")
    }
  } else {
    reSt <- NULL
  }
    
  mfArgs <- list(formula = asOneFormula(formula(reSt), fixed),
                 data = newdata, na.action = na.action,
                 drop.unused.levels = TRUE)
  dataMix <- do.call(model.frame, mfArgs)
  origOrder <- row.names(dataMix)	# preserve the original order
  whichRows <- match(origOrder, row.names(newdata))
  
  if (maxQ > 0) {
    ## sort the model.frame by groups and get the matrices and parameters
    ## used in the estimation procedures
    grps <- getGroups(newdata,
                      eval(substitute(~ 1 | GRPS,
                                      list(GRPS = groups[[2]]))))
    ## ordering data by groups
    if (inherits(grps, "factor")) {	# single level
      grps <- grps[whichRows, drop = TRUE]
      oGrps <- data.frame(grps)
      ## checking if there are missing groups
      if (any(naGrps <- is.na(grps))) {
        grps[naGrps] <- levels(grps)[1L]	# input with existing level
      }
      ord <- order(grps)     #"order" treats a single named argument peculiarly
      grps <- data.frame(grps)
      row.names(grps) <- origOrder
      names(grps) <- names(oGrps) <- as.character(deparse((groups[[2L]])))
    } else {
      grps <- oGrps <-
        do.call(data.frame, ## FIXME?  better  lapply(*, drop)   ??
                lapply(grps[whichRows, ], function(x) x[drop = TRUE]))
      ## checking for missing groups
      if (any(naGrps <- is.na(grps))) {
        ## need to input missing groups
        for(i in names(grps)) {
          grps[naGrps[, i], i] <- levels(grps[,i])[1L]
        }
        naGrps <- t(apply(naGrps, 1, cumsum)) # propagating NAs
      }
      ord <- do.call(order, grps)
      ## making group levels unique
      grps[, 1] <- grps[, 1][drop = TRUE]
      for(i in 2:ncol(grps)) {
        grps[, i] <-
          as.factor(paste(as.character(grps[, i-1]),
                          as.character(grps[, i  ]), sep = "/"))
      }
    }
    naGrps <- cbind(FALSE, naGrps)[ord, , drop = FALSE]
    grps <- grps[ord, , drop = FALSE]
    dataMix <- dataMix[ord, ,drop = FALSE]
  }
  ## making sure factor levels are the same as in contrasts
  contr <- object$contrasts
  for(i in names(dataMix)) {
    if (inherits(dataMix[,i], "factor") && !is.null(contr[[i]])) {
      levs <- levels(dataMix[,i])
      levsC <- dimnames(contr[[i]])[[1L]]
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(sprintf(ngettext(sum(wch),
                              "level %s not allowed for %s",
                              "levels %s not allowed for %s"),
                     paste(levs[wch], collapse = ",")),
             domain = NA)
      }
      ## if (length(levs) < length(levsC)) {
      ##   if (inherits(dataMix[,i], "ordered")) {
      ##     dataMix[,i] <- ordered(as.character(dataMix[,i]), levels = levsC)
      ##   } else {
      ##     dataMix[,i] <- factor(as.character(dataMix[,i]), levels = levsC)
      ##   }
      ## }
      attr(dataMix[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
    }
  }
  if (maxQ > 0) {
    revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order
    Z <- model.matrix(reSt, dataMix)
    ncols <- attr(Z, "ncols")
    Names(lmeSt$reStruct) <- attr(Z, "nams")
  }
  N <- nrow(dataMix)
  X <- if (length(all.vars(fixed)) > 0) {
    model.matrix(fixed, model.frame(delete.response(Terms), dataMix))
  } else if(attr(terms(fixed), "intercept")) {
    array(1, c(N, 1), list(row.names(dataMix), "(Intercept)"))
  } else {
    array(, c(N, 0))
  }
  ## This is used when population level fits are requested
  if (maxQ == 0) {
    
    if(psim == 0) fixf <- fixef(object)
    
    if(psim == 1){
      fixf <- MASS::mvrnorm(n = 1, mu = fixef(object), Sigma = vcov(object))
    }
    ## only population predictions
    val <-  if(ncol(X)) c(X %*% fixf) else rep(0, nrow(X))
    attr(val, "label") <- "Predicted values"
    return(val)
  }
    
  ncols <- c(ncols, dim(X)[2L], 1)
  ## creating the condensed linear model
  attr(lmeSt, "conLin") <-
    list(Xy = array(c(Z, X, double(N)), c(N, sum(ncols)),
                    list(row.names(dataMix), c(colnames(Z), colnames(X), "resp"))),
         dims = nlraa_MEdims(grps, ncols))
  ## Getting the appropriate BLUPs of the random effects
  re <- object$coefficients$random[1:maxQ]
  j <- 1
  for(i in names(re)) {
    ugrps <- unique(as.character(grps[, i]))
    val <- array(NA, c(length(ugrps), ncol(re[[i]])),
                 list(ugrps, dimnames(re[[i]])[[2L]]))
    mGrps <- match(ugrps, dimnames(re[[i]])[[1L]])
    mGrps <- mGrps[!is.na(mGrps)]
    re[[i]] <- re[[i]][mGrps, , drop = FALSE]
    val[dimnames(re[[i]])[[1L]], ] <- re[[i]]
    
    if(psim == 3){
      reDelta <- as.matrix(object$modelStruct$reStruct[[j]]) * sigma(object)^2
      
      val <- MASS::mvrnorm(n = nrow(re[[i]]), 
                           mu = rep(0, ncol(re[[i]])), 
                           Sigma = reDelta, empirical = TRUE)
      dimnames(val) <- dimnames(re[[i]])
      j <- j + 1
    }
    re[[i]] <- val
  }
  
  ##---------- FEM insert: 2020-05-30 ---------------##
  ## Returns fitted values
  if(psim == 0){
    fix <- fixef(object)
  }
  
  ## Simualte from vector of fixed effects only
  if(psim == 1){
    fix <- MASS::mvrnorm(n = 1, mu = fixef(object), Sigma = vcov(object))
  }
    
  if(psim == 2 || psim == 3){
    fix <- MASS::mvrnorm(n = 1, mu = fixef(object), Sigma = vcov(object))
    
    if(is.null(object$modelStruct$corStruct)){
      if(is.null(args$newdata) || is.null(object$modelStruct$varStruct)){
        rsds.std <- stats::rnorm(N, 0, 1)
        rsds <- rsds.std * attr(object[["residuals"]], "std") ## This last term is 'sigma'        
      }else{
        rsds.std <- stats::rnorm(nrow(newdata), 0, 1)
        rsds <- rsds.std * predict_varFunc(object, newdata = newdata)
      }
    }else{
      ## For details on this see simulate_gnls
      var.cov.err <- var_cov(object, sparse = TRUE, data = newdata)
      chol.var.cov.err <- Matrix::chol(var.cov.err)
      rsds <- Matrix::as.matrix(chol.var.cov.err %*% rnorm(nrow(chol.var.cov.err)))
    }
  }
 
  attr(lmeSt, "lmeFit") <- list(beta = fix, b = re) ## I edited this line "fix" instead of "fixef(object)"
  val <- fitted(lmeSt, level = 0:maxQ)
  val[as.logical(naGrps)] <- NA			# setting missing groups to NA
  ## putting back in original order and extracting levels
  val <- val[revOrder, level + 1L]		# predictions

  if(psim == 2 || psim == 3){
      ## The reordering can be tricky, I'm assuming the previous step put them in the 
      ## original order and that rsds are also in the original order
      val <- val + rsds
  } 
    
  if (maxQ > 1) {                      # making groups unique
    for(i in 2:maxQ)
      oGrps[, i] <-
        as.factor(paste(as.character(oGrps[,i-1]),
                        as.character(oGrps[,i  ]), sep = "/"))
  }
  if (nlev == 1) {
    grps <- as.character(oGrps[, level])
    if (asList) {
      val <- split(val, ordered(grps, levels = unique(grps)))
    } else {
      names(val) <- grps
    }
    lab <- "Predicted values"
    if (!is.null(aux <- attr(object, "units")$y)) {
      lab <- paste(lab, aux)
    }
    attr(val, "label") <- lab
    val
  } else {
    data.frame(oGrps, predict = val)
  }
}

nlraa_MEdims <- function(groups, ncols)
{
  ## define constants used in matrix decompositions and log-lik calculations
  ## first need some local functions
  glengths <-
    ## returns the group lengths from a vector of last rows in the group
    function(lstrow) diff(c(0, lstrow))
  offsets <-
    ## converts total number of columns(N), columns per level(ncols), and
    ## a list of group lengths to offsets in C arrays
    function(N, ncols, lstrow, triangle = FALSE)
    {
      pop <- function(x) x[-length(x)]
      cstart <- c(0, cumsum(N * ncols))
      for (i in seq_along(lstrow)) {
        lstrow[[i]] <- cstart[i] +
          if (triangle) {
            lstrow[[i]] - ncols[i]        # storage offsets style
          } else {
            pop(c(0, lstrow[[i]]))        # decomposition style
          }
      }
      lstrow
    }
  Q <- ncol(groups)                     # number of levels
  N <- nrow(groups)                     # number of observations
  ## 'isLast' indicates if the row is the last row in the group at that level.
  ## this version propagates changes from outer groups to inner groups
  ## isLast <- (array(unlist(lapply(c(rev(as.list(groups)),
  ##                                list(X = rep(0, N), y = rep(0, N))),
  ##                               function(x) c(0 != diff(codes(x)), TRUE))),
  ##                 c(N, Q+2), list(NULL, c(rev(names(groups)), "X", "y")))
  ##            %*% (row(diag(Q+2)) >= col(diag(Q+2)))) != 0
  ## this version does not propagate changes from outer to inner.
  isLast <- array(FALSE, dim(groups) + c(0, 2),
                  list(NULL, c(rev(names(groups)), "X", "y")))
  for(i in 1:Q) {
    isLast[, Q + 1L - i] <- c(0 != diff(as.integer(groups[[i]])), TRUE)
  }
  isLast[N,  ] <- TRUE
  lastRow <- apply(isLast, 2, function(x) seq_along(x)[x])
  if(!is.list(lastRow)) {
    nm <- names(lastRow)
    lastRow <- as.list(lastRow)
    names(lastRow) <- nm
  }
  
  isLast <- t(isLast)
  strSizes <- cumsum(ncols * isLast) * isLast # required storage sizes
  lastStr <- apply(t(strSizes), 2, function(x) x[x != 0])
  if(!is.list(lastStr)) {
    nm <- names(lastStr)
    lastStr <- as.list(lastStr)
    names(lastStr) <- nm
  }
  strRows <- max(lastStr[[length(lastStr)]])
  lastBlock <- vector("list", Q)
  names(lastBlock) <- rownames(strSizes)[1:Q]
  for(i in 1:Q) lastBlock[[i]] <- c(strSizes[i, -N], strRows)
  maxStr <- do.call(pmax, lastBlock)
  for(i in 1:Q) lastBlock[[i]] <- maxStr[as.logical(lastBlock[[i]])]
  lastBlock <- c(lastBlock, list(X = strRows, y = strRows))
  list(N = N,                   # total number of rows in data
       ZXrows = N,              # no. of rows in array
       ZXcols = sum(ncols),     # no. of columns in array
       Q = Q,                   # no. of levels of random effects
       StrRows = strRows,       # no. of rows required for storage
       qvec = ncols * c(rep(1, Q), 0, 0), # lengths of random effects
       # no. of groups at each level
       ### This looks wrong: ")" at wrong place: unlist(*, N, N) !!
       ngrps = c(unlist(lapply(lastRow, length), N, N)),
       ###?ok ngrps = c(lengths(lastRow), N, N),# no. of groups at each level
       DmOff = c(0, cumsum(ncols^2))[1:(Q+2)],# offsets into DmHalf array by level
       ncol = ncols,            # no. of columns decomposed per level
       nrot = rev(c(0, cumsum(rev(ncols))))[-1L],# no. of columns rotated per level
       
       ZXoff = offsets(N, ncols, lastRow), # offsets into ZXy
       ZXlen = lapply(lastRow, glengths), # lengths of ZXy groups
       # storage array offsets
       SToff = offsets(strRows, ncols, lastStr, triangle = TRUE),
       # decomposition offsets
       DecOff = offsets(strRows, ncols, lastBlock),
       # decomposition lengths
       DecLen = lapply(lastBlock, glengths)
  )
}