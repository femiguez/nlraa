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

#' Simulate nlme based on \sQuote{nlme::predict.nlme} function 
#' 
#' @title Simulate fitted values from an object of class \sQuote{nlme}
#' @name simulate_nlme_one
#' @param object object of class \sQuote{gnls}
#' @param psim parameter simulation level, 0: for fitted values, 1: for simulation from 
#' fixed parameters (assuming a fixed vcov matrix), 2: for simulation considering the 
#' uncertainty in vcov (not implemented yet)
#' @param na.action default \sQuote{na.fail}
#' @param naPattern missing value pattern
#' @details It uses function \sQuote{MASS::mvrnorm} to generate new values for the coefficients
#' of the model using the Variance-Covariance matrix \sQuote{vcov}

simulate_nlme_one <- function(object, psim = 1, level = Q, asList = FALSE, na.action = na.fail,
           naPattern = NULL, ...){
    ##
    ## method for predict() designed for objects inheriting from class nlme
    ##
   
  if(!inherits(object, "nlme")) stop("Only for objects of class 'nlme'")
  
  Q <- object$dims$Q
    
  newdata <- eval(object$call$data)
    
  maxQ <- max(level)			# maximum level for predictions
  nlev <- length(level)
  newdata <- as.data.frame(newdata)
    
  if(maxQ > 0){			# predictions with random effects
    whichQ <- Q - (maxQ-1):0
    reSt <- object$modelStruct$reStruct[whichQ]
    ##    nlmeSt <- nlmeStruct(reStruct = reSt)
    groups <- nlme::getGroupsFormula(reSt)
    if (any(is.na(match(all.vars(groups), names(newdata))))) {
      ## groups cannot be evaluated in newdata
      stop("cannot evaluate groups for desired levels on 'newdata'")
    }
  }else{
    reSt <- NULL
  }
    
  mfArgs <- list(formula = nlme::asOneFormula(
    formula(object),
    object$call$fixed, formula(reSt), naPattern,
      omit = c(names(object$plist), "pi",
               deparse(nlme::getResponseFormula(object)[[2]]))),
      data = newdata, na.action = na.action,
      drop.unused.levels = TRUE)
    dataMix <- do.call(model.frame, mfArgs)
    origOrder <- row.names(dataMix)	# preserve the original order
    whichRows <- match(origOrder, row.names(newdata))
    
    if(maxQ > 0){
      ## sort the model.frame by groups and get the matrices and parameters
      ## used in the estimation procedures
      grps <- nlme::getGroups(newdata, eval(substitute(~ 1 | GRPS, list(GRPS = groups[[2]]))))
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
      }else{
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
      ## if (match(0, level, nomatch = 0)) {
      ##   naGrps <- cbind(FALSE, naGrps)
      ## }
      ## naGrps <- as.matrix(naGrps)[ord, , drop = FALSE]
      naGrps <- cbind(FALSE, naGrps)[ord, , drop = FALSE]
      grps <- grps[ord, , drop = FALSE]
      dataMix <- dataMix[ord, ,drop = FALSE]
      revOrder <- match(origOrder, row.names(dataMix)) # putting in orig. order
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
        attr(dataMix[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
      }
    }
    
    N <- nrow(dataMix)
    ##
    ## evaluating the naPattern expression, if any
    ##
    naPat <- if(is.null(naPattern)) rep(TRUE, N)
    else
      as.logical(eval(stats::asOneSidedFormula(naPattern)[[2]], dataMix))
    ##
    ## Getting  the plist for the new data frame
    ##
    plist <- object$plist
    fixed <- eval(object$call$fixed)
    if (!is.list(fixed))
      fixed <- list(fixed)
    
    fixed <- do.call(c, lapply(fixed, function(fix.i) {
      if (is.name(fix.i[[2]]))
        list(fix.i)
      else
        ## multiple parameters on left hand side
        eval(parse(text = paste0("list(", paste(paste(all.vars(fix.i[[2]]),
                                                      deparse (fix.i[[3]]),
                                                      sep = "~"),
                                                collapse = ","), ")")))
    }))
    
    fnames <- unlist(lapply(fixed, function(el) deparse(el[[2]])))
    names(fixed) <- fnames

    ## This section was included by FEM    
    if(psim == 0){
      fix <- fixef(object)  
    }
    
    if(psim == 1){
      fix <- MASS::mvrnorm(n = 1, mu = fixef(object), Sigma = vcov(object))
    }
    
    if(psim == 2){
      ## This method would need to sample Sigma before feeding it into mvrnorm
      stop("not implemented yet")
    }

    for(nm in fnames){
      if (!is.logical(plist[[nm]]$fixed)) {
        oSform <- stats::asOneSidedFormula(fixed[[nm]][[3]])
        plist[[nm]]$fixed <- model.matrix(oSform, model.frame(oSform, dataMix))
      }
    }
    
    if(maxQ > 0){
      grpsRev <- lapply(rev(grps), as.character)
      ranForm <- formula(reSt)[whichQ]
      namGrp <- names(ranForm)
      rnames <- lapply(ranForm, function(el)
        unlist(lapply(el, function(el1) deparse(el1[[2]]))))
      
      for(i in seq_along(ranForm)) {
        names(ranForm[[i]]) <- rnames[[i]]
      }
      ran <- ranef(object)
      ran <- if(is.data.frame(ran)) list(ran) else rev(ran)
      ##    rn <- lapply(ran[whichQ], names)
      ran <- lapply(ran, t)
      ranVec <- unlist(ran)
      
      for(nm in names(plist)) {
        for(i in namGrp) {
          if (!is.logical(plist[[nm]]$random[[i]])) {
            wch <- which(!is.na(match(rnames[[i]], nm)))
            plist[[nm]]$random[[i]] <-
              if (length(wch) == 1) {         # only one formula for nm
                oSform <- stats::asOneSidedFormula(ranForm[[i]][[nm]][[3]])
                model.matrix(oSform, model.frame(oSform, dataMix))
              } else {                        # multiple formulae
                lapply(ranForm[[i]][wch], function(el) {
                  if (el[[3]] == "1") {
                    TRUE
                  } else {
                    oSform <- stats::asOneSidedFormula(el[[3]])
                    model.matrix(oSform, model.frame(oSform, dataMix))
                  } })
              }
          }
        }
      }
    }else{
      namGrp <- ""
      grpsRev <- ranVec <- ran <- NULL
    }
    
    val <- vector("list", nlev)
    modForm <- nlme::getCovariateFormula(object)[[2]]
    omap <- object$map
    
    for(i in 1:nlev) {
      val[[i]] <- eval(modForm,
                       data.frame(dataMix,
                                  nlme:::getParsNlme(plist, omap$fmap, omap$rmapRel,
                                                     omap$bmap, grpsRev, fix, ranVec, ran,
                                                     level[i], N)))[naPat]
    }
    
    names(val) <- c("fixed", rev(namGrp))[level + 1]
    val <- as.data.frame(val)
    
    if(maxQ > 0){
      val <- val[revOrder, , drop = FALSE]
      if (any(naGrps)) {
        val[naGrps] <- NA
      }
    }
    ## putting back in original order
    
    if(maxQ > 1){                      # making groups unique
      for(i in 2:maxQ)
        oGrps[, i] <-
          as.factor(paste(as.character(oGrps[,i-1]),
                          as.character(oGrps[,i  ]), sep = "/"))
    }
    
    if(length(level) == 1){
      val <- val[,1] ## ?? not in predict.lme()
      if(level > 0){ # otherwise 'oGrps' are typically undefined
        grps <- as.character(oGrps[, level])
        if(asList){
          val <- split(val, ordered(grps, levels = unique(grps)))
        }else{
          names(val) <- grps
        }
      }
      lab <- "Simulated values"
      
      if(!is.null(aux <- attr(object, "units")$y)){
        lab <- paste(lab, aux)
      }
      
      attr(val, "label") <- lab
      return(val)
    }else{
      return(data.frame(oGrps, predict = val))
    }
}
