#' 
#' @title Average predictions from several (non)linear models based on IC weights
#' @name predict_nls
#' @rdname predict_nls
#' @description Computes weights based on AIC, AICc, or BIC and it generates weighted predictions by
#' the relative value of the IC values
#' @param ... \sQuote{nls} or \sQuote{lm} objects (\sQuote{glm} and \sQuote{gam} objects inherit \sQuote{lm}). 
#' @param criteria either \sQuote{AIC}, \sQuote{AICc} or \sQuote{BIC}.
#' @param interval either \sQuote{none}, \sQuote{confidence} or \sQuote{prediction}.
#' @param level probability level for the interval (default 0.95)
#' @param nsim number of simulations to perform for intervals. Default 1000.
#' @param resid.type either \sQuote{none}, \dQuote{resample}, \dQuote{normal} or \dQuote{wild}.
#' @param newdata new data frame for predictions
#' @param weights vector of weights of the same length as the number of models. It should sum up to one and 
#' it will override the information-criteria based weights. The weights should match the order of the models.
#' @return numeric vector of the same length as the fitted object when interval is equal to \sQuote{none}. Otherwise,
#' a data.frame with columns named (for a 0.95 level) \sQuote{Estimate}, \sQuote{Est.Error}, \sQuote{Q2.5} and \sQuote{Q97.5}
#' @note all the objects should be fitted to the same data. Weights are
#' based on the chosen IC value (exp(-0.5 * delta IC)). 
#' For models of class \code{\link[mgcv]{gam}} there is very limited support.
#' @seealso \code{\link[stats]{predict.lm}}, \code{\link[stats]{predict.nls}}, \code{\link[mgcv]{predict.gam}}, \code{\link{simulate_nls}}, \code{\link{simulate_gam}}
#' @export
#' @examples
#' \donttest{
#' ## Example
#' require(ggplot2)
#' require(mgcv)
#' data(barley, package = "nlraa")
#' 
#' fm.L <- lm(yield ~ NF, data = barley)
#' fm.Q <- lm(yield ~ NF + I(NF^2), data = barley)
#' fm.A <- nls(yield ~ SSasymp(NF, Asym, R0, lrc), data = barley)
#' fm.LP <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
#' fm.QP <- nls(yield ~ SSquadp3(NF, a, b, c), data = barley)
#' fm.BL <- nls(yield ~ SSblin(NF, a, b, xs, c), data = barley)
#' fm.G <- gam(yield ~ s(NF, k = 6), data = barley)
#' 
#' ## Print the table with weights
#' IC_tab(fm.L, fm.Q, fm.A, fm.LP, fm.QP, fm.BL, fm.G)
#' 
#' ## Each model prediction is weighted according to their AIC values
#' prd <- predict_nls(fm.L, fm.Q, fm.A, fm.LP, fm.QP, fm.BL, fm.G)
#' 
#' ggplot(data = barley, aes(x = NF, y = yield)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fm.L), color = "Linear")) +
#'   geom_line(aes(y = fitted(fm.Q), color = "Quadratic")) +
#'   geom_line(aes(y = fitted(fm.A), color = "Asymptotic")) +  
#'   geom_line(aes(y = fitted(fm.LP), color = "Linear-plateau")) + 
#'   geom_line(aes(y = fitted(fm.QP), color = "Quadratic-plateau")) + 
#'   geom_line(aes(y = fitted(fm.BL), color = "Bi-linear")) + 
#'   geom_line(aes(y = fitted(fm.G), color = "GAM")) + 
#'   geom_line(aes(y = prd, color = "Avg. Model"), linewidth = 1.2)
#' }

predict_nls <- function(..., criteria = c("AIC", "AICc", "BIC"), 
                        interval = c("none", "confidence", "prediction"),
                        level = 0.95, nsim = 1e3,
                        resid.type = c("none", "resample", "normal", "wild"),
                        newdata = NULL, weights){
  
  ## all objects should be of class 'nls' or inherit 'lm' (but this includes 'gam' and 'glm')
  nls.objs <- list(...)
  criteria <- match.arg(criteria)
  interval <- match.arg(interval)
  resid.type <- match.arg(resid.type)
  nls.nms <- get_mnames(match.call())
  
  lobjs <- length(nls.objs)
  wtab <- data.frame(model = character(lobjs), IC = NA)
  
  nr <- stats::nobs(nls.objs[[1]])
  if(!is.null(newdata)) nr <- nrow(newdata)
  
  if(interval == "none"){
    prd.mat <- matrix(nrow = nr, ncol = lobjs)
  }else{
    prd.mat <- matrix(nrow = nr, ncol = lobjs)
    prd.mat.se <- matrix(nrow = nr, ncol = lobjs)
    prd.mat.lwr <- matrix(nrow = nr, ncol = lobjs)
    prd.mat.upr <- matrix(nrow = nr, ncol = lobjs)
  } 
  
  data.name <- as.character(nls.objs[[1]]$call$data)
  
  for(i in seq_len(lobjs)){
    nls.obj <- nls.objs[[i]]
    if(!inherits(nls.obj, c("nls","lm"))) stop("All objects should be of class 'nls' or 'lm'")
    if(data.name != as.character(nls.obj$call$data)) stop("All models should be fitted to the same data")
    
    wtab$model[i] <- nls.nms[i]
    
    if(criteria == "AIC") wtab$IC[i] <- stats::AIC(nls.obj)
    if(criteria == "AICc") wtab$IC[i] <- AICc(nls.obj)
    if(criteria == "BIC") wtab$IC[i] <- stats::BIC(nls.obj)
  }
  
  ## Check weights
  if(!missing(weights)){
    if(length(weights) != lobjs)
      stop("'weights' should be a vector of length equal to the number of models", call. = FALSE)
    if(isFALSE(all.equal(sum(weights), 1)))
      stop("'sum of 'weights' should be equal to 1.", call. = FALSE)
    if(any(weights < 0))
      stop("all weights should be greater than zero", call. = FALSE)
    if(any(weights > 1))
      stop("all weights should be greater than zero", call. = FALSE)
  }
  
  ## Predictions
  if(interval == "none"){
    for(i in seq_len(lobjs)){
      nls.obj <- nls.objs[[i]]
      if(!is.null(newdata)){
        if(inherits(nls.obj, "gam")){
          prd.mat[,i] <- predict(nls.obj, newdata = newdata, type = "response")  
        }else{
          prd.mat[,i] <- predict(nls.obj, newdata = newdata)    
        }
      }else{
        if(inherits(nls.obj, "gam")){
          prd.mat[,i] <- predict(nls.obj, type = "response")  
        }else{
          prd.mat[,i] <- predict(nls.obj)    
        }
      }
    }
  }
  
  if(interval == "confidence" || interval == "prediction"){
    
    psim <- ifelse(interval == "confidence", 1, 2)

    lb <- (1 - level)/2
    ub <- 1 - lb 
    
    for(i in seq_len(lobjs)){
      nls.obj <- nls.objs[[i]]
      
      if(inherits(nls.obj, "lm") && !inherits(nls.obj, "glm")) 
        tmp.sim <- simulate_lm(nls.obj, psim = psim, nsim = nsim, 
                               resid.type = resid.type, newdata = newdata) 
      
      if(inherits(nls.obj, "gam")) 
        tmp.sim <- simulate_gam(nls.obj, psim = psim, nsim = nsim, 
                                resid.type = resid.type, newdata = newdata) 
      
      if(inherits(nls.obj, "nls")) 
        tmp.sim <- simulate_nls(nls.obj, psim = psim, nsim = nsim, 
                                resid.type = resid.type, newdata = newdata)

      prd.mat[,i] <- apply(tmp.sim, 1, quantile, probs = 0.5)
      prd.mat.se[,i] <- apply(tmp.sim, 1, sd)
      prd.mat.lwr[,i] <- apply(tmp.sim, 1, quantile, probs = lb)
      prd.mat.upr[,i] <- apply(tmp.sim, 1, quantile, probs = ub)
    }
  }

  wtab$dIC <- wtab$IC - min(wtab$IC)
  
  if(missing(weights)){
    wtab$weight <- exp(-0.5 * wtab$dIC) / sum(exp(-0.5 * wtab$dIC))  
  }else{
    wtab$weight <- weights
  }
  
  
  if(interval == "none"){
    ans <- rowSums(sweep(prd.mat, 2, wtab$weight, "*"))
  }else{ 
    prd <- rowSums(sweep(prd.mat, 2, wtab$weight, "*"))
    se <- rowSums(sweep(prd.mat.se, 2, wtab$weight, "*"))
    lwr <- rowSums(sweep(prd.mat.lwr, 2, wtab$weight, "*"))
    upr <- rowSums(sweep(prd.mat.upr, 2, wtab$weight, "*"))
    ans <- cbind(prd, se, lwr, upr)
    colnames(ans) <- c("Estimate", "Est.Error", 
                       paste0("Q", (1 - level)/2 * 100),
                       paste0("Q", (1 - (1 - level)/2)*100))
  }
  
  return(ans)
}


#' @title Information Criteria Table
#' @name IC_tab
#' @description Information criteria table with weights
#' @param ... model fit objects fitted to the same data
#' @param criteria either \sQuote{AIC}, \sQuote{AICc} or \sQuote{BIC}.
#' @param sort whether to sort by weights (default to TRUE)
#' @seealso \code{\link[bbmle]{ICtab}}
#' @note The delta and weights are calculated based on the \sQuote{criteria}
#' @export
#' 
IC_tab <- function(..., criteria = c("AIC","AICc","BIC"), sort = TRUE){
  
  objs <- list(...)
  criteria <- match.arg(criteria)
  nms <- get_mnames(match.call())
  lobjs <- length(objs)
  
  ictab <- data.frame(model = character(lobjs), df = NA, AIC = NA, AICc = NA, BIC = NA)  
  
  if(inherits(objs[[1]], c("lmerMod", "merMod"))){
    data.name <- as.character(deparse(objs[[1]]@call$data))
  }else{
    data.name <- as.character(deparse(objs[[1]]$call$data))  
  }
  
  for(i in seq_len(lobjs)){
    obj <- objs[[i]]
  
    if(inherits(obj, c("lmerMod", "merMod"))){
      if(data.name != as.character(deparse(obj@call$data))) 
        stop("All models should be fitted to the same data")
    }else{
      if(data.name != as.character(deparse(obj$call$data))) 
        stop("All models should be fitted to the same data")  
    }  
    
    ictab$model[i] <- nms[i]
    ictab$df[i] <- attr(stats::logLik(obj), "df")
    ictab$AIC[i] <- stats::AIC(obj)
    ictab$AICc[i] <- AICc(obj) ## This works for any object??
    ictab$BIC[i] <- stats::BIC(obj)
  }
  
  ## Calculating weights
  ## http://brianomeara.info/aic.html
  if(criteria == "AIC"){
    ictab$dAIC <- ictab$AIC - min(ictab$AIC)
    ictab$weight <- exp(-0.5 * ictab$dAIC) / sum(exp(-0.5 * ictab$dAIC))    
  } 
  if(criteria == "AICc"){
    ictab$dAICc <- ictab$AICc - min(ictab$AICc)
    ictab$weight <- exp(-0.5 * ictab$dAICc) / sum(exp(-0.5 * ictab$dAICc))    
  }
  if(criteria == "BIC"){
    ictab$dBIC <- ictab$BIC - min(ictab$BIC)
    ictab$weight <- exp(-0.5 * ictab$dBIC) / sum(exp(-0.5 * ictab$dBIC))    
  } 
  
  if(sort) ictab <- ictab[order(ictab$weight, decreasing = TRUE),]
  ictab
}

## Internal function to calculate small sample "corrected" AIC
AICc <- function(x){
  n <- stats::nobs(x)
  ## k <- length(coef(x)) + 1 ## Plus one is for sigma
  k <- attr(logLik(x), "df")
  cf <- (2 * k * (k + 1))/(n - k - 1)
  ans <- stats::AIC(x) + cf
  return(ans)
}

get_mnames <- function(x){
  mnames <- as.character(x)[-1]
  if (anyDuplicated(mnames)) 
        stop("model names must be distinct")
  mnames
}

#' @name predict2_gam
#' @rdname predict_nls
#' @description predict function for objects of class \code{\link[mgcv]{gam}}
#' @export
predict2_gam <- predict_nls

#' @title Modified prediciton function based on predict.gam
#' @name predict_gam
#' @description Largely based on predict.gam, but with some minor modifications to make it compatible
#' with \code{\link{predict_nls}}
#' @param object object of class \sQuote{gam} or as returned by function \sQuote{gamm}
#' @param newdata see \code{\link[mgcv]{predict.gam}}
#' @param type see \code{\link[mgcv]{predict.gam}}
#' @param se.fit see \code{\link[mgcv]{predict.gam}}. Notice that the default is changed to TRUE.
#' @param terms see \code{\link[mgcv]{predict.gam}}
#' @param exclude see \code{\link[mgcv]{predict.gam}}
#' @param block.size see \code{\link[mgcv]{predict.gam}}
#' @param newdata.guaranteed see \code{\link[mgcv]{predict.gam}}
#' @param na.action see \code{\link[mgcv]{predict.gam}}
#' @param unconditional see \code{\link[mgcv]{predict.gam}}
#' @param iterms.type see \code{\link[mgcv]{predict.gam}}
#' @param interval either \sQuote{none}, \sQuote{confidence} or \sQuote{prediction}.
#' @param level probability level for the interval (default 0.95)
#' @param tvalue t-value statistic used for constructing the intervals
#' @param ... additional arguments to be passed to \code{\link[mgcv]{predict.gam}}.
#' @return numeric vector of the same length as the fitted object when interval is equal to \sQuote{none}. 
#' Otherwise, a data.frame with columns named (for a 0.95 level) 
#' \sQuote{Estimate}, \sQuote{Est.Error}, \sQuote{Q2.5} and \sQuote{Q97.5}
#' @note this is a very simple wrapper for \code{\link[mgcv]{predict.gam}}.
#' @seealso \code{\link[stats]{predict.lm}}, \code{\link[stats]{predict.nls}}, \code{\link[mgcv]{predict.gam}}, \code{\link{simulate_nls}}, \code{\link{simulate_gam}}
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' require(mgcv)
#' data(barley)
#' 
#' fm.G <- gam(yield ~ s(NF, k = 6), data = barley)
#' 
#' ## confidence and prediction intervals
#' cis <- predict_gam(fm.G, interval = "conf")
#' pis <- predict_gam(fm.G, interval = "pred")
#' 
#' barleyA.ci <- cbind(barley, cis)
#' barleyA.pi <- cbind(barley, pis)
#' 
#' ggplot() + 
#'   geom_point(data = barleyA.ci, aes(x = NF, y = yield)) + 
#'   geom_line(data = barleyA.ci, aes(x = NF, y = Estimate)) + 
#'   geom_ribbon(data = barleyA.ci, aes(x = NF, ymin = Q2.5, ymax = Q97.5), 
#'               color = "red", alpha = 0.3) + 
#'   geom_ribbon(data = barleyA.pi, aes(x = NF, ymin = Q2.5, ymax = Q97.5), 
#'               color = "blue", alpha = 0.3) + 
#'   ggtitle("95% confidence and prediction bands")
#'   
#' }
predict_gam <- function(object, newdata=NULL, type="link", se.fit=TRUE, terms=NULL,
                        exclude=NULL, block.size=NULL, newdata.guaranteed=FALSE,
                        na.action=na.pass, unconditional=FALSE, iterms.type=NULL,
                        interval = c("none", "confidence", "prediction"),
                        level = 0.95, tvalue=NULL, ...){
  
  if(is.list(object) && inherits(object[[1]], "lme") && inherits(object[[2]], "gam")){
    object <- object[[2]]
  }
  
  if(!inherits(object, "gam")) stop("object should be of class 'gam'")
  
  interval <- match.arg(interval)
  
  if(missing(newdata)) newdata <- eval(object$call$data)

  prd <- mgcv::predict.gam(object, newdata=newdata, type=type, se.fit=se.fit,
                           terms=terms, exclude=exclude, block.size=block.size, 
                           newdata.guaranteed=newdata.guaranteed, na.action=na.action,
                           unconditional = unconditional, iterms.type = iterms.type, ...)
  
  if(interval == "none"){
    ans <- prd$fit
  }
  
  if(interval == "confidence"){
    if(missing(tvalue)){
      tvalue <- stats::qt(level, df = object$df.residual)
    }
    lwr <- prd$fit - tvalue * prd$se.fit
    upr <- prd$fit + tvalue * prd$se.fit
    ans <- data.frame(Estimate = prd$fit, Est.Error = prd$se.fit, lwr = lwr, upr = upr)
    names(ans) <- c("Estimate", "Est.Error", paste0("Q", 100*(1 - level)/2), paste0("Q", 100 * (1 - (1 - level)/2)))
  }
  
  if(interval == "prediction"){
    if(missing(tvalue)){
      tvalue <- stats::qt(level, df = object$df.residual)
    }
    lwr <- prd$fit - tvalue * sqrt(prd$se.fit^2 + sigma(object)^2)
    upr <- prd$fit + tvalue * sqrt(prd$se.fit^2 + sigma(object)^2)
    ans <- data.frame(Estimate = prd$fit, Est.Error = prd$se.fit, lwr = lwr, upr = upr)
    names(ans) <- c("Estimate", "Est.Error", paste0("Q", 100*(1 - level)/2), paste0("Q", 100 * (1 - (1 - level)/2)))
  }
  return(ans)
}

#' The method used in this function is described in Battes and Watts (2007)
#' Nonlinear Regression Analysis and Its Applications (see pages 58-59). It is 
#' known as the Delta method.
#' 
#' This method is approximate and it works better when the distribution of the 
#' parameter estimates are normally distributed. This assumption can be evaluated
#' by using bootstrap.
#' 
#' The method currently works well for any nonlinear function, but if predictions
#' are needed on new data, then it is required that a selfStart function is used.
#' 
#' 
#' @title Prediction Bands for Nonlinear Regression
#' @name predict2_nls
#' @param object object of class \sQuote{nls}
#' @param newdata data frame with values for the predictor
#' @param interval either \sQuote{none}, \sQuote{confidence} or \sQuote{prediction}
#' @param level probability level (default is 0.95)
#' @return a data frame with Estimate, Est.Error, lower interval bound and upper 
#' interval bound. For example, if the level = 0.95, 
#' the lower bound would be named Q2.5 and the upper bound would be name Q97.5
#' @export
#' @seealso \code{\link{predict.nls}} and \code{\link{predict_nls}}
#' @examples 
#' \donttest{
#' require(ggplot2)
#' require(nlme)
#' data(Soybean)
#' 
#' SoyF <- subset(Soybean, Variety == "F" & Year == 1988)
#' fm1 <- nls(weight ~ SSlogis(Time, Asym, xmid, scal), data = SoyF)
#' ## The SSlogis also supplies analytical derivatives
#' ## therefore the predict function returns the gradient too
#' prd1 <- predict(fm1, newdata = SoyF)
#' 
#' ## Gradient
#' head(attr(prd1, "gradient"))
#' ## Prediction method using gradient
#' prds <- predict2_nls(fm1, interval = "conf")
#' SoyFA <- cbind(SoyF, prds)
#' ggplot(data = SoyFA, aes(x = Time, y = weight)) + 
#'    geom_point() + 
#'    geom_line(aes(y = Estimate)) + 
#'    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) +
#'    ggtitle("95% Confidence Bands")
#'   
#' ## This is equivalent 
#' fm2 <- nls(weight ~ Asym/(1 + exp((xmid - Time)/scal)), data = SoyF,
#'            start = c(Asym = 20, xmid = 56, scal = 8))
#'            
#' ## Prediction interval
#' prdi <- predict2_nls(fm1, interval = "pred")
#' SoyFA.PI <- cbind(SoyF, prdi) 
#' ## Make prediction interval plot
#' ggplot(data = SoyFA.PI, aes(x = Time, y = weight)) + 
#'    geom_point() + 
#'    geom_line(aes(y = Estimate)) + 
#'    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) + 
#'    ggtitle("95% Prediction Band")
#'    
#' ## For these data we should be using gnls instead with an increasing variance
#' fmg1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), 
#'              data = SoyF, weights = varPower())
#'              
#' IC_tab(fm1, fmg1)
#' prdg <- predict_gnls(fmg1, interval = "pred")
#' SoyFA.GPI <- cbind(SoyF, prdg) 
#' 
#' ## These prediction bands are not perfect, but they could be smoothed
#' ## to eliminate the ragged appearance
#'  ggplot(data = SoyFA.GPI, aes(x = Time, y = weight)) + 
#'    geom_point() + 
#'    geom_line(aes(y = Estimate)) + 
#'    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) + 
#'    ggtitle("95% Prediction Band. NLS model which \n accomodates an increasing variance")
#' }
predict2_nls <- function(object, newdata = NULL, 
                         interval = c("none", "confidence", "prediction"), 
                         level = 0.95){
  
  if(!inherits(object, "nls"))
    stop("This function is only for objects of class 'nls'")
  
  interval <- match.arg(interval)
  npar <- length(coef(object))
  degf <- length(fitted(object)) - npar
  
  prd <- predict(object, newdata = newdata) ## predictions plus gradient
  se.prd <- NA
  Lwr <- NA
  Upr <- NA
  if(interval != "none"){
    grd <- attr(prd, "gradient") ## retrieves the gradient  
    if(is.null(grd)){
      if(is.null(newdata)){
        grd <- object$m$gradient()
        colnames(grd) <- names(coef(object))
      }else{
        ## Need to approximate the gradient numerically
        form <- object$m$formula()
        grd0 <- with(newdata, numericDeriv(form[[3]], names(coef(object))))
        grd <- attr(grd0, "gradient")
      }
    }
    Rmat <- object$m$Rmat() ## R matrix
    vtR <- grd %*% solve(Rmat) ## gradient times R^-1
    vtR.L <- sqrt(apply(vtR^2, 1, sum)) ## length of vector
    se.prd <- sigma(object) * vtR.L 
    ci.width <- se.prd * sqrt(npar * stats::qf(level, npar, degf))
    Lwr <- as.vector(prd) - ci.width
    Upr <- as.vector(prd) + ci.width
  }
  
  if(interval == "prediction"){
    ci.width <- sqrt(se.prd^2 + sigma(object)^2) * sqrt(npar * stats::qf(level, npar, degf))
    Lwr <- as.vector(prd) - ci.width
    Upr <- as.vector(prd) + ci.width
  }
  
  ans <- data.frame(Estimate = as.vector(prd), Est.Error = se.prd, 
                    Lwr = Lwr, Upr = Upr)
  names(ans) <- c("Estimate", "Est.Error", 
                  paste0("Q", 100 * (1 - level)/2), 
                  paste0("Q", 100 - 100 * (1 - level)/2))
  return(ans)
}