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
#' @return numeric vector of the same length as the fitted object.
#' @note all the \code{\link[stats]{nls}} or \code{\link[stats]{lm}} 
#' objects should be fitted to the same data. Weights are
#' based on the chosen IC value (exp(-0.5 * IC)). 
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
#' fm.BL <- nls(yield ~ SSblin(NF, a, b, xs, c), data = barley)
#' fm.G <- gam(yield ~ NF + s(NF^2, k = 3), data = barley)
#' 
#' ## Print the table with weights
#' IC_tab(fm.L, fm.Q, fm.A, fm.LP, fm.BL, fm.G)
#' 
#' ## Each model prediction is weighted according to their AIC values
#' prd <- predict_nls(fm.L, fm.Q, fm.A, fm.LP, fm.BL, fm.G)
#' 
#' ggplot(data = barley, aes(x = NF, y = yield)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fm.L), color = "Linear")) +
#'   geom_line(aes(y = fitted(fm.Q), color = "Quadratic")) +
#'   geom_line(aes(y = fitted(fm.A), color = "Asymptotic")) +  
#'   geom_line(aes(y = fitted(fm.LP), color = "Linear-plateau")) + 
#'   geom_line(aes(y = fitted(fm.BL), color = "Bi-linear")) + 
#'   geom_line(aes(y = fitted(fm.G), color = "GAM")) + 
#'   geom_line(aes(y = prd, color = "Avg. Model"), size = 1.2)
#' }

predict_nls <- function(..., criteria = c("AIC", "AICc", "BIC"), 
                        interval = c("none", "confidence", "prediction"),
                        level = 0.95, nsim = 1e3,
                        resid.type = c("none", "resample", "normal", "wild"),
                        newdata = NULL){
  
  ## all objects should be of class 'nls' or inherit 'lm'
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
  
  ## Predictions
  if(interval == "none"){
    for(i in seq_len(lobjs)){
      nls.obj <- nls.objs[[i]]
      if(!is.null(newdata)){
        prd.mat[,i] <- predict(nls.obj, newdata = newdata)  
      }else{
        prd.mat[,i] <- predict(nls.obj)  
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
      
      if(inherits(nls.obj, "glm")) 
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
  wtab$weight <- exp(-0.5 * wtab$dIC) / sum(exp(-0.5 * wtab$dIC))
  
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
  data.name <- as.character(objs[[1]]$call$data)
  
  for(i in seq_len(lobjs)){
    obj <- objs[[i]]
    if(data.name != as.character(obj$call$data)) 
      stop("All models should be fitted to the same data")
    
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


#' @name predict_gam
#' @rdname predict_nls
#' @description predict function for objects of class \code{\link[mgcv]{gam}}
#' @export
predict_gam <- predict_nls
