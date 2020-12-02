#' 
#' @title Average predictions from several (non)linear models based on IC weights
#' @name predict_nls_mm
#' @rdname predict_nls_mm
#' @description Computes weights based on AIC, AICc, or BIC and it generates weighted predictions by
#' the relative value of the IC values
#' @param ... nls or lm objects. 
#' @param criteria either \sQuote{AIC}, \sQuote{AICc} or \sQuote{BIC}.
#' @return numeric vector of the same length as the fitted object.
#' @note all the nls or lm objects should be fitted to the same data. The weights are
#' based on the inverse of the IC value.
#' @export
#' @examples
#' \donttest{
#' ## Example
#' require(ggplot2)
#' data(barley, package = "nlraa")
#' 
#' fm.L <- lm(yield ~ NF, data = barley)
#' fm.Q <- lm(yield ~ NF + I(NF^2), data = barley)
#' fm.A <- nls(yield ~ SSasymp(NF, Asym, R0, lrc), data = barley)
#' fm.LP <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
#' fm.BL <- nls(yield ~ SSblin(NF, a, b, xs, c), data = barley)
#'
#' ## Print the table with weights
#' AIC_tab(fm.L, fm.Q, fm.A, fm.LP, fm.BL)
#' 
#' ## Each model prediction is weighted by the inverse of their AIC values
#' prd <- predict_nls_mm(fm.L, fm.Q, fm.A, fm.LP, fm.BL)
#' 
#' ggplot(data = barley, aes(x = NF, y = yield)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fm.L), color = "Linear")) +
#'   geom_line(aes(y = fitted(fm.Q), color = "Quadratic")) +
#'   geom_line(aes(y = fitted(fm.A), color = "Asymptotic")) +  
#'   geom_line(aes(y = fitted(fm.LP), color = "Linear-plateau")) + 
#'   geom_line(aes(y = fitted(fm.BL), color = "Bi-linear")) + 
#'   geom_line(aes(y = prd, color = "Avg. Model"), size = 1.2)
#' }

predict_nls_mm <- function(..., criteria = c("AIC", "AICc", "BIC")){
  
  ## all objects should be of class 'nls'
  ## 1. Create table of AIC based weights
  nls.objs <- list(...)
  
  criteria <- match.arg(criteria)

  nls.nms <- get_mnames(match.call())
  
  lobjs <- length(nls.objs)
  
  wtab <- data.frame(model = character(lobjs), IC = NA)
  
  lprd <- length(fitted(nls.objs[[1]]))
  prd.mat <- matrix(nrow = lprd, ncol = lobjs)
  
  data.name <- as.character(nls.objs[[1]]$call$data)
  
  for(i in 1:lobjs){
    nls.obj <- nls.objs[[i]]
    
    if(!inherits(nls.obj, c("nls","lm"))) stop("All objects should be of class 'nls' or 'lm'")
    if(data.name != as.character(nls.obj$call$data)) stop("All models should be fitted to the same data")
    
    wtab$model[i] <- nls.nms[i]
    
    if(criteria == "AIC") wtab$IC[i] <- stats::AIC(nls.obj)
    if(criteria == "AICc") wtab$IC[i] <- AICc_nls(nls.obj)
    if(criteria == "BIC") wtab$IC[i] <- stats::BIC(nls.obj)
  
    ## Predictions
    prd.mat[,i] <- predict(nls.obj)
  }

  wtab$dIC <- wtab$IC - min(wtab$IC)
  wtab$weight <- exp(-0.5 * wtab$dIC) / sum(exp(-0.5 * wtab$dIC))
  prd <- rowSums(sweep(prd.mat, 2, wtab$weight, "*"))
  return(prd)
}


#' @rdname predict_nls_mm
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
  
  for(i in 1:lobjs){
    obj <- objs[[i]]
    if(data.name != as.character(obj$call$data)) 
      stop("All models should be fitted to the same data")
    
    ictab$model[i] <- nms[i]
    ictab$df[i] <- attr(stats::logLik(obj), "df")
    ictab$AIC[i] <- stats::AIC(obj)
    ictab$AICc[i] <- AICc_nls(obj)
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
AICc_nls <- function(x){
  n <- stats::nobs(x)
  k <- length(coef(x)) + 1 ## Plus one is for sigma
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
