#' 
#' @title Average predictions from several (non)linear models based on IC weights
#' @name predict_nls_mm
#' @rdname predict_nls_mm
#' @description Computes weights based on AIC or BIC and it generates weighted predictions by
#' the relative value of the AIC or BIC values
#' @param ... nls or lm objects. 
#' @param criteria either \sQuote{AIC} or \sQuote{BIC}
#' @param print.table wether to print the table with IC and weights
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
#' IC_tab(fm.L, fm.Q, fm.A, fm.LP, fm.BL)
#' 
#' ## Each model prediction is weighted by the inverse of their AIC values
#' prd <- predict_nls_mm(fm.L, fm.Q, fm.A, fm.LP, fm.BL, print.table = TRUE)
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

predict_nls_mm <- function(..., criteria = c("AIC","BIC"), print.table = FALSE){
  
  criteria <- match.arg(criteria)
  ## all objects should be of class 'nls'
  ## 1. Create table of AIC based weights
  nls.objs <- list(...)

  if(!requireNamespace("bbmle", quietly = TRUE)){
    warning("bbmle is required for this function")
    return(NULL)
  }
  
  nls.nms <- bbmle:::get.mnames(match.call())
  
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
    
    if(criteria == "AIC") wtab$IC[i] <- AIC(nls.obj)
    if(criteria == "BIC") wtab$IC[i] <- BIC(nls.obj)
    
    ## Predictions
    prd.mat[,i] <- predict(nls.obj)
  }

  wtab$weight <- 1/wtab$IC / sum(1/wtab$IC)
  
  if(print.table){
    names(wtab) <- c("model", criteria, "weights")
    print(wtab)
  } 
  
  prd <- rowSums(prd.mat * wtab$weight)
  return(prd)
}


#' @rdname predict_nls_mm
#' @description Information criteria table with weights
#' @param ... model fit objects fitted to the same data
#' @export
#' 
IC_tab <- function(...){
  
  objs <- list(...)
  
  if(!requireNamespace("bbmle", quietly = TRUE)){
    warning("bbmle is required for this function")
    return(NULL)
  }
  
  nms <- bbmle:::get.mnames(match.call())
  
  lobjs <- length(objs)
  
  ictab <- data.frame(model = character(lobjs), AIC = NA, BIC = NA)  
  
  data.name <- as.character(objs[[1]]$call$data)
  
  for(i in 1:lobjs){
    obj <- objs[[i]]
    if(data.name != as.character(obj$call$data)) 
      stop("All models should be fitted to the same data")
    
    ictab$model[i] <- nms[i]
    
    ictab$AIC[i] <- AIC(obj)
    ictab$BIC[i] <- BIC(obj)
  }
  
  ictab$AIC.weight <- 1/ictab$AIC / sum(1/ictab$AIC)
  ictab$BIC.weight <- 1/ictab$BIC / sum(1/ictab$BIC)
  
  ictab
}
