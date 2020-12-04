#' 
#' @title Average predictions from several (non)linear models based on IC weights
#' @name predict_nlme
#' @rdname predict_nlme
#' @description Computes weights based on AIC, AICc, or BIC and it generates weighted predictions by
#' the relative value of the IC values
#' @param ... nlme, lme, gls or gnls objects. 
#' @param criteria either \sQuote{AIC}, \sQuote{AICc} or \sQuote{BIC}.
#' @param interval either \sQuote{none}, \sQuote{confidence} or \sQuote{prediction}.
#' @param level probability level for the interval (default 0.95)
#' @param nsim number of simulations to perform for intervals. Default 1000.
#' @param plevel parameter level prediction to be passed to prediciton functions.
#' @param newdata new data frame for predictions
#' @return numeric vector of the same length as the fitted object.
#' @note all the nls or lm objects should be fitted to the same data. The weights are
#' based on the inverse of the IC value.
#' @seealso \code{\link{predict.lme}} \code{\link{predict.gnls}}
#' @export
#' @examples
#' \donttest{
#' ## Example
#' require(ggplot2)
#' require(nlme)
#' data(Orange)
#' 
#' ## All models should be fitted using Maximum Likelihood
#' fm.L <- nlme(circumference ~ SSlogis(age, Asym, xmid, scal), 
#'                 random = pdDiag(Asym + xmid + scal ~ 1), 
#'                 method = "ML", data = Orange)
#' fm.G <- nlme(circumference ~ SSgompertz(age, Asym, b2, b3), 
#'                 random = pdDiag(Asym + b2 + b3 ~ 1), 
#'                 method = "ML", data = Orange)
#' fm.F <- nlme(circumference ~ SSfpl(age, A, B, xmid, scal), 
#'                 random = pdDiag(A + B + xmid + scal ~ 1), 
#'                 method = "ML", data = Orange)
#' fm.B <- nlme(circumference ~ SSbg4rp(age, w.max, lt.e, ldtm, ldtb), 
#'                 random = pdDiag(w.max + lt.e + ldtm + ldtb ~ 1), 
#'                 method = "ML", data = Orange)
#'
#' ## Print the table with weights
#' IC_tab(fm.L, fm.G, fm.F, fm.B)
#' 
#' ## Each model prediction is weighted according to their AIC values
#' prd <- predict_nlme(fm.L, fm.G, fm.F, fm.B)
#' 
#' ggplot(data = Orange, aes(x = age, y = circumference)) + 
#'   geom_point() + 
#'   geom_line(aes(y = predict(fm.L, level = 0), color = "Logistic")) +
#'   geom_line(aes(y = predict(fm.G, level = 0), color = "Gompertz")) +
#'   geom_line(aes(y = predict(fm.F, level = 0), color = "4P-Logistic")) +  
#'   geom_line(aes(y = predict(fm.B, level = 0), color = "Beta")) +
#'   geom_line(aes(y = prd, color = "Avg. Model"), size = 1.2)
#' }

predict_nlme <- function(..., criteria = c("AIC", "AICc", "BIC"), 
                         interval = c("none", "confidence", "prediction"),
                         level = 0.95, nsim = 1e3, plevel = 0,
                         newdata = NULL){
  
  objs <- list(...)
  criteria <- match.arg(criteria)
  interval <- match.arg(interval)
  nms <- get_mnames(match.call())
  
  lobjs <- length(objs)
  wtab <- data.frame(model = character(lobjs), IC = NA)
  
  nr <- stats::nobs(objs[[1]])
  if(!is.null(newdata)) nr <- nrow(newdata)
  
  if(interval == "none"){
    prd.mat <- matrix(nrow = nr, ncol = lobjs)
  }else{
    prd.mat <- matrix(nrow = nr, ncol = lobjs)
    prd.mat.se <- matrix(nrow = nr, ncol = lobjs)
    prd.mat.lwr <- matrix(nrow = nr, ncol = lobjs)
    prd.mat.upr <- matrix(nrow = nr, ncol = lobjs)
  } 
  
  data.name <- as.character(objs[[1]]$call$data)
  
  for(i in seq_len(lobjs)){
    obj <- objs[[i]]
    if(!inherits(obj, c("nlme","lme","gnls","gls"))) stop("All objects should be of class 'nlme', 'lme', 'gnls' or 'gls'")
    if(data.name != as.character(obj$call$data)) stop("All models should be fitted to the same data")
    
    wtab$model[i] <- nms[i]
    
    if(criteria == "AIC") wtab$IC[i] <- stats::AIC(obj)
    if(criteria == "AICc") wtab$IC[i] <- AICc(obj)
    if(criteria == "BIC") wtab$IC[i] <- stats::BIC(obj)
  }
  
  ## Predictions
  if(interval == "none"){
    for(i in seq_len(lobjs)){
      obj <- objs[[i]]
      
      if(inherits(obj, "lme")){
        if(is.null(newdata)){
          prd.mat[,i] <- predict(obj, level = plevel)  
        }else{
          prd.mat[,i] <- predict(obj, newdata = newdata, level = plevel)
        }
      }
        
      if(inherits(obj, "gls")){
        if(is.null(newdata)){
          prd.mat[,i] <- predict(obj)
        }else{
          prd.mat[,i] <- predict(obj, newdata = newdata)  
        }
      }
    }
  }
  
  if(interval == "confidence" || interval == "prediction"){
    
    psim <- ifelse(interval == "confidence", 1, 2)
    lb <- (1 - level)/2
    ub <- 1 - lb 
    
    for(i in seq_len(lobjs)){
      obj <- objs[[i]]
      
      if(inherits(obj, "lme")){
        if(inherits(obj, "nlme")){
          tmp.sim <- simulate_nlme(obj, psim = psim, nsim = nsim, 
                                   level = plevel, newdata = newdata)           
        }else{
          tmp.sim <- simulate_lme(obj, psim = psim, nsim = nsim, 
                                   level = plevel, newdata = newdata)           
        }
      } 
        
      if(inherits(obj, "gls")){
        if(inherits(obj, "gnls")){
          tmp.sim <- simulate_nlme(obj, psim = psim, nsim = nsim, newdata = newdata)  
        }else{
          tmp.sim <- simulate_lme(obj, psim = psim, nsim = nsim, newdata = newdata)  
        }
      } 
        
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


