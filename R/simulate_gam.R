#'
#' This function is probably redundant. Simply using \code{\link{simulate}}
#' generates data from the correct distribution for objects which inherit 
#' class \code{\link[stats]{lm}}. The difference is that I'm trying to add the 
#' uncertainty in the parameter estimates.
#' 
#' @title Simulate responses from a generalized additive linear model \code{\link[mgcv]{gam}}
#' @name simulate_gam
#' @description By sampling from the vector of coefficients it is possible to simulate
#' data from a \sQuote{gam} model. 
#' @param object object of class \code{\link[mgcv]{gam}} or \code{\link[stats]{glm}}.
#' @param psim parameter simulation level (an integer, 0, 1, 2 or 3).
#' @param nsim number of simulations to perform
#' @param resid.type type of residual to use. \sQuote{none}, \sQuote{resample}, \sQuote{normal} or \sQuote{wild}.
#' @param value either \sQuote{matrix} or \sQuote{data.frame}
#' @param ... additional arguments (none used at the moment)
#' @note psim = 3 is not implemented at the moment.
#' @return matrix or data.frame with responses
#' @details 
#' These are the options that control the parameter simulation level
#' \describe{
#'   \item{psim = 0}{returns the fitted values}
#'   \item{psim = 1}{simulates from a beta vector (mean response)}
#'   \item{psim = 2}{simulates observations according to the residual type (similar to observed data)}
#'   \item{psim = 3}{simulates a beta vector, considers uncertainty in the variance covariance matrix of beta and adds residuals (prediction)}
#'  }
#' The residual type (resid.type) controls how the residuals are generated. 
#' They are either resampled, simulated from a normal distribution or \sQuote{wild} where the
#' Rademacher distribution is used (\url{https://en.wikipedia.org/wiki/Rademacher_distribution}).
#' Resampled and normal both assume iid, but \sQuote{normal} makes the stronger assumption of normality.
#' \sQuote{wild} does not assume constant variance, but it assumes symmetry.
#' @seealso \code{\link{predict}}, \code{\link[mgcv]{predict.gam}}, \code{\link{simulate}} and \code{\link{simulate_lm}}.
#' @references Generalized Additive Models. An Introduction with R. Second Edition. (2017) Simon N. Wood. CRC Press. 
#' @note The purpose of this function is to make it compatible with other functions in this
#' package. It has some limitations compared to the functions in the \sQuote{see also} section.
#' @export
#' @examples
#' \donttest{
#' require(ggplot2)
#' require(mgcv)
#' ## These count data are from GAM book by Simon Wood (pg. 132) - see reference
#' y <- c(12, 14, 33, 50, 67, 74, 123, 141, 165, 204, 253, 246, 240)
#' t <- 1:13
#' dat <- data.frame(y = y, t = t)
#' fit <- gam(y ~ t + I(t^2), family = poisson, data = dat)
#' sims <- simulate_gam(fit, nsim = 100, value = "data.frame")
#' 
#' ggplot(data = sims) + 
#'   geom_line(aes(x = t, y = sim.y, group = ii), 
#'             color = "gray", alpha = 0.5) + 
#'   geom_point(aes(x = t, y = y)) 
#' }
#' 

simulate_gam <- function(object, nsim = 1, psim = 1, 
                         resid.type = c("none", "resample", "normal", "wild"),
                         value = c("matrix", "data.frame"), ...){
  
  if(!inherits(object, "glm")) stop("Object should be of class 'glm' or 'gam'")
  
  resid.type <- match.arg(resid.type)
  value <- match.arg(value)
  
  xargs <- list(...)
  if(is.null(xargs$newdata)){
    newdata <- NULL 
  }else{
    newdata <- xargs$newdata
  } 
  
  if(is.null(xargs$exclude)){
    exclude <- NULL 
  }else{
    exclude <- xargs$exclude
  } 

  if(!is.null(newdata) && psim > 1)
    stop("'newdata' is not compatible with 'psim > 1'")
  
  nr <- ifelse(is.null(newdata), stats::nobs(object), nrow(newdata))
  ans.mat <- matrix(nrow = nr, ncol = nsim)  
  
  if(!is.null(newdata) && value == "data.frame")
    stop("'newdata' is not compatible with 'value = data.frame'")
  
  ## They say that replicate is faster than this,
  ## but not when storage is pre-allocated
  for(i in 1:nsim){
    ans.mat[,i] <- simulate_gam_one(object, psim = psim, resid.type = resid.type, 
                                   newdata = newdata, exclude = exclude)  
  }
  
  if(value == "matrix"){
    colnames(ans.mat) <- paste0("sim_",1:nsim)
    return(ans.mat)
  }else{
    dfr <- eval(getCall(object)$data)
    if(is.null(dfr)) stop("data argument should be supplied")
    ## I think this is guaranteed to exist
    ## but it can be weird
    ans.dat <- data.frame(ii = as.factor(rep(1:nsim, each = nrow(dfr))),
                          dfr,
                          sim.y = c(ans.mat), 
                          row.names = 1:c(nsim * nr))
    return(ans.dat)
  }
}

simulate_gam_one <- function(object, psim = 1, 
                             resid.type = c("none", "resample", "normal", "wild"),
                             newdata = NULL, exclude = NULL){
  
  resid.type <- match.arg(resid.type)
  
  n <- stats::nobs(object)
  
  if(!is.null(newdata)){
    ## See section 7.2.6 in Wood's GAM book. pg. 339
    ## Check that names are correct, Is this needed?
    if(!identical(names(newdata), attr(object$terms, "term.labels")))
      stop("names in 'newdata' do not correspond to 'term.labels'")
    X <- predict(object, newdata = newdata, type = "lpmatrix", exclude = exclude)
  }else{
    X <- predict(object, type = "lpmatrix", exclude = exclude)
  }
  
  if(resid.type == "resample"){
    rsd0 <- residuals(object, type = "scaled.pearson") ## Are these ~ N(0, 1)?
    rsd.sd <- sqrt(object$family$variance(fitted(object)))
    rsds <- sample(rsd0, size = n, replace = TRUE) * rsd.sd     
  }
  
  if(resid.type == "normal"){
    rsd.sd <- sqrt(object$family$variance(fitted(object)))
    rsds <- rnorm(n = n, sd = rsd.sd)
  }
  
  if(resid.type == "wild"){
    rsd0 <- residuals(object, type = "response")
    rsds <- sample(c(-1, 1), size = n, replace = TRUE) * rsd0
  }
  
  ## Nothing is simulated, fitted values
  if(psim == 0){
    betav <- stats::coef(object)
    ans <- object$family$linkinv(X %*% betav)
  }
  
  ## Simulate from coefficient vector only (beta)
  if(psim == 1){
    betav <- MASS::mvrnorm(mu = coef(object), Sigma = vcov(object))  
    ans <- object$family$linkinv(X %*% betav)
  }
  
  ## Simulate using default simulate.lm method
  ## I think this works because simulate automatically works in the 
  ## response scale
  if(psim == 2 && resid.type == "none"){
    ans <- simulate(object, nsim = 1)[,1]
  }
  
  if(psim == 2 && resid.type != "none"){
    betav <- MASS::mvrnorm(mu = coef(object), Sigma = vcov(object))  
    ans <- object$family$linkinv(X %*% betav) + rsds
    ## not sure for how many families only a discrete response makes sense
    if(object$family$family == "poisson") ans <- round(ans)
  }
  
  if(psim == 3){
    stop("Not implemented yet")
  }
  
  return(ans)
}