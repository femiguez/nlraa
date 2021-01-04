#'
#' Simulate responses from a linear model \code{\link[stats]{lm}}
#' 
#' @title Simulate responses from a linear model \code{\link[stats]{lm}}
#' @name simulate_lm
#' @description The function \code{\link[stats]{simulate}} does not consider the 
#' uncertainty in the estimation of the model parameters. This function will attempt 
#' to do this.
#' @param object object of class \code{\link[stats]{lm}}
#' @param psim parameter simulation level (an integer, 0, 1, 2, 3 or 4).
#' @param nsim number of simulations to perform
#' @param resid.type type of residual to include (none, resample, normal or wild)
#' @param value either \sQuote{matrix} or \sQuote{data.frame}
#' @param data the data argument might be needed when using this function inside user defined functions.
#' At least it is expected to be safer.
#' @param ... additional arguments (none used at the moment)
#' @return matrix or data.frame with responses
#' @details 
#' These are the options that control the parameter simulation level
#' \describe{
#'   \item{psim = 0}{returns the fitted values}
#'   \item{psim = 1}{simulates a beta vector (mean response)}
#'   \item{psim = 2}{simulates a beta vector and adds resampled residuals (similar to observed data)}
#'   \item{psim = 3}{simulates a beta vector, considers uncertainty in the variance covariance matrix of beta and adds residuals (prediction)}
#'   \item{psim = 4}{only adds residuals according to resid.type (similar to simulate.lm)}
#'  }
#' The residual type (resid.type) controls how the residuals are generated. 
#' They are either resampled, simulated from a normal distribution or \sQuote{wild} where the
#' Rademacher distribution is used (\url{https://en.wikipedia.org/wiki/Rademacher_distribution}).
#' Resampled and normal both assume iid, but \sQuote{normal} makes the stronger assumption of normality.
#' When psim = 2 and resid.type = none, \code{\link{simulate}} is used instead.
#' \sQuote{wild} does not assume constant variance, but it assumes symmetry.
#' @references See
#' \dQuote{Inference Based on the Wild Bootstrap} James G. MacKinnon
#' \url{https://www.math.kth.se/matstat/gru/sf2930/papers/wild.bootstrap.pdf}
#' \dQuote{Bootstrap in Nonstationary Autoregression} Zuzana Praskova 
#' \url{https://dml.cz/bitstream/handle/10338.dmlcz/135473/Kybernetika_38-2002-4_1.pdf}
#' \dQuote{Jackknife, Bootstrap and other Resampling Methods in Regression Analysis} C. F. J. Wu.
#' The Annals of Statistics. 1986. Vol 14. 1261-1295.
#' @export
#' @examples
#' \donttest{
#' require(ggplot2)
#' data(Orange)
#' fit <- lm(circumference ~ age, data = Orange)
#' sims <- simulate_lm(fit, nsim = 100, value = "data.frame")
#' 
#' ggplot(data = sims) + 
#'   geom_line(aes(x = age, y = sim.y, group = ii), 
#'             color = "gray", alpha = 0.5) + 
#'   geom_point(aes(x = age, y = circumference)) 
#' }
#' 
#' 

simulate_lm <- function(object, psim = 1, nsim = 1, 
                        resid.type = c("none", "resample","normal","wild"),
                        value = c("matrix","data.frame"), 
                        data = NULL, ...){

  if(!inherits(object, "lm"))
    stop("only for objects which inherit class 'lm' ")
  
  value <- match.arg(value)
  resid.type <- match.arg(resid.type)

  xargs <- list(...)

  if(!is.null(xargs$newdata) && psim > 1)
    stop("'newdata' is not compatible with 'psim > 1'")  
  
  if(is.null(xargs$newdata) && is.null(data)){
    newdata <- try(eval(getCall(object)$data), silent = TRUE)
    if(is.null(newdata)) stop("data not found. If you are using simulate_lm inside another function, please pass the data")
  }else{
    if(!is.null(xargs$newdata) && !is.null(data))
      stop("either supply 'data' or 'newdata'")
    
    if(!is.null(data)){
      newdata <- data
    }else{
      newdata <- xargs$newdata
    }    
  }
  
  nr <- ifelse(is.null(newdata), stats::nobs(object), nrow(newdata))
  ans.mat <- matrix(nrow = nr, ncol = nsim)  
  
  ## They say that replicate is faster than this,
  ## but not when storage is pre-allocated
  for(i in 1:nsim){
      ans.mat[,i] <- simulate_lm_one(object, psim = psim, resid.type = resid.type, 
                                     newdata = newdata)  
  }
  
  if(value == "matrix"){
    colnames(ans.mat) <- paste0("sim_",1:nsim)
    return(ans.mat)
  }else{
    if(is.null(xargs$newdata)){
      dfr <- eval(getCall(object)$data)
      if(is.null(dfr)) stop("'data' argument should be supplied")      
    }else{
      dfr <- newdata
    }
    ## I think this is guaranteed to exist
    ## but it can be weird
    ans.dat <- data.frame(ii = as.factor(rep(1:nsim, each = nrow(dfr))),
                          dfr,
                          sim.y = c(ans.mat), 
                          row.names = 1:c(nsim * nr))
    return(ans.dat)
  }
}

simulate_lm_one <- function(object, psim = 1, 
                            resid.type = c("none", "resample","normal","wild"), 
                            newdata = NULL){
  
  resid.type <- match.arg(resid.type)
  
  X <- stats::model.matrix(object, data = newdata)
  n <- stats::nobs(object)
  rsd0 <- stats::resid(object)
  
  if(is.null(newdata)) stop("'newdata' is needed for this function")
  
  if(resid.type %in% c("none", "resample")){
    rsds <- sample(rsd0, size = n, replace = TRUE)      
  }
  if(resid.type == "normal"){
    rsds <- stats::rnorm(n = n, mean = 0, sd = stats::sigma(object))
  }
  if(resid.type == "wild"){
    rsds <- sample(c(-1, 1), size = n, replace = TRUE) * rsd0
  }
  
  ## Nothing is simulated, fitted values
  if(psim == 0){
    betav <- stats::coef(object)
    ans <- X %*% betav
  }
  
  ## Simulate from coefficient vector only (beta)
  if(psim == 1){
    betav <- MASS::mvrnorm(mu = coef(object), Sigma = vcov(object))  
    ans <- X %*% betav
  }
  
  ## This is the 'default' simulate method, which also seems to work for 'glm' objects
  if(psim == 2 && resid.type == "none"){
    ans <- simulate(object, nsim = 1)[,1]
  }
    
  if(psim == 2 && resid.type != "none"){
    betav <- MASS::mvrnorm(mu = coef(object), Sigma = vcov(object))  
    ans <- X %*% betav + rsds
  }
  
  ## Simulate from betav with noise on the VC matrix
  if(psim == 3){
    rsgms <- stats::var(rsds)/sigma(object) ## Ratio of sigmas
    Sigma0 <- rsgms * stats::vcov(object) ## Is this guaranteed to be positive-definite? (I think so)
    betav <- MASS::mvrnorm(mu = coef(object), Sigma = Sigma0)  
    ans <- X %*% betav + rsds
  }
  
  ## What simulate.lm does is just simulate residuals (I think)
  ## This does not work for objects of class 'glm' or 'gam'
  if(psim == 4){
    if(inherits(object, "glm")) stop("'glm' object are not supported with this option")
    ans <- stats::fitted(object) + rsds 
  }
  
  return(ans)
} 