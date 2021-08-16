#' 
#' The equation for this function is:
#' 
#' \deqn{f(x) = asym2 + (asym1 - asym2)/(1 + exp(iscal * (log(x) - log(xmid))))^theta}
#' 
#' @title self start for five-parameter logistic function
#' @name SSlogis5
#' @rdname SSlogis5
#' @description Self starter for a five-parameter logistic function.
#' @param x input vector (x) 
#' @param asym1 asymptotic value for low values of x
#' @param asym2 asymptotic value for high values of x
#' @param iscal steepness of transition from asym1 to asym2 (inverse of the scale)
#' @param xmid value of x at which y = (asym1 + asym2)/2 (only when theta = 1)
#' @param theta asymmetry parameter, if it is equal to 1, this is the four parameter logistic
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @details This is known as the Richards' function or the log-logistic and it is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506).
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- seq(0, 2000, 100)
#' y <- logis5(x, 35, 10, 800, 5, 2) + rnorm(length(x), 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSlogis5(x, asym1, asym2, xmid, iscal, theta), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

logis5Init <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 5){
    stop("Too few distinct input values to fit a logis5")
  }
  asym1 <- xy[1,"y"]
  asym2 <- xy[nrow(xy),"y"]
  xmid <- NLSstClosestX(xy, mean(c(asym1, asym2)))
  iscal <- 1/(max(xy[,"x"], na.rm = TRUE) - xmid)
  theta <- 1
  
  objfun <- function(cfs){
    pred <- logis5(xy[,"x"], asym1=cfs[1], asym2=cfs[2], xmid=cfs[3], iscal=cfs[4], theta = cfs[5])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(asym1, asym2, xmid, iscal, theta)
  op <- try(stats::optim(cfs, objfun), silent = TRUE)
  
  if(!inherits(op, "try-error")){
    asym1 <- op$par[1]
    asym2 <- op$par[2]
    xmid <- op$par[3]
    iscal <- op$par[4]
    theta <- op$par[5]
  }
  
  value <- c(asym1, asym2, xmid, iscal, theta)
  names(value) <- mCall[c("asym1","asym2","xmid","iscal","theta")]
  value
  
}

#' @rdname SSlogis5
#' @return logis5: vector of the same length as x (time) using the 5-parameter logistic
#' @examples 
#' x <- seq(0, 2000)
#' y <- logis5(x, 30, 10, 800, 5, 2)
#' plot(x, y)
#' @export
#' 
logis5 <- function(x, asym1, asym2, xmid, iscal, theta){
  
  if(any(x < 0))
    stop("Input (x) should be positive for this equation", call. = FALSE)
  
  .lxmid <- suppressWarnings(log(xmid))
  .expre1 <- 1 + exp(iscal * (log(x) - .lxmid))
  .expre2 <- .expre1^theta
  .expre3 <- (asym1 - asym2) / .expre2
  .value <- asym2 + .expre3 
  .value <- ifelse(is.nan(.value),0,.value)
  
  ## Derivative with respect to asym1
  ## deriv(~asym2 + (asym1 - asym2)/((1 + exp(iscal * (log(x) - log(xmid))))^theta),"asym1")
  .expr8 <- (1 + exp(iscal * (log(x) - .lxmid)))^theta
  .exp1 <- ifelse(is.nan(.expr8), 0, 1/.expr8)
   
  ## Derivative with respect to asym2
  ## deriv(~asym2 + (asym1 - asym2)/((1 + exp(iscal * (log(x) - log(xmid))))^theta),"asym2")          
  .exp2 <- ifelse(is.nan(.expr8), 0, 1 - 1/.expr8)
  
  ## Derivative with respect to xmid
  ## deriv(~asym2 + (asym1 - asym2)/((1 + exp(iscal * (log(x) - log(xmid))))^theta),"xmid")          
  .expr1 <- asym1 - asym2
  .expr6 <- exp(iscal * (log(x) - .lxmid))
  .expr7 <- 1 + .expr6
  .exp3 <- .expr1 * (.expr7^(theta - 1) * (theta * (.expr6 * (iscal * (1/xmid)))))/.expr8^2
  .exp3 <- ifelse(is.nan(.exp3), 0, .exp3)
  
  ## Derivative with respect to iscal
  ## deriv(~asym2 + (asym1 - asym2)/((1 + exp(iscal * (log(x) - log(xmid))))^theta),"iscal") 
  .expr4 <- log(x) - .lxmid
  .exp4 <- -(.expr1 * (.expr7^(theta - 1) * (theta * (.expr6 * .expr4)))/.expr8^2)
  .exp4 <- ifelse(is.nan(.exp4), 0, .exp4)
  
  ## Derivative with respect to theta
  ## deriv(~asym2 + (asym1 - asym2)/((1 + exp(iscal * (log(x) - log(xmid))))^theta),"theta") 
  .exp5 <- -(.expr1 * (.expr8 * log(.expr7))/.expr8^2)
  .exp5 <- ifelse(is.nan(.exp5),0,.exp5)
  
  .actualArgs <- as.list(match.call()[c("asym1", "asym2", "xmid", "iscal", "theta")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 5L), list(NULL, c("asym1", "asym2", "xmid","iscal","theta")))
    .grad[, "asym1"] <- .exp1
    .grad[, "asym2"] <- .exp2
    .grad[, "xmid"] <- .exp3
    .grad[, "iscal"] <- .exp4
    .grad[, "theta"] <- .exp5
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSlogis5
#' @export
SSlogis5 <- selfStart(logis5, initial = logis5Init, c("asym1", "asym2", "xmid", "iscal","theta"))
