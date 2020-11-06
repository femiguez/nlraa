#' 
#' @title self start for Ricker Function
#' @name SSricker
#' @rdname SSricker
#' @description Self starter for Ricker function with parameters a and b
#' @param time input vector (x) which is normally \sQuote{time}, the smallest value should be close to zero.
#' @param a which is related to the initial growth slope
#' @param b which is related to the slowing down or decline
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506) and originally in Ricker, W. E. (1954) Stock and Recruitment Journal of the Fisheries Research Board of Canada, 11(5): 559–623. (doi:10.1139/f54-039).
#' The equation is: \eqn{a * time * exp(-b * time)}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- 30 * x * exp(-0.3 * x) + rnorm(30, 0, 0.25)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSricker(x, a, b), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

rickerInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a ricker")
  }
  ## Use the log transform
  xy <- subset(xy, xy[,"y"] > 0 & xy[,"x"] > 0)
  
  ry <- log(xy[,"y"]) - log(xy[,"x"])
  fit <- stats::lm(ry ~ xy[,"x"])
  a <- exp(coef(fit)[1])
  b <- -coef(fit)[2]
  value <- c(a, b)
  names(value) <- mCall[c("a","b")]
  value
  
}

#' @rdname SSricker
#' @return ricker: vector of the same length as x (time) using the ricker function
#' @export
ricker <- function(time, a, b){
  
  .value <- a * time * exp(-b * time)
  
  ## Derivative with respect to a
  ## .exp1 <- deriv(~ a * time * exp(-b * time), "a")
  .exp1 <- time * exp(-b * time)
  
  ## Derivative with respect to b
  ## .exp2 <- deriv(~ a * time * exp(-b * time), "b")
  .exp2 <- -(a * time * (exp(-b * time) * time))
  
  .actualArgs <- as.list(match.call()[c("a", "b")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a", "b")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSricker
#' @export
SSricker <- selfStart(ricker, initial = rickerInit, c("a", "b"))

