#' 
#' J. GOUDRIAAN, J. L. MONTEITH, A Mathematical Function for Crop Growth Based 
#' on Light Interception and Leaf Area Expansion, Annals of Botany, Volume 66, 
#' Issue 6, December 1990, Pages 695–701,
#'  \url{https://doi.org/10.1093/oxfordjournals.aob.a088084}
#'  
#' The equation is: \deqn{ (cm/rm) * log(1 + exp(rm * (t - tb)))}
#' 
#' @title self start for the exponential-linear growth equation
#' @name SSexplin
#' @rdname SSexplin
#' @description Self starter for an exponential-linear growth equation
#' @param t input vector (time) 
#' @param cm parameter related to the maximum growth during the linear phase
#' @param rm parameter related to the maximum growth during the exponential phase
#' @param tb time at which switch happens
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506).  
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(12345)
#' x <- seq(1,100, by = 5)
#' y <- explin(x, 20, 0.14, 30) + rnorm(length(x), 0, 5)
#' y <- abs(y)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSexplin(x, cm, rm, tb), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

explinInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["t"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a explin function.")
  }
  
  ## First phase is exponential and the second phase is linear
  cm <- coef(lm(y ~ x, data = xy))[2]
  y2 <- xy[,"y"]/max(xy[,"y"])
  rm <- coef(lm(y2 ~ xy[,"x"]))[2]
  tb <- floor(max(xy[,"x"])/2)

  value <- c(cm, rm, tb)
  names(value) <- mCall[c("cm","rm","tb")]
  value
  
}

#' @rdname SSexplin
#' @return explin: vector of the same length as x using a explin function
#' @export
#' 
explin <- function(t, cm, rm, tb){
  
  .expre1 <- cm / rm
  .expre2 <- rm * (t - tb)
  .expre3 <- log(1 + exp(.expre2))
  .value <- .expre1 * .expre3
  
  ## Derivative with respect to cm
  ## deriv(~ (cm / rm) * log(1 + exp(rm * (t - tb))), "cm")
  .expr6 <- suppressWarnings(log(1 + exp(rm * (t - tb))))
  .exprr6 <- suppressWarnings(1/rm * .expr6)
  .exp1 <- ifelse(is.nan(.expr6),0,.exprr6)
  
  ## Derivative with respect to rm
  ## deriv(~ (cm / rm) * log(1 + exp(rm * (t - tb))), "rm")
  .expr1 <- cm/rm
  .expr2 <- t - tb
  .expr4 <- exp(rm * .expr2)
  .expr5 <- 1 + .expr4
  .expr6 <- suppressWarnings(log(1 + exp(rm * (t - tb))))
  .exprr6 <- suppressWarnings(1/rm * .expr6)
  .exprrr6 <- suppressWarnings(.expr1 * (.expr4 * .expr2/.expr5) - cm/rm^2 * .expr6)
  .exp2 <- ifelse(is.nan(.exprrr6), 0, .exprrr6)
  
  ## Derivative with respect to tb
  ## deriv(~ (cm / rm) * log(1 + exp(rm * (t - tb))), "tb")
  .exp3 <- -(.expr1 * (.expr4 * rm/.expr5))
  
  .actualArgs <- as.list(match.call()[c("cm", "rm", "tb")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("cm", "rm", "tb")))
    .grad[, "cm"] <- .exp1
    .grad[, "rm"] <- .exp2
    .grad[, "tb"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSexplin
#' @export
SSexplin <- selfStart(explin, initial = explinInit, c("cm", "rm", "tb"))

