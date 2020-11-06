#' 
#' @title self start for a bell-shaped curve
#' @name SSbell
#' @rdname SSbell
#' @description Self starter for a type of bell-shaped curve
#' @param x input vector 
#' @param ymax maximum value of y
#' @param a parameter controlling the spread (associated with a quadratic term)
#' @param b parameter controlling the spread (associated with a cubic term)
#' @param xc centering parameter
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506). One example application is Hammer et al. (2009) (doi:10.2135/cropsci2008.03.0152).
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:20
#' y <- bell(x, 8, -0.0314, 0.000317, 13) + rnorm(length(x), 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSbell(x, ymax, a, b, xc), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

bellInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a bell-shaped curve.")
  }
  ymax <- max(xy[,"y"])
  xc <- NLSstClosestX(xy, ymax)
  wygz <- which(xy[,"y"] > 0)
  xy1 <- xy[wygz,]
  sy <- log(xy1[,"y"]/ymax)
  sx <- xy1[,"x"] - xc
  lm.cfs <- coef(stats::lm(sy ~ I(sx^2) + I(sx^3) - 1))
  a <- lm.cfs[1]
  b <- lm.cfs[2]
  value <- c(ymax, a, b, xc)
  names(value) <- mCall[c("ymax","a","b","xc")]
  value
  
}

#' @rdname SSbell
#' @return bell: vector of the same length as x using a bell-shaped curve
#' @export
#' 
bell <- function(x, ymax, a, b, xc){
  
  .expre1 <- a * (x - xc)^2
  .expre2 <- b * (x - xc)^3
  .value <- ymax * exp(.expre1 + .expre2)
    
  ## Derivative with respect to ymax
  ## deriv(~ ymax * exp(a * (x - xc)^2 + b * (x - xc)^3), "ymax")
  .expr1 <- x - xc
  .expr7 <- exp(a * .expr1^2 + b * .expr1^3)
  .exp1 <- .expr7
  
  ## Derivative with respect to a
  ## deriv(~ ymax * exp(a * (x - xc)^2 + b * (x - xc)^3), "a")
  .expr2 <- .expr1^2
  .exp2 <- ymax * (.expr7 * .expr2)
  
  ## Derivative with respect to b
  ## deriv(~ ymax * exp(a * (x - xc)^2 + b * (x - xc)^3), "b")
  .expr4 <- .expr1^3
  .exp3 <- ymax * (.expr7 * .expr4)
  
  ## Derivative with respect to xc
  ## deriv(~ ymax * exp(a * (x - xc)^2 + b * (x - xc)^3), "xc")
  .exp4 <- -(ymax * (.expr7 * (b * (3 * .expr2) + a * (2 * .expr1))))
  
  .actualArgs <- as.list(match.call()[c("ymax", "a", "b", "xc")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("ymax", "a", "b","xc")))
    .grad[, "ymax"] <- .exp1
    .grad[, "a"] <- .exp2
    .grad[, "b"] <- .exp3
    .grad[, "xc"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSbell
#' @export
SSbell <- selfStart(bell, initial = bellInit, c("ymax", "a", "b", "xc"))
