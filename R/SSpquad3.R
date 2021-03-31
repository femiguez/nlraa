#' 
#' @title self start for plateau-quadratic function
#' @name SSpquad3
#' @rdname SSpquad3
#' @description Self starter for plateau-quadratic function with (three) parameters a (intercept), b (slope), c (quadratic)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param c quadratic term
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details Reference for nonlinear regression Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506). 
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' require(minpack.lm)
#' set.seed(123)
#' x <- 1:30
#' y <- pquad3(x, 20.5, 0.36, -0.012) + rnorm(30, 0, 0.3)
#' dat <- data.frame(x = x, y = y)
#' fit <- nlsLM(y ~ SSpquad3(x, a, b, c), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

pquad3Init <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a plateau-quadratic-3.")
  }
  ## Guess for a, b and c is to fit a quadratic linear regression to all the data
  half.xy <- xy[floor(nrow(xy)/2):nrow(xy),]
  fit <- lm(half.xy[,"y"] ~ half.xy[,"x"] + I(half.xy[,"x"]^2))
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  c <- coef(fit)[3]
  
  value <- c(a, b, c)
  names(value) <- mCall[c("a","b","c")]
  value
}

#' @rdname SSpquad3
#' @return quadp: vector of the same length as x using the quadratic-plateau function
#' @export
pquad3 <- function(x, a, b, c){
  
  .xs <- -0.5 * b/c
  
  .value <- (x >= .xs) * (a + b * x + c * x^2) + (x < .xs) * (a + (-b^2)/(4 * c))
  
  ## Derivative with respect to a
  .exp1 <- 1 
  
  ## Derivative with respect to b
  ## .exp2 <- deriv(~ a  + b * x + c * x^2, "b")
  .exp2 <- ifelse(x < .xs, x, .xs)
  
  ## Derivative with respect to c
  ## .exp3 <- deriv(~ a  + b * x + c * x^2, "c")
  .exp3 <- ifelse(x < .xs, x^2, .xs^2)
  
  .actualArgs <- as.list(match.call()[c("a","b","c")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a","b","c")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    .grad[, "c"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSpquad3
#' @export
SSpquad3 <- selfStart(pquad3, initial = pquad3Init, c("a","b","c"))

