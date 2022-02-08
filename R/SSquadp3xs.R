#' 
#' The equation is, for a response (y) and a predictor (x): \cr
#'   \eqn{y ~ (x <= xs) * (a + b * x + (-0.5 * b/xs) * x^2) + (x > xs) * (a + (b^2)/(-2 * b/xs))} \cr
#'   
#' where the quadratic term is (c) is -0.5*b/xs \cr
#' and the asymptote is (a + (b^2)/(4 * c)).
#' 
#' This model does not estimate the quadratic parameter \sQuote{c} directly. 
#' If this is required, the model \sQuote{SSquadp3} should be used instead.
#' 
#' @title self start for quadratic-plateau function (xs)
#' @name SSquadp3xs
#' @rdname SSquadp3xs
#' @description Self starter for quadratic plateau function with (three) parameters a (intercept), b (slope), xs (break-point)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param xs break-point
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- quadp3xs(x, 5, 1.7, 20) + rnorm(30, 0, 0.6)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSquadp3xs(x, a, b, xs), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

quadp3xsInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a quadratic-platueau-3-xs.")
  }
  ## Guess for a, b and xs is to fit a quadratic linear regression to all the data
  fit <- lm(xy[,"y"] ~ xy[,"x"] + I(xy[,"x"]^2))
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  c <- coef(fit)[3]
  xs <- -0.5 * b/c
  ## If I fix a and b maybe I can try to optimze xs only
  value <- c(a, b, xs)
  names(value) <- mCall[c("a","b","xs")]
  value
}

#' @rdname SSquadp3xs
#' @return quadp3xs: vector of the same length as x using the quadratic-plateau function
#' @export
quadp3xs <- function(x, a, b, xs){
  
  .value <- (x <= xs) * (a + b * x + (-0.5 * b/xs) * x^2) + (x > xs) * (a + (-b^2)/(4 * -0.5 * b/xs))
  
  ## Derivative with respect to a
  .exp1 <- 1 
  
  ## Derivative with respect to b
  ## .exp2 <- deriv(~ a + b * x + (-0.5 * b/xs) * x^2, "b")
  ## .exp2b <- deriv(~ a + (-b^2)/(4 * -0.5 * b/xs), "b")
  .expr2 <- -b^2
  .expr4 <- 4 * -0.5
  .expr6 <- .expr4 * b/xs
  .exp2 <- ifelse(x <= xs, x - 0.5/xs * x^2, -(2 * b/.expr6 + .expr2 * (.expr4/xs)/.expr6^2))
  
  ## Derivative with respect to xs
  ## .exp2 <- deriv(~ (a + b * x + (-0.5 * b/xs) * x^2), "xs")
  ## .exp2b <- deriv(~ a + (-b^2)/(4 * -0.5 * b/xs), "xs")
  .expr4 <- -0.5 * b
  .expr6 <- x^2
  .expr2 <- -b^2
  .expr5 <- 4 * -0.5 * b
  .expr7 <- .expr5/xs
  .exp3 <- ifelse(x <= xs, -(.expr4/xs^2 * .expr6), .expr2 * (.expr5/xs^2)/.expr7^2)
  
  .actualArgs <- as.list(match.call()[c("a","b","xs")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a","b","xs")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    .grad[, "xs"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSquadp3xs
#' @export
SSquadp3xs <- selfStart(quadp3xs, initial = quadp3xsInit, c("a","b","xs"))

