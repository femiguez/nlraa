#' 
#' The equation is, for a response (y) and a predictor (x): \cr
#'   \eqn{y ~ (x <= xs) * (a + b * x + c * x^2) + (x >= xs) * (a + (-b^2)/(4 * c))} \cr
#'   
#' where the break-point (xs) is -b/c \cr
#' and the asymptote is (a + (-b^2)/(4 * c))
#' 
#' @title self start for quadratic-plateau function
#' @name SSquadp3
#' @rdname SSquadp3
#' @description Self starter for quadratic plateau function with (three) parameters a (intercept), b (slope), c (quadratic)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param c quadratic term
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- quadp3(x, 5, 1.7, -0.04) + rnorm(30, 0, 0.6)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSquadp3(x, a, b, c), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

quadp3Init <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a quadratic-platueau-3.")
  }
  ## Guess for a, b and c is to fit a quadratic linear regression to all the data
  fit <- lm(xy[,"y"] ~ xy[,"x"] + I(xy[,"x"]^2))
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  c <- coef(fit)[3]
  ## If I fix a and b maybe I can try to optimze xs only
  value <- c(a, b, c)
  names(value) <- mCall[c("a","b","c")]
  value
}

#' @rdname SSquadp3
#' @return quadp: vector of the same length as x using the quadratic-plateau function
#' @export
quadp3 <- function(x, a, b, c){
  
  .xs <- -0.5 * b/c
  
  .value <- (x <= .xs) * (a + b * x + c * x^2) + (x >= .xs) * (a + (-b^2)/(4 * c))
  
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

#' @rdname SSquadp3
#' @export
SSquadp3 <- selfStart(quadp3, initial = quadp3Init, c("a","b","c"))

