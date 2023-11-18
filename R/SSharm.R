#' 
#' @title self start for a harmonic regression model
#' @name SSharm1
#' @rdname SSharm
#' @description Self starter for a harmonic regression
#' @param x input vector 
#' @param b0 intercept of the harmonic regression
#' @param b1 slope of the harmonic regression 
#' @param cos1 coefficient associated with the cosine of the harmonic regression
#' @param sin1 coefficient associated with the sine of the harmonic regression
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details Harmonic regression is actually a type of linear regression. Just adding it for convenience.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- seq(0, 3, length.out = 100)
#' y <- harm1(x, 0, 0, 0.05, 0) + rnorm(length(x), 0, 0.002)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSharm1(x, b0, b1, cos1, sin1), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

harm1Init <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a harmonic curve.")
  }
  
  cfs <- coef(lm(xy[,"y"] ~ xy[,"x"] + I(cos(2 * pi * xy[,"x"])) + 
                   I(sin(2 * pi * xy[,"x"]))))
  
  b0 <- cfs[1]
  b1 <- cfs[2]
  cos1 <- cfs[3]
  sin1 <- cfs[4]
  
  value <- c(b0, b1, cos1, sin1)
  names(value) <- mCall[c("b0", "b1", "cos1", "sin1")]
  value
  
}

#' @rdname SSharm
#' @return harm1: vector of the same length as x using a harmonic regression
#' @export
#' 
harm1 <- function(x, b0, b1, cos1, sin1){
  
  .expre1 <- b0 + (b1 * x) 
  .expre2 <- cos1 * cos( 2 * pi * x) 
  .expre3 <- sin1 * sin( 2 * pi * x)

  .value <- .expre1 + .expre2 + .expre3
  
  ## Derivative with respect to b0
  .exp1 <- 1
  
  ## Derivative with respect to b1
  .exp2 <- x
  
  ## Derivative with respect to cos1
  ## deriv(~ b0 + (b1 * x) + (cos1 * cos( 2 * pi * x)) + (sin1 * sin( 2 * pi * x)) , "cos1")
  .exp3 <- cos(2 * pi * x)
  
  ## Derivative with respect to sin1
  ## deriv(~ b0 + (b1 * x) + (cos1 * cos( 2 * pi * x)) + (sin1 * sin( 2 * pi * x)) , "sin1")
  .exp4 <- sin(2 * pi * x)
  
  .actualArgs <- as.list(match.call()[c("b0", "b1", "cos1", "sin1")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("b0", "b1", "cos1", "sin1")))
    .grad[, "b0"] <- .exp1
    .grad[, "b1"] <- .exp2
    .grad[, "cos1"] <- .exp3
    .grad[, "sin1"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSharm
#' @export
SSharm1 <- selfStart(harm1, initial = harm1Init, c("b0", "b1", "cos1", "sin1"))
