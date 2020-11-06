#' 
#' @title self start for quadratic-plateau function
#' @name SSquadp
#' @rdname SSquadp
#' @description Self starter for quadratic plateau function with parameters a (intercept), b (slope), c (quadratic), xs (break-point)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param c quadratic term
#' @param xs break point of transition between quadratic and plateau 
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details Reference for nonlinear regression Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506). 
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- quadp(x, 5, 1.7, -0.04, 20) + rnorm(30, 0, 0.6)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSquadp(x, a, b, c, xs), data = dat, algorithm = "port")
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

quadpInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a quadratic-platueau.")
  }
  ## Dumb guess for a and b is to fit a quadratic linear regression to all the data
  fit <- lm(xy[,"y"] ~ xy[,"x"] + I(xy[,"x"]^2))
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  c <- coef(fit)[3]
  ## If I fix a and b maybe I can try to optimze xs only
  objfun <- function(xs, a, b, c){
    pred <- quadp(xy[,"x"], a, b, c, xs)
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  op.xs <- try(stats::optimize(objfun, range(xy[,"x"]), 
                               a = a, b = b,
                               c = c), silent = TRUE)

  if(class(op.xs) != "try-error"){
    xs <- op.xs$minimum
  }else{
    ## If it fails I use the mean
    xs <- mean(xy[,"x"])
  }
  value <- c(a, b, c, xs)
  names(value) <- mCall[c("a","b","c","xs")]
  value
}

#' @rdname SSquadp
#' @return quadp: vector of the same length as x using the quadratic-plateau function
#' @export
quadp <- function(x, a, b, c, xs){
  
  .asym <- a + b * xs + c * xs^2
  
  .value <- (x < xs) * (a + b * x + c * x^2) + (x >= xs) * .asym
  
  ## Derivative with respect to a
  .exp1 <- 1 # ifelse(x < xs, 1, 1) 
  
  ## Derivative with respect to b
  ## .exp2 <- deriv(~ a  + b * x + c * x^2, "b")
  .exp2 <- ifelse(x < xs, x, xs)
  
  ## Derivative with respect to c
  ## .exp3 <- deriv(~ a  + b * x + c * x^2, "c")
  .exp3 <- ifelse(x < xs, x^2, xs^2)
  
  ## Derivative with respect to xs
  ## .exp4 <- deriv(~ a  + b * xs + c * xs^2, "xs")
  .exp4 <- ifelse(x < xs, 0, b + 2 * c * xs)
  
  .actualArgs <- as.list(match.call()[c("a","b","c","xs")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a","b","c","xs")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    .grad[, "c"] <- .exp3
    .grad[, "xs"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
   }
  .value
}

#' @rdname SSquadp
#' @export
SSquadp <- selfStart(quadp, initial = quadpInit, c("a","b","c","xs"))

