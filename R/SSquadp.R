#' 
#' @title self start for Quadratic Plateua Function
#' @name SSquadp
#' @rdname SSquadp
#' @description Self starter for Quadratic Plateau function with parameters a (intercept), b (slope), c (quadratic), xs (break-point)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param c quadratic term
#' @param xs break point of transition between quadratic and plateau 
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506) 
#' @export
#' @examples 
#' \dontrun{
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

quadpInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 5){
    stop("Too few distinct input values to fit a quadratic-platueau")
  }
  
  ## Dumb guess for a and b is to fit a linear regression to all the data
  fit <- lm(xy[,"y"] ~ xy[,"x"] + I(xy[,"x"]^2))
  ## If I fix a and b maybe I can try to optimze xs only
  objfun <- function(xs, a, b, c){
    pred <- quadp(xy[,"x"], a, b, c, xs)
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  op.xs <- try(optimize(objfun, range(xy[,"x"]), 
                        a = coef(fit)[1], b = coef(fit)[2],
                        c = coef(fit)[3]), silent = TRUE)
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  c <- coef(fit)[3]
  if(class(op.xs) != "try-error"){
    xs <- op.xs$minimum
  }else{
    ## If everything fails I use the mean
    xs <- mean(xy[,"x"])
  }
  ## Guess xs
  value <- c(a, b, c, xs)
  names(value) <- mCall[c("a","b","c","xs")]
  value
}

#' @rdname SSquadp
#' @return vector of the same length as x using the quadratic-plateau function
#' @export
quadp <- function(x, a, b, c, xs){
  
  .asym <- a + b * xs + c * xs^2
  
  .value <- (x < xs) * (a + b * x + c * x^2) + (x >= xs) * .asym
  
  ## Derivative with respect to a
  ## .exp1 <- deriv(~ a * time * exp(-b * time), "a")
  .exp1 <- (x <= xs) 
  
  ## Derivative with respect to b
  ## .exp2 <- deriv(~ a * time * exp(-b * time), "b")
  .exp2 <- (x <= xs) * b
  
  ## Derivative with respect to xs
  .exp3 <- (x >= xs) * .asym
  
  .actualArgs <- as.list(match.call()[c("a","b","c","xs")])
  
  ##  Gradient
  ## if (all(unlist(lapply(.actualArgs, is.name)))) {
  ##  .grad <- array(0, c(length(.value), 2L), list(NULL, c("a", "b")))
  ##  .grad[, "a"] <- .exp1
  ##  .grad[, "b"] <- .exp2
  ##  dimnames(.grad) <- list(NULL, .actualArgs)
  ##  attr(.value, "gradient") <- .grad
  ## }
  .value
}

#' @rdname SSquadp
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @export
SSquadp <- selfStart(quadp, initial = quadpInit, c("a","b","c","xs"))

