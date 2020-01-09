#' 
#' The equation is: \deqn{ a * x ^ c / (1 + b * x ^ d)}
#' 
#' @title self start for a rational curve
#' @name SSratio
#' @rdname SSratio
#' @description Self starter for a rational curve
#' @param x input vector 
#' @param a parameter related to the maximum value of the response (numerator)
#' @param b power exponent for numerator
#' @param c parameter related to the maximum value of the response (denominator)
#' @param d power exponent for denominator
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506).  
#' One example application is in Bril et al. (1994) \url{https://edepot.wur.nl/333930} - pages 19 and 21. 
#' The parameters are difficult to interpret, but the function is very flexible. I have not tested this, 
#' but it might be beneficial to re-scale x and y to the (0,1) range if this function is hard to fit.
#' \url{https://en.wikipedia.org/wiki/Rational_function}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' require(minpack.lm)
#' set.seed(1234)
#' x <- 1:100
#' y <- ratio(x, 1, 0.5, 1, 1.5) + rnorm(length(x), 0, 0.025)
#' dat <- data.frame(x = x, y = y)
#' fit <- nlsLM(y ~ SSratio(x, a, b, c, d), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

ratioInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a rational function.")
  }
  
  objfun <- function(cfs){
    pred <- ratio(xy[,"x"], a=cfs[1], b=cfs[2], c=cfs[3], d=cfs[4])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  
  cfs <- c(1,1,1,1)
  op <- try(stats::optim(cfs, objfun), silent = TRUE)
  
  if(class(op) != "try-error"){
    a <- op$par[1]
    b <- op$par[2]
    c <- op$par[3]
    d <- op$par[4]
  }else{
    op <- try(stats::optim(cfs, objfun, method = "SANN"), silent = TRUE)
    if(class(op) != "try-error"){
      a <- op$par[1]
      b <- op$par[2]
      c <- op$par[3]
      d <- op$par[4]
      warning('Used method = "SANN" in optim')
    }else{
      a <- 1
      b <- 1
      c <- 1
      d <- 1
      warning("Could not find suitable starting values")
    }
  }
  value <- c(a, b, c, d)
  names(value) <- mCall[c("a","b","c","d")]
  value
  
}

#' @rdname SSratio
#' @return ratio: vector of the same length as x using a rational function
#' @export
#' 
ratio <- function(x, a, b, c, d){
  
  .expre1 <- a * x^c
  .expre2 <- 1 +  b * x^d
  .value <- .expre1/.expre2
  
  ## Derivative with respect to a
  ## deriv(~ (a * x^c)/(1 + b * x^d), "a")
  .expr1 <- x^c
  .expr5 <- 1 + b * x^d
  .exp1 <- .expr1/.expr5
  .exp1 <- ifelse(is.nan(.exp1),0,.exp1)
  
  ## Derivative with respect to b
  ## deriv(~ (a * x^c)/(1 + b * x^d), "b")
  .expr2 <- a * x^c
  .expr3 <- x^d
  .exp2 <- -(.expr2 * .expr3/.expr5^2)
  .exp2 <- ifelse(is.nan(.expr2), 0, .expr2)
  
  ## Derivative with respect to c
  ## deriv(~ (a * x^c)/(1 + b * x^d), "c")
  .lx <- suppressWarnings(log(x))
  .exp3 <- a * (.expr1 * .lx)/.expr5
  .exp3 <- ifelse(is.nan(.expr3), 0, .expr3)
  
  ## Derivative with respect to d
  .exp4 <- -(.expr2 * (b * (.expr3 * .lx))/.expr5^2)
  .exp4 <- ifelse(is.nan(.exp4), 0, .exp4)
  
  .actualArgs <- as.list(match.call()[c("a", "b", "c", "d")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
     .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "b", "c","d")))
     .grad[, "a"] <- .exp1
     .grad[, "b"] <- .exp2
     .grad[, "c"] <- .exp3
     .grad[, "d"] <- .exp4
     dimnames(.grad) <- list(NULL, .actualArgs)
     attr(.value, "gradient") <- .grad
   }
  .value
}

#' @rdname SSratio
#' @export
SSratio <- selfStart(ratio, initial = ratioInit, c("a", "b", "c", "d"))

