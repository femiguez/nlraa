#' 
#' @title self start for plateau-quadratic function
#' @name SSpquad
#' @rdname SSpquad
#' @description Self starter for plateau-quadratic function with parameters a (plateau), xs (break-point), b (slope), c (quadratic)
#' @param x input vector 
#' @param a the plateau value
#' @param xs break-point of transition between plateau and quadratic
#' @param b the slope (linear term)
#' @param c quadratic term
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details Reference for nonlinear regression Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506). 
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(12345)
#' x <- 1:40
#' y <- pquad(x, 5, 20, 1.7, -0.04) + rnorm(40, 0, 0.6)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSpquad(x, a, xs, b, c), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' confint(fit)
#' }
NULL

pquadInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a plateau-quadratic")
  }
  ## Dumb guess for a and b is to fit a quadratic linear regression to 
  ## Second half of the data
  xy1 <- xy[1:floor(nrow(xy)/2),]
  xy2 <- xy[floor(nrow(xy)/2):nrow(xy),]
  xy2$x2 <- xy2[,"x"] - min(xy2[,"x"])
  fit2 <- stats::lm(xy2[,"y"] ~ xy2[,"x2"] + I(xy2[,"x2"]^2))
  a <- coef(fit2)[1]
  b <- coef(fit2)[2]
  c <- coef(fit2)[3]
  ## If I fix a and b maybe I can try to optimze xs only
  objfun <- function(cfs){
    pred <- pquad(xy[,"x"], a=cfs[1], xs=cfs[2], b=cfs[3], c=cfs[4])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  op <- try(stats::optim(c(a, mean(xy[,"x"]),b,c), objfun,
                         method = "L-BFGS-B",
                         upper = c(Inf, max(xy[,"x"]), Inf, Inf),
                         lower = c(-Inf, min(xy[,"x"]), -Inf, -Inf)), silent = TRUE)

  if(!inherits(op, "try-error")){
    a <- op$par[1]
    xs <- op$par[2]
    b <- op$par[3]
    c <- op$par[4]
  }else{
    ## If it fails I use...
    a <- mean(xy1[,"y"])
    xs <- mean(xy[,"x"])
    b <- b
    c <- c
  }
  value <- c(a, xs, b, c)
  names(value) <- mCall[c("a","xs","b","c")]
  value
}

#' @rdname SSpquad
#' @return pquad: vector of the same length as x using the plateau-quadratic function
#' @export
pquad <- function(x, a, xs, b, c){
  
  .value <- (x < xs) * a + (x >= xs) * (a + b * (x - xs) + c * (x - xs)^2)
  
  ## Derivative with respect to a
  .exp1 <- 1 # ifelse(x < xs, 1, 1) 
  
  ## Derivative with respect to xs
  ## .exp2 <- deriv(~ a  + b * (x - xs) + c * (x - xs)^2, "xs")
  .exp2 <- ifelse(x < xs, 0, -b + -(c * (2 * (x - xs))))
  
  ## Derivative with respect to b
  ## .exp3 <- deriv(~ a  + b * (x - xs) + c * (x - xs)^2, "b")
  .exp3 <- ifelse(x < xs, 0, x - xs)
  
  ## Derivative with respect to c
  ## .exp4 <- deriv(~ a  + b * (x - xs) + c * (x - xs)^2, "c")
  .exp4 <- ifelse(x < xs, 0, (x - xs)^2)
  
  .actualArgs <- as.list(match.call()[c("a","xs","b","c")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a","xs","b","c")))
    .grad[, "a"] <- .exp1
    .grad[, "xs"] <- .exp2
    .grad[, "b"] <- .exp3
    .grad[, "c"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSpquad
#' @export
SSpquad <- selfStart(pquad, initial = pquadInit, c("a","xs","b","c"))

