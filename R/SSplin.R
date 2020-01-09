#' 
#' @title self start for plateau-linear function
#' @name SSplin
#' @rdname SSplin
#' @description Self starter for plateau-linear function with parameters a (plateau), xs (break-point), b (slope) 
#' @param x input vector 
#' @param a the initial plateau
#' @param xs break-point of transition between plateau and linear 
#' @param b the slope
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details Initial plateau with a second linear phase. When \eqn{x < xs: y = a} and when \eqn{x >= xs: y = a + b * (x - xs)}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- plin(x, 10, 20, 1) + rnorm(30, 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSplin(x, a, xs, b), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Confidence intervals
#' confint(fit)
#' }
#' 
NULL

plinInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a plateau-linear")
  }
  ## Split the data in half as an initial guess
  xy1 <- xy[1:floor(nrow(xy)/2),]
  xy2 <- xy[floor(nrow(xy)/2):nrow(xy),]
  fit2 <- stats::lm(xy2[,"y"] ~ xy2[,"x"])
  ## Atomic bomb approach to kill a mosquito
  objfun <- function(cfs){
    pred <- plin(xy[,"x"], a=cfs[1], xs=cfs[2], b=cfs[3])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(mean(xy1[,"y"]), mean(xy[,"x"]),coef(fit2)[2])
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-U",
                         upper = c(Inf, max(xy[,"x"]), Inf),
                         lower = c(-Inf, min(xy[,"x"])), -Inf), silent = TRUE)
  
  if(class(op) != "try-error"){
    a <- op$par[1]
    xs <- op$par[2]
    b <- op$par[3]
  }else{
    ## If it fails I use the mean
    a <- mean(xy1[,"y"])
    b <- coef(fit2)[2]
    xs <- mean(xy[,"x"])
  }
  
  value <- c(a, xs, b)
  names(value) <- mCall[c("a","xs","b")]
  value
}

#' @rdname SSplin
#' @return plin: vector of the same length as x using the plateau-linear function
#' @export
plin <- function(x, a, xs, b){
  
  .value <- (x < xs) * a + (x >= xs) * (a + b * (x - xs))
  
  ## Derivative with respect to a when (x < xs)
  .exp1 <- 1 ## ifelse(x < xs, 1, 1)
  
  ## Derivative with respect to xs
  ## .exp2 <- deriv(~ a + b * (x - xs),"xs")
  .exp2 <- ifelse(x < xs, 0, -1 * b)
  
  ## Derivative with respect to b
  ## .exp3 <- deriv(~ a + b * (x - xs),"b")
  .exp3 <- ifelse(x < xs, 0, x - xs)
  
  .actualArgs <- as.list(match.call()[c("a","xs","b")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a","xs","b")))
    .grad[, "a"] <- .exp1
    .grad[, "xs"] <- .exp2
    .grad[, "b"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSplin
#' @export
SSplin <- selfStart(plin, initial = plinInit, c("a","xs","b"))

