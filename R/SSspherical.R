#' 
#'  This equation was found in a publication by Dobermann et al. (2011) \doi{doi:10.2134/agronj2010.0179}
#' 
#' @title self start for spherical function
#' @name SSspherical
#' @rdname SSspherical
#' @description Self starter for a spherical function with parameters a (intercept), b (see below), xs (break-point)
#' @param x input vector 
#' @param a the intercept
#' @param b the difference between the intercept and the asymptote, so that a + b = asymptote
#' @param xs break-point of transition between nonlinear and plateau 
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is nonlinear when \eqn{x < xs} and flat (\eqn{asymptote = a + b}) when \eqn{x >= xs}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- seq(0, 400, length.out = 50)
#' y <- spherical(x, 2, 5, 200) + rnorm(length(x), sd = 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSspherical(x, a, b, xs), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Confidence intervals
#' confint(fit)
#' }
#' 
NULL

sphericalInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a spherical")
  }
  
  ## Guessing a from an intercept model
  fit1 <- stats::lm(xy[,"y"] ~ xy[,"x"] + I(xy[,"x"]^2))
  ## Atomic bomb approach to kill a mosquito
  objfun <- function(cfs){
    pred <- spherical(xy[,"x"], a=cfs[1], b=cfs[2], xs=cfs[3])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(coef(fit1)[1], max(xy[,"y"]) - coef(fit1)[1], mean(xy[,"x"]))
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-B",
                         upper = c(Inf, Inf, max(xy[,"x"])),
                         lower = c(-Inf, -Inf, min(xy[,"x"]))), silent = TRUE)
  
  if(!inherits(op, "try-error")){
    a <- op$par[1]
    b <- op$par[2]
    xs <- op$par[3]
  }else{
    ## If it fails I use the mean for the breakpoint
    a <- coef(fit1)[1]
    b <- max(xy[,"y"]) - coef(fit1)[1]
    xs <- mean(xy[,"x"])
  }
  
  value <- c(a, b, xs)
  names(value) <- mCall[c("a", "b", "xs")]
  value
}

#' @rdname SSspherical
#' @return spherical: vector of the same length as x using the spherical function
#' @export
spherical <- function(x, a, b, xs){
  
  .sp1 <- a + b * ((3 * x)/(2 * xs) - 0.5 * (x/xs)^3)
  
  .value <- ((xs - x) < 0) * (a + b) + ((xs - x) >= 0) * .sp1
  
  ## Derivative with respect to a 
  .exp1 <- 1
  
  ## Derivative with respect to b
  .exp2 <- ifelse((xs - x) < 0, 1, 3 * x/(2 * xs) - 0.5 * (x/xs)^3) 
  
  ## Derivative with respect to xs
  ## deriv(~a + b * ((3 * x)/(2 * xs) - 0.5 * (x/xs)^3), "xs")
  .expr1 <- 3 * x
  .expr2 <- 2 * xs
  .expr4 <- x/xs
  .expr3 <- -(b * (.expr1 * 2/.expr2^2 - 0.5 * (3 * (x/xs^2 * .expr4^2))))
  .exp3 <- ifelse((xs - x) < 0, 0, .expr3)
  
  .actualArgs <- as.list(match.call()[c("a", "b", "xs")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a", "b", "xs")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    .grad[, "xs"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSspherical
#' @export
SSspherical <- selfStart(spherical, initial = sphericalInit, c("a", "b", "xs"))


set.seed(123)
x <- seq(0, 400, length.out = 50)
y <- spherical(x, 2, 5, 200) + rnorm(length(x), sd = 0.5)
dat <- data.frame(x = x, y = y)
