#' 
#' @title self start for linear-plateau function
#' @name SSlinp
#' @rdname SSlinp
#' @description Self starter for linear-plateau function with parameters a (intercept), b (slope), xs (break-point)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param xs break-point of transition between linear and plateau 
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is linear when \eqn{x < xs: (a + b * x)} and flat (\eqn{asymptote = a + b * xs}) when \eqn{x >= xs}.
#' @seealso package \CRANpkg{segmented}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- linp(x, 0, 1, 20) + rnorm(30, 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSlinp(x, a, b, xs), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Confidence intervals
#' confint(fit)
#' }
#' 
NULL

linpInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a linear-plateau")
  }

  ## Dumb guess for a and b is to fit a linear regression to half the data
  xy1 <- xy[1:floor(nrow(xy)/2),]
  fit1 <- stats::lm(xy1[,"y"] ~ xy1[,"x"])
  ## Atomic bomb approach to kill a mosquito
  objfun <- function(cfs){
    pred <- linp(xy[,"x"], a=cfs[1], b=cfs[2], xs=cfs[3])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(coef(fit1),mean(xy[,"x"]))
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-B",
                         upper = c(Inf, Inf, max(xy[,"x"])),
                         lower = c(-Inf, -Inf, min(xy[,"x"]))), silent = TRUE)

  if(class(op) != "try-error"){
    a <- op$par[1]
    b <- op$par[2]
    xs <- op$par[3]
  }else{
    ## If it fails I use the mean for the breakpoint
    a <- coef(fit1)[1]
    b <- coef(fit1)[2]
    xs <- mean(xy[,"x"])
  }
  
  value <- c(a, b, xs)
  names(value) <- mCall[c("a","b","xs")]
  value
}

#' @rdname SSlinp
#' @return linp: vector of the same length as x using the linear-plateau function
#' @export
linp <- function(x, a, b, xs){
  
  .asym <- a + b * xs
  
  .value <- (x < xs) * (a + b * x) + (x >= xs) * .asym
  
  ## Derivative with respect to a when (x < xs)
  ## .exp1 <- deriv(~ a +  b * x + c * x^2, "a")
  .exp1 <- 1 ## ifelse(x < xs, 1, 1)
  
  ## Derivative with respect to b
  ## if x < xs: .exp2 <- deriv(~ a +  b * x + c * x^2, "b")
  ## if x >= xs: .exp2 <- deriv(~ a +  b * xs + c * x^2, "b")
  .exp2 <- ifelse(x < xs, x, xs)
  
  ## Derivative with respect to xs
  ## .exp3 <- deriv(~ a +  b * xs, "xs")
  .exp3 <- ifelse(x < xs, 0, b)
  
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

#' @rdname SSlinp
#' @export
SSlinp <- selfStart(linp, initial = linpInit, c("a","b","xs"))

