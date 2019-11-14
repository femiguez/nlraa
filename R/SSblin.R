#' 
#' @title self start for Bilinear Function
#' @name SSblin
#' @rdname SSblin
#' @description Self starter for Bilinear function with parameters a (intercept), b (slope), xs (break-point)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param xs break point of transition between linear and plateau 
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506) 
#' @export
#' @examples 
#' \dontrun{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- blin(x, 0, 1, 20) + rnorm(30, 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSblin(x, a, b, xs), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

blinInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a blinear")
  }

  ## Dumb guess for a and b is to fit a linear regression to all the data
  fit <- lm(xy[,"y"] ~ xy[,"x"])
  ## Atomic bomb approach to kill a mosquito
  objfun <- function(cfs){
    pred <- blin(xy[,"x"], a=cfs[1], b=cfs[2], xs=cfs[3])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(mean(xy[,"x"]),coef(fit))
  op <- try(optim(cfs, objfun, method = "L-BFGS-U",
                  upper = c(Inf, Inf, max(xy[,"x"])),
                  lower = c(-Inf, -Inf, min(xy[,"x"]))), silent = TRUE)
  
  if(class(op) != "try-error"){
    a <- op$par[1]
    b <- op$par[2]
    xs <- op$par[3]
  }else{
    ## If everything fails I use the mean
    a <- coef(fit)[1]
    b <- coef(fit)[2]
    xs <- mean(xy[,"x"])
  }
  
  value <- c(a, b, xs)
  names(value) <- mCall[c("a","b","xs")]
  value
}

#' @rdname SSblin
#' @return vector of the same length as x using the bilinear function
#' @export
blin <- function(x, a, b, xs){
  
  .asym <- a + b * xs
  
  .value <- (x < xs) * (a + b * x) + (x >= xs) * .asym
  
  ## Derivative with respect to a when (x < xs)
  ## .exp1 <- deriv(~ a +  b * x, "a")
  .exp1 <- ifelse(x <= xs, 1, 0)
  
  ## Derivative with respect to b
  ## .exp2 <- deriv(~ a +  b * x, "b")
  .exp2 <- ifelse(x <= xs, x, xs)
  
  ## Derivative with respect to xs
  ## .exp3 <- deriv(~ a +  b * xs, "xs")
  .exp3 <- ifelse(x <= xs, 0, b)
  
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

#' @rdname SSblin
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @export
SSblin <- selfStart(blin, initial = blinInit, c("a","b","xs"))

