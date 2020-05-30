#' For details see https://en.wikipedia.org/wiki/Hill_equation_(biochemistry)
#' 
#' @title self start for Hill Function
#' @name SShill
#' @rdname SShill
#' @description Self starter for Hill function with parameters Ka, n and a
#' @param x input vector (x). Concentration of substrate in the original Hill model.
#' @param Ka parameter representing the concentration at which half of maximum y is attained
#' @param n parameter which controls the curvature 
#' @param a parameter which controls the maximum value of the response (asymptote)
#' @details The form of the equations are: \cr
#' hill1: \deqn{1 / (1 + (Ka/x))}. \cr 
#' hill2: \deqn{1 / (1 + (Ka/x)^n)}. \cr
#' hill3: \deqn{a / (1 + (Ka/x)^n)}. \cr
#' @note Zero values are not allowed.
#' @examples 
#' \donttest{
#' require(ggplot2)
#' ## Example for hill1
#' set.seed(1234)
#' x <- 1:20
#' y <- hill1(x, 10) + rnorm(20, sd = 0.03)
#' dat1 <- data.frame(x = x, y = y)
#' fit1 <- nls(y ~ SShill1(x, Ka), data = dat1)
#' 
#' ## Example for hill2
#' y <- hill2(x, 10, 1.5) + rnorm(20, sd = 0.03)
#' dat2 <- data.frame(x = x, y = y)
#' fit2 <- nls(y ~ SShill2(x, Ka, n), data = dat2)
#' 
#' ## Example for hill3
#' y <- hill3(x, 10, 1.5, 5) + rnorm(20, sd = 0.03)
#' dat3 <- data.frame(x = x, y = y)
#' fit3 <- nls(y ~ SShill3(x, Ka, n, a), data = dat3)
#' 
#' ggplot(data = dat3, aes(x, y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit3)))
#' }
NULL

## define hill1
hill1Init <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a hill1")
  }
  
  xy2 <- xy[,"x" > 0 & "y" > 0]
  y1 <- log(xy2[,"y"]/(1 - xy2[,"y"])) ## When y == 1 it causes Inf 
  y2 <- ifelse(is.finite(y1), y1, NA)
  cfs <- try(coef(lm(y2 ~ log(xy2[,"x"]), na.action = "na.omit")), silent = TRUE)
  
  if(inherits(cfs, "try-error")){
    Ka <- mean(xy2[,"x"], na.rm = TRUE)
  }else{
    Ka <- exp(-cfs[1])
  } 

  value <- c(Ka)
  names(value) <- mCall[c("Ka")]
  value
  
}

#' @rdname SShill
#' @return hill1: vector of the same length as x (time) using the Hill 1 function
#' @export
hill1 <- function(x, Ka){
  
  if(any(identical(x, 0))) stop("zero x is not allowed")
  
  .value <- 1 / (1 + (Ka/x))
  
  ## Derivative with respect to Ka
  ## deriv(~1 / (1 + (Ka/x)),"Ka")
  .expr2 <- 1 + Ka/x
  .expi1 <- -(1/x/.expr2^2)

  .actualArgs <- as.list(match.call()[c("Ka")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("Ka")))
    .grad[, "Ka"] <- .expi1
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SShill
#' @export
SShill1 <- selfStart(hill1, initial = hill1Init, c("Ka"))

## define hill2
hill2Init <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a hill2")
  }
  
  xy2 <- xy[,"x" > 0 & "y" > 0]
  y1 <- log(xy2[,"y"]/(1 - xy2[,"y"]))
  y2 <- ifelse(is.finite(y1), y1, NA)
  cfs <- try(coef(lm(y2 ~ log(xy2[,"x"]), na.action = "na.omit")), silent = TRUE)
  
  if(inherits(cfs, "try-error")){
    n <- 1
    Ka <- mean(xy2[,"x"], na.rm = TRUE)
  }else{
    n <- cfs[2]
    Ka <- exp(-cfs[1]/n)
  } 
  
  value <- c(Ka, n)
  names(value) <- mCall[c("Ka","n")]
  value
}

#' @rdname SShill
#' @return hill2: vector of the same length as x (time) using the Hill 2 function
#' @export
hill2 <- function(x, Ka, n){
  
  if(any(identical(x, 0))) stop("zero x is not allowed")
  
  .value <- 1 / (1 + (Ka/x)^n)
  
  ## Derivative with respect to Ka
  ## deriv(~1 / (1 + (Ka/x)^n),"Ka")
  .expr1 <- Ka/x
  .expr3 <- 1 + .expr1^n
  .expi1 <- -(.expr1^(n - 1) * (n * (1/x))/.expr3^2)
  
  ## Derivative with respect to n
  ## deriv(~1 / (1 + (Ka/x)^n),"n")
  .expr2 <- .expr1^n
  .expi2 <- -(.expr2 * log(.expr1)/.expr3^2)
  
  .actualArgs <- as.list(match.call()[c("Ka","n")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("Ka","n")))
    .grad[, "Ka"] <- .expi1
    .grad[, "n"] <- .expi2
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SShill
#' @export
SShill2 <- selfStart(hill2, initial = hill2Init, c("Ka","n"))

## define hill3
hill3Init <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a hill3")
  }
  
  xy2 <- xy[,"x" > 0]
  y0 <- xy2[,"y"]/max(xy2[,"y"])
  y1 <- log(y0/(1 - y0))
  y2 <- ifelse(is.finite(y1), y1, NA)
  cfs <- try(coef(lm(y2 ~ log(xy2[,"x"]), na.action = "na.omit")), silent = TRUE)
  
  if(inherits(cfs, "try-error")){
    a <- max(xy[,"y"])
    n <- 1
    Ka <- mean(xy2[,"x"], na.rm = TRUE)
  }else{
    a <- max(xy[,"y"])
    n <- cfs[2]
    Ka <- exp(-cfs[1]/n)
  } 
  
  value <- c(Ka, n, a)
  names(value) <- mCall[c("Ka","n","a")]
  value
}

#' @rdname SShill
#' @return hill3: vector of the same length as x (time) using the Hill 3 function
#' @export
hill3 <- function(x, Ka, n, a){
  
  if(any(identical(x, 0))) stop("zero x is not allowed")
  
  .value <- a / (1 + (Ka/x)^n)
  
  ## Derivative with respect to Ka
  ## deriv(~a / (1 + (Ka/x)^n),"Ka")
  .expr1 <- Ka/x
  .expr3 <- 1 + .expr1^n
  .expi1 <- -(a * (.expr1^(n - 1) * (n * (1/x)))/.expr3^2)
  
  ## Derivative with respect to n
  ## deriv(~1 / (1 + (Ka/x)^n),"n")
  .expr2 <- .expr1^n
  .expi2 <- -(a * (.expr2 * log(.expr1))/.expr3^2)
  
  ## Derivative with respect to a
  ## deriv(~1 / (1 + (Ka/x)^n),"a")
  .expi3 <- 1/.expr3
  
  .actualArgs <- as.list(match.call()[c("Ka","n","a")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("Ka","n","a")))
    .grad[, "Ka"] <- .expi1
    .grad[, "n"] <- .expi2
    .grad[, "a"] <- .expi3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SShill
#' @export
SShill3 <- selfStart(hill3, initial = hill3Init, c("Ka","n","a"))