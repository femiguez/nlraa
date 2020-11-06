#' Collatz GJ , Ribas-Carbo M Berry JA (1992) Coupled Photosynthesis-Stomatal Conductance Model for Leaves of C4 Plants. Functional Plant Biology 19, 519-538.
#' https://doi.org/10.1071/PP9920519
#' 
#' @title self start for Collatz temperature response
#' @name SStemp3
#' @rdname SStemp
#' @description Self starter for Collatz temperature response function 
#' @param x input vector (x) which is normally \sQuote{temperature}.
#' @param t.m medium temperature
#' @param t.l low temparature
#' @param t.h high temperature
#' @export
#' @examples 
#' \donttest{
#' ## A temperature response function
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:50
#' y <- temp3(x, 25, 13, 36) + rnorm(length(x), sd = 0.05)
#' dat1 <- data.frame(x = x, y = y)
#' fit1 <- nls(y ~ SStemp3(x, t.m, t.l, t.h), data = dat1)
#' 
#' ggplot(data = dat1, aes(x, y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit1)))
#' }
NULL

temp3Init <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a temp3 function.")
  }
  
  meany <- mean(xy[,"y"], na.rm = TRUE)
  ## t.m is actually the medium temparature, usually 25
  mdn.temp <- NLSstClosestX(xy, meany)
  min.temp <- min(xy[,"x"]) 
  max.temp <- max(xy[,"x"])
  
  ## Eductaed guesses
  t.m <- mdn.temp
  t.l <- (mdn.temp - min.temp)/2
  t.h <- mdn.temp + (max.temp - mdn.temp)/2
  
  value <- c(t.m, t.l, t.h)
  names(value) <- mCall[c("t.m","t.l","t.h")]
  value
}

#' @rdname SStemp
#' @return temp3: vector of the same length as x using a temp function
#' @export
#' 
temp3 <- function(x, t.m, t.l, t.h){
  
  Q10 <- 2
  num <- Q10 ^ ((x - t.m)/10)
  den1 <- 1 + exp(0.3 * (t.l - x))
  den2 <- 1 + exp(0.3 * (x - t.h))
  
  .value <- num / (den1 * den2)
  
  ## Derivative wrt t.m
  ## deriv(~(2^((x - t.m)/10))/((1 + exp(0.3 * (t.l - x)))*(1 + exp(0.3 * (x - t.h)))), c("t.m","t.l","t.h"))
  .expr3 <- 2^((x - t.m)/10)
  .expr12 <- (1 + exp(0.3 * (t.l - x))) * (1 + exp(0.3 * (x - t.h)))
  .expi1 <- -(.expr3 * (log(2) * (1/10))/.expr12)
  ## Derivative wrt t.l
  .expr6 <- exp(0.3 * (t.l - x))
  .expr11 <- 1 + exp(0.3 * (x - t.h))
  .expi2 <- -(.expr3 * (.expr6 * 0.3 * .expr11)/.expr12^2)
  ## Derivatie wrt t.h
  .expr7 <- 1 + exp(0.3 * (t.l - x))
  .expr10 <- exp(0.3 * (x - t.h))
  .expi3 <- .expr3 * (.expr7 * (.expr10 * 0.3))/.expr12^2
  
  .actualArgs <- as.list(match.call()[c("t.m","t.l","t.h")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("t.m","t.l","t.h")))
    .grad[, "t.m"] <- .expi1
    .grad[, "t.l"] <- .expi2
    .grad[, "t.h"] <- .expi3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SStemp
#' @export
SStemp3 <- selfStart(temp3, initial = temp3Init, c("t.m", "t.l", "t.h"))