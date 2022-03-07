#' For details see the publication by Yin et al. (2003) \dQuote{A Flexible Sigmoid Function of Determinate Growth}.
#' 
#' @title self start for Beta Growth Function
#' @name SSbgf
#' @rdname SSbgf
#' @description Self starter for Beta Growth function with parameters w.max, t.m and t.e
#' @param time input vector (x) which is normally \sQuote{time}, the smallest value should be close to zero.
#' @param w.max value of weight or mass at its peak
#' @param t.m time at which half of the maximum weight or mass has been reached.
#' @param t.e time at which the weight or mass reaches its peak.
#' @details The form of the equation is: \deqn{w.max * (1 + (t.e - time)/(t.e - t.m)) * (time/t.e)^(t.e / (t.e - t.m))}.
#' Given this function weight is expected to decay and reach zero again at \eqn{2*t.e - t.m}.
#' @export
#' @examples 
#' \donttest{
#' ## See extended example in vignette 'nlraa-AgronJ-paper'
#' x <- seq(0, 17, by = 0.25)
#' y <- bgf(x, 5, 15, 7)
#' plot(x, y)
#' }
NULL

bgfInit <- function(mCall, LHS, data, ...){

  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a bgf")
  }

  w.max <- max(xy[,"y"])
  t.e <- NLSstClosestX(xy, w.max)
  t.m <- t.e / 2
  value <- c(w.max, t.e, t.m)
  names(value) <- mCall[c("w.max","t.e","t.m")]
  value

}

#' @rdname SSbgf
#' @return bgf: vector of the same length as x (time) using the beta growth function
#' @export
bgf <- function(time, w.max, t.e, t.m){

  .expr1 <- t.e / (t.e - t.m)
  .expr2 <- (time/t.e)^.expr1
  .expr3 <- (1 + (t.e - time)/(t.e - t.m))
  .value <- w.max * .expr3 * .expr2

  ## Derivative with respect to w.max
  ## deriv(~ w.max * (1 + (t.e - time)/(t.e - t.m)) * (time/t.e)^(t.e / (t.e - t.m)),"w.max")
  .expr2 <- t.e - t.m
  .expr4 <- 1 + (t.e - time)/.expr2
  .expr8 <- (time/t.e)^(t.e/.expr2)
  .expi1 <- .expr4 * .expr8
  .expi1 <- ifelse(is.nan(.expi1),0,.expi1)
  
  ## Derivative with respect to t.e
  .expr1 <- t.e - time
  .expr5 <- w.max * (1 + .expr1/.expr2)
  .expr6 <- time/t.e
  .lexpr6 <- suppressWarnings(log(.expr6))
  .expr7 <- t.e/.expr2
  .expr8 <- .expr6^.expr7
  .expr10 <- 1/.expr2
  .expr11 <- .expr2^2
  .expi2 <- w.max * (.expr10 - .expr1/.expr11) * .expr8 + .expr5 * (.expr8 * (.lexpr6 * (.expr10 - t.e/.expr11)) - .expr6^(.expr7 - 1) * (.expr7 * (time/t.e^2)))
  .expi2 <- ifelse(is.nan(.expi2),0,.expi2)

  ## Derivative with respect to t.m
  ## deriv(~ w.max * (1 + (t.e - time)/(t.e - t.m)) * (time/t.e)^(t.e / (t.e - t.m)),"t.m")
  .expr10 <- .expr2^2
  .expi3 <- w.max * (.expr1/.expr10) * .expr8 + .expr5 * (.expr8 * (.lexpr6 * (t.e/.expr10)))
  .expi3 <- ifelse(is.nan(.expi3),0,.expi3)
  
  .actualArgs <- as.list(match.call()[c("w.max", "t.e", "t.m")])

##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("w.max", "t.e", "t.m")))
    .grad[, "w.max"] <- .expi1
    .grad[, "t.e"] <- .expi2
    .grad[, "t.m"] <- .expi3 
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
    .value
}

#' @rdname SSbgf
#' @export
SSbgf <- selfStart(bgf, initial = bgfInit, c("w.max", "t.e", "t.m"))

#' Beta growth initial growth
#' @rdname SSbgf
#' @return bgf2: a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @param w.b weight or biomass at initial time
#' @param t.b initial time offset
#' @export
bgf2 <- function(time, w.max, w.b, t.e, t.m, t.b){

  .expr1 <- (t.e - t.b) / (t.e - t.m) 
  .expr11 <- (time - t.b) 
  .expr2 <- .expr11/(t.e-t.b)
  .expr3 <- .expr2 ^ (.expr1) 
  .expr4 <- 1 + (t.e - time)/(t.e - t.m)
  .value <- w.b + (w.max - w.b) * .expr4 * .expr3

  .value[is.nan(.value)] <- 0
  .value
}

bgf_maximum_growth_rate <- function(coefs){
  
  ## This is equation 9 in the paper by Yin et al. (2003)
  
  if(!missing(coefs)){
    if(length(coefs) != 3)
      stop("Length of coefs should be equal to 3", call. = FALSE)
    w.max <- coefs[1]
    t.e <- coefs[2]
    t.m <- coefs[3]
  }else{
     stop(" 'coefs' is required for this function", call. = FALSE)
  }
  
  num <- 2 * t.e - t.m
  den <- t.e * (t.e - t.m)
  prt3 <- (t.m / t.e)^ (t.m/ (t.e - t.m))
  
  ans <- as.vector(num/den * prt3 * w.max)
  
  ans
}