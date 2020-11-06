#' Response function:  \deqn{y = (asym - a2) / (1 + exp((xmid - time)/scal))) + a2}. 
#' 
#' \itemize{
#'  \item asym: upper asymptote
#'  \item xmid: time when y is midway between w and a 
#'  \item scal: controls the spread
#'  \item a2: lower asymptote 
#'  }
#'  
#'  The four parameter logistic \code{\link{SSfpl}} is essentially equivalent to this function,
#'  but here the interpretation of the parameters is different and this is the form used in 
#'  Oddi et. al. (2019) (see vignette). 
#' 
#' @title self start for Declining Logistic Function
#' @name SSdlf
#' @rdname SSdlf
#' @description Self starter for declining logistic function with parameters asym, a2, xmid and scal
#' @param time input vector (x) which is normally \sQuote{time}, the smalles value should be close to zero.
#' @param asym value of weight or mass at its peak (maximum)
#' @param a2 value of weight or mass at its trough (minimum)
#' @param xmid time at which half of the maximum weight or mass has bean reached.
#' @param scal scale parameter which controls the spread also interpreted in terms of time to go from xmid to approx. 0.63 asym
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @export
#' @examples 
#' \donttest{
#' ## Extended example in the vignette 'nlraa-Oddi-LFMC'
#' x <- seq(0, 17, by = 0.25)
#' y <- dlf(x, 2, 10, 8, 1)
#' plot(x, y, type = "l")
#' }
NULL

dlfInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a dlf.")
  }
  z <- xy[["y"]]
  x1 <- xy[["x"]]
  asym <- max(xy[,"y"])
  a2 <- min(xy[,"y"])
  xmid <- NLSstClosestX(xy, (asym + a2)/2)
  scal <- -xmid/2 ## This is a crude estimate for scal
  ## This is an alternative method of getting a good guess for scal
  ## It uses the idea that this decreasing logistic is really an inverted logistic
  zz <- z - a2 ## subtract a2 to have zero as the lowest value
  ## Use SSlogis to find out a guess for scal, but be silent about it
  scal2 <- try(coef(nls(zz ~ SSlogis(-x1, asym, xmid, scal)))[3], silent = TRUE)
  ## If the previous line returns an error we are silent about it and maintain the crude 
  ## guess for scal
  if(class(scal2) != "try-error")  scal <- -scal2
  value <- c(asym, a2, xmid, scal)
  names(value) <- mCall[c("asym","a2","xmid","scal")]
  value
  
}

#' @rdname SSdlf
#' @return dlf: vector of the same length as x (time) using the declining logistic function
#' @export
dlf <- function(time, asym, a2, xmid, scal){
  
  .expr1 <- (xmid - time)/scal
  .expr2 <- 1 + exp(.expr1)
  .value <- (asym - a2) / .expr2 + a2
  
  ## Derivative with respect to asym
  .expi1 <- 1/.expr2
  
  ## Derivative with respect to a2
  .expi2 <- 1 - 1/.expr2
  
  ## Derivative with respect to xmid
  .expr3 <- asym - a2
  .expr4 <- exp((xmid - time)/scal)
  .expr5 <- 1 + .expr4
  .expi3 <- -(.expr3 * (.expr4 * (1/scal))/.expr5^2)
  
  ## Derivative with respect to scal
  .eexpr1 <- asym - a2
  .eexpr2 <- xmid - time
  .eexpr4 <- exp(.eexpr2/scal)
  .eexpr5 <- 1 + .eexpr4
  .expi4 <- .eexpr1 * (.eexpr4 * (.eexpr2/scal^2))/.eexpr5^2
  
  .actualArgs <- as.list(match.call()[c("asym", "a2", "xmid", "scal")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("asym", "a2", "xmid","scal")))
    .grad[, "asym"] <- .expi1
    .grad[, "a2"] <- .expi2
    .grad[, "xmid"] <- .expi3 
    .grad[, "scal"] <- .expi4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSdlf
#' @export
SSdlf <- selfStart(dlf, initial = dlfInit, c("asym","a2","xmid","scal"))

## How to calculate partial derivatives:
##  deriv(~ (asym - a2)/(1 + exp((xmid - time)/scal)) + a2, "asym")
##  deriv(~ (asym - a2)/(1 + exp((xmid - time)/scal)) + a2, "a2")
##  deriv(~ (asym - a2)/(1 + exp((xmid - time)/scal)) + a2, "xmid")
##  deriv(~ (asym - a2)/(1 + exp((xmid - time)/scal)) + a2, "scal")
