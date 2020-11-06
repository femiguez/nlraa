#' 
#' @title self start for an exponential function
#' @name SSexpf
#' @rdname SSexpf
#' @description Self starter for a simple exponential function
#' @param x input vector (x) 
#' @param a represents the value at x = 0
#' @param c represents the exponential rate
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506). 
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:15
#' y <- expf(x, 10, -0.3) + rnorm(15, 0, 0.2)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSexpf(x, a, c), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

expfInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit an exponential")
  }
  
  if(any(xy[,"y"] < 0)) stop("negative values in y are not allowed.")
  
  ## On the log scale for 'y'
  fit <- try(stats::lm(log(xy[,"y"]) ~ xy[,"x"], na.action = "na.omit"), silent = TRUE)
  
  if(class(fit) == "try-error"){
    ## I don't see any reason why 'fit' should fail..., but in that case...
    a <- xy[1,"y"] ## First observation in the sorted data
    c <- (xy[nrow(xy),"y"] - xy[1,"y"])/(xy[nrow(xy),"x"] - xy[1,"x"]) ## Average slope
  }else{
    a <- exp(coef(fit)[1])
    c <- coef(fit)[2]
  }
  
  value <- c(a, c)
  names(value) <- mCall[c("a","c")]
  value
  
}

#' @rdname SSexpf
#' @return expf: vector of the same length as x using the profd function
#' @export
expf <- function(x, a, c){
  
  .value <- a * exp(c * x)
  
  ## Derivative with respect to a, b, c, d
  ## deriv(~ a * exp(c * x), c("a"))
  .exp1 <- exp(c * x)
  
  ## Derivative with respect to c
  ## deriv(~ a * exp(c * x), c("c"))
  .exp2 <- a * (.exp1 * x)
  
  .actualArgs <- as.list(match.call()[c("a","c")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 2L), list(NULL, c("a", "c")))
    .grad[, "a"] <- .exp1
    .grad[, "c"] <- .exp2
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSexpf
#' @export
SSexpf <- selfStart(expf, initial = expfInit, c("a", "c"))

