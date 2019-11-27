#' 
#' @title self start for exponential unction
#' @name SSexpf
#' @rdname SSexpf
#' @description Self starter for a simple exponential function
#' @param x input vector (x) 
#' @param a represents the value at x = 0
#' @param c represents the exponential rate
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506) and originally Johnson et al. (2010) Annals of Botany 106: 735–749, 2010. doi:10.1093/aob/mcq183
#' @export
#' @examples 
#' \dontrun{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:10
#' y <- expf(x, 10, -0.1) + rnorm(10, 0, 0.2)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSexpf(x, a, c), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

expfInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit an exponential")
  }
  
  ## On the log scale
  fit <- try(lm(log(xy[,"y"]) ~ xy[,"x"]), silent = TRUE)
  
  if(class(fit) == "try-error"){
    a <- 1
    c <- 0.25
  }else{
    a <- exp(coef(fit)[1])
    c <- coef(fit)[2]
  }
  
  value <- c(a, c)
  names(value) <- mCall[c("a","c")]
  value
  
}

#' @rdname SSexpf
#' @return vector of the same length as x using the profd function
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
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @export
SSexpf <- selfStart(expf, initial = expfInit, c("a", "c"))

