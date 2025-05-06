#' 
#' @title self start for a generalized Michaelis-Menten function
#' @name SSgmicmen
#' @rdname SSgmicmen
#' @description Self starter for a Michealis-Menten function with parameters  
#' @param x input vector 
#' @param ymax the maximum value
#' @param c parameter controlling the shape of the function
#' @param k parameter controlling the time for maximum value
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details Source: A generalized Michaelis-Menten equation for the analysis of growth. Lopez et al. J. Anim. Sci. 2000. 78:1816-1828.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 1:30
#' y <- gmicmen(x, 1, 10, 0.7, 10) + rnorm(30, 0, 0.01)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSgmicmen(x, y0, ymax, c, k), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Confidence intervals
#' confint(fit)
#' 
#' ## Different value for c
#' y <- gmicmen(x, 10, 5, 10) + rnorm(30, 0, 0.01)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSgmicmen(x, ymax, c, k), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Confidence intervals
#' confint(fit)
#' }
#' 
NULL

gmicmenInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a generalize Michaelis-Menten")
  }

  ## brute-force
  objfun <- function(cfs){
    pred <- gmicmen(xy[,"x"], y0 = cfs[1], ymax = cfs[2], c = cfs[3], k = cfs[4])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(min(xy[,"y"]), max(xy[,"y"]), 1, mean(xy[,"x"]))
  op <- try(stats::optim(cfs, objfun), silent = TRUE)
  
  if(!inherits(op, "try-error")){
    y0 <- op$par[1]
    ymax <- op$par[2]
    c <- op$par[3]
    k <- op$par[4]
  }else{
    ## If it fails...
    y0 <- min(xy[,"y"])
    ymax <- max(xy[,"y"])
    c <- 1
    k <- mean(xy[,"x"])
  }
  
  value <- c(y0, ymax, c, k)
  names(value) <- mCall[c("y0", "ymax", "c", "k")]
  value
}

#' @rdname SSgmicmen
#' @return plin: vector of the same length as x using the plateau-linear function
#' @export
gmicmen <- function(x, y0, ymax, c, k){
  
  .v1 <- y0 * k^c
  .v2 <- ymax * x^c
  .value <- (.v1 + .v2) / (k^c + x^c) 

  ## Derivative with respect to y0
  ## deriv(~ (y0 * k^c + ymax * x^c) / (k^c + x^c), "y0")
  .expr1 <- k^c
  .expr3 <- x^c
  .expr6 <- .expr1 + .expr3
  .exp1 <- .expr1/.expr6
    
  ## Derivative with respect to ymax
  ## deriv(~ (y0 * k^c + ymax * x^c) / (k^c + x^c), "ymax")
  .expr1 <- k^c
  .expr3 <- x^c
  .expr6 <- .expr1 + .expr3
  .exp2 <- .expr3/.expr6 
  
  ## Derivative with respect to c
  ## deriv(~ (y0 * k^c + ymax * x^c) / (k^c + x^c), "c")
  .expr5 <- y0 * .expr1 + ymax * .expr3
  .expr9 <- .expr1 * suppressWarnings(log(k))
  .expr12 <- .expr3 * suppressWarnings(log(x))
  .exp3 <- (y0 * .expr9 + ymax * .expr12)/.expr6 - .expr5 * 
    (.expr9 + .expr12)/.expr6^2
  
  ## Derivative with respect to k
  ## deriv(~ (y0 * k^c + ymax * x^c) / (k^c + x^c), "k")
  .expr10 <- k^(c - 1) * c
  .exp4 <- y0 * .expr10/.expr6 - .expr5 * .expr10/.expr6^2
  
  .actualArgs <- as.list(match.call()[c("y0", "ymax", "c", "k")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("y0", "ymax", "c", "k")))
    .grad[, "y0"] <- .exp1
    .grad[, "ymax"] <- .exp2
    .grad[, "c"] <- .exp3
    .grad[, "k"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSgmicmen
#' @export
SSgmicmen <- selfStart(gmicmen, initial = gmicmenInit, c("y0", "ymax", "c", "k"))

