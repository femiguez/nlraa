#' 
#' @title self start for modified hyperbola (photosynthesis)
#' @name SSmoh
#' @rdname SSmoh
#' @description Self starter for modified Hyperbola with parameters: asymp, xmin and k
#' @param x input vector (x) which is normally a controlling variable such as nitrogen
#' @param asym asymptotic value when x tends to infinity
#' @param xmin value of x for which y equals zero
#' @param k curvature parameter
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506).
#' See Table S3 (Eq. 3.8)
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- seq(3, 30)
#' y <- moh(x, 35, 3, 0.83) + rnorm(length(x), 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSmoh(x, asym, xmin, k), data = dat)
#' ## Visualize observed and simulated
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Testing predict function
#' prd <- predict_nls(fit, interval = "confidence")
#' datA <- cbind(dat, prd)
#' ## Plotting
#' ggplot(data = datA, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit))) + 
#'   geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), 
#'   fill = "purple", alpha = 0.3)
#' }
NULL

mohInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a moh")
  }
  
  asym <- max(xy[,"y"])
  xmin <- NLSstClosestX(xy, 0)
  k <- 1
    
  objfun <- function(cfs){
    pred <- moh(xy[,"x"], asym=cfs[1], xmin=cfs[2], k=cfs[3])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(asym, xmin, k)
  op <- try(stats::optim(cfs, objfun), silent = TRUE)
  
  if(class(op) != "try-error"){
    asym <- op$par[1]
    xmin <- op$par[2]
    k <- op$par[3]
  }
  
  value <- c(asym, xmin, k)
  names(value) <- mCall[c("asym", "xmin", "k")]
  value
}

#' @rdname SSmoh
#' @return moh: vector of the same length as x (time) using the modified hyperbola
#' @examples 
#' x <- seq(0, 20)
#' y <- moh(x, 30, 3, 0.9)
#' plot(x, y)
#' @export
#' 
moh <- function(x, asym, xmin, k){
  
  .value <- asym * (x - xmin)/(x + k)
  
  ## Derivative with respect to asym
  ## deriv(~ asym * (x - xmin)/(x + k),"asym")
  .expr1 <- x - xmin
  .expr3 <- x + k
  .exp1 <- .expr1 / .expr3
  
  ## Derivative with respect to xmin
  ## deriv(~ asym * (x - xmin)/(x + k),"xmin")
  .exp2 <- -(asym/.expr3)
  
  ## Derivative with respect to theta
  ## deriv(~ asym * (x - xmin)/(x + k),"k")
  .expr2 <- asym * (x - xmin)
  .exp3 <- -(.expr2/.expr3^2)
  
  .actualArgs <- as.list(match.call()[c("asym", "xmin", "k")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("asym", "xmin", "k")))
    .grad[, "asym"] <- .exp1
    .grad[, "xmin"] <- .exp2
    .grad[, "k"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSmoh
#' @export
SSmoh <- selfStart(moh, initial = mohInit, c("asym", "xmin", "k"))
