#'
#' @title self start for profile decay function
#' @name SSprofd
#' @rdname SSprofd
#' @description Self starter for a decay of a variable within a canopy (e.g. nitrogen concentration).
#' @param x input vector (x) 
#' @param a represents the maximum value
#' @param b represents the minimum value
#' @param c represents the rate of decay
#' @param d represents an empirical coefficient which provides flexibility
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506) and originally in Johnson et al. (2010) Annals of Botany 106: 735–749, 2010. (doi:10.1093/aob/mcq183).
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:10
#' y <- profd(x, 0.3, 0.05, 0.5, 4) + rnorm(10, 0, 0.01)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSprofd(x, a, b, c, d), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## profiling
#' ## It does not work at a lower alphamax level
#' ## Use this if you want to look at all four parameters
#' ## par(mfrow=c(2,2))
#' plot(profile(fit, alphamax = 0.016))
#' ## Reset graphical parameter as appropriate: par(mfrow=c(1,1))
#' ## parameter 'd' is not well constrained
#' confint(fit, level = 0.9)
#' }
#' 
NULL

profdInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a profd.")
  }
  ## Use the log transform
  a <- max(xy[,"y"])
  b <- min(xy[,"y"])
  
  ry <- (xy[,"y"] - b)/a 
  dat <- data.frame(ry = ry, x = xy[,"x"])
  dat <- subset(dat, dat[,"ry"] > 0 & dat[,"x"] > 0)
  ## Initial crude guess for c 
  fit <- try(nls(ry ~ 1 - (1 - exp(-c * x))^d, data = dat, 
                 start = list(c = 0.5, d = 1),
                 algorithm = "port",
                 lower = c(0,0)), silent = TRUE)
  
  if(class(fit) == "try-error"){
    c <- 0.5
    d <- 1
  }else{
    c <- coef(fit)[1]
    d <- coef(fit)[2]
  }
  
  value <- c(a, b, c, d)
  names(value) <- mCall[c("a","b","c","d")]
  value
}

#' @rdname SSprofd
#' @return profd: vector of the same length as x using the profd function
#' @export
profd <- function(x, a, b, c, d){
  
  .value <- a  - (a - b) * (1 - exp(-c * x))^d
  
  ## Derivative with respect to a, b, c, d
  ## deriv(~ a - (a - b) * (1 - exp(-c * x))^d, c("a","b","c","d"))
  .expr6 <- (1 - exp(-c * x))^d
  .exp1 <- 1 - .expr6
  .exp1 <- ifelse(is.nan(.exp1), 0, .exp1)
  
  ## Derivative with respect to b
  .exp2 <- .expr6
  .exp2 <- ifelse(is.nan(.exp2), 0, .exp2)
  
  ## Derivative with respect to c
  .expr1 <- a - b
  .expr4 <- exp(-c * x)
  .expr5 <- 1 - .expr4
  .exp3 <- -(.expr1 * (.expr5^(d - 1) * (d * (.expr4 * x))))
  .exp3 <- ifelse(is.nan(.exp3), 0, .exp3)
  
  ## Derivative with respect to d
  .lexpr5 <- suppressWarnings(log(.expr5))
  .exp4 <- -(.expr1 * (.expr6 * .lexpr5))
  .exp4 <- ifelse(is.nan(.exp4),0,.exp4)
  
  .actualArgs <- as.list(match.call()[c("a", "b", "c", "d")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a", "b", "c", "d")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    .grad[, "c"] <- .exp3
    .grad[, "d"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSprofd
#' @export
SSprofd <- selfStart(profd, initial = profdInit, c("a", "b", "c", "d"))

