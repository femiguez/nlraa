#' 
#' @title self start for an asymmetric Gaussian bell-shaped curve
#' @name SSagauss
#' @rdname SSagauss
#' @description Self starter for a type of bell-shaped curve
#' @param x input vector 
#' @param eta maximum value of y
#' @param beta parameter controlling the lower values
#' @param delta break-point separating the first and second half-bell curve
#' @param sigma1 scale for the first half
#' @param sigma2 scale for the second half
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in (doi:10.3390/rs12050827). 
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:30
#' y <- agauss(x, 10, 2, 10, 2, 6) + rnorm(length(x), 0, 0.02)
#' dat <- data.frame(x = x, y = y)
#' fit <- minpack.lm::nlsLM(y ~ SSagauss(x, eta, beta, delta, sigma1, sigma2), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

agaussInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 5){
    stop("Too few distinct input values to fit an assymmetric Gaussian curve.")
  }
  eta <- max(xy[,"y"])
  beta <- min(xy[,"y"])
  delta <- NLSstClosestX(xy, eta)
  ## Better guess for sigma1
  sigma1 <- abs((min(xy[,"x"]) - delta)/2) 
  sigma2 <- abs((max(xy[,"x"]) - delta)/2) 
  value <- c(eta, beta, delta, sigma1, sigma2)
  names(value) <- mCall[c("eta", "beta", "delta", "sigma1", "sigma2")]
  value
  
}

#' @rdname SSagauss
#' @return agauss: vector of the same length as x using an asymmetric bell-shaped Gaussian curve
#' @export
#' 
agauss <- function(x, eta, beta, delta, sigma1, sigma2){
  
  .expre0 <- (eta - beta)
  .expre1 <-  exp(-((x - delta)^2) / (2 * sigma1^2)) 
  .expre2 <- exp(-((x - delta)^2) / (2 * sigma2^2)) 
  .value <- ifelse(x < delta, 
                   beta + .expre0 * .expre1, 
                   beta + .expre0 * .expre2)
  
  ## Derivative with respect to beta
  ## deriv(~ beta + (eta - beta) * exp(-((x - delta)^2) / (2 * sigma1^2)) , "beta")
  .exp1 <- ifelse(x < delta,
                  1 - exp(-(x - delta)^2/(2 * sigma1^2)),
                  1 - exp(-(x - delta)^2/(2 * sigma2^2)))
  
  ## Derivative with respect to eta
  .exp2 <- ifelse(x < delta,
                  exp(-(x - delta)^2/(2 * sigma1^2)),
                  exp(-(x - delta)^2/(2 * sigma2^2)))
  
  ## Derivative with respect to delta
  .expr1 <- eta - beta
  .expr2 <- x - delta
  .expr4 <- 2 * sigma1^2
  .expr41 <- 2 * sigma2^2
  .expr6 <- exp(-.expr2^2/.expr4)
  .expr61 <- exp(-.expr2^2/.expr41)
  .expr7 <- -(.expr1 * (.expr6 * (2 * .expr2/.expr4)))
  .expr8 <- -(.expr1 * (.expr61 * (2 * .expr2/.expr41)))
  .exp3 <- ifelse(x < delta, 
                  .expr7,
                  .expr8)
  
  ## Derivative with respect to sigma1 and sigma2
  .expr1 <- eta - beta
  .expr3 <- (x - delta)^2
  .expr6 <- 2 * sigma1^2
  .expr8 <- exp(-.expr3/.expr6)
  .exp4 <- ifelse(x < delta, 
                  .expr1 * (.expr8 * (.expr3 * (2 * (2 * 
                                                       sigma1))/.expr6^2)),
                  0)
  .exp5 <- ifelse(x < delta, 0,
                  .expr1 * (.expr8 * (.expr3 * (2 * (2 * 
                                                       sigma1))/.expr6^2)))
  
  .actualArgs <- as.list(match.call()[c("eta", "beta", "delta", "sigma1", "sigma2")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 5L), list(NULL, c("eta", "beta", "delta", "sigma1", "sigma2")))
    .grad[, "eta"] <- .exp1
    .grad[, "beta"] <- .exp2
    .grad[, "delta"] <- .exp3
    .grad[, "sigma1"] <- .exp4
    .grad[, "sigma2"] <- .exp5
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSagauss
#' @export
SSagauss <- selfStart(agauss, initial = agaussInit, c("eta", "beta", "delta", "sigma1", "sigma2"))
