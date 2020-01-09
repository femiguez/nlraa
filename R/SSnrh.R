#' 
#' @title self start for non-rectangualr hyperbola (photosynthesis)
#' @name SSnrh
#' @rdname SSnrh
#' @description Self starter for Non-rectangual Hyperbola with parameters: asymptote, quantum efficiency, curvature and dark respiration
#' @param x input vector (x) which is normally light intensity (PPFD, Photosynthetic Photon Flux Density).
#' @param asym asymptotic value for photosynthesis
#' @param phi quantum efficiency (mol CO2 per mol of photons) or initial slope of the light response
#' @param theta curvature parameter for smooth transition between limitations
#' @param rd dark respiration or value of CO2 uptake at zero light levels
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506).
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- seq(0, 2000, 100)
#' y <- nrh(x, 35, 0.04, 0.83, 2) + rnorm(length(x), 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSnrh(x, asym, phi, theta, rd), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

nrhInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a nrh")
  }
  asym <- max(xy[,"y"])
  xy1 <- xy[1:floor(nrow(xy)/2),]
  lm.cfs <- coef(stats::lm(xy1[,"y"] ~ xy1[,"x"]))
  phi <- lm.cfs[1]
  rd <- lm.cfs[2]
  theta <- 0.8
  
  objfun <- function(cfs){
    pred <- nrh(xy[,"x"], asym=cfs[1], phi=cfs[2], theta=cfs[3], rd=cfs[4])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(asym, phi, theta, rd)
  op <- try(stats::optim(cfs, objfun), silent = TRUE)
  
  if(class(op) != "try-error"){
    asym <- op$par[1]
    phi <- op$par[2]
    theta <- op$par[3]
    rd <- op$par[4]
  }
  
  value <- c(asym, phi, theta, rd)
  names(value) <- mCall[c("asym","phi","theta","rd")]
  value
  
}

#' @rdname SSnrh
#' @return nrh: vector of the same length as x (time) using the non-rectangular hyperbola
#' @examples 
#' x <- seq(0, 2000)
#' y <- nrh(x, 30, 0.04, 0.85, 2)
#' plot(x, y)
#' @export
#' 
nrh <- function(x, asym, phi, theta, rd){
  
  .expre1 <- 4 * theta * phi * x * asym
  .expre2 <- (phi * x + asym)^2
  .expre25 <- suppressWarnings(sqrt(.expre2 - .expre1))
  .expre3 <- phi * x + asym - .expre25
  .expre4 <- .expre3 / (2 * theta)
  .value <- .expre4 - rd
  .value <- ifelse(is.nan(.value),0,.value)
  
  ## Derivative with respect to asym
  ## deriv(~ (phi * x + asym - sqrt((phi * x + asym)^2 - (4 * theta * phi * x * asym)))/(2 * theta) - rd,"asym")
  .expr2 <- phi * x + asym
  .expr6 <- 4 * theta * phi * x
  .expr8 <- .expr2^2 - .expr6 * asym
  .isqrtExpr8 <- suppressWarnings(1/sqrt(.expr8))
  .expr11 <- 2 * theta
  .exp1 <- (1 - 0.5 * ((2 * .expr2 - .expr6) * .isqrtExpr8))/.expr11
  .exp1 <- ifelse(is.nan(.exp1), 0, .exp1)
  
  ## Derivative with respect to phi
  ## deriv(~ (phi * x + asym - sqrt((phi * x + asym)^2 - (4 * theta * phi * x * asym)))/(2 * theta) - rd,"phi")
  .expr4 <- 4 * theta
  .expr11 <- 2 * theta
  .exp2 <- (x - 0.5 * ((2 * (x * .expr2) - .expr4 * x * asym) * .isqrtExpr8))/.expr11
  .exp2 <- ifelse(is.nan(.exp2), 0, .exp2)
  
  ## Derivative with respect to theta
  ## deriv(~ (phi * x + asym - sqrt((phi * x + asym)^2 - (4 * theta * phi * x * asym)))/(2 * theta) - rd,"theta")
  .expr95 <- suppressWarnings(sqrt(.expr8))
  .expr10 <- .expr2 - .expr95
  .exp3 <-  0.5 * (4 * phi * x * asym * .isqrtExpr8)/.expr11 - .expr10 * 2/.expr11^2
  .exp3 <- ifelse(is.nan(.exp3), 0, .exp3)
  ## Derivative with respect to rd
  ## deriv(~ (phi * x + asym - sqrt((phi * x + asym)^2 - (4 * theta * phi * x * asym)))/(2 * theta) - rd,"rd")
  .exp4 <- -1
  
  .actualArgs <- as.list(match.call()[c("asym", "phi", "theta", "rd")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("asym", "phi", "theta","rd")))
    .grad[, "asym"] <- .exp1
    .grad[, "phi"] <- .exp2
    .grad[, "theta"] <- .exp3
    .grad[, "rd"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSnrh
#' @export
SSnrh <- selfStart(nrh, initial = nrhInit, c("asym", "phi", "theta", "rd"))
