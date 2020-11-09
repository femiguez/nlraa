#' For details see the publication by Yin et al. (1995) \dQuote{A nonlinear model for crop development as a function of temperature }.
#' Agricultural and Forest Meteorology 77 (1995) 1-16 
#'
#' @title self start for Beta 5-parameter function
#' @name SSbeta5
#' @rdname SSbeta5
#' @description Self starter for Beta 5-parameter function 
#' @param temp input vector which is normally \sQuote{temperature}
#' @param mu mu parameter (see equation) 
#' @param tb base (low) temperature at which no expansion accurs
#' @param a parameter describing the initial increasing phase
#' @param tc critical (high) temperature at which no expasion occurs
#' @param b parameter describing the decreasing phase
#' @details The form of the equation is: \deqn{exp(mu) * (temp - tb)^a * (tc - temp)^b}.
#' @export
#' @examples 
#' \donttest{
#' require(minpack.lm)
#' require(ggplot2)
#' ## Temperature response example
#' data(maizeleafext)
#' ## Fit model
#' fit <- nlsLM(rate ~ SSbeta5(temp, mu, tb, a, tc, b), data = maizeleafext)
#' ## Visualize
#' ndat <- data.frame(temp = 0:45)
#' ndat$rate <- predict(fit, newdata = ndat)
#' ggplot() + 
#'    geom_point(data = maizeleafext, aes(temp, rate)) + 
#'    geom_line(data = ndat, aes(x = temp, y = rate))
#' }
NULL

beta5Init <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["temp"]], LHS, data)
  if(nrow(xy) < 5){
    stop("Too few distinct input values to fit a beta5")
  }
  
  ## If we guess values of 10 and 35 for tb and tc
  ## exclude zero values
  xy <- xy[xy$y > 0, ]
  logy <- log(xy[,"y"])
  tb <- min(xy[,"x"])
  tc <- max(xy[,"x"])
  x1 <- suppressWarnings(log(xy[,"x"] - tb))
  x2 <- suppressWarnings(log(tc - xy[,"x"]))
  dat <- data.frame(logy = logy, x1 = x1, x2 = x2)
  dat <- subset(dat, is.finite(logy) & is.finite(x1) & is.finite(x2))
  fit <- lm(logy ~ x1 + x2, data = dat, na.action = na.omit)
  mu <- coef(fit)[1]
  a <- coef(fit)[2]
  b <- abs(coef(fit)[3])
  
  value <- c(mu, tb, a, tc, b)
  names(value) <- mCall[c("mu","tb","a","tc","b")]
  value
}

#' @rdname SSbeta5
#' @return beta5: vector of the same length as x (temp) using the beta5 function
#' @export
#' 
beta5 <- function(temp, mu, tb, a, tc, b){
  
  .value <- exp(mu) * ((temp - tb)^a) * ((tc - temp)^b)
  .value <- ifelse(is.nan(.value), 0, .value)
  .value <- ifelse(!is.finite(.value), 0, .value)
  
  ## Get derivatives
  ## deriv(~exp(mu) * (temp - tb)^a * (tc - temp)^b, "mu")
  ## wrt mu
  .exp1 <- exp(mu) * (temp - tb)^a * (tc - temp)^b
  .exp1 <- ifelse(is.nan(.exp1), 0, .exp1)
  .exp1 <- ifelse(!is.finite(.exp1), 0, .exp1)  
  ## wrt tb
  .expr1 <- exp(mu)
  .expr2 <- temp - tb
  .expr6 <- (tc - temp)^b
  .exp2 <- -(.expr1 * (.expr2^(a - 1) * a) * .expr6)
  .exp2 <- ifelse(is.nan(.exp2), 0, .exp2)
  .exp2 <- ifelse(!is.finite(.exp2), 0, .exp2)
  ## wrt a
  .expr3 <- .expr2^a
  .lexpr2 <- suppressWarnings(log(.expr2))
  .exp3 <- .expr1 * (.expr3 * .lexpr2) * .expr6
  .exp3 <- ifelse(is.nan(.exp3), 0, .exp3)
  .exp3 <- ifelse(!is.finite(.exp3), 0, .exp3)
  ## wrt tc
  .expr4 <- exp(mu) * (temp - tb)^a
  .expr5 <- tc - temp
  .exp4 <- .expr4 * (.expr5^(b - 1) * b)
  .exp4 <- ifelse(is.nan(.exp4), 0, .exp4)
  .exp4 <- ifelse(!is.finite(.exp4), 0, .exp4)
  ## wrt b
  .lexpr5 <- suppressWarnings(log(.expr5))
  .exp5 <- .expr4 * (.expr6 * .lexpr5)
  .exp5 <- ifelse(is.nan(.exp5), 0, .exp5)
  .exp5 <- ifelse(!is.finite(.exp5), 0, .exp5)
  
  .actualArgs <- as.list(match.call()[c("mu", "tb", "a", "tc", "b")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 5L), list(NULL, c("mu", "tb", "a", "tc", "b")))
    .grad[, "mu"] <- .exp1
    .grad[, "tb"] <- .exp2
    .grad[, "a"] <- .exp3 
    .grad[, "tc"] <- .exp4 
    .grad[, "b"] <- .exp5 
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSbeta5
#' @export
SSbeta5 <- selfStart(beta5, initial = beta5Init, c("mu", "tb", "a", "tc", "b"))

## The data below is from Tollenar, but not sure how to use it yet

# maizedev0 <- c(
# trt = 1, day = 35, night = 25, rate = 0.524, rate.se = 0.101,
# trt = 2, day = 32, night = 19, rate = 0.473, rate.se = 0.097,
# trt = 3, day = 31, night = 25, rate = 0.521, rate.se = 0.052,
# trt = 4, day = 29, night = 11, rate = 0.383, rate.se = 0.073,
# trt = 5, day = 28, night = 13, rate = 0.307, rate.se = 0.073,
# trt = 6, day = 26, night = 15, rate = 0.307, rate.se = 0.057,
# trt = 7, day = 25, night = 15, rate = 0.344, rate.se = 0.048,
# trt = 8, day = 19, night = 16, rate = 0.256, rate.se = 0.048,
# trt = 9, day = 16, night = 5, rate = 0.155, rate.se =  0.020,
# trt = 10, day = 15, night = 25, rate = 0.296, rate.se = 0.051,
# trt = 11, day = 15, night = 5, rate = 0.158, rate.se = 0.033,
# trt = 12, day = 13, night = 28, rate = 0.307, rate.se = 0.046,
# trt = 13, day = 32, night = 16, rate = 0.437, rate.se = 0.054,
# trt = 14, day = 29, night = 18, rate = 0.391, rate.se = 0.096,
# trt = 15, day = 24, night = 8, rate = 0.272, rate.se = 0.070, 
# trt = 16, day = 22, night = 9, rate = 0.255, rate.se = 0.045)
# 
# w1 <- which(names(maizedev) == "trt")
# w2 <- which(names(maizedev) == "day")
# w3 <- which(names(maizedev) == "night")
# w4 <- which(names(maizedev) == "rate")
# w5 <- which(names(maizedev) == "rate.se")
# 
# maizedev <- data.frame(trt = maizedev0[w1],
#                        day = maizedev0[w2],
#                        night = maizedev[w3],
#                        rate = maizedev0[w4],
#                        rate.se = maizedev0[w5])
# 
# maizedev$temp <- (maizedev$day + maizedev$night) / 2
# 
# ggplot(data = maizedev, aes(temp, rate)) + geom_point()


