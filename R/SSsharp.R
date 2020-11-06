#' For details see Schoolfield, R. M., Sharpe, P. J. & Magnuson, C. E. Non-linear regression of biological temperature-dependent rate models based on absolute reaction-rate theory. Journal of Theoretical Biology 88, 719-731 (1981)
#' 
#' @title self start for temperature response
#' @name SSsharp
#' @rdname SSsharp
#' @description Self starter for temperature response function 
#' @param temp input vector (x) which is normally \sQuote{temperature}.
#' @param r_tref rate at the standardised temperature, tref
#' @param e activation energy (eV)
#' @param el low temperature de-activation energy (eV)
#' @param tl temperature at which the enzyme is half active and half suppressed due to low temperatures
#' @param eh high temperature de-activation energy (eV)
#' @param th temperature at which enzyme is half active and half suppressed dut to high temperatures
#' @param tref standardisation temperature in degrees centigrade. Temperature at which rates are not inactivated by either high or low temperatures. Typically, 25 degrees.
#' @note I do not recommend using this function.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' require(minpack.lm)
#' 
#' temp <- 0:45
#' rate <- sharp(temp, 1, 0.03, 1.44, 28, 19, 44) + rnorm(length(temp), 0, 0.05)
#' dat <- data.frame(temp = temp, rate = rate)
#' ## Fit model
#' fit <- nlsLM(rate ~ SSsharp(temp, r_tref, e, el, tl, eh, th, tref = 20), data = dat)
#' ## Visualize
#' ggplot(data = dat, aes(temp, rate)) + geom_point() + geom_line(aes(y = fitted(fit)))
#'
#' }
NULL

sharpInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["temp"]], LHS, data)
  if(nrow(xy) < 6){
    stop("Too few distinct input values to fit a sharp function.")
  }
  
  ## Is r_tref close to the maximum
  r_tref <- max(xy[,"y"])
  opt.temp <- NLSstClosestX(xy, r_tref)
  min.temp <- min(xy[,"x"]) 
  max.temp <- max(xy[,"x"])
  ## Simple guesses
  e <- 0.1
  el <- 1
  tl <-  min.temp + (opt.temp - min.temp)/2
  eh <- min.temp
  th <- max.temp
  
  value <- c(r_tref, e, el, tl, eh, th)
  names(value) <- mCall[c("r_tref","e","el","tl","eh","th")]
  value
  
}

#' @rdname SSsharp
#' @return sharp: vector of the same length as x using a sharp function
#' @export
#' 

sharp <- function(temp, r_tref, e, el, tl, eh, th, tref = 25){
  
  ## Set up the model
  tref <- 273.15 + tref
  k <- 8.62e-05
  boltzmann.term <- r_tref * exp(e/k * (1/tref - 1/(temp + 273.15)))
  inactivation.term <- 1/(1 + exp(-el/k * (1/(tl + 273.15) - 1/(temp + 273.15))) + exp(eh/k * (1/(th + 273.15) - 1/(temp + 273.15))))
  
  .value <- boltzmann.term * inactivation.term

  ## Derivatives with respect to r_tref 1/6
  ## deriv(~ (r_tref * exp(e/8.62e-05 * (1/(tref + 273.15) - 1/(temp + 273.15)))) * 1/(1 + exp(-el/k * (1/(tl + 273.15) - 1/(temp + 273.15))) + exp(eh/k * (1/(th + 273.15) - 1/(temp + 273.15)))),"r_tref")
  .expr5 <- 1/(temp + 273.15)
  .expr6 <- 1/(tref + 273.15) - .expr5
  .expr8 <- exp(e/k * .expr6)
  .expr25 <- 1 + exp(-el/k * (1/(tl + 273.15) - .expr5)) + exp(eh/k * (1/(th + 273.15) - .expr5))
  .expi1 <- .expr8/.expr25
    
  ## Derivatives with respect to e 2/6
  ## deriv(~ (r_tref * exp(e/8.62e-05 * (1/(tref + 273.15) - 1/(temp + 273.15)))) * 1/(1 + exp(-el/k * (1/(tl + 273.15) - 1/(temp + 273.15))) + exp(eh/k * (1/(th + 273.15) - 1/(temp + 273.15)))),"e")
  .expr5 <- 1/(temp + 273.15)
  .expr6 <- 1/(tref + 273.15) - .expr5
  .expr8 <- exp(e/k * .expr6)
  .expr25 <- 1 + exp(-el/k * (1/(tl + 273.15) - .expr5)) + exp(eh/k * (1/(th + 273.15) - .expr5))
  .expi2 <- r_tref * (.expr8 * (1/k * .expr6))/.expr25
  
  ## Derivatives with respect to el 3/6
  .expr10 <- r_tref * exp(e/k * (1/(tref + 273.15) - .expr5)) * 1
  .expr15 <- 1/(tl + 273.15) - .expr5
  .expr17 <- exp(-el/k * .expr15)
  .expi3 <- .expr10 * (.expr17 * (1/k * .expr15))/.expr25^2
  
  ## Derivative with respect to tl 4/6
  .expr12 <- -el/k
  .expr13 <- tl + 273.15
  .expr17 <- exp(.expr12 * (1/.expr13 - .expr5))
  .expi4 <- .expr10 * (.expr17 * (.expr12 * (1/.expr13^2)))/.expr25^2
  
  ## Derivative with respect to eh 5/6
  .expr22 <- 1/(th + 273.15) - .expr5
  .expr24 <- exp(eh/k * .expr22)
  .expi5 <- -(.expr10 * (.expr24 * (1/k * .expr22))/.expr25^2)
  
  ## Derivative with respect to th 6/6
  .expr19 <- eh/k
  .expr20 <- th + 273.15
  .expr24 <- exp(.expr19 * (1/.expr20 - .expr5))
  .expi6 <- .expr10 * (.expr24 * (.expr19 * (1/.expr20^2)))/.expr25^2
  
  .actualArgs <- as.list(match.call()[c("r_tref", "e", "el","tl","eh","th")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 6L), list(NULL, c("r_tref", "e", "el","tl","eh","th")))
    .grad[, "r_tref"] <- .expi1
    .grad[, "e"] <- .expi2
    .grad[, "el"] <- .expi3
    .grad[, "tl"] <- .expi4
    .grad[, "eh"] <- .expi5
    .grad[, "th"] <- .expi6
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSsharp
#' @export
SSsharp <- selfStart(sharp, initial = sharpInit, c("r_tref", "e", "el","tl","eh","th"))