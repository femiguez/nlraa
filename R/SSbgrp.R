#' For details see the publication by Yin et al. (2003) \dQuote{A Flexible Sigmoid Function of Determinate Growth}.
#' This is a reparameterization of the beta growth function with guaranteed constraints, so it is expected to 
#' behave numerically better than \code{\link{SSbgf}}.
#' 
#' @title self start for the reparameterized Beta growth function
#' @name SSbgrp
#' @rdname SSbgrp
#' @description Self starter for Beta Growth function with parameters w.max, lt.m and ldt
#' @param time input vector (x) which is normally \sQuote{time}, the smallest value should be close to zero.
#' @param w.max value of weight or mass at its peak
#' @param lt.e log of the time at which the maximum weight or mass has been reached.
#' @param ldt log of the difference between time at which the weight or mass reaches its peak and half its peak (\eqn{log(t.e - t.m)}).
#' @details The form of the equation is: \deqn{w.max * (1 + (exp(lt.e) - time)/exp(ldt)) * (time/exp(lt.e))^(exp(lt.e) / exp(ldt))}.
#' Given this function weight is expected to decay and reach zero again at \eqn{2*ldt}. This is a reparameterized version 
#' of the Beta-Growth function in which the parameters are unconstrained, but they are expressed in the log-scale.
#' @note In a few tests it seems that zero values of \sQuote{time} can cause the error message \sQuote{NA/NaN/Inf in foreign function call (arg 1)}, so it might be better to remove them before running this function.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' x <- 1:30
#' y <- bgrp(x, 20, log(25), log(5)) + rnorm(30, 0, 1)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSbgrp(x, w.max, lt.e, ldt), data = dat)
#' ## We are able to recover the original values
#' exp(coef(fit)[2:3])
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

bgrpInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit a bgrp.")
  }
  
  w.max <- max(xy[,"y"])
  lt.e <- log(NLSstClosestX(xy, w.max))
  ldt <- lt.e / 2
  value <- c(w.max, lt.e, ldt)
  names(value) <- mCall[c("w.max","lt.e","ldt")]
  value
}

#' @rdname SSbgrp
#' @return bgrp: vector of the same length as x (time) using the beta growth function (reparameterized).
#' @export
bgrp <- function(time, w.max, lt.e, ldt){
  
  .expr1 <- exp(lt.e) / exp(ldt)
  .expr2 <- (time/exp(lt.e))^.expr1
  .expr3 <- (1 + (exp(lt.e) - time)/exp(ldt))
  .value <- w.max * .expr3 * .expr2
  
  ## Derivative with respect to lt.e
  ## deriv(~w.max * (1 + (exp(lt.e) - time)/exp(ldt)) * (time/exp(lt.e))^(exp(lt.e) / exp(ldt)),"lt.e")
  .exp1 <- exp(lt.e)
  .exp3 <- exp(ldt)
  .exp6 <- w.max * (1 + (.exp1 - time)/.exp3)
  .exp7 <- time/.exp1
  .exp8 <- .exp1/.exp3
  .exp9 <- .exp7^.exp8
  .exp10 <- w.max * .exp8 * .exp9 + .exp6 * (.exp9 * (log(.exp7) * .exp8) - .exp7^(.exp8 - 1) * (.exp8 * (time * .exp1/.exp1^2)))
  
  ## Derivative with respect to ldt
  .ex2 <- .exp1 - time
  .ex6 <- w.max * (1 + .ex2/.exp3)
  .ex7 <- time/.exp1
  .ex9 <- .ex7^(.exp1/.exp3)
  .ex13 <- .exp3^2
  .ex14 <- -(.ex6 * (.ex9 * (log(.ex7) * (.exp1 * .exp3/.ex13))) + w.max * (.ex2 * .exp3/.ex13) * .ex9)
  
  .actualArgs <- as.list(match.call()[c("w.max", "lt.e", "ldt")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("w.max", 
                                                          "lt.e", "ldt")))
    .grad[, "w.max"] <- .expr3 * .expr2
    .grad[, "lt.e"] <- .exp10
    .grad[, "ldt"] <- .ex14 
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSbgrp
#' @export
SSbgrp <- selfStart(bgrp, initial = bgrpInit, c("w.max", "lt.e", "ldt"))


