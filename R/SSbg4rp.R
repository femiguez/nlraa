#' Beta Growth Function
#' 
#' For details see the publication by Yin et al. (2003) "A Flexible Sigmoid Function of Determinate Growth".
#' This is a reparameterization of the beta growth function (4 parameters) with guaranteed constraints, so it is expected to 
#' behave numerically better than SSbgf4
#' 
#' Reparameterizing the four parameter beta growth
#' ldtm = log(t.e - t.m)
#' ldtb = log(t.m - t.b)
#' t.e = exp(lt.e)
#' t.m = exp(lt.e) - exp(ldtm)
#' t.b = (exp(lt.e) - exp(ldtm)) - exp(ldtb)
#' 
#' @title self start for the reparameterized Beta growth function with four parameters
#' @name SSbg4rp
#' @rdname SSbg4rp
#' @description Self starter for Beta Growth function with parameters w.max, lt.m, ldt
#' @param time input vector (x) which is normally 'time', the smallest value should be close to zero.
#' @param w.max value of weight or mass at its peak
#' @param lt.e log of the time at which the maximum weight or mass has been reached.
#' @param ldtm log of the difference between time at which the weight or mass reaches its peak and half its peak.
#' @param ldtb log of the difference between time at which the weight or mass reaches its peak and when it starts growing
#' @details The form of the equation is: \deqn{w.max * (1 + (exp(lt.e) - time)/exp(ldtm)) * ((time - (exp(lt.e) - exp(ldtb)))/exp(ldtb))^(exp(ldtb)/exp(ldtm))}.
#' This is a reparameterized version of the Beta-Growth function in which the parameters are unconstrained, but they are expressed in the log-scale.
#' This version is not fully unconstrained since ldtb > ldtm is not enforced.
#' @export
#' @examples 
#' \dontrun{
#' set.seed(1234)
#' x <- 1:100
#' y <- bg4rp(x, 20, log(70), log(30), log(20)) + rnorm(100, 0, 1)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSbg4rp(x, w.max, lt.e, ldtm, ldtb), data = dat)
#' ## We are able to recover the original values
#' exp(coef(fit)[2:4])
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

bg4rpInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a bg4rp.")
  }
  w.max <- max(xy[,"y"])
  lt.e <- log(NLSstClosestX(xy, w.max))
  ## Let's assume that t.b is the minimum of x
  t.b <- min(xy[,"x"])
  tetb <- exp(lt.e) - t.b
  ldtm <- log(tetb / 2)
  ldtb <- ldtm
  value <- c(w.max, lt.e, ldtm, ldtb)
  names(value) <- mCall[c("w.max","lt.e","ldtm","ldtb")]
  value
}

#' @rdname SSbg4rp
#' @return bg4rp: vector of the same length as x (time) using the beta growth function with four parameters
#' @export
bg4rp <- function(time, w.max, lt.e, ldtm, ldtb){
  
  ## Reparameterizing the four parameter beta growth
  ## ldtm = log(t.e - t.m)
  ## ldtb = log(t.m - t.b)
  ## t.e = exp(lt.e)
  ## t.m = exp(lt.e) - exp(ldtm)
  ## t.b = (exp(lt.e) - exp(ldtm)) - exp(ldtb)
  .exp0 <- (exp(lt.e) - exp(ldtm)) - exp(ldtb) ## This is t.b
  .exp1 <- exp(lt.e) - .exp0 ## t.e - t.b
  .exp2 <- .exp1 / exp(ldtm) ## This is (t.e - t.b)/(t.e - t.m)
  .exp3 <- time - .exp0 ## time - t.b
  .exp4 <- .exp3 / .exp1 ## (time - t.b)/(t.e - t.b)
  .exp5 <- .exp4^.exp2
  .exp6 <- 1 + (exp(lt.e) - time)/exp(ldtm)
  .value <- w.max * .exp6 * .exp5
  
  ## This function returns zero when time is less than t.b
  ## .value <- ifelse(time < t.b, 0, .value)  
  .value[is.nan(.value)] <- 0
  .value[.value < 0] <- 0
  
  ## The gradient is problematic
  # ## Derivative with respect to lt.e
  # ## deriv(~w.max * (1 + (exp(lt.e) - time)/exp(ldtm)) * ((time - (exp(lt.e) - exp(ldtm)) - exp(ldtb))/(exp(lt.e) - (exp(lt.e) - exp(ldtm)) - exp(ldtb)))^((exp(lt.e) - (exp(lt.e) - exp(ldtm)) - exp(ldtb))/exp(ldtm)),"lt.e")
  # .expr1 <- exp(lt.e)
  # .expr3 <- exp(ldtm)
  # .expr6 <- w.max * (1 + (.expr1 - time)/.expr3)
  # .expr7 <- .expr1 - .expr3
  # .expr9 <- exp(ldtb)
  # .expr10 <- time - .expr7 - .expr9
  # .expr12 <- .expr1 - .expr7 - .expr9
  # .expr13 <- .expr10/.expr12
  # .lexpr13 <- suppressWarnings(log(.expr13))
  # .expr14 <- .expr12/.expr3
  # .expr15 <- .expr13^.expr14
  # .expr21 <- .expr1 - .expr1
  # .expi1 <- w.max * (.expr1/.expr3) * .expr15 + .expr6 * (.expr15 * (.lexpr13 * (.expr21/.expr3)) - .expr13^(.expr14 - 1) * (.expr14 * (.expr1/.expr12 + .expr10 * .expr21/.expr12^2)))
  # .expi1 <- ifelse(!is.finite(.expi1), 0, .expi1)
  # ## Derivative with respect to ldtm
  # .expr1 <- exp(lt.e)
  # .expr2 <- .expr1 - time
  # .expr3 <- exp(ldtm)
  # .expr6 <- w.max * (1 + .expr2/.expr3)
  # .expr7 <- .expr1 - .expr3
  # .expr9 <- exp(ldtb)
  # .expr12 <- .expr1 - .expr7 - .expr9
  # .expr13 <- .expr10/.expr12
  # .expr14 <- .expr12/.expr3
  # .expr15 <- .expr13^.expr14
  # .expr29 <- .expr3^2
  # .expi2 <- .expr6 * (.expr13^(.expr14 - 1) * (.expr14 * (.expr3/.expr12 - .expr10 * .expr3/.expr12^2)) + .expr15 * (.lexpr13 * (.expr3/.expr3 - .expr12 * .expr3/.expr29))) - w.max * (.expr2 * .expr3/.expr29) * .expr15
  # .expi2 <- ifelse(!is.finite(.expi2), 0, .expi2)
  # ## Derivative with respect to ldtb
  # .expr1 <- exp(lt.e)
  # .expr3 <- exp(ldtm)
  # .expr6 <- w.max * (1 + (.expr1 - time)/.expr3)
  # .expr7 <- .expr1 - .expr3
  # .expr9 <- exp(ldtb)
  # .expr12 <- .expr1 - .expr7 - .expr9
  # .expr13 <- .expr10/.expr12
  # .expr14 <- .expr12/.expr3
  # .expr15 <- .expr13^.expr14
  # .expi3 <- -(.expr6 * (.expr15 * (.lexpr13 * (.expr9/.expr3)) + .expr13^(.expr14 - 1) * (.expr14 * (.expr9/.expr12 - .expr10 * .expr9/.expr12^2))))
  # .expi3 <- ifelse(!is.finite(.expi3), 0, .expi3)
  # 
  # .actualArgs <- as.list(match.call()[c("w.max", "lt.e", "ldtm", "ldtb")])
  # 
  # ##  Gradient
  # if (all(unlist(lapply(.actualArgs, is.name)))) {
  #    .grad <- array(0, c(length(.value), 4L), list(NULL, c("w.max", "lt.e", "ldtm","ldtb")))
  #    .grad[, "w.max"] <- .exp6 * .exp5
  #    .grad[, "lt.e"] <- .expi1
  #    .grad[, "ldtm"] <- .expi2
  #    .grad[, "ldtb"] <- .expi3
  #    dimnames(.grad) <- list(NULL, .actualArgs)
  #    attr(.value, "gradient") <- .grad
  #  }
  .value
}

#' @rdname SSbg4rp
#' @export
SSbg4rp <- selfStart(bg4rp, initial = bg4rpInit, c("w.max", "lt.e", "ldtm", "ldtb"))


