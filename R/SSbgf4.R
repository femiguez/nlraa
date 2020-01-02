#' Beta Growth Function
#' 
#' For details see the publication by Yin et al. (2003) "A Flexible Sigmoid Function of Determinate Growth".
#' This is a difficult function to fit because the linear constrains are not explicitly introduced 
#' in the optimization process. 
#' 
#' @title self start for Beta growth function with four parameters
#' @name SSbgf4
#' @rdname SSbgf4
#' @description Self starter for Beta Growth function with parameters w.max, t.e, t.m and t.b
#' @param time input vector (x) which is normally 'time'.
#' @param w.max value of weight or mass at its peak
#' @param t.e time at which the weight or mass reaches its peak.
#' @param t.m time at which half of the maximum weight or mass has bean reached.
#' @param t.b time at which biomass growth starts
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @details Given this function weight is expected to decay and reach zero again at \eqn{2*t.e - t.m}. This is equation
#' 11 (pg. 368) in the Yin paper, but with errata correction (https://academic.oup.com/aob/article/91/6/753/211049) and with w.b equal to zero.
#' @export
#' @examples 
#' \dontrun{
#' data(sm)
#' ## Let's just pick one crop
#' sm2 <- subset(sm, Crop == "M")
#' ## default 'nls' fails so we are using 'minpack.lm'
#' require(minpack.lm)
#' fit <- nlsLM(Yield ~ SSbgf4(DOY, w.max, t.e, t.m, t.b), data = sm2)
#' plot(Yield ~ DOY, data = sm2)
#' lines(sm2$DOY,fitted(fit))
#' ## For this particular problem it is easier to 'fix' t.b and w.b
#' fit0 <- nls(Yield ~ bgf2(DOY, w.max, w.b = 0, t.e, t.m, t.b = 141), 
#'            data = sm2, start = list(w.max = 16, t.e= 240, t.m = 200))
#' }
NULL

bgf4Init <- function(mCall, LHS, data){

  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 5){
    stop("Too few distinct input values to fit a bgf4")
  }
  w.max <- max(xy[,"y"])
  t.e <- NLSstClosestX(xy, w.max)
  t.b <- min(xy[,"x"])
  t.m <- (t.e + t.b) / 2
  value <- c(w.max, t.e, t.m, t.b)
  names(value) <- mCall[c("w.max","t.e","t.m","t.b")]
  value

}

#' @rdname SSbgf4
#' @return bgf4: vector of the same length as x (time) using the beta growth function with four parameters
#' @examples 
#' x <- seq(0, 17, by = 0.25)
#' y <- bgf4(x, 20, 15, 10, 2)
#' plot(x, y)
#' @export
#' 
bgf4 <- function(time, w.max, t.e, t.m, t.b){

  if(t.m > t.e) stop("t.m should be smaller than t.e")
  if(t.b > t.m) stop("t.b should be smaller than t.m")
  if(t.b > t.e) stop("t.b should be smaller than t.e")
  
  .expre1 <- (t.e - t.b)/(t.e - t.m)
  .expre11 <- (time - t.b)
  .expre2 <- .expre11/(t.e - t.b)
  .expre3 <- .expre2^(.expre1)
  .expre4 <- 1 + (t.e - time)/(t.e - t.m)
  .value <- w.max * .expre4 * .expre3

  ## This function returns zero when time is less than t.b
  .value <- ifelse(time < t.b, 0, .value)  
  .value[is.nan(.value)] <- 0
  .value[.value < 0] <- 0

#   err.val <- 0
#   ## Derivative with respect to w.max
#   ## deriv(~w.max * (1 + (t.e - time)/(t.e - t.m)) * ((time - t.b)/(t.e - t.b))^((t.m - t.b)/(t.e - t.m)),"w.max")
#   .exp01 <- .expre4 * .expre3
#   .exp1 <- ifelse(is.nan(.exp01),err.val,.exp01)
#   
#   ## Derivative with respect to t.e
#   ## deriv(~w.max * (1 + (t.e - time)/(t.e - t.m)) * ((time - t.b)/(t.e - t.b))^((t.m - t.b)/(t.e - t.m)),"t.e")
#   .expr1 <- t.e - time
#   .expr2 <- t.e - t.m
#   .expr5 <- w.max * (1 + .expr1/.expr2)
#   .expr6 <- time - t.b
#   .expr7 <- t.e - t.b
#   .expr8 <- .expr6/.expr7
#   .lexpr8 <- suppressWarnings(log(.expr8))
#   .expr9 <- t.m - t.b
#   .expr10 <- .expr9/.expr2
#   .expr11 <- .expr8^.expr10
#   .expr14 <- .expr2^2
#   .exp02 <- w.max * (1/.expr2 - .expr1/.expr14) * .expr11 - 
#     .expr5 * (.expr11 * (.lexpr8 * (.expr9/.expr14)) + 
#                 .expr8^(.expr10 - 1) * (.expr10 * (.expr6/.expr7^2)))
#   .exp2 <- ifelse(is.nan(.exp02),err.val,.exp02)
#     
#   ## Derivative with respect to t.m
#   ## deriv(~w.max * (1 + (t.e - time)/(t.e - t.m)) * ((time - t.b)/(t.e - t.b))^((t.m - t.b)/(t.e - t.m)),"t.m")
#   .expr1 <- t.e - time
#   .expr2 <- t.e - t.m
#   .expr5 <- w.max * (1 + .expr1/.expr2)
#   .expr8 <- (time - t.b)/(t.e - t.b)
#   .lexpr8 <- suppressWarnings(log(.expr8))
#   .expr9 <- t.m - t.b
#   .expr11 <- .expr8^(.expr9/.expr2)
#   .expr13 <- .expr2^2
#   .exp03 <- w.max * (.expr1/.expr13) * .expr11 + .expr5 * 
#     (.expr11 * (.lexpr8 * (1/.expr2 + .expr9/.expr13)))
#   .exp3 <- ifelse(is.nan(.exp03),err.val,.exp03)
#     
#   ## Derivative with respect to t.b
#   ## deriv(~w.max * (1 + (t.e - time)/(t.e - t.m)) * ((time - t.b)/(t.e - t.b))^((t.m - t.b)/(t.e - t.m)),"t.b")
#   .expr2 <- t.e - t.m
#   .expr5 <- w.max * (1 + (t.e - time)/.expr2)
#   .expr6 <- time - t.b
#   .expr7 <- t.e - t.b
#   .expr8 <- .expr6/.expr7
#   .lexpr8 <- suppressWarnings(log(.expr8))
#   .expr10 <- (t.m - t.b)/.expr2
#   .expr11 <- .expr8^.expr10
#   .exp04 <- -(.expr5 * (.expr11 * (.lexpr8 * (1/.expr2)) + 
#                          .expr8^(.expr10 - 1) * (.expr10 * (1/.expr7 - .expr6/.expr7^2))))
#   .exp4 <- ifelse(is.nan(.exp04),err.val,.exp04)
#     
#   .actualArgs <- as.list(match.call()[c("w.max", "t.e", "t.m", "t.b")])
# 
# ##  Gradient
#   if (all(unlist(lapply(.actualArgs, is.name)))) {
#     .grad <- array(0, c(length(.value), 4L), list(NULL, c("w.max", "t.e", "t.m","t.b")))
#     .grad[, "w.max"] <- .exp1
#     .grad[, "t.e"] <- .exp2
#     .grad[, "t.m"] <- .exp3
#     .grad[, "t.b"] <- .exp4
#     dimnames(.grad) <- list(NULL, .actualArgs)
#     attr(.value, "gradient") <- .grad
#   }
    .value
}

#' @rdname SSbgf4
#' @export
SSbgf4 <- selfStart(bgf4, initial = bgf4Init, c("w.max", "t.e", "t.m", "t.b"))


