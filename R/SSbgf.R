#' Beta Growth Function
#' 
#' For details see the publication by Yin et al. (2003) "A Flexible Sigmoid Function of Determinate Growth"
#' 
#' @title self start for Beta Growth Function
#' @name SSbgf
#' @rdname SSbgf
#' @description Self starter for Beta Growth function with parameters w.max, t.m and t.e
#' @param time input vector (x) which is normally 'time', the smallest value should be close to zero.
#' @param w.max value of weight or mass at its peak
#' @param t.m time at which half of the maximum weight or mass has bean reached.
#' @param t.e time at which the weight or mass reaches its peak.
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @details Given this function weight is expected to decay and reach zero again at 2*t.e - t.m
#' @export
#' @examples 
#' \dontrun{
#' data(sm)
#' Examples from old vignette
#' }
NULL

bgfInit <- function(mCall, LHS, data){

  xy <- sortedXyData(mCall[["time"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a bgf")
  }
  w.max <- max(xy[,"y"])
  t.e <- NLSstClosestX(xy, w.max)
  t.m <- t.e / 2
  value <- c(w.max, t.e, t.m)
  names(value) <- mCall[c("w.max","t.e","t.m")]
  value

}

#' @rdname SSbgf
#' @return vector of the same length as x (time) using the beta growth function
#' @examples 
#' x <- seq(0, 17, by = 0.25)
#' y <- bgf(x, 5, 10, 3)
#' plot(x, y)
#' @export
bgf <- function(time, w.max, t.e, t.m){

  .expr1 <- t.e / (t.e - t.m)
  .expr2 <- (time/t.e)^.expr1
  .expr3 <- (1 + (t.e - time)/(t.e - t.m))
  .value <- w.max * .expr3 * .expr2

  ## Derivative with respect to t.e
  .exp1 <- ((time/t.e)^(t.e/(t.e - t.m))) * ((t.e-time)/(t.e-t.m) + 1)
  .exp2 <- (log(time/t.e)*((1/(t.e-t.m) - (t.e/(t.e-t.m)^2) - (1/(t.e - t.m)))))*w.max
  .exp3 <- (time/t.e)^(t.e/(t.e-t.m))
  .exp4 <- w.max * ((1/(t.e-t.m)) - ((t.e - time)/(t.e-t.m)^2))
  .exp5 <- .exp1 * .exp2 + .exp3 * .exp4 

  ## Derivative with respect to t.m
  .ex1 <- t.e * (time/t.e)^((t.e/(t.e - t.m))) * log(time/t.e) * ((t.e - time)/(t.e - t.m) + 1) * w.max
  .ex2 <- (t.e - time) * w.max * (time/t.e)^(t.e/(t.e-t.m))
  .ex3 <- (t.e - t.m)^2
  .ex4 <- .ex1 / .ex3 + .ex2 / .ex3
  
  .actualArgs <- as.list(match.call()[c("w.max", "t.e", "t.m")])

##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("w.max", 
                                                          "t.e", "t.m")))
    .grad[, "w.max"] <- .expr3 * .expr2
    .grad[, "t.e"] <- .exp5
    .grad[, "t.m"] <- .ex4 
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
    .value
}

#' @rdname SSbgf
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @examples 
#' y <- bgf(1:10, 5, 3, 10)
#' @export
SSbgf <- selfStart(bgf, initial = bgfInit, c("w.max", "t.e", "t.m"))

## Beta growth initial growth

#' @rdname SSbgf
#' @return a numeric vector of the same length as x (time) containing parameter estimates for equation specified
#' @param w.b weight or biomass at initial time
#' @param t.b initial time offset
#' @export
bgf2 <- function(time, w.max, w.b, t.e, t.m, t.b){

  .expr1 <- (t.e - t.b) / (t.e - t.m) 
  .expr11 <- (time - t.b) 
  .expr2 <- .expr11/(t.e-t.b)
  .expr3 <- .expr2 ^ (.expr1) 
  .expr4 <- 1 + (t.e - time)/(t.e - t.m)
  .value <- w.b + (w.max - w.b) * .expr4 * .expr3

  .value[is.nan(.value)] <- 0
  .value
}
