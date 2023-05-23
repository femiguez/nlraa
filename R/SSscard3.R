#' 
#' @title self start for cardinal temperature response
#' @name SScard3
#' @rdname SScard3
#' @description Self starter for cardinal temperature response function 
#' @param x input vector (x) which is normally \sQuote{temperature}.
#' @param tb base temperature
#' @param to optimum temperature
#' @param tm maximum temperature
#' @author Caio dos Santos and Fernando Miguez
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506) 
#' @export
#' @examples 
#' \donttest{
#' ## A temperature response function
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:50
#' y <- card3(x, 13, 25, 36) + rnorm(length(x), sd = 0.05)
#' dat1 <- data.frame(x = x, y = y)
#' fit1 <- nls(y ~ SScard3(x, tb, to, tm), data = dat1)
#' 
#' ggplot(data = dat1, aes(x, y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit1)))
#' }
NULL

## Function initializer
card3Init<- function(mCall, LHS, data, ...) {
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit this model.")
  }
  
  ymax <- max(xy[,"y"])
  to <- NLSstClosestX(xy, ymax)
  tb <- mean(c(xy[1, 'x'], to))
  tm <- mean(c(xy[nrow(xy), 'x'], to))
  
  objfun <- function(cfs){
    pred <- card3(xy[,"x"], 
                  tb = cfs[1],
                  to = cfs[2],
                  tm = cfs[3]
    )
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  
  cfs <- c(tb, to, tm)
  
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-B",
                         upper = c(max(xy[,"x"]), max(xy[,"x"]), max(xy[,"x"])),
                         lower = c(min(xy[,"x"]), min(xy[,"x"]), min(xy[,"x"]))), silent = TRUE)
  
  if(!inherits(op, "try-error")){
    tb <- op$par[1]
    to <- op$par[2]
    tm <- op$par[3]
  }
  
  value <- c(tb, to, tm)
  names(value) <- mCall[c('tb', 'to', 'tm')]
  return(value)
}

## the function itself
#' @rdname SScard3
#' @return card3: vector of the same length as x using a card3 function
#' @export
#' 
card3 <- function(x,  tb, to, tm){
  
  a <- ifelse(x<tb, 0,
              ifelse(x<to, (x-tb)*(1/(to-tb)),
                     ifelse(x<tm, 1- ((x-to) * (1/(tm-to)) ), 0)
              ))
  return(a)
}

## Self starting function
#' @rdname SScard3
#' @export
SScard3 <- selfStart(card3, 
                     initial = card3Init, 
                     c('tb', 'to', 'tm'))

#' 
#' An example application can be found in (doi:10.1016/j.envsoft.2014.04.009)
#' 
#' @title self start for smooth cardinal temperature response
#' @name SSscard3
#' @rdname SSscard3
#' @description Self starter for smooth cardinal temperature response function 
#' @param x input vector (x) which is normally \sQuote{temperature}.
#' @param tb base temperature
#' @param to optimum temperature
#' @param tm maximum temperature
#' @param curve curvature (default is 2)
#' @author Caio dos Santos and Fernando Miguez
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506) - Equation 5.1 in Table 1.
#' @export
#' @examples 
#' \donttest{
#' ## A temperature response function
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:50
#' y <- scard3(x, 13, 25, 36) + rnorm(length(x), sd = 0.05)
#' dat1 <- data.frame(x = x, y = y)
#' fit1 <- nls(y ~ SSscard3(x, tb, to, tm), data = dat1)
#' 
#' ggplot(data = dat1, aes(x, y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit1)))
#' }
NULL

## Function to sef-start the previous one
scard3Init<- function(mCall, LHS, data, ...) {
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit this model.")
  }
  
  ymax <- max(xy[,"y"])
  to <- NLSstClosestX(xy, ymax)
  tb <- mean(c(xy[1, 'x'], to))
  tm <- mean(c(xy[nrow(xy), 'x'], to))
  
  crv <- mCall$curve
  
  objfun <- function(cfs){
    pred <- scard3(xy[,"x"], 
                   tb = cfs[1],
                   to = cfs[2],
                   tm = cfs[3],
                   curve = crv)
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  
  cfs <- c(tb, to, tm)
  
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-B",
                         upper = c(max(xy[,"x"]), max(xy[,"x"]), max(xy[,"x"])),
                         lower = c(min(xy[,"x"]), min(xy[,"x"]), min(xy[,"x"]))), silent = TRUE)
  
  
  if(!inherits(op, "try-error")){
    tb <- op$par[1]
    to <- op$par[2]
    tm <- op$par[3]
  }
  
  value <- c(tb, to, tm)
  names(value) <- mCall[c('tb', 'to', 'tm')]
  return(value)
}


## The function itself
#' @rdname SSscard3
#' @return scard3: vector of the same length as x using a scard3 function
#' @export
#' 
scard3 <- function(x, tb, to, tm, curve = 2) {
  .expr1 <- (tm - x) / (tm - to)
  .expr2 <- (x - tb) / (to - tb)
  .expr3 <- (to - tb) / (tm - to)
  .expr4 <- (.expr1 * .expr2) ^ .expr3
  .expr5 <- .expr4 ^ curve
  
  .expr6 <- ifelse(x < tb, 0, 
                   ifelse(x > tm, 0, 
                          .expr5 ))
  return(.expr6)
}

## Defining the self starting function
#' @rdname SSscard3
#' @export
SSscard3 <- selfStart(scard3, 
                           initial = scard3Init, 
                           c('tb', 'to', 'tm'))