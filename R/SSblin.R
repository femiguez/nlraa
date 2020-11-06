#' 
#' @title self start for a bilinear Function
#' @name SSblin
#' @rdname SSblin
#' @description Self starter for a bilinear function with parameters a (intercept), b (first slope), xs (break-point), c (second slope)
#' @param x input vector 
#' @param a the intercept
#' @param b the first-phase slope
#' @param xs break-point of transition between first-phase linear and second-phase linear
#' @param c the second-phase slope 
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This is a special case with just two parts but a more general approach is to consider a segmented 
#' function with several breakpoints and linear segments. Splines would be even more general. Also this
#' model assumes that there is a break-point that needs to be estimated.
#' @seealso package \pkg{segmented}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:30
#' y <- blin(x, 0, 0.75, 15, 1.75) + rnorm(30, 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSblin(x, a, b, xs, c), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Minimal example
#' ## This is probably about the smallest dataset you 
#' ## should use with this function
#' dat2 <- data.frame(x = 1:5, y = c(1.1, 1.9, 3.1, 2, 0.9))
#' fit2 <- nls(y ~ SSblin(x, a, b, xs, c), data = dat2)
#' ggplot(data = dat2, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit2)))
#' 
#' }
NULL

blinInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a blinear")
  }
  ## Dumb guess for a and b is to fit a linear regression to the first
  ## half and another linear regression to the second half
  xy1 <- xy[1:(floor(nrow(xy)/2)),]
  xy2 <- xy[floor(nrow(xy)/2):nrow(xy),]
  fit1 <- stats::lm(xy1[,"y"] ~ xy1[,"x"])
  fit2 <- stats::lm(xy2[,"y"] ~ xy2[,"x"])

  objfun <- function(cfs){
    pred <- blin(xy[,"x"], a=cfs[1], b=cfs[2], xs=cfs[3], c = cfs[4])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(coef(fit1),mean(xy[,"x"]),coef(fit2)[2])
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-B",
                         upper = c(Inf, Inf, max(xy[,"x"]),Inf),
                         lower = c(-Inf, -Inf, min(xy[,"x"])),-Inf), silent = TRUE)
  
  if(class(op) != "try-error"){
    a <- op$par[1]
    b <- op$par[2]
    xs <- op$par[3]
    c <- op$par[4]
  }else{
    ## If it fails I use the mean
    a <- coef(fit1)[1]
    b <- coef(fit1)[2]
    xs <- mean(xy[,"x"])
    c <- coef(fit2)[2]
  }
  
  value <- c(a, b, xs, c)
  names(value) <- mCall[c("a","b","xs","c")]
  value
}

#' @rdname SSblin
#' @return blin: vector of the same length as x using the bilinear function
#' @export
blin <- function(x, a, b, xs, c){
  
  .a2 <- a + b * xs ## This is the second intercept
  
  .value <- (x < xs) * (a + b * x) + (x >= xs) * (.a2 + c * (x - xs))
  
  ## Derivative with respect to a when (x < xs)
  ## exp1 <- deriv(~ a +  b * x, "a")
  ## exp1 <- deriv(~((a + b * xs) + c * (x - xs)), "a")
  .exp1 <- 1 ## ifelse(x <= xs, 1, 1)
  
  ## Derivative with respect to b
  ## exp2 <- deriv(~ a +  b * x, "b")
  ## exp2 <- deriv(~((a + b * xs) + c * (x - xs)), "a")
  .exp2 <- ifelse(x < xs, x, xs)
  
  ## Derivative with respect to xs
  ## exp3 <- deriv(~ a +  b * xs, "xs")
  ## exp3 <- deriv(~ (a + b*xs) +  c * (x - xs), "xs")
  .exp3 <- ifelse(x < xs, 0, b - c)
  
  ## Derivative with respect to c
  ## exp4 <- deriv(~ a +  b * xs, "c")
  ## exp4 <- deriv(~ (a + b*xs) +  c * (x - xs), "c")
  .exp4 <- ifelse(x < xs, 0, x - xs)
  
  .actualArgs <- as.list(match.call()[c("a","b","xs","c")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 4L), list(NULL, c("a","b","xs","c")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    .grad[, "xs"] <- .exp3
    .grad[, "c"] <- .exp4
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
   }
  .value
}

#' @rdname SSblin
#' @export
SSblin <- selfStart(blin, initial = blinInit, c("a","b","xs","c"))

