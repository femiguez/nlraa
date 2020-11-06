#' 
#' @title self start for a trilinear Function
#' @name SStrlin
#' @rdname SStrlin
#' @description Self starter for a tri-linear function with parameters a (intercept), b (first slope), xs1 (first break-point), c (second slope), xs2 (second break-point) and d (third slope)
#' @param x input vector 
#' @param a the intercept
#' @param b the first-phase slope
#' @param xs1 first break-point of transition between first-phase linear and second-phase linear
#' @param c the second-phase slope 
#' @param xs2 second break-point of transition between second-phase linear and third-phase linear
#' @param d the third-phase slope 
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This is a special case with just three parts (and two break points) but a more general approach 
#' is to consider a segmented function with several breakpoints and linear segments. 
#' Splines would be even more general. Also this model assumes that there are two break-points that needs 
#' to be estimated. The guess for the initial values splits the dataset in half, so it this will work
#' when one break-point is in the first half and the second is in the second half.
#' @seealso package \pkg{segmented}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:30
#' y <- trlin(x, 0.5, 2, 10, 0.1, 20, 1.75) + rnorm(30, 0, 0.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SStrlin(x, a, b, xs1, c, xs2, d), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' ## Minimal example
#' ## This is probably about the smallest dataset you 
#' ## should use with this function
#' dat2 <- data.frame(x = 1:8, y = c(1.1, 1.9, 3.1, 2.5, 1.4, 0.9, 2.2, 2.9))
#' fit2 <- nls(y ~ SStrlin(x, a, b, xs1, c, xs2, d), data = dat2)
#' ## expangin for plotting
#' ndat <- data.frame(x = seq(1, 8, by = 0.1))
#' ndat$prd <- predict(fit2, newdata = ndat)
#' ggplot() + 
#'   geom_point(data = dat2, aes(x = x, y = y)) + 
#'   geom_line(data = ndat, aes(x = x, y = prd))
#' 
#' }
NULL

trlinInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)

  if(nrow(xy) < 5){
    stop("Too few distinct input values to fit a trlinear")
  }
  ## Splitting this into two bilinear problems
  xy1 <- xy[1:(ceiling(nrow(xy)/2)),]
  xy2 <- xy[floor(nrow(xy)/2):nrow(xy),]
  cfs1 <- getInitial(y ~ SSblin(x, a, b, xs, c), data = xy1)
  cfs2 <- getInitial(y ~ SSblin(x, a, b, xs, c), data = xy2)
  
  a <- cfs1[1]
  b <- cfs1[2]
  xs1 <- cfs1[3]
  c <- mean(c(cfs1[4], cfs2[2]))
  xs2 <- cfs2[3]
  d <- cfs2[4]
  
  value <- c(a, b, xs1, c, xs2, d)
  names(value) <- mCall[c("a","b","xs1","c","xs2","d")]
  value
}

#' @rdname SStrlin
#' @return trlin: vector of the same length as x using the tri-linear function
#' @export
trlin <- function(x, a, b, xs1, c, xs2, d){
  
  .a2 <- a + b * xs1 ## This is the second intercept
  
  .a3 <- a + b * xs1 + c * (xs2 - xs1) ## This is the third intercept
  
  .value <- (x < xs1) * (a + b * x) + 
    (x >= xs1) * (x < xs2) * (.a2 + c * (x - xs1)) +
    (x >= xs2) * (.a3 + d * (x - xs2))
  
  ## Derivative with respect to a 
  ## if x > xs2 := exp1 <- deriv(~(a + b * xs1 + c * (xs2 - xs1)) + d * (x - xs2), "a")
  ## if x > xs1 := exp2 <- deriv(~(a + b * xs1) + c * (x - xs1), "a")
  ## if x < xs1 : = exp2 <- deriv(~(a + b * x), "b")
  .exp1 <- 1 
  
  ## Derivative with respect to b
  ## if x > xs2 := exp1 <- deriv(~(a + b * xs1 + c * (xs2 - xs1)) + d * (x - xs2), "b")
  ## if x > xs1 := exp2 <- deriv(~(a + b * xs1) + c * (x - xs1), "b")
  ## if x < xs1 : = exp2 <- deriv(~(a + b * x), "b")
  .exp2 <- ifelse(x > xs1, xs1, x)
  
  ## Derivative with respect to xs1
  ## if x > xs2 := exp1 <- deriv(~(a + b * xs1 + c * (xs2 - xs1)) + d * (x - xs2), "xs1")
  ## if x > xs1 := exp3 <- deriv(~(a + b * xs1) + c * (x - xs1), "xs1")
  ## if x < xs1 := exp3 <- 0
  .exp3 <- ifelse(x > xs1, b - c, 0)
  
  ## Derivative with respect to c
  ## if x > xs2 := exp1 <- deriv(~(a + b * xs1 + c * (xs2 - xs1)) + d * (x - xs2), "c")
  ## if x > xs1 := exp4 <- deriv(~(a + b * xs1) + c * (x - xs1), "c")
  ## if x < xs1 := exp4 <- deriv(~(a + b * xs1) + c * (x - xs1), "c")
  .exp4 <- ifelse(x > xs2, xs2 - xs1, ifelse(x > xs1, x - xs1, 0))
  
  ## Derivative with respect to xs2
  ## if x > xs2 := exp1 <- deriv(~(a + b * xs1 + c * (xs2 - xs1)) + d * (x - xs2), "xs2")
  ## if x > xs1 := exp5 <- deriv(~(a + b * xs1) + c * (x - xs1), "xs2")
  ## if x < xs1 := exp5 <- deriv(~(a + b * xs1) + c * (x - xs1), "xs2")
  .exp5 <- ifelse(x > xs2, c - d, 0)
  
  ## Derivative with respect to d
  ## if x > xs2 := exp1 <- deriv(~(a + b * xs1 + c * (xs2 - xs1)) + d * (x - xs2), "d")
  ## if x > xs1 := exp6 <- deriv(~(a + b * xs1) + c * (x - xs1), "xs2")
  ## if x < xs1 := exp6 <- deriv(~(a + b * xs1) + c * (x - xs1), "c")
  .exp6 <- ifelse(x > xs2, x - xs2, 0)
  
  .actualArgs <- as.list(match.call()[c("a","b","xs1","c","xs2","d")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 6L), list(NULL, c("a","b","xs1","c","xs2","d")))
    .grad[, "a"] <- .exp1
    .grad[, "b"] <- .exp2
    .grad[, "xs1"] <- .exp3
    .grad[, "c"] <- .exp4
    .grad[, "xs2"] <- .exp5
    .grad[, "d"] <- .exp6
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SStrlin
#' @export
SStrlin <- selfStart(trlin, initial = trlinInit, c("a","b","xs1","c","xs2","d"))

