#' 
#' @title self start for plateau-exponential function
#' @name SSpexpf
#' @rdname SSpexpf
#' @description Self starter for an plateau-exponential function
#' @param x input vector (x) 
#' @param a represents the value for the plateau
#' @param xs represents the breakpoint at which the plateau ends
#' @param c represents the exponential rate
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details The equation is: \eqn{for x < xs: y = a and x >= xs: a * exp(c * (x-xs))}.
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(1234)
#' x <- 1:30
#' y <- pexpf(x, 20, 15, -0.2) + rnorm(30, 0, 1)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSpexpf(x, a, xs, c), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

pexpfInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit an plateau-exponential.")
  }
  
  if(any(xy[,"y"] < 0)) stop("negative values are not allowed.")
  ## On the log scale
  xy2 <- xy[floor(nrow(xy)/2):nrow(xy),]
  ## Fit to second half of the data
  fit <- try(stats::lm(log(xy2[,"y"]) ~ xy2[,"x"]), silent = TRUE)
  
  if(class(fit) == "try-error"){
    ## I don't see any reason why 'fit' should fail..., but in that case...
    c <- (xy2[nrow(xy2),"y"] - xy2[1,"y"])/(xy2[nrow(xy2),"x"] - xy2[1,"x"]) ## Average slope
  }else{
    c <- coef(fit)[2]
  }
  
  objfun <- function(cfs){
    pred <- pexpf(xy[,"x"], a=cfs[1], xs=cfs[2], c=cfs[3])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  a <- xy2[1,"y"] ## First observation in the sorted data
  cfs <- c(a,mean(xy[,"x"]),c)
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-B",
                         upper = c(Inf, max(xy[,"x"]),Inf),
                         lower = c(-Inf, min(xy[,"x"]),-Inf)), silent = TRUE)
 
  if(class(op) != "try-error"){
    a <- op$par[1]
    xs <- op$par[2]
    c <- op$par[3]
  }else{
    ## If it fails I use the mean for the breakpoint
    xs <- mean(xy[,"x"])
  }
  
  value <- c(a, xs, c)
  names(value) <- mCall[c("a","xs","c")]
  value
  
}

#' @rdname SSpexpf
#' @return pexpf: vector of the same length as x using the pexpf function
#' @export
pexpf <- function(x, a, xs, c){
  
  .value <- (x < xs) * a + (x >= xs) * (a * exp(c *(x - xs)))
  
  ## Derivative with respect to a, xs, c
  .exp1 <- ifelse(x < xs, 1, exp(c *(x - xs))) 
  
  ## Derivative with respect to xs
  ## deriv(~ a * exp(c * (x - xs)), c("xs"))
  .exp2 <- ifelse(x < xs, 0, -(a * exp(c * (x - xs)) * c))
  
  ## Derivative with respect to c
  ## deriv(~ a * exp(c * (x - xs)), c("c"))
  .exp3 <- ifelse(x < xs, 0, a * exp(c * (x - xs)) * (x - xs))
  
  .actualArgs <- as.list(match.call()[c("a","xs","c")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a", "xs", "c")))
    .grad[, "a"] <- .exp1
    .grad[, "xs"] <- .exp2
    .grad[, "c"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSpexpf
#' @export
SSpexpf <- selfStart(pexpf, initial = pexpfInit, c("a", "xs", "c"))

