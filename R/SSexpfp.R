#' 
#' @title self start for an exponential-plateau function
#' @name SSexpfp
#' @rdname SSexpfp
#' @description Self starter for an exponential-plateau function
#' @param x input vector (x) 
#' @param a represents the value at x = 0
#' @param c represents the exponential rate
#' @param xs represents the breakpoint at which the plateau starts
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @details This function is described in Archontoulis and Miguez (2015) - (doi:10.2134/agronj2012.0506). 
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(12345)
#' x <- 1:30
#' y <- expfp(x, 10, 0.1, 15) + rnorm(30, 0, 1.5)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSexpfp(x, a, c, xs), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit)))
#' }
NULL

expfpInit <- function(mCall, LHS, data){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 3){
    stop("Too few distinct input values to fit an exponential-plateau.")
  }
  
  if(any(xy[,"y"] < 0)) stop("negative values in y are not allowed.")
  ## On the log scale
  xy1 <- xy[1:floor(nrow(xy)/2),]
  ## Fit to half the data
  fit <- try(stats::lm(log(xy1[,"y"]) ~ xy1[,"x"]), silent = TRUE)
  
  if(class(fit) == "try-error"){
    ## I don't see any reason why 'fit' should fail..., but in that case...
    a <- xy1[1,"y"] ## First observation in the sorted data
    c <- (xy1[nrow(xy1),"y"] - xy1[1,"y"])/(xy1[nrow(xy1),"x"] - xy1[1,"x"]) ## Average slope
  }else{
    a <- exp(coef(fit)[1])
    c <- coef(fit)[2]
  }
  
  objfun <- function(cfs){
    pred <- expfp(xy[,"x"], a=cfs[1], c=cfs[2], xs=cfs[3])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  cfs <- c(a,c,mean(xy[,"x"]))
  op <- try(stats::optim(cfs, objfun, method = "L-BFGS-B",
                         upper = c(Inf, Inf, max(xy[,"x"])),
                         lower = c(-Inf, -Inf, min(xy[,"x"]))), silent = TRUE)

  if(class(op) != "try-error"){
    a <- op$par[1]
    c <- op$par[2]
    xs <- op$par[3]
  }else{
    ## If it fails I use the mean for the breakpoint
    xs <- mean(xy[,"x"])
  }
  
  value <- c(a, c, xs)
  names(value) <- mCall[c("a","c","xs")]
  value
  
}

#' @rdname SSexpfp
#' @return expfp: vector of the same length as x using the expfp function
#' @export
expfp <- function(x, a, c, xs){
  
  .value <- (x < xs) * a * exp(c * x) + (x >= xs) * (a * exp(c * xs))
  
  ## Derivative with respect to a, c, xs
  ## deriv(~ a * exp(c * x), c("a"))
  .exp1 <- ifelse(x < xs, exp(c * x), exp(c * xs)) 
  
  ## Derivative with respect to c
  ## deriv(~ a * exp(c * x), c("c"))
  .exp2 <- ifelse(x < xs, a * (exp(c * x) * x), a * (exp(c * xs) * xs)) 
  
  ## Derivative with respect to xs
  ## deriv(~ a * exp(c * xs), c("xs"))
  .exp3 <- ifelse(x < xs, 0, a * (exp(c * xs) * c))
  
  .actualArgs <- as.list(match.call()[c("a","c","xs")])
  
  ##  Gradient
  if (all(unlist(lapply(.actualArgs, is.name)))) {
    .grad <- array(0, c(length(.value), 3L), list(NULL, c("a", "c", "xs")))
    .grad[, "a"] <- .exp1
    .grad[, "c"] <- .exp2
    .grad[, "xs"] <- .exp3
    dimnames(.grad) <- list(NULL, .actualArgs)
    attr(.value, "gradient") <- .grad
  }
  .value
}

#' @rdname SSexpfp
#' @export
SSexpfp <- selfStart(expfp, initial = expfpInit, c("a", "c", "xs"))

