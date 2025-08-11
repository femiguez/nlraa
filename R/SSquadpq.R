#' 
#' The equation is, for a response (y) and a predictor (x): \cr
#'   \deqn{y ~ (x <= xs) * (a + b * x + (-0.5 * b/xs) * x^2) + 
#'   (x > xs \text{\&} x <= (xs + dxs)) * (a + (-b^2)/(4 * -0.5 * b/xs)) + 
#'   (x > (xs + dxs)) * ((a + (-b^2)/(4 * -0.5 * b/xs)) + ((-0.5 * b/xs) * (x - (xs + dxs))^2))} \cr
#'   
#'  This is a somewhat complicated equation. The interpretation of the parameters are simple.
#'  The model is parameterized in terms of the log of dxs (or ldxs). The parameter is ensured
#'  to be positive by taking the exponential.
#'  
#' @title self start for quadratic-plateau-quadratic (QPQ) function 
#' @name SSquadpq
#' @rdname SSquadpq
#' @description Self starter for quadratic plateau function with (four) parameters a (intercept), b (slope), xs (break-point), ldxs (log of the difference between break-point and second break-point)
#' @param x input vector 
#' @param a the intercept
#' @param b the slope
#' @param xs first break-point
#' @param ldxs log of the difference between break-point and second break-point
#' @return a numeric vector of the same length as x containing parameter estimates for equation specified
#' @export
#' @examples 
#' \donttest{
#' require(ggplot2)
#' set.seed(123)
#' x <- 0:25
#' y <- quadpq(x, 1, 0.5, 10, 1.5) + rnorm(length(x), 0, 0.3)
#' dat <- data.frame(x = x, y = y)
#' fit <- nls(y ~ SSquadpq(x, a, b, xs, dxs), data = dat)
#' ## plot
#' ggplot(data = dat, aes(x = x, y = y)) + 
#'   geom_point() + 
#'   geom_line(aes(y = fitted(fit))) + 
#'   geom_vline(aes(xintercept = coef(fit)[3]), linetype = 2) +
#'   geom_vline(aes(xintercept = coef(fit)[3] + exp(coef(fit)[4])), linetype = 3)
#'      
#' }
NULL

quadpqInit <- function(mCall, LHS, data, ...){
  
  xy <- sortedXyData(mCall[["x"]], LHS, data)
  if(nrow(xy) < 4){
    stop("Too few distinct input values to fit a quadratic-platueau-quadratic.")
  }
  ## Guess for a, b and xs is to fit a quadratic linear regression to all the data
  fit <- lm(xy[,"y"] ~ xy[,"x"] + I(xy[,"x"]^2))
  a <- coef(fit)[1]
  b <- coef(fit)[2]
  c <- coef(fit)[3]
  xs <- -0.5 * b/c
  ldxs <- log(max(xy[, "x"], na.rm = TRUE) * 0.25)
  
  objfun <- function(cfs){
    pred <- quadpq(xy[,"x"], a=cfs[1], b=cfs[2], xs=cfs[3], ldxs = cfs[4])
    ans <- sum((xy[,"y"] - pred)^2)
    ans
  }
  
  res <- list()
  j <- 1
  for(i in 1:5){
    cfs <- c(a, b, max(xy[,"x"])/i, log(0.3 * max(xy[,"x"])/i))
    op <- try(stats::optim(cfs, objfun), silent = TRUE)
    if(!inherits(op, "try-error")){
      res[[j]] <- op
      j <- j + 1
    }
  }
  
  if(all(sapply(res, function(x) inherits(x, 'try-error')))){
    ## If all fail
    a <- coef(fit)[1]
    b <- coef(fit)[2]
    c <- coef(fit)[3]
    xs <- -0.5 * b/c
    ldxs <- log(max(xy[, "x"], na.rm = TRUE) * 0.25)   
  }else{
    vals <- sapply(res, function(x) x$value)
    if(all(is.na(vals))){
      a <- coef(fit)[1]
      b <- coef(fit)[2]
      c <- coef(fit)[3]
      xs <- -0.5 * b/c
      ldxs <- log(max(xy[, "x"], na.rm = TRUE) * 0.25)   
    }else{
      wminvals <- which.min(vals)
      op <- res[[wminvals]]      
      a <- op$par[1]
      b <- op$par[2]
      xs <- op$par[3]
      ldxs <- op$par[4] 
      if(ldxs > log(diff(range(xy[, "x"])))){
        ldxs <- log(max(xy[, "x"], na.rm = TRUE) * 0.25)   
      }
    }
  }
  
  value <- c(a, b, xs, ldxs)
  names(value) <- mCall[c("a", "b", "xs", "ldxs")]
  value
}

#' @rdname SSquadpq
#' @return quadpq: vector of the same length as x using the quadratic-plateau-quadratic function
#' @export
quadpq <- function(x, a, b, xs, ldxs){
  
  dxs <- exp(ldxs)
  .value <- (x <= xs) * (a + b * x + (-0.5 * b/xs) * x^2) + 
    (x > xs & x <= (xs + dxs)) * (a + (-b^2)/(4 * -0.5 * b/xs)) + 
    (x > (xs + dxs)) * ((a + (-b^2)/(4 * -0.5 * b/xs)) + ((-0.5 * b/xs) * (x - (xs + dxs))^2))
 
  ## Derivative with respect to a
  ## The derivative is always 1, regardless of the phase?
  .exp1 <- 1
  
  ## Derivative with respect to b
  ## .exp2 <- deriv(~ a + b * x + (-0.5 * b/xs) * x^2, "b")
  ## .exp2b <- deriv(~ a + (-b^2)/(4 * -0.5 * b/xs), "b")
  ## .exp2c <- deriv(~ (a + (-b^2)/(4 * -0.5 * b/xs)) + ((-0.5 * b/xs) * (x - (xs + dxs))^2), "b")
  
  .expr2 <- -b^2
  .expr4 <- 4 * -0.5
  .expr6 <- .expr4 * b/xs
  .expr61 <- (x - (xs + dxs))^2
  .expr62 <- -(0.5/xs * .expr61 + (2 * b/.expr61 + .expr2 * (.expr4/xs)/.expr61^2)) 
  .exp2 <- ifelse(x < xs, x - 0.5/xs * x^2,
                  ifelse(x > (x + dxs), .expr62, -(2 * b/.expr6 + .expr2 * (.expr4/xs)/.expr6^2)))
  ## print(.exp2)
  ## Derivative with respect to xs
  ## .exp2 <- deriv(~ (a + b * x + (-0.5 * b/xs) * x^2), "xs")
  ## .exp2b <- deriv(~ a + (-b^2)/(4 * -0.5 * b/xs), "xs")
  .expr4 <- -0.5 * b
  .expr6 <- x^2
  .expr2 <- -b^2
  .expr5 <- 4 * -0.5 * b
  .expr7 <- .expr5/xs
  .expr9 <- .expr4
  .expr10 <- .expr9/xs
  .expr12 <- x - (xs + dxs)
  .expr13 <- .expr12^2
  .expr16 <- xs^2
  .expr20 <- .expr2 * (.expr5/.expr16)/.expr6^2 - (.expr10 * (2 * .expr12) + .expr9/.expr16 * .expr13)
  .exp3 <- ifelse(x <= xs, -(.expr4/xs^2 * .expr6), 
                  ifelse(x > (x + dxs), .expr20, .expr2 * (.expr5/xs^2)/.expr7^2))

  ## print(.exp3)
  ## .exp2d <- deriv(~ (a + (-b^2)/(4 * -0.5 * b/xs)) + ((-0.5 * b/xs) * (x - (xs + dxs))^2), "dxs")
  .expr3 <- -0.5
  .expr10 <- .expr3 * b/xs
  .expr12 <- x - (xs + dxs)
  .exp4 <- ifelse(x < (xs + dxs), 0, -(.expr10 * (2 * .expr12)))
  ## print(.exp4)
  
  .actualArgs <- as.list(match.call()[c("a", "b", "xs", "ldxs")])
  
  if(FALSE){
    cat("Gradients... a:", round(.exp1, 3), 
        "_ b:", round(.exp2, 3),
        "_ xs:", round(.exp3, 3),
        "_ ldxs:", round(.exp4, 3), "\n")
  }
  
  ## Try to only return a gradient if the values are good?
  ## Something wrong with the gradient, need to fix it
  ##  Gradient
  # if (all(unlist(lapply(.actualArgs, is.name)))) {
  #   .grad <- array(0, c(length(.value), 4L), list(NULL, c("a","b","xs","ldxs")))
  #   .grad[, "a"] <- .exp1
  #   .grad[, "b"] <- .exp2
  #   .grad[, "xs"] <- .exp3
  #   .grad[, "ldxs"] <- .exp4
  #   dimnames(.grad) <- list(NULL, .actualArgs)
  #   attr(.value, "gradient") <- .grad
  # }
  .value
}

#' @rdname SSquadpq
#' @export
SSquadpq <- selfStart(quadpq, initial = quadpqInit, c("a", "b", "xs", "ldxs"))

