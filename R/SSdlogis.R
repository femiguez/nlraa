#' ### This equation is not ready for production
#' 
#' #' For details see \sQuote{Double logistic curve in regression modeling}
#' #' 
#' #' @title self start for a double logistic function
#' #' @name SSdlogis
#' #' @rdname SSdlogis
#' #' @description Self starter for Hill function with parameters Ka, n and a
#' #' @param x input vector (x). Concentration of substrate in the original Hill model.
#' #' @param Ka parameter representing the concentration at which half of maximum y is attained
#' #' @param n parameter which controls the curvature 
#' #' @param a parameter which controls the maximum value of the response (asymptote)
#' #' @details The form of the equations are: \cr
#' #' hill1: \deqn{1 / (1 + (Ka/x))}. \cr 
#' #' hill2: \deqn{1 / (1 + (Ka/x)^n)}. \cr
#' #' hill3: \deqn{a / (1 + (Ka/x)^n)}. \cr
#' #' @note Zero values are not allowed.
#' #' @examples 
#' #' \donttest{
#' #' require(ggplot2)
#' #' ## Example for hill1
#' #' set.seed(1234)
#' #' x <- 1:20
#' #' y <- hill1(x, 10) + rnorm(20, sd = 0.03)
#' #' dat1 <- data.frame(x = x, y = y)
#' #' fit1 <- nls(y ~ SShill1(x, Ka), data = dat1)
#' #' 
#' #' ## Example for hill2
#' #' y <- hill2(x, 10, 1.5) + rnorm(20, sd = 0.03)
#' #' dat2 <- data.frame(x = x, y = y)
#' #' fit2 <- nls(y ~ SShill2(x, Ka, n), data = dat2)
#' #' 
#' #' ## Example for hill3
#' #' y <- hill3(x, 10, 1.5, 5) + rnorm(20, sd = 0.03)
#' #' dat3 <- data.frame(x = x, y = y)
#' #' fit3 <- nls(y ~ SShill3(x, Ka, n, a), data = dat3)
#' #' 
#' #' ggplot(data = dat3, aes(x, y)) + 
#' #'   geom_point() + 
#' #'   geom_line(aes(y = fitted(fit3)))
#' #' }
#' NULL
#' 
#' dlogis7Init <- function(mCall, LHS, data, ...){
#'   
#'   xy <- sortedXyData(mCall[["x"]], LHS, data)
#'   if(nrow(xy) < 7){
#'     stop("Too few distinct input values to fit a double logistic 7.")
#'   }
#' 
#'   objfun <- function(cfs){
#'     pred <- dlogis7(xy[,"x"], y.min=cfs[1], y.mid=cfs[2], y.max=cfs[3], xmid1=cfs[4], scal1=cfs[5], xmid2=cfs[6], scal2=cfs[7])
#'     ans <- sum((xy[,"y"] - pred)^2)
#'     ans
#'   }
#'   
#'   res <- list()
#'   j <- 1
#'   xmid <- max(xy[,"x"], na.rm = TRUE)/2
#'   for(i in 1:5){
#'     cfs <- c(min(xy[,"y"], na.rm = TRUE), ## Guess for y.min
#'              max(xy[,"y"], na.rm = TRUE)/2, ## Guess for y.mid
#'              max(xy[,"y"], na.rm = TRUE), ## Guess for y.max
#'              xmid/i, ## Guess for xmid1
#'              (xmid/4)/i, ## Guess for scal1
#'              (xmid + xmid/i), ## Guess for xmid2
#'              (xmid/4)/i) ## Guess for scal2
#'     op <- try(stats::optim(cfs, objfun), silent = TRUE)
#'     if(!inherits(op, "try-error")){
#'       res[[j]] <- op
#'       j <- j + 1
#'     }
#'   }
#'   
#'   if(all(sapply(res, \(x) inherits(x, 'try-error')))){
#'     ## If all fail
#'     y.min <- min(xy[,"y"], na.rm = TRUE) ## Guess for y.min
#'     y.mid <- max(xy[,"y"], na.rm = TRUE)/2 ## Guess for y.mid
#'     y.max <- max(xy[,"y"], na.rm = TRUE) ## Guess for y.max
#'     xmid1 <- xmid / 2 ## Guess for xmid1
#'     scal1 <- xmid / 4 ## Guess for scal1
#'     xmid2 <- 2 * xmid  ## Guess for xmid2
#'     scal2 <- xmid / 4 ## Guess for scal2
#'   }else{
#'     vals <- sapply(res, \(x) x$value)
#'     if(all(is.na(vals))){
#'       y.min <- min(xy[,"y"], na.rm = TRUE) ## Guess for y.min
#'       y.mid <- max(xy[,"y"], na.rm = TRUE)/2 ## Guess for y.mid
#'       y.max <- max(xy[,"y"], na.rm = TRUE) ## Guess for y.max
#'       xmid1 <- xmid / 2 ## Guess for xmid1
#'       scal1 <- xmid / 4 ## Guess for scal1
#'       xmid2 <- 2 * xmid  ## Guess for xmid2
#'       scal2 <- xmid / 4 ## Guess for scal2
#'     }else{
#'       wminvals <- which.min(vals)
#'       op <- res[[wminvals]]      
#'       y.min <- op$par[1]
#'       y.mid <- op$par[2]
#'       y.max <- op$par[3]
#'       xmid1 <- op$par[4] 
#'       scal1 <- op$par[5] 
#'       xmid2 <- op$par[6]
#'       scal2 <- op$par[7] 
#'     }
#'   }
#'   
#'   # y.min <- min(xy[,"y"], na.rm = TRUE) ## Guess for y.min
#'   # y.mid <- max(xy[,"y"], na.rm = TRUE)/2 ## Guess for y.mid
#'   # y.max <- max(xy[,"y"], na.rm = TRUE) ## Guess for y.max
#'   # xmid1 <- xmid / 2 ## Guess for xmid1
#'   # scal1 <- xmid / 4 ## Guess for scal1
#'   # xmid2 <- 2 * xmid  ## Guess for xmid2
#'   # scal2 <- xmid / 4 ## Guess for scal2
#'   
#'   value <- c(y.min, y.mid, y.max, xmid1, scal1, xmid2, scal2)
#'   names(value) <- mCall[c("y.min", "y.mid", "y.max", "xmid1", "scal1", "xmid2", "scal2")]
#'   value
#' }
#' 
#' 
#' #' @rdname SSdlogis
#' #' @return dlogis7: vector of the same length as x using the quadratic-plateau-quadratic function
#' #' @export
#' dlogis7 <- function(x, y.min, y.mid, y.max, xmid1, scal1, xmid2, scal2){
#'   
#'   .value <- y.min + (y.mid - y.min)/(1 + exp((xmid1 - x)/scal1)) + (y.max - y.mid)/(1 + exp((xmid2 - x)/scal2))
#'   
#'   .actualArgs <- as.list(match.call()[c("y.min", "y.mid", "y.max", "xmid1", "scal1", "xmid2", "scal2")])
#' 
#'   .value
#' }
#' 
#' #' @rdname SSdlogis
#' #' @export
#' SSdlogis7 <- selfStart(dlogis7, initial = dlogis7Init, c("y.min", "y.mid", "y.max", "xmid1", "scal1", "xmid2", "scal2"))
#' 
#' 
#' #### Simpler double logistic ----
#' 
#' dlogis4Init <- function(mCall, LHS, data, ...){
#'   
#'   xy <- sortedXyData(mCall[["x"]], LHS, data)
#'   if(nrow(xy) < 5){
#'     stop("Too few distinct input values to fit a double logistic 4.")
#'   }
#'   
#'   objfun <- function(cfs){
#'     pred <- dlogis4(xy[,"x"], yt=cfs[1], yd=cfs[2], xmid=cfs[3], scal=cfs[4])
#'     ans <- sum((xy[,"y"] - pred)^2)
#'     ans
#'   }
#'   
#'   res <- list()
#'   j <- 1
#'   ymin <- min(xy[,"y"], na.rm = TRUE) 
#'   ymax <- max(xy[,"y"], na.rm = TRUE)
#'   for(i in 1:5){
#'     cfs <- c((min(xy[,"y"], na.rm = TRUE) + max(xy[,"y"], na.rm = TRUE))/2, ## Guess for yt
#'              (max(xy[,"y"], na.rm = TRUE) - min(xy[,"y"], na.rm = TRUE))/2, ## Guess for yd
#'              max(xy[,"x"], na.rm = TRUE)/i, ## Guess for xmid
#'              (max(xy[,"x"], na.rm = TRUE)/4)/i) ## Guess for scal
#'     op <- try(stats::optim(cfs, objfun), silent = TRUE)
#'     if(!inherits(op, "try-error")){
#'       res[[j]] <- op
#'       j <- j + 1
#'     }
#'   }
#'   
#'   if(all(sapply(res, \(x) inherits(x, 'try-error')))){
#'     ## If all fail
#'     yt <- (min(xy[,"y"], na.rm = TRUE) + max(xy[,"y"], na.rm = TRUE))/2 ## Guess for yt
#'     yd <- (max(xy[,"y"], na.rm = TRUE) - min(xy[,"y"], na.rm = TRUE))/2 ## Guess for yd
#'     xmid <- max(xy[,"x"], na.rm = TRUE)/2 ## Guess for xmid
#'     scal <-  max(xy[,"x"], na.rm = TRUE)/4 ## Guess for scal
#'   }else{
#'     vals <- sapply(res, \(x) x$value)
#'     if(all(is.na(vals))){
#'       yt <- (min(xy[,"y"], na.rm = TRUE) + max(xy[,"y"], na.rm = TRUE))/2 ## Guess for yt
#'       yd <- (max(xy[,"y"], na.rm = TRUE) - min(xy[,"y"], na.rm = TRUE))/2 ## Guess for yd
#'       xmid <- max(xy[,"x"], na.rm = TRUE)/2 ## Guess for xmid
#'       scal <-  max(xy[,"x"], na.rm = TRUE)/4 ## Guess for scal
#'     }else{
#'       wminvals <- which.min(vals)
#'       op <- res[[wminvals]]      
#'       yt <- op$par[1]
#'       yd <- op$par[2]
#'       xmid <- op$par[3] 
#'       scal <- op$par[4] 
#'     }
#'   }
#'   
#'   value <- c(yt, yd, xmid, scal)
#'   names(value) <- mCall[c("yt", "yd", "xmid", "scal")]
#'   value
#' }
#' 
#' 
#' #' @rdname SSdlogis
#' #' @return dlogis4: vector of the same length as x using the quadratic-plateau-quadratic function
#' #' @export
#' dlogis4 <- function(x, yt, yd, xmid, scal){
#'   
#'   .value <- yt + 2 * yd * ((1 / (1 + exp(-1 * abs((x - xmid)/scal)))) - 0.5) * sign((x - xmid)/scal)
#' 
#'   .actualArgs <- as.list(match.call()[c("yt", "yd", "xmid", "scal")])
#'   
#'   .value
#' }
#' 
#' #' @rdname SSdlogis
#' #' @export
#' SSdlogis4 <- selfStart(dlogis4, initial = dlogis4Init, c("yt", "yd", "xmid", "scal"))
