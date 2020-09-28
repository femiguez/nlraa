#' Simulate multiple samples from a nonlinear model
#' 
#' @title Simulate samples from a nonlinear mixed model from fixed effects
#' @name simulate_nlme
#' @param object object of class \code{\link[nlme]{gnls}} or \code{\link[nlme]{nlme}}
#' @param nsim number of samples, default 1
#' @param psim simulation level for vector of fixed parameters for \code{\link{simulate_nlme_one}}
#' @param value whether to return a matrix (default) or an augmented data frame
#' @param ... additional arguments to be passed to either \code{\link{simulate_gnls}} or \code{\link{simulate_nlme_one}}
#' @details The details can be found in either \code{\link{simulate_gnls}} or \code{\link{simulate_nlme_one}}.
#' This function is very simple and it only sets up a matrix and a loop in order to simulate several instances of 
#' model outputs.
#' @return It returns a matrix with simulated values from the original object
#' with number of rows equal to the number of rows of \code{\link{fitted}} and number
#' of columns equal to the number of simulated samples (\sQuote{nsim}). In the case of 'data.frame'
#' it returns an augmented data.frame, which can potentially be a very large object, but which
#' makes furhter plotting more convenient.  
#' @export
#' @examples 
#' \donttest{
#' require(car)
#' require(nlme)
#' data(barley)
#' barley2 <- subset(barley, year < 1974)
#' fit.lp.gnls2 <- gnls(yield ~ SSlinp(NF, a, b, xs), data = barley2)
#' barley2$year.f <- as.factor(barley2$year)
#' cfs <- coef(fit.lp.gnls2)
#' fit.lp.gnls3 <- update(fit.lp.gnls2, 
#'                       params = list(a + b + xs ~ year.f),
#'                       start = c(cfs[1], 0, 0, 0, 
#'                                 cfs[2], 0, 0, 0,
#'                                 cfs[3], 0, 0, 0))
#'                                 
#' sims <- simulate_nlme(fit.lp.gnls3, nsim = 3)
#' }
#' 

simulate_nlme <- function(object, 
                          nsim = 1, 
                          psim = 1, 
                          value = c("matrix", "data.frame"),...){
  
  ## Error checking
  if(!inherits(object, c("gnls","nlme"))) stop("object should be of class 'gnls' or 'nlme'")
  
  value <- match.arg(value)
  
  sim.mat <- matrix(ncol = nsim, nrow = length(fitted(object)))
  
  ## First example for the gnls case
  for(i in seq_len(nsim)){
      if(inherits(object, "gnls")){
        sim.mat[,i] <- as.vector(simulate_gnls(object, psim = psim, ...))
      }
      if(inherits(object, "nlme")){
        sim.mat[,i] <- as.vector(simulate_nlme_one(object, psim = psim, ...))
      }
  }

  if(value == "matrix"){
    colnames(sim.mat) <- paste0("sim_",1:nsim)
    return(sim.mat)  
  }else{
    dat <- eval(object$call$data)
    adat <- data.frame(ii = as.factor(rep(1:nsim, each = nrow(dat))),
                      dat,
                      sim.y = c(sim.mat),
                      row.names = 1:c(nsim * nrow(dat)))   
    return(adat)
  }
}

