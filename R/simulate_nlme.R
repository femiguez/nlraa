#' Simulate samples from a nonlinear model
#' 
#' @title Simulate samples from a nonlinear mixed model from fixed effects
#' @name simulate_nlme
#' @param object object of class \sQuote{gnls} or \sQuote{nlme}
#' @param f function to be used for the simulation
#' @param nsim number of samples, default 1
#' @param ... additional arguments (none used at the moment)
#' @details details to come
#' @return It returns a matrix with simulated values from the original object
#' with number of rows equal to the number of rows of \sQuote{fitted} and number
#' of columns equal to the number of simulated samples (\sQuote{nsim}).
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
#' ## This will take a few seconds                               
#' fit.lp.gnls.Bt3 <- boot_nlme(fit.lp.gnls3) 
#' confint(fit.lp.gnls.Bt3, type = "perc")
#' }
#' 

simulate_nlme <- function(object, 
                          f = NULL, 
                          R = 1, 
                          psim = 1, ...){
  
  ## Error checking
  if(!inherits(object, c("gnls","nlme"))) stop("object should be of class 'gnls' or 'nlme'")
  
  sim.mat <- matrix(ncol = R, nrow = length(fitted(object)))
  
  ## First example for the gnls case
  for(i in seq_len(R)){
      if(inherits(object, "gnls")){
        sim.mat[,i] <- as.vector(simulate_gnls(object, psim = psim, ...))
      }
      if(inherits(x, "nlme")){
        sim.mat[,i] <- as.vector(simulate_nlme_one(object, psim = psim, ...))
      }
  }

  return(sim.mat)
}

