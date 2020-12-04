#'
#' @title Summarize a matrix of simulations by their mean (median), sd (mad), and quantiles
#' @name summary_simulate
#' @description Utility function to summarize the output from \sQuote{simulate} 
#' functions in this package
#' @param object nobs x nsim matrix where nobs are the number of observations in the
#' dataset and nsim are the number of simulations
#' @param probs the percentiles to be computed by the quantile function
#' @param robust 	If FALSE (the default) the mean is used as the measure of central tendency 
#' and the standard deviation as the measure of variability. 
#' If TRUE, the median and the median absolute deviation (MAD) are applied instead. 
#' @param ... additional arguments to be passed. (none used at the moment)
#' @export
#' @examples 
#' \donttest{
#' data(barley, package = "nlraa")
#' fit <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
#' sim <- simulate_nls(fit, nsim = 100)
#' sims <- summary_simulate(sim)
#' }
#'

summary_simulate <- function(object, probs = c(0.025, 0.975), robust = FALSE, ...){
  
  if(!inherits(object, "matrix")) stop("'object' should be a matrix")
  
  mat <- matrix(nrow = nrow(object), ncol = 4)
  
  lwr <- probs[1]
  upr <- probs[2]
  lwr.lbl <- paste0("Q", lwr * 1e2)
  upr.lbl <- paste0("Q", upr * 1e2)
  colnames(mat) <- c("Estimate", "Est.Error", lwr.lbl, upr.lbl)
  
  mat[,1] <- apply(object, 1, mean)
  mat[,2] <- apply(object, 1, stats::sd)
  mat[,3] <- apply(object, 1, stats::quantile, probs = lwr)
  mat[,4] <- apply(object, 1, stats::quantile, probs = upr)

  if(robust){
    mat[,1] <- apply(object, 1, stats::median)
    mat[,2] <- apply(object, 1, stats::mad)
  } 
  
  return(mat)
}