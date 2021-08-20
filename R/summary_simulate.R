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
#' @param data the original data.frame used to fit the model. A data.frame will be
#' returned instead of a matrix in this case.
#' @param by opionally aggregate the results by some factor in the data.frame. It 
#' will be coarced to a formula. This should either be a character or a formula (starting with \sQuote{~}).
#' The aggregation follows the \sQuote{robust} argument above.
#' @param ... additional arguments to be passed. (none used at the moment)
#' @return By default it returns a matrix unless the \sQuote{data} argument is present and then
#' it will return a data.frame
#' @export
#' @examples 
#' \donttest{
#' data(barley, package = "nlraa")
#' fit <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
#' sim <- simulate_nls(fit, nsim = 100)
#' sims <- summary_simulate(sim)
#' 
#' ## If we want to combine the data.frame
#' simd <- summary_simulate(sim, data = barley)
#' ## If we also want to aggregate by nitrogen rate
#' simda <- summary_simulate(sim, data = barley, by = "NF")
#'
#' }
#'

summary_simulate <- function(object, probs = c(0.025, 0.975), robust = FALSE, 
                             data, by, ...){
  
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
  
  if(missing(data)){
    ans <- mat
  }else{
    
    if(nrow(data) != nrow(mat))
      stop("number of rows in the data argument should match the number of rows in the matrix")
    
    dat <- cbind(data, mat)
    
    if(missing(by)){
      ans <- dat
    }else{
      if(!robust){
        ## This new code allows me to use formula instead of a character
        ## Testing for a formula
        ## https://stackoverflow.com/questions/36361158/how-to-test-if-an-object-is-a-formula-in-base-r
        if(is.call(by) && by[[1]] == quote(`~`)){
          by <- as.character(by)[[2]]
        }
        agf0 <- paste0(c("cbind(", c("Estimate,", "Est.Error,", lwr.lbl, ",", upr.lbl, ")")), collapse = " ")
        agf <- as.formula(paste0(c(agf0, paste0("~", by)), collapse = " "))
        ans <- stats::aggregate(formula = agf, data = dat, FUN = mean)
      }else{
        ans <- stats::aggregate(formula = agf, data = dat, FUN = median)
      } 
    }
  } 
  
  return(ans)
}