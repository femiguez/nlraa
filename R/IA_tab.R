#' This function returns several indexes that might be useful for interpretation
#' 
#' For obbjects of class \sQuote{lm} or \sQuote{nls} \cr
#' bias: mean(obs - sim) \cr
#' intercept: intercept of the model obs ~ beta_0 + beta_1 * sim + error \cr
#' slope: slope of the model obs ~ beta_0 + beta_1 * sim + error \cr
#' RSS (deviance): residual sum of squares of the previous model \cr
#' MSE (RSS / n): mean squared error; where n is the number of observations \cr
#' RMSE: squared root of the previous index \cr
#' R2.1: R-squared extracted from an \sQuote{lm} object \cr
#' R2.2: R-squared computed as the correlation between observed and simulated to the power of 2. \cr
#' ME: model efficiency \cr
#' NME: Normalized model efficiency \cr
#' Corr: correlation between observed and simulated \cr
#' ConCorr: concordance correlation  \cr
#' 
#' For objects of class \sQuote{gls}, \sQuote{gnls}, \sQuote{lme} or \sQuote{nlme} there
#' are additional metrics such as:
#' Pseudo_R2: See reference by Nakagawa and Schielzeth
#' Normalized_Pseudo_R2: See reference by Nakagawa and Schielzeth
#' NS_R2_marginal: See reference by Nakagawa and Schielzeth
#' NS_R2_conditional: See reference by Nakagawa and Schielzeth
#' 
#'
#'  https://en.wikipedia.org/wiki/Coefficient_of_determination \cr
#'  https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient \cr
#'  https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
#'  
#' @title Indexes of Agreement Table
#' @name IA_tab
#' @rdname IA_tab
#' @description Indexes of agreement
#' @param obs vector with observed data
#' @param sim vector with simulated data (should be the same length as observed)
#' @param object alternative to the previous two arguments. An object of class \sQuote{lm}, \sQuote{nls} or \sQuote{lme}
#' @param null.object optional object which represents the \sQuote{null} model. It is an intercept-only model
#' by default.
#' @seealso \code{\link{IC_tab}}
#' @references Nakagawa and Schielzeth Methods in Ecology and Evolution (doi: 10.1111/j.2041-210x.2012.00261.x)
#' @note Not sure that all the formulas are correct for Mixed Models yet.
#' @export
#' @examples 
#' \donttest{
#' require(nlme)
#' require(ggplot2)
#' ## Fit a simple model and then compute IAs
#' data(swpg)
#' #' ## Linear model
#' fit0 <- lm(lfgr ~ ftsw + I(ftsw^2), data = swpg)
#' ias0 <- IA_tab(object = fit0)
#' ias0$IA_tab
#' ## Nonlinear model
#' fit1 <- nls(lfgr ~ SSblin(ftsw, a, b, xs, c), data = swpg)
#' ias1 <- IA_tab(object = fit1)
#' ias1$IA_tab
#' plot(ias1)
#' ## Linear Mixed Models
#' data(barley, package = "nlraa")
#' fit2 <- lme(yield ~ NF + I(NF^2), random = ~ 1 | year, data = barley)
#' ias2 <- IA_tab(object = fit2)
#' ias2$IA_tab
#' ## The warning is becase no null model was provided
#' ## This means the model is compared to an intercept-only model
#' ## Nonlinear Mixed Model
#' barleyG <- groupedData(yield ~ NF | year, data = barley)
#' fit3L <- nlsLMList(yield ~ SSquadp3(NF, a, b, c), data = barleyG)
#' fit3 <- nlme(fit3L, random = pdDiag(a + b ~ 1))
#' ias3 <- IA_tab(object = fit3)
#' ias3$IA_tab
#' plot(ias3)
#' ## Plotting model
#' prds <- predict_nlme(fit3, interval = "conf", plevel = 0)
#' barleyGA <- cbind(barleyG, prds)
#' ggplot(data = barleyGA, aes(x = NF, y = yield)) + 
#'    geom_point() + 
#'    geom_line(aes(y = Estimate)) + 
#'    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), 
#'                fill = "purple", alpha = 0.2)
#' }
#' 

IA_tab <- function(obs, sim, object, null.object){
  
  if(missing(object)){
    if(length(obs) != length(sim))
      stop("obs and sim vectors should be of equal length")
  }
  
  if(!missing(obs) && !inherits(obs, "numeric"))
    stop("object obs should be of class 'numeric'")
  
  if(!missing(sim) && !inherits(sim, "numeric"))
    stop("object sim should be of class 'numeric'")
  
  if(!missing(object)){
    sim <- fitted(object)
    ## How do I retrieve data from different R objects?
    if(inherits(object, c("lm", "nls"))){
      resp.var <- all.vars(formula(object))[1]
      if(is.null(object$call$data))
        stop("data argument should be used when fitting lm or nls models")
      obs <- eval(getCall(object)$data)[[resp.var]]
    }
    if(inherits(object, c("gls", "gnls", "lme", "nlme")))
      obs <- nlme::getResponse(object)
  }
  ## Bias
  bias <- mean(obs - sim)
  ## correlation and correlation squared
  corr <- cor(obs, sim)
  R2.1 <- cor(obs, sim)^2
  ## Two different types of R-squared
  if(inherits(object, c("nls", "lm"))){
    R2.2 <- 1 - deviance(object)/deviance(lm(obs ~ 1))
  }
  
  if(inherits(object, c("gls", "gnls", "lme", "nlme"))){
    if(missing(null.object)){
      warning("No null model was provided. Assuming an intercept only model.")
      ll.nm <- logLik(lm(obs ~ 1))[[1]]
    }else{
      ll.nm <- logLik(null.object)[[1]]
    }
    ## Pseudo-R-squared for mixed models?
    Pseudo_R2 <- 1 - exp((2/nobs(object)) * (ll.nm - logLik(object)[[1]]))
    Normalized_Pseudo_R2 <- Pseudo_R2 / (1 - exp(ll.nm) ^ (2 / nobs(object)))
    ## This is an atempt at computing the Nakagawa - Schielzeth R-squared
    ## doi: 10.1111/j.2041-210x.2012.00261.x
    NS.R2m.num <- var(fitted(object))
    NS.R2m.random <- 0
    if(inherits(object, c("lme", "nlme")))
      NS.R2m.random <- sum(diag(nlraa::var_cov(object, type = "random")))
    NS.R2m.den <- NS.R2m.num + NS.R2m.random + sigma(object)^2
    NS.R2.marginal <- NS.R2m.num / NS.R2m.den
    ## Conditional
    NS.R2c.num <- NS.R2m.num + NS.R2m.random
    NS.R2c.den <- NS.R2c.num + sigma(object)^2
    NS.R2.conditional <- NS.R2c.num / NS.R2c.den
  }
  ## Nash-Sutclife Model Efficiency
  ## https://en.wikipedia.org/wiki/Nash%E2%80%93Sutcliffe_model_efficiency_coefficient
  mod.eff <- 1 - (sum((obs - sim)^2))/sum((obs - mean(obs))^2) 
  norm.mod.eff <- 1 / (2 - mod.eff)
  
  ## Concordance correlation
  ## https://en.wikipedia.org/wiki/Concordance_correlation_coefficient
  concor.num <- 2 * cor(obs, sim) * sd(obs) * sd(sim)
  concor.den <- var(obs) + var(sim) + (mean(obs) - mean(sim))^2
  concor <- concor.num/concor.den
  
  ## Things to add:
  fit <- lm(obs ~ sim)
  intercept <- as.vector(coef(fit)[1])
  slope <- as.vector(coef(fit)[2])
  RSS <- deviance(fit)
  MSE <- RSS / length(obs)
  RMSE <- sqrt(MSE)
  
  if(missing(object) || inherits(object, c("lm", "nls"))){
    ia_tab <- data.frame(bias = bias, 
                         intercept = intercept,
                         slope = slope,
                         RSS = RSS, MSE = MSE, RMSE = RMSE,
                         R2.1 = R2.1, R2.2 = R2.2, 
                         ME = mod.eff, NME = mod.eff, Corr = corr,
                         ConCorr = concor)
  }else{
    ia_tab <- data.frame(bias = bias, 
                         intercept = intercept,
                         slope = slope,
                         RSS = RSS, MSE = MSE, RMSE = RMSE,
                         R2 = R2.1, 
                         Pseudo_R2 = Pseudo_R2,
                         Normalized_Pseudo_R2 = Normalized_Pseudo_R2,
                         NS_R2_marginal = NS.R2.marginal,
                         NS_R2_conditional = NS.R2.conditional,
                         ME = mod.eff, NME = mod.eff, Corr = corr,
                         ConCorr = concor)
  }

  lst <- list(IA_tab = ia_tab, obs = obs, sim = sim)
  ans <- structure(lst, class = "IA_tab")
  ans
}

#' @rdname IA_tab
#' @description plotting function for a IA_tab, it requires \sQuote{ggplot2}
#' @param x object of class \sQuote{IA_tab}.
#' @param y not used at the moment
#' @param ... additional plotting arguments (none use at the moment).
#' @param type either \dQuote{OvsS} (observed vs. simulated) or \dQuote{RvsS} (residuals vs. simulated). 
#' @export 

plot.IA_tab <- function(x, y, ..., type = c("OvsS", "RvsS")){
  
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    warning("ggplot2 is required for this plotting function")
    return(NULL)
  }
  
  type <- match.arg(type)
  
  if(type == "OvsS"){
    gp1 <- ggplot2::ggplot(mapping = ggplot2::aes(x = x$sim, y = x$obs)) + 
      ggplot2::geom_point() + 
      ggplot2::geom_smooth(method = "lm", se = FALSE) + 
      ggplot2::ylab("Observed") + ggplot2::xlab("Simulated")    
  }else{
    residu <- x$obs - x$sim
    gp1 <- ggplot2::ggplot(mapping = ggplot2::aes(x = x$sim, y = residu)) + 
      ggplot2::geom_point() + 
      ggplot2::geom_hline(yintercept = 0) + 
      ggplot2::ylab("Residual (Obs - Sim)") + ggplot2::xlab("Simulated")    
  }

  gp1
}

