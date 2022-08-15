#' This function returns several indexes that might be useful for interpretation
#' 
#' For objects of class \sQuote{lm} and \sQuote{nls} \cr
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
#' 
#'
#'  \url{https://en.wikipedia.org/wiki/Coefficient_of_determination} \cr
#'  \url{https://en.wikipedia.org/wiki/Nash-Sutcliffe_model_efficiency_coefficient} \cr
#'  \url{https://en.wikipedia.org/wiki/Concordance_correlation_coefficient}
#'  
#' @title Indexes of Agreement Table
#' @name IA_tab
#' @rdname IA_tab
#' @description Indexes of agreement
#' @param obs vector with observed data
#' @param sim vector with simulated data (should be the same length as observed)
#' @param object alternative to the previous two arguments. An object of class \sQuote{lm}, \sQuote{nls}, \sQuote{lme} or \sQuote{nlme}
#' @param null.object optional object which represents the \sQuote{null} model. It is an intercept-only model
#' by default. (Not used at the moment).
#' @seealso \code{\link{IC_tab}}
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
#' ias0
#' ## Nonlinear model
#' fit1 <- nls(lfgr ~ SSblin(ftsw, a, b, xs, c), data = swpg)
#' ias1 <- IA_tab(object = fit1)
#' ias1
#' plot(ias1)
#' ## Linear Mixed Models
#' data(barley, package = "nlraa")
#' fit2 <- lme(yield ~ NF + I(NF^2), random = ~ 1 | year, data = barley)
#' ias2 <- IA_tab(object = fit2)
#' ias2
#' ## Nonlinear Mixed Model
#' barleyG <- groupedData(yield ~ NF | year, data = barley)
#' fit3L <- nlsLMList(yield ~ SSquadp3(NF, a, b, c), data = barleyG)
#' fit3 <- nlme(fit3L, random = pdDiag(a + b ~ 1))
#' ias3 <- IA_tab(object = fit3)
#' ias3
#' plot(ias3)
#' ## Plotting model
#' prds <- predict_nlme(fit3, interval = "conf", plevel = 0)
#' barleyGA <- cbind(barleyG, prds)
#' ggplot(data = barleyGA, aes(x = NF, y = yield)) + 
#'    geom_point() + 
#'    geom_line(aes(y = Estimate)) + 
#'    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), 
#'                fill = "purple", alpha = 0.2)
#' ## R2M for model 2
#' R2M(fit2)
#' ## R2M for model 3
#' R2M(fit3)
#' 
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
    
    if(inherits(object, c("lmerMod", "merMod"))){
      rsp.nm <- gsub("\\s+", "", strsplit(deparse(formula(object)), "~")[[1]][1])
      obs <- getData(object)[[rsp.nm]]
    }
      
  }
  ## Bias
  bias <- mean(obs - sim)
  ## correlation and correlation squared
  corr <- cor(obs, sim)
  R2.1 <- cor(obs, sim)^2
  ## Two different types of R-squared
  if(inherits(object, c("nls", "lm", "gls", "gnls"))){
    R2.2 <- nlraa::R2M(object)
  }
  
  if(inherits(object, c("lme", "nlme"))){
    R2 <- nlraa::R2M(object)
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
  
  if(missing(object) || inherits(object, c("lm", "nls", "gls", "gnls"))){
    ia_tab <- data.frame(bias = bias, 
                         intercept = intercept,
                         slope = slope,
                         RSS = RSS, MSE = MSE, RMSE = RMSE,
                         R2.1 = R2.1, R2.2 = R2.2$R2, 
                         ME = mod.eff, NME = norm.mod.eff, Corr = corr,
                         ConCorr = concor)
  }else{
    if(inherits(object, c("lmerMod", "merMod"))){
      warning("R2 for 'lmerMod' or 'merMod' objects are not available.
              There are some other packages which can compute this.
              For example, MuMIn.")
      ia_tab <- data.frame(bias = bias, 
                           intercept = intercept,
                           slope = slope,
                           RSS = RSS, MSE = MSE, RMSE = RMSE,
                           R2 = NA, 
                           R2.marginal = NA,
                           R2.conditional = NA,
                           ME = mod.eff, NME = norm.mod.eff, Corr = corr,
                           ConCorr = concor) 
    }else{
      ia_tab <- data.frame(bias = bias, 
                           intercept = intercept,
                           slope = slope,
                           RSS = RSS, MSE = MSE, RMSE = RMSE,
                           R2 = R2.1, 
                           R2.marginal = R2$R2.marginal,
                           R2.conditional = R2$R2.conditional,
                           ME = mod.eff, NME = norm.mod.eff, Corr = corr,
                           ConCorr = concor)      
    }
  }

  lst <- list(IA_tab = ia_tab, obs = obs, sim = sim)
  ans <- structure(lst, class = "IA_tab")
  ans
}


#' @rdname IA_tab
#' @description printing function for IA_tab
#' @param x object of class \sQuote{IA_tab}.
#' @param ... additional plotting arguments (none use at the moment).
#' @param digits number of digits for rounding (default is 2)
#' @export
#' 
print.IA_tab <- function(x, ..., digits = 2){
  return(print(x$IA_tab, digits = digits))
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

#' I have read some papers about computing an R-squared for (non)linear (mixed) models
#' and I am not sure that it makes sense at all. However, here they are and
#' the method is general enough that it extends to nonlinear mixed models. What do
#' these numbers mean and why would you want to compute them are good questions to 
#' ponder... \cr
#' 
#' Recommended reading: \cr
#' Nakagawa and Schielzeth Methods in Ecology and Evolution \doi{10.1111/j.2041-210x.2012.00261.x} \cr
#' 
#' \url{https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/} \cr
#' 
#' Spiess, AN., Neumeyer, N. An evaluation of R2 as an inadequate measure for nonlinear models in 
#' pharmacological and biochemical research: a Monte Carlo approach. BMC Pharmacol 10, 6 (2010). 
#' \doi{10.1186/1471-2210-10-6} \cr
#' 
#' \url{https://stat.ethz.ch/pipermail/r-sig-mixed-models/2010q1/003363.html} \cr
#' 
#' \url{https://blog.minitab.com/en/adventures-in-statistics-2/why-is-there-no-r-squared-for-nonlinear-regression} \cr
#' 
#' \url{https://stats.stackexchange.com/questions/111150/calculating-r2-in-mixed-models-using-nakagawa-schielzeths-2013-r2glmm-me/225334#225334} \cr
#'  
#' Other R pacakges which calculate some version of an R-squared: performance, rcompanion, MuMIn
#' 
#' @title R-squared for nonlinear mixed models
#' @name R2M
#' @rdname R2M
#' @description R-squared \sQuote{modified} for nonlinear (mixed) models
#' @param x object of class \sQuote{lm}, \sQuote{nls}, \sQuote{gls}, \sQuote{gnls}, 
#' \sQuote{lme} or \sQuote{nlme} .
#' @param ... additional arguments (none use at the moment).
#' @note The references here strongly discourage the use of R-squared in anything
#' but linear models.
#' @seealso \code{\link{IA_tab}}
#' @return it returns a list with the following structure: \cr
#' for an object of class \sQuote{lm}, \sQuote{nls}, \sQuote{gls} or \sQuote{gnls}, \cr
#' R2: R-squared \cr
#' var.comp: variance components with var.fixed and var.resid \cr
#' var.perc: variance components percents (should add up to 100) \cr
#' for an object of class \sQuote{lme} or \sQuote{nlme} in addition it also returns: \cr
#' R2.marginal: marginal R2 which only includes the fixed effects \cr
#' R2.conditional: conditional R2 which includes both the fixed and random effects \cr
#' var.random: the variance contribution of the random effects
#' @export 
#' @examples 
#' \donttest{
#' require(nlme)
#' data(barley, package = "nlraa")
#' fit2 <- lme(yield ~ NF + I(NF^2), random = ~ 1 | year, data = barley)
#' R2M(fit2)
#' ## Nonlinear Mixed Model
#' barleyG <- groupedData(yield ~ NF | year, data = barley)
#' fit3L <- nlsLMList(yield ~ SSquadp3(NF, a, b, c), data = barleyG)
#' fit3 <- nlme(fit3L, random = pdDiag(a + b ~ 1))
#' R2M(fit3)
#' }

R2M <- function(x, ...){
  UseMethod("R2M")
}

#' @rdname R2M
#' @export
R2M.nls <- function(x, ...){
  
  var.fixed <- stats::var(fitted(x))
  var.resid <- stats::var(residuals(x))
  
  var.comp <- c(var.fixed = var.fixed, var.resid = var.resid)
  var.perc <- round((var.comp / sum(var.comp)) * 100, 1)
  
  ans <- list(R2 = var.fixed / (var.fixed + var.resid), 
              var.comp = var.comp,
              var.perc = var.perc)
  ans
}

#' @rdname R2M
#' @export
R2M.lm <- R2M.nls

#' @rdname R2M
#' @export
R2M.gls <- R2M.nls

#' @rdname R2M
#' @export
R2M.gnls <- R2M.nls

#' @rdname R2M
#' @export
R2M.lme <- function(x, ...){
  
  ## Marginal
  Q <- x$dims$Q

  prd0 <- predict(x, level = 0) ## This only includes the fixed effects
  prdQ <- predict(x, level = Q) ## This is fixed effects plus random effects
  ## Computing fixed, random and residual
  var.fixed <- stats::var(prd0)
  var.random <- stats::var(prd0 - prdQ)
  var.resid <- stats::var(residuals(x))

  var.comp <- c(var.fixed = var.fixed,
                var.random = var.random, 
                var.resid = var.resid)
  
  var.perc <- round((var.comp / sum(var.comp)) * 100, 1)  
  
  R2.marginal <- var.fixed / (var.fixed + var.random + var.resid)
  R2.conditional <- (var.fixed + var.random) / (var.fixed + var.random + var.resid)
  
  ans <- list(R2.marginal = R2.marginal,
              R2.conditional = R2.conditional,
              var.comp = var.comp,
              var.perc = var.perc)
  ans  
  
}

#' @rdname R2M
#' @export
R2M.nlme <- R2M.lme




