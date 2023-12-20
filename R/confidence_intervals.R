#' Confidence interval methods for (non)-linear models 
#' 
#' @title Confidence interval methods for (non)-linear models 
#' @name confidence_intervals
#' @param x object of class \code{\link[stats]{lm}}, \code{\link[stats]{nls}}, \code{\link[nlme]{nlme}}, \code{\link[nlme]{gls}} or \code{\link[nlme]{gnls}}
#' @param method method or methods to use. Possible options are: \sQuote{Wald}, \sQuote{profile}, \sQuote{bootstrap}, \sQuote{all}
#' @param parm optional character string to select the parameter
#' @param level probability level
#' @param verbose logical (default FALSE) whether to print messages
#' @param ... additional arguments to be passed to function \code{\link[boot]{boot}}
#' @export
#' @examples 
#' \donttest{
#' require(car)
#' data(barley, package = "nlraa")
#' ## Fit a linear model (quadratic)
#' fit.lm <- lm(yield ~ NF + I(NF^2), data = barley)
#' 
#' cfs.int <- confidence_intervals(fit.lm, method = c("wald", "bootstrap"))
#' 
#' fit.nls <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
#' 
#' cfs.int2 <- confidence_intervals(fit.nls, method = c("wald", "profile", "bootstrap"))
#' 
#' }
#' 


confidence_intervals <- function(x, 
                                 method = c("wald", "profile", "bootstrap", "simple-bayes", "all"), 
                                 parm,
                                 level = 0.95,
                                 verbose = FALSE,
                                 ...){
  
  method <- match.arg(method, several.ok = TRUE)
  
  if(length(method) == 5) method <- "wald"
  
  ## Should test for objects which are not
  ## acceptable
  
  xargs <- list(...)
  
  ### Collect bootstrap arguments
  if(!is.null(xargs$R)){
    R <- xargs$R
  }else{
    R <- 999
  } 
  
  j <- 0
  k <- 0
  
  if(inherits(x, "lm") || inherits(x, "gls") || inherits(x, "gnls") || inherits(x, "lme")){
    
    len.cfs <- length(coef(x))
    if(any(method == "all")){
      if(inherits(x, "lm") || inherits(x, "gls") || inherits(x, "gnls") || inherits(x, "lme"))
        nrows <- 2 * len.cfs
      if(inherits(x, "nls"))
        nrows <- 3 * len.cfs
    }else{
      nrows <- length(method) * len.cfs
    }
    ans.dat <- data.frame(method = rep(NA, nrows), parm = NA, level = level,
                          lower = NA, estimate = NA, upper = NA)
    
    if(any(method == "profile"))
      stop("method = profile not available for objects of class 'lm', 'lme', 'gls' or 'gnls'")
    
    if(any(method == "wald") || any(method == "all")){
      if(missing(parm)){
        if(inherits(x, "lme")){
          ans <- nlme::intervals(x, which = "fixed")$fixed 
        }else{
          ans <- confint(x, level = level)  
        }
      }else{
        if(inherits(x, "lme")){
          ans <- nlme::intervals(x, which = "fixed")$fixed 
          which.parm <- which(row.names(ans) == parm)
          ans <- ans[which.parm, ,drop = FALSE]
        }else{
          ans <- confint(x, parm = parm, level = level)  
        }
      }

      if(inherits(x, "lme")){
        for(i in seq_len(dim(ans)[1])){
          ans.dat[i, "method"] <- "wald"
          ans.dat[i, "parm"] <- names(fixef(x))[i]
          ans.dat[i, "lower"] <- ans[i, 1]
          ans.dat[i, "upper"] <- ans[i, 3]
          ans.dat[i, "estimate"] <- fixef(x)[i]
        }
      }else{
        for(i in seq_len(dim(ans)[1])){
          ans.dat[i, "method"] <- "wald"
          ans.dat[i, "parm"] <- row.names(ans)[i]
          ans.dat[i, "lower"] <- ans[i, 1]
          ans.dat[i, "upper"] <- ans[i, 2]
          ans.dat[i, "estimate"] <- coef(x)[i]
        } 
      }
      j <- len.cfs
    }
    
    if(any(method == "bootstrap") || any(method == "all")){
      
      if(verbose){
        if(inherits(x, "lm"))
          ans.bt <- try(boot_lm(x, R = R), silent = TRUE)  
        if(inherits(x, "gls"))
          ans.bt <- try(boot_gls(x, R = R), silent = TRUE)  
        if(inherits(x, "gnls"))
          ans.bt <- try(boot_gnls(x, R = R), silent = TRUE) 
        if(inherits(x, "lme") && !inherits(x, "nlme"))
          ans.bt <- try(boot_lme(x, R = R), silent = TRUE) 
        if(inherits(x, "nlme"))
          ans.bt <- try(boot_nlme(x, R = R), silent = TRUE) 
      }else{
        if(inherits(x, "lm"))
          ans.bt <- try(suppressMessages(boot_lm(x, R = R)), silent = TRUE)
        if(inherits(x, "gls"))
          ans.bt <- try(suppressMessages(boot_gls(x, R = R)), silent = TRUE)
        if(inherits(x, "gnls"))
          ans.bt <- try(suppressMessages(boot_gnls(x, R = R)), silent = TRUE)
        if(inherits(x, "lme") && !inherits(x, "nlme"))
          ans.bt <- try(suppressMessages(boot_lme(x, R = R)), silent = TRUE)
        if(inherits(x, "nlme"))
          ans.bt <- try(suppressMessages(boot_nlme(x, R = R)), silent = TRUE)
      }
      
      if(inherits(ans.bt, 'try-error'))
        warning('bootstrap method failed')
      
      if(!inherits(ans.bt, 'try-error')){
        
        if(missing(parm)){
          ans <- suppressWarnings(car:::confint.boot(ans.bt, level = level))
        }else{
          ans <- suppressWarnings(car:::confint.boot(ans.bt, parm = parm, level = level))  
        }
        
        if(inherits(x, "lme")){
          for(i in seq_len(dim(ans)[1])){
            ans.dat[i + j, "method"] <- "bootstrap"
            ans.dat[i + j, "parm"] <- names(fixef(x))[i]
            ans.dat[i + j, "lower"] <- ans[i, 1]
            ans.dat[i + j, "upper"] <- ans[i, 2]
            ans.dat[i + j, "estimate"] <- fixef(x)[i]
          } 
        }else{
          for(i in seq_len(dim(ans)[1])){
            ans.dat[i + j, "method"] <- "bootstrap"
            ans.dat[i + j, "parm"] <- names(coef(x))[i]
            ans.dat[i + j, "lower"] <- ans[i, 1]
            ans.dat[i + j, "upper"] <- ans[i, 2]
            ans.dat[i + j, "estimate"] <- coef(x)[i]
          }          
        }
      }
    }
  }
  
  if(inherits(x, "nls")){
  
    len.cfs <- length(coef(x))
    nrows <- length(method) * len.cfs
    ans.dat <- data.frame(method = rep(NA, nrows), parm = NA, level = level,
                          lower = NA, estimate = NA, upper = NA)
    
    if(any(method == "wald") || any(method == "all")){
      require(nlstools)
      ans <- nlstools::confint2(x)
      for(i in seq_len(dim(ans)[1])){
        ans.dat[i, "method"] <- "wald"
        ans.dat[i, "parm"] <- row.names(ans)[i]
        ans.dat[i, "lower"] <- ans[i, 1]
        ans.dat[i, "upper"] <- ans[i, 2]
        ans.dat[i, "estimate"] <- coef(x)[i]
      }
      j <- len.cfs
    }
    
    if(any(method == "profile") || any(method == "all")){
      
      if(verbose){
        ans <- try(confint(x, parm = parm, level = level), silent = TRUE)
      }else{
        ans <- try(suppressMessages(confint(x, parm = parm, level = level)), silent = TRUE)
      }
      
      if(inherits(ans, 'try-error'))
        warning('profile method failed')
      
      if(!inherits(ans, 'try-error')){
        for(i in seq_len(dim(ans)[1])){
          ans.dat[i + j, "method"] <- "profile"
          ans.dat[i + j, "parm"] <- row.names(ans)[i]
          ans.dat[i + j, "lower"] <- ans[i, 1]
          ans.dat[i + j, "upper"] <- ans[i, 2]
          ans.dat[i + j, "estimate"] <- coef(x)[i]
        }
      }
      k <- len.cfs
    }
    
    if(any(method == "bootstrap") || any(method == "all")){
      
      require(car)
      
      if(verbose){
        ans.bt <- try(boot_nls(x, R = R), silent = TRUE)  
      }else{
        ans.bt <- try(suppressMessages(boot_nls(x, R = R)), silent = TRUE)
      }
      
      if(inherits(ans.bt, 'try-error'))
        warning('bootstrap method failed')
      
      if(!inherits(ans.bt, 'try-error')){
        
        if(missing(parm)){
          ans <- suppressWarnings(car:::confint.boot(ans.bt, level = level))
        }else{
          ans <- suppressWarnings(car:::confint.boot(ans.bt, parm = parm, level = level))  
        }
        
        for(i in seq_len(dim(ans)[1])){
          ans.dat[i + j + k, "method"] <- "bootstrap"
          ans.dat[i + j + k, "parm"] <- names(coef(x))[i]
          ans.dat[i + j + k, "lower"] <- ans[i, 1]
          ans.dat[i + j + k, "upper"] <- ans[i, 2]
          ans.dat[i + j + k, "estimate"] <- coef(x)[i]
        }
      }
    }
  }
  return(ans.dat)
}


