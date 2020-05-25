#' Bootstraping for linear models 
#' 
#' @title Bootstrapping for linear models 
#' @name boot_lm
#' @param object object of class \code{\link[stats]{lm}}
#' @param f function to be applied (and bootstrapped), default coef 
#' @param R number of bootstrap samples, default 999
#' @param psim simulation level for \code{\link{simulate_lm}}
#' @param ... additional arguments to be passed to function \code{\link[boot]{boot}}
#' @details This function is inspired by \code{\link[car]{Boot}}, but has more flexibility
#' in terms of how data are simulated. This function makes multiple copies 
#' of the original data, so it can be very hungry in terms of memory use, but
#' I do not believe this to be a big problem given the models we typically fit.
#' @export
#' @examples 
#' \donttest{
#' require(car)
#' data(barley, package = "nlraa")
#' ## Fit a linear model (quadratic)
#' fit.lm <- lm(yield ~ NF + I(NF^2), data = barley)
#' 
#' ## Bootstrap coefficients by default
#' fit.lm.bt <- boot_lm(fit.lm)
#' ## Compute confidence intervals
#' confint(fit.lm.bt, type = "perc")
#' ## Visualize
#' hist(fit.lm.bt, 1, ci = "perc", main = "Intercept")
#' hist(fit.lm.bt, 2, ci = "perc", main = "NF term")
#' hist(fit.lm.bt, 3, ci = "perc", main = "I(NF^2) term")
#' }
#' 

boot_lm <- function(object, 
                    f = NULL, 
                    R = 999, 
                    psim = 2, 
                    resid.type = c("resample","normal","wild"),
                    verbose = FALSE, ...){
  ## I chose to write it in this way instead of UseMethod
  ## because I feel it is more efficient and results in less code
  ## Maybe I'm wrong and will change this in the future
  ## Error checking
  if(!inherits(object, "lm")) stop("object should be of class 'lm'")
  
  if(psim < 2) stop("psim should be 2 or higher")
  
  resid.type <- match.arg(resid.type)
  
  ## extract the data
  dat <- eval(getCall(object)$data)

  if(missing(f)) f <- coef
    
  f0 <- f(object) ## To model the behavior of the orignial function
  NA.j <- 0
  
  boot_fun_resid <- function(data, fn, model, psim,...){
    
    ## Copy that will be later modified
    bdat <- data
    ## I need a hack to make the first iteration be t0 for boot
    if(identical(get(".k.boot.lm", envir = nlraa.lm.env), 0L)){
      ## This makes the assumption that the first time
      ## boot calls 'boot_fun_resid' it will generate 't0' 
      new.y <- fitted(model) + resid(model)
      assign(".k.boot.lm", 1L, envir = nlraa.lm.env)
    }else{
      ## Fitted values, which are simulated values if psim = 1
      new.y <- simulate_lm(model, psim = psim, resid.type = resid.type, ...)
    }
    
    resp.var <- all.vars(formula(model))[1]
    bdat[[resp.var]] <- new.y
    
    assign(".bdat.lm", bdat, envir = nlraa.lm.env)
    umod <- tryCatch(update(model, data = get(".bdat.lm", envir = nlraa.lm.env)), error = function(e){invisible(e)})
    
    ## It should be rare for an lm model to fail to converge though
    ## Trying to catch any condition in which the model above does not converge
    if(inherits(umod, "error") || any(is.na(fn(umod))) || is.null(fn(umod)) || any(is.nan(fn(umod)))){
      out <- rep(NA, length(f0))
      NA.j <<- NA.j + 1
    }else{
      out <- fn(umod)
    }
    out
  }
  
  ans <- boot::boot(data = dat, 
                    sim = "parametric",
                    statistic = boot_fun_resid, 
                    R = R, fn = f,
                    model = object,
                    psim = psim, ...)
  
  ## There is no reason why a lm should not converge, but just in case
  if(verbose) cat("Number of times model fit did not converge",NA.j,"out of",R,"\n")
  
  assign(".k.boot.lm", 0L, envir = nlraa.lm.env)
  assign(".bdat.lm", NULL, envir = nlraa.lm.env)
  return(ans)
}


### Just initializing variables in the environment
### Not sure if this is even necessary
nlraa.lm.env <- new.env(parent = emptyenv())
assign('.bdat.lm', NA, nlraa.lm.env)
assign('.k.boot.lm', 0L, nlraa.lm.env)