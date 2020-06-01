#' Bootstraping tools for linear mixed-models using a consistent interface
#' 
#' @title Bootstraping for linear mixed models
#' @name boot_lme
#' @param object object of class \code{\link[nlme]{lme}} or \code{\link[nlme]{gnls}}
#' @param f function to be applied (and bootstrapped), default coef (gls) or fixef (lme)
#' @param R number of bootstrap samples, default 999
#' @param psim simulation level for vector of fixed parameters either for \code{\link{simulate_gls}} or \code{\link{simulate_lme}}
#' @param cores number of cores to use for parallel computation
#' @param ... additional arguments to be passed to function \code{\link[boot]{boot}}
#' @details This function is inspired by \code{\link[car]{Boot}}, which does not
#' seem to work with \sQuote{gls} or \sQuote{lme} objects. This function makes multiple copies 
#' of the original data, so it can be very hungry in terms of memory use, but
#' I do not believe this to be a big problem given the models we typically fit.
#' @export
#' @examples 
#' \donttest{
#' require(nlme)
#' require(car)
#' data(Orange)
#'
#' fm1 <- lme(circumference ~ age, random = ~ 1 | Tree, data = Orange)
#' fm1.bt <- boot_lme(fm1, R = 50)
#' 
#' hist(fm1.bt)
#' 
#' }
#' 

boot_lme <- function(object, 
                      f = NULL, 
                      R = 999, 
                      psim = 1, 
                      cores = 1L,
                      ...){
  ## I chose to write it in this way instead of UseMethod
  ## because I feel it is more efficient and results in less code
  ## Maybe I'm wrong and will change this in the future
  ## Error checking
  if(!inherits(object, c("gls","lme"))) stop("object should be of class 'gls' or 'lme'")
  
  ## extract the original data
  ## This is needed for 'boot'
  ## dat <- eval(object$call$data) -- old version
  dat <- nlme::getData(object)
  
  if(missing(f)){
    if(inherits(object, "gls")) f <- coef
    if(inherits(object, "lme")) f <- fixef
  }
  
  f0 <- f(object) ## To model the behavior of the orignial function
  
  boot_fun_resid <- function(data, indices, fn, model, psim, ...){
    
    ## Copy for bootstrap
    bdat <- data
    ## I need a hack to make the first iteration be t0 for boot
    if(identical(get(".k.boot.lme", envir = nlraa.lme.env), 0L)){
      ## This makes the assumption that the first time
      ## boot calls 'boot_fun_resid' it will generate 't0' 
      fttd <- fitted(model)
      assign(".k.boot.lme", 1L, envir = nlraa.lme.env)
    }else{
      ## Fitted values, which are simulated values if psim = 1
      ## The newdata argument is really only needed to be able to parallelize the code under Windows
      if(inherits(model, "gls")) fttd <- simulate_gls(model, psim = psim, newdata = data)
      if(inherits(model, "lme")) fttd <- simulate_lme_one(model, psim = psim, newdata = data)
    }
    
    rsds.std <- residuals(model, type = "pearson") ## Extract 'pearson' residuals
    
    ## Different extraction of sigma depending on the object type
    if(inherits(model, "gls")) rsds.sigma <- attr(rsds.std, "std")
    if(inherits(model, "lme")) rsds.sigma <- attr(model[["residuals"]], "std")
    
    new.y <- as.vector(fttd + rsds.std[indices] * rsds.sigma)
    resp.var <- all.vars(formula(model))[1]
    bdat[[resp.var]] <- new.y
    
    assign(".bdat.lme", bdat, envir = nlraa.lme.env)
    umod <- tryCatch(update(model, data = get(".bdat.lme", envir = nlraa.lme.env)), error = function(e){invisible(e)})
    ## umod <- update(model, data = get(".bdat", envir = nlraa.env))
    
    ## Trying to catch any condition in which the model above does not converge
    if(inherits(umod, "error") || any(is.na(fn(umod))) || is.null(fn(umod)) || any(is.nan(fn(umod)))){
      out <- rep(NA, length(f0))
    }else{
      out <- fn(umod)
    }
    out
  }
  
  prll <- "no" ## default
  clst <- NULL
  if(.Platform$OS.type == "unix" && cores > 1L){
    prll <- "multicore"
  }
  
  if(.Platform$OS.type == "windows" && cores > 1L){
    prll <- "snow"
    clst <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(clst))
    ## Potentially I need to extract the function name
    ## The specific environment depends on the object
    ## Since this is only for gnls or nlme objects it should work...right?
    vr.lst <- c("nlraa.lme.env", class(object)[1], deparse(object$call$data), "fixef")
    
    SSfun <- grep("^SS", as.character(getCovariateFormula(object)[[2]]), value = TRUE)
    if(length(SSfun) > 0) vr.lst <- c(vr.lst, SSfun)
    
    parallel::clusterExport(clst, 
                            varlist = vr.lst) ## This exports everything that might be needed?
  }
  
  ## Override defaults
  args <- list(...)
  if(!is.null(args$parallel)){prll <- args$parallel; parallel <- NULL}
  if(!is.null(args$ncpus)){cores <- args$ncpus; ncpus <- NULL}
  if(!is.null(args$cl)){clst <- args$cl; cl <- NULL}
  
  ans <- boot::boot(data = dat, 
                    stype = "i",
                    statistic = boot_fun_resid, 
                    R = R, fn = f,
                    model = object,
                    psim = psim,
                    parallel = prll,
                    ncpus = cores,
                    cl = clst,
                    ...)
  
  if(sum(is.na(ans$t[,1])) > 0){
    cat("Number of times model fit did not converge",
        sum(is.na(ans$t[,1])),
        "out of",R,"\n")
  }
  
  assign(".bdat.lme", NA, envir = nlraa.lme.env)
  assign(".k.boot.lme", 0L, envir = nlraa.lme.env)
  return(ans)
}

#' Create an nlraa environment for bootstrapping lme
#' 
#' @title Initialize nlraa.env for boot_lme
#' @description Environment which stores indecies and data for bootstraping mostly
#' @noRd
#' 
nlraa.lme.env <- new.env(parent = emptyenv())
assign('.bdat.lme', NA, nlraa.lme.env)
assign('.k.boot.lme', 0L, nlraa.lme.env)