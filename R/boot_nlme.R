#' Bootstraping tools for nonlinear models using a consistent interface
#' 
#' @title Bootstraping for generalized nonlinear models and nonlinear mixed models
#' @name boot_nlme
#' @rdname boot_nlme
#' @param object object of class \code{\link[nlme]{nlme}} or \code{\link[nlme]{gnls}}
#' @param f function to be applied (and bootstrapped), default coef (gnls) or fixef (nlme)
#' @param R number of bootstrap samples, default 999
#' @param psim simulation level for vector of fixed parameters either for \code{\link{simulate_gnls}} or \code{\link{simulate_nlme_one}}
#' @param cores number of cores to use for parallel computation
#' @param data optional data argument (useful/needed when data are not in an available environment).
#' @param ... additional arguments to be passed to function \code{\link[boot]{boot}}
#' @details This function is inspired by \code{\link[car]{Boot}}, which does not
#' seem to work with 'gnls' or 'nlme' objects. This function makes multiple copies 
#' of the original data, so it can be very hungry in terms of memory use, but
#' I do not believe this to be a big problem given the models we typically fit.
#' @export
#' @examples 
#' \donttest{
#' require(car)
#' require(nlme)
#' data(barley, package = "nlraa")
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
#' fit.lp.gnls.Bt3 <- boot_nlme(fit.lp.gnls3, R = 300) 
#' confint(fit.lp.gnls.Bt3, type = "perc")
#' }
#' 

boot_nlme <- function(object, 
                      f = NULL, 
                      R = 999, 
                      psim = 1, 
                      cores = 1L,
                      data = NULL,
                      ...){
  ## I chose to write it in this way instead of UseMethod
  ## because I feel it is more efficient and results in less code
  ## Maybe I'm wrong and will change this in the future
  ## Error checking
  if(!inherits(object, c("gnls","nlme"))) stop("object should be of class 'gnls' or 'nlme'")

  ## extract the original data
  ## This is needed for 'boot'
  dat <- try(nlme::getData(object), silent = TRUE)
  if(inherits(dat, "try-error") && missing(data)){
    stop("'data' argument is required. It is likely you are using boot_nlme inside another function. 
          The data are not in an available environment.", call. = FALSE)    
  } 
  
  if(!missing(data)){
    dat <- data
  }

  if(missing(f)){
    if(inherits(object, "gnls")) f <- coef
    if(inherits(object, "nlme")) f <- fixef
  }
  
  f0 <- f(object) ## To model the behavior of the orignial function
  
  boot_fun_resid <- function(data, indices, fn, model, psim, ...){
  
    ## Copy for bootstrap
    bdat <- data
    ## I need a hack to make the first iteration be t0 for boot
    if(identical(get(".k.boot", envir = nlraa.env), 0L)){
      ## This makes the assumption that the first time
      ## boot calls 'boot_fun_resid' it will generate 't0' 
      fttd <- fitted(model)
      assign(".k.boot", 1L, envir = nlraa.env)
    }else{
      ## Fitted values, which are simulated values if psim = 1
      ## The newdata argument is really only needed to be able to parallelize the code under Windows
      if(inherits(model, "gnls")) fttd <- simulate_gnls(model, psim = psim, newdata = data)
      if(inherits(model, "nlme")) fttd <- simulate_nlme_one(model, psim = psim, newdata = data)
    }
    
    rsds.std <- residuals(model, type = "pearson") ## Extract 'pearson' residuals
   
    ## Different extraction of sigma depending on the object type
    if(inherits(model, "gnls")) rsds.sigma <- attr(rsds.std, "std")
    if(inherits(model, "nlme")) rsds.sigma <- attr(model[["residuals"]], "std")

    new.y <- as.vector(fttd + rsds.std[indices] * rsds.sigma)
    resp.var <- all.vars(formula(model))[1]
    bdat[[resp.var]] <- new.y
    
    assign(".bdat", bdat, envir = nlraa.env)
    umod <- tryCatch(update(model, data = get(".bdat", envir = nlraa.env)), error = function(e){invisible(e)})
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
    vr.lst <- c("nlraa.env", class(object)[1], deparse(object$call$data), "fixef")

    SSfun <- grep("^SS", as.character(getCovariateFormula(object)[[2]]), value = TRUE)
    if(length(SSfun) > 0) vr.lst <- c(vr.lst, SSfun)
    
    parallel::clusterExport(clst, 
                            varlist = vr.lst) ## This exports everything that might be needed?
    ## It seems like in this case the start argument is also needed
    if(inherits(object, "gnls")) object <- update(object, start = coef(object))
    if(inherits(object, "nlme")) object <- update(object, start = fixef(object)) ## Don't know if this is enough for nlme
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
  
  cat("Number of times model fit did not converge",
      sum(is.na(ans$t[,1])),
      "out of",R,"\n")
  
  assign(".bdat", NA, envir = nlraa.env)
  assign(".k.boot", 0L, envir = nlraa.env)
  return(ans)
}

#' @rdname boot_nlme
#' @description bootstrap function for objects of class \code{\link[nlme]{gnls}}
#' @export
boot_gnls <- boot_nlme

#' Create an nlraa environment for bootstrapping
#' 
#' @title Environment to store options and data for nlraa
#' @description Environment which stores indecies and data for bootstraping mostly
#' @export
#' 
nlraa.env <- new.env(parent = emptyenv())
assign('.bdat', NA, nlraa.env)
assign('.k.boot', 0L, nlraa.env)