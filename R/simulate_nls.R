#' This function is based on \code{\link[nlme]{predict.gnls}} function 
#' 
#' @title Simulate fitted values from an object of class \code{\link[stats]{nls}}
#' @description Simulate values from an object of class nls.  
#' @name simulate_nls
#' @param object object of class \code{\link[stats]{nls}}
#' @param nsim number of simulations to perform
#' @param psim parameter simulation level, 0: for fitted values, 1: for simulation from 
#' fixed parameters (assuming a fixed vcov matrix), 2: simulation from sampling
#' both from the parameters and the residuals, 3: for simulation considering the 
#' uncertainty in the residual standard error only (sigma) and fixing the
#' parameter estimates at their original value; this will result in simulations 
#' similar to the observed values.
#' @param resid.type either \sQuote{none}, \dQuote{resample}, \dQuote{normal} or \dQuote{wild}.
#' @param value either \sQuote{matrix} or \sQuote{data.frame}
#' @param data the data argument is needed when using this function inside user defined functions.
#' @param ... additional arguments (it is possible to supply a newdata this way)
#' @return It returns a vector with simulated values with length equal to the number of rows 
#' in the original data
#' @details It uses function \code{\link[MASS]{mvrnorm}} to generate new values for the coefficients
#' of the model using the Variance-Covariance matrix \code{\link{vcov}}. This variance-covariance matrix 
#' refers to the one for the parameters \sQuote{beta}, not the one for the residuals.
#' @note The default behavior is that simulations are perfomed for the mean function only.
#' When \sQuote{psim = 2} this function will silently choose \sQuote{resample} as the 
#' \sQuote{resid.type}. This is not ideal design for this function, but I made this choice for 
#' compatibility with other types of simulation originating from \code{\link[stats]{glm}} and
#' \code{\link[mgcv]{gam}}.
#' @seealso \code{\link[nlme]{predict.gnls}}, \code{\link{predict_nls}}
#' @export
#' @examples 
#' \donttest{
#' data(barley, package = "nlraa")
#' 
#' fit <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
#' 
#' sim <- simulate_nls(fit, nsim = 100)
#' }

simulate_nls <- function(object, 
                         nsim = 1, 
                         psim = 1, 
                         resid.type = c("none", "resample", "normal", "wild"),
                         value = c("matrix", "data.frame"), 
                         data = NULL, ...){
  
  ## Error checking
  if(!inherits(object, "nls")) stop("object should be of class 'nls' ")
  
  value <- match.arg(value)
  
  resid.type <- match.arg(resid.type)
  
  if(resid.type == "none" && psim == 2) resid.type <- "resample"
  
  if(is.null(list(...)$newdata)){
    sim.mat <- matrix(ncol = nsim, nrow = stats::nobs(object))
  }else{
    sim.mat <- matrix(ncol = nsim, nrow = nrow(list(...)$newdata))  
  } 
  
  for(i in seq_len(nsim)){
      sim.mat[,i] <- as.vector(simulate_nls_one(object, psim = psim, resid.type = resid.type, data = data, ...))
  }
  
  if(value == "matrix"){
    colnames(sim.mat) <- paste0("sim_",1:nsim)
    return(sim.mat)  
  }else{
    xargs <- list(...)
    if(is.null(xargs$newdata)){
      dat <- eval(object$call$data)
      if(is.null(dat)) stop("'data' argument should be supplied")      
    }else{
      dat <- xargs$newdata
    }
    adat <- data.frame(ii = as.factor(rep(1:nsim, each = nrow(dat))),
                       dat,
                       sim.y = c(sim.mat),
                       row.names = 1:c(nsim * nrow(dat)))   
    return(adat)
  }
}

simulate_nls_one <- function(object, 
                             psim = 1, 
                             resid.type = c("none", "resample","normal","wild"),
                             na.action = na.fail, naPattern = NULL, 
                             data = NULL, ...){
  ##
  ## method for predict() designed for objects inheriting from class gnls
  ##
  if(!inherits(object, "nls")) stop("This function is only for 'nls' objects")
  
  mCall <- object$call
  
  resid.type <- match.arg(resid.type)
  if(resid.type == "none") resid.type <- "resample"
  
  ## Is this more robust?
  args <- list(...)
  if(!is.null(args$newdata)){
    ndata <- args$newdata
    if(psim > 1) stop("'newdata' is not compatible with psim > 1")
  }else{
    if(is.null(data)){
      ndata <- eval(object$data)      
      if(is.null(ndata)) 
        stop("'data' argument is required. It is likely you are using simulate_nls inside another function")
    }else{
      ndata <- data
    } 
  } 
  
  mfArgs <- list(formula =
                   nlme::asOneFormula(formula(object),
                                      mCall$params, naPattern,
                                      omit = c(names(object$m$getPars()), "pi",
                                               deparse(getResponseFormula_nls(object)[[2]]))),
                 data = ndata, na.action = na.action,
                 drop.unused.levels = TRUE)
  
  dataMod <- do.call("model.frame", mfArgs)
  
  ## making sure factor levels are the same as in contrasts
  contr <- object$contrasts
  
  for(i in names(dataMod)) {
    if (inherits(dataMod[,i], "factor") &&
        !is.null(contr[[i]]) && is.matrix(contr[[i]]) ) {
      levs <- levels(dataMod[,i])
      levsC <- dimnames(contr[[i]])[[1]]
      if (any(wch <- is.na(match(levs, levsC)))) {
        stop(sprintf(ngettext(sum(wch),
                              "level %s not allowed for %s",
                              "levels %s not allowed for %s"),
                     paste(levs[wch], collapse = ",")), domain = NA)
      }
      attr(dataMod[,i], "contrasts") <- contr[[i]][levs, , drop = FALSE]
    }
  }
  N <- nrow(dataMod)
  ##
  ## evaluating the naPattern expression, if any
  ##
  naPat <- if(is.null(naPattern)){
    rep(TRUE, N)
  }else{
    as.logical(eval(stats::asOneSidedFormula(naPattern)[[2]], dataMod))
  } 
  ##
  ## Getting  the plist for the new data frame
  ##
  
  ## Create structures similar to gnls
  plist <- vector("list", length(object$m$getPars()))
  for(i in 1:length(object$m$getPars())) plist[[i]] <- TRUE
  names(plist) <- names(object$m$getPars())
  
  pmap <- vector("list", length(object$m$getPars()))
  for(i in 1:length(object$m$getPars())) pmap[[i]] <- i
  names(pmap) <- names(object$m$getPars())
  
  pnames <- names(plist)
  ## This should always result in null with objects of class 'nls'
  if(is.null(params <- eval(object$call$params))){
    params <- list(formula(paste0(paste(pnames, collapse = "+"), "~ 1")))
  }else{ 
    if (!is.list(params)) params <- list(params)
  }
  
  params <- unlist(lapply(params, function(pp){
    if(is.name(pp[[2]])){
      list(pp)
    }else{
      ## multiple parameters on left hand side
      eval(parse(text = paste("list(",
                              paste(paste(all.vars(pp[[2]]), deparse(pp[[3]]), sep = "~"),
                                    collapse = ","),
                              ")")))
    }
  }), recursive=FALSE)
  
  names(params) <- pnames
  
  ##------ This section is by FEM 2020-04-11 ------##
  if(psim == 0){
    prs <- coef(object)  
  }
  
  if(psim == 1){
    ## Only sample from the vector of coefficiencts
    prs <- MASS::mvrnorm(n = 1, mu = coef(object), Sigma = vcov(object))
  }
  
  if(psim == 2){
    prs <- MASS::mvrnorm(n = 1, mu = coef(object), Sigma = vcov(object))
    ## Simply add residuals
    n <- stats::nobs(object)
    ## Change made May 6 2021. Many resources scale the residuals before re-sampling
    ## I always knew this, but it seemed like a vey minor issue
    ## Before this date the residuals were not scaled
    ## See, for example, MASS Venables and Ripley page 225-226
    rsd0 <- scale(stats::resid(object), scale = FALSE)
    if(resid.type == "resample"){
      rsds <- sample(rsd0, size = n, replace = TRUE)      
    }
    if(resid.type == "normal"){
      rsds <- stats::rnorm(n = n, mean = 0, sd = stats::sigma(object))
    }
    if(resid.type == "wild"){
      rsds <- sample(c(-1, 1), size = n, replace = TRUE) * rsd0
    }
  }
  
  if(psim == 3){
    prs <- coef(object)
    ## Simply add residuals
    n <- stats::nobs(object)
    rsd0 <- stats::resid(object)
    if(resid.type == "resample"){
      rsds <- sample(rsd0, size = n, replace = TRUE)      
    }
    if(resid.type == "normal"){
      rsds <- stats::rnorm(n = n, mean = 0, sd = stats::sigma(object))
    }
    if(resid.type == "wild"){
      rsds <- sample(c(-1, 1), size = n, replace = TRUE) * rsd0
    }
  }
  ##------ End FEM section ----------------------##
  
  for(nm in pnames) {
    if (!is.logical(plist[[nm]])) {
      form1s <- stats::asOneSidedFormula(params[[nm]][[3]])
      plist[[nm]] <- model.matrix(form1s, model.frame(form1s, dataMod))
    }
  }
  
  modForm <- nlme::getCovariateFormula(object)[[2]]
  val <- eval(modForm, data.frame(dataMod,
                                  nlraa_getParsNls(plist, pmap, prs, N)))[naPat]
  
  ## ------------ FEM 04-14-2020 ---------------- ##
  ## If the simulation level is 2 we also add residuals
  if(psim >= 2) val <- val + rsds
  ## --------------------------------------------##
  
  names(val) <- row.names(ndata)
  lab <- "Simulated values"
  if (!is.null(aux <- attr(object, "units")$y)) {
    lab <- paste(lab, aux)
  }
  attr(val, "label") <- lab
  
  return(val)
}

##' This function is similar to nlme::getResponseFormula
##' @name getResponseFormula_nls
##' @param object object of class 'nls'
##' @noRd
getResponseFormula_nls <- function(object){
  
  form <- stats::formula(object)
  if (!(inherits(form, "formula") && (length(form) == 3))) {
    stop("'form' must be a two-sided formula")
  }
  eval(parse(text = paste("~", deparse(form[[2]]))))
}

##' @name nlraa_getParsNls
##' @param plist parameter list
##' @param pmap parameter map
##' @param beta beta coefficients
##' @param N I presume the number of observations
##' @noRd
nlraa_getParsNls <- function(plist, pmap, beta, N)
{
  pars <- array(0, c(N, length(plist)), list(NULL, names(plist)))
  for (nm in names(plist)) {
    pars[, nm] <-
      if (is.logical(p <- plist[[nm]]))
        beta[pmap[[nm]]]
    else
      p %*% beta[pmap[[nm]]]
  }
  pars
}
