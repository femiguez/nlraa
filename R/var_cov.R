#' Variance Covariance matrix for (non)linear mixed models 
#' 
#' @title Variance Covariance matrix of for g(n)ls and (n)lme models
#' @name var_cov
#' @description Extracts the variance covariance matrix (residuals, random or all)
#' @param object object which inherits class \code{\link[stats]{lm}}, \code{\link[nlme]{gls}} or \code{\link[nlme]{lme}}
#' @param type \dQuote{residual} for the variance-covariance for the residuals, \dQuote{random}
#' for the variance-covariance of the random effects or \dQuote{all} for the sum of both.
#' @param aug whether to augment the matrix of the random effects to the dimensions of the data
#' @param sparse whether to return a sparse matrix (default is FALSE)
#' @param data optional passing of \sQuote{data}, probably needed when using this function inside other functions.
#' @note See Chapter 5 of Pinheiro and Bates. This returns potentially a very large 
#' matrix of N x N, where N is the number of rows in the data.frame. 
#' The function \code{\link[nlme]{getVarCov}} only works well for  
#' \code{\link[nlme]{lme}} objects. \cr
#' The equivalence is more or less: \cr
#' getVarCov type = \dQuote{random.effects} equivalent to var_cov type = \dQuote{random}. \cr
#' getVarCov type = \dQuote{conditional} equivalent to var_cov type = \dQuote{residual}. \cr
#' getVarCov type = \dQuote{marginal} equivalent to var_cov type = \dQuote{all}. \cr
#' The difference is that getVarCov has an argument that specifies the individual 
#' for which the matrix is being retrieved and var_cov returns the full matrix only.
#' @return It returns a \code{\link[base]{matrix}} or a sparse matrix \code{\link[Matrix]{Matrix}}.
#' @seealso \code{\link[nlme]{getVarCov}}
#' @export
#' @examples 
#' \donttest{
#' require(graphics)
#' require(nlme)
#' data(ChickWeight)
#' ## First a linear model
#' flm <- lm(weight ~ Time, data = ChickWeight)
#' vlm <- var_cov(flm)
#' ## First model with no modeling of the Variance-Covariance
#' fit0 <- gls(weight ~ Time, data = ChickWeight)
#' v0 <- var_cov(fit0)
#' ## Only modeling the diagonal (weights)
#' fit1 <- gls(weight ~ Time, data = ChickWeight, weights = varPower())
#' v1 <- var_cov(fit1)
#' ## Only the correlation structure is defined and there are no groups
#' fit2 <- gls(weight ~ Time, data = ChickWeight, correlation = corAR1())
#' v2 <- var_cov(fit2)
#' ## The correlation structure is defined and there are groups present
#' fit3 <- gls(weight ~ Time, data = ChickWeight, correlation = corCAR1(form = ~ Time | Chick))
#' v3 <- var_cov(fit3)
#' ## There are both weights and correlations
#' fit4 <- gls(weight ~ Time, data = ChickWeight, 
#'             weights = varPower(),
#'             correlation = corCAR1(form = ~ Time | Chick))
#' v4 <- var_cov(fit4)
#' ## Tip: you can visualize these matrices using
#' image(log(v4[,ncol(v4):1]))
#' }

var_cov <- function(object, type = c("residual","random","all", "conditional", "marginal"), 
                    aug = FALSE, sparse = FALSE, data = NULL){
  
  type <- match.arg(type)
  
  if(type == "conditional") type <- "residual"
  if(type == "marginal") type <- "all"
  
  if(type == "random" && inherits(object, c("lm","nls","gls")))
     stop("The variance-covariance of the random effects is only available for \n
          objects which inherit class 'lme' ")
  
  if(isTRUE(sparse)){
    if(!requireNamespace("Matrix", quietly = TRUE)){
      warning("Matrix package is required for this option")
      return(NULL)
    }
  }
     
  if(inherits(object, c("lm","nls"))){
    ans <- diag(nrow = length(fitted(object))) * sigma(object)^2  
    if(sparse) ans <- Matrix::Matrix(ans, sparse = TRUE)
  }
  
  if(inherits(object, c("gls","lme"))){
    if(type == "residual"){
      ans <- var_cov_lme_resid(object, sparse = sparse, data = data)  
    }
    if(type == "random"){
      ans <- var_cov_lme_ranef(object, aug = aug, sparse = sparse, data = data)  
    }
    if(type == "all"){
      ans <- var_cov_lme_resid(object, sparse = sparse, data = data) + var_cov_lme_ranef(object, aug = TRUE, sparse = sparse, data = data)  
    }
  }

  return(ans)
}

var_cov_lme_resid <- function(object, sparse = FALSE, data = data){
  
  if(!inherits(object, c("gls", "lme"))) 
    stop("Only for objects which inherit the 'gls' or 'lme' class")
  
  if(inherits(object, "gls")){
    sgms <- attr(residuals(object), "std")  
  }else{
    sgms <- attr(object[["residuals"]], "std")
  }

  ## If there is no correlation structure this should work
  if(is.null(object$modelStruct$corStruct)){
    Lambda <- diag(sgms^2)     
  }else{
    if(is.null(object$groups)){
      ## When groups are not present corMatrix returns a matrix of appropriate dimensions
      Lambda <- t(sgms * corMatrix(object$modelStruct$corStruct)) * sgms
    }else{
      ## When groups are present corMatrix returns a list and it needs to be converted to
      ## a block diagonal matrix. However, the order of the elements will not necessarily match
      ## how they apprear in the dataset and this is a problem
      corrMat <- corMatrix(object$modelStruct$corStruct)
      
      if(is.list(corrMat)){
        ## If it is not a list it should return a matrix of appropriate dimensions
        ## This extracts the group name
        grp.nm <- deparse(nlme::getGroupsFormula(object$modelStruct$corStruct)[[2]])
        ## I need to extract the groups in the original order in which they appear
        ## in the dataset
        ## Is unique guranteed to return the levels in the order in which they appear?
        ## This is needed because of an evaulation issue when this function is nested
        ## within other functions
        if(is.null(data)){
          gdat <- nlme::getData(object)
        }else{
          gdat <- data
        } 
        ogrpo <- unique(gdat[[grp.nm]])
        ## This reorders the list in the order in which they appear in the dataset
        corrMat <- corrMat[ogrpo]
      }
      
      Lambda <- t(sgms * as.matrix(Matrix::bdiag(corrMat))) * sgms
    }
  }

  ## package Matrix is 'suggested'
  if(sparse){
    Lambda <- Matrix::Matrix(Lambda, sparse = TRUE)
  }
  
  return(Lambda)
}

## This works for objects of class gls, lme and nlme (1 level)
var_cov_lme_ranef <- function(object, aug = FALSE, sparse = FALSE, data = NULL){
  
  if(!inherits(object, c("gls", "lme"))) 
    stop("Only for objects which inherit the 'gls' or 'lme' class")

  if(inherits(object, "nlme") && aug)
    stop("Don't know how to augment nlme random effects.")
  
  ## Number of levels for the random effects
  lreg <- length(names(object$modelStruct$reStruct))
  
  ## If there is only one level and we do not augment
  if(lreg == 1L && aug == FALSE){
    ## This should work for lme and nlme
    ## But not for every object as reStruct[[1]] might not
    ## be coercible to a matrix
      ans <- as.matrix(object$modelStruct$reStruct[[1]]) * sigma(object)^2    
  }else{
    if(!inherits(object, "nlme")){
      ## Does this work for 'lme' object with more than one group?
      V <- mgcv::extract.lme.cov(object, data = data) ## This is V + Z %*% Vr %*% t(Z)
      ## This is clearly potentially dangerous as it might return
      ## negative values, need to test it
      ans <- V - var_cov_lme_resid(object, data = data)      
    }else{
      ans <- vector("list", lreg)
      names(ans) <- names(object$modelStruct$reStruct)
      for(i in 1:lreg){
        tm <- as.matrix(object$modelStruct$reStruct[[i]]) * sigma(object)^2 
        if(sparse) tm <- Matrix::Matrix(tm, sparse = TRUE)
        ans[[i]] <- tm
      }
    }
  }

  if(sparse && !is.list(ans)) ans <- Matrix::Matrix(ans, sparse = TRUE)
  
  return(ans)
}