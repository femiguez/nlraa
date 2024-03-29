% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_gls.R
\name{simulate_gls}
\alias{simulate_gls}
\title{Simulate fitted values from an object of class \code{\link[nlme]{gls}}}
\usage{
simulate_gls(
  object,
  psim = 1,
  na.action = na.fail,
  naPattern = NULL,
  data = NULL,
  ...
)
}
\arguments{
\item{object}{object of class \code{\link[nlme]{gls}}}

\item{psim}{parameter simulation level, 0: for fitted values, 1: for simulation from 
fixed parameters (assuming a fixed vcov matrix), 2: for simulation considering the 
uncertainty in the residual standard error (sigma), this returns data which
will appear similar to the observed values}

\item{na.action}{default \sQuote{na.fail}. See \code{\link[nlme]{predict.gls}}}

\item{naPattern}{missing value pattern. See \code{\link[nlme]{predict.gls}}}

\item{data}{the data argument is needed when using this function inside user defined functions.
It should be identical to the data used to fit the model.}

\item{...}{additional arguments (it is possible to supply a newdata this way)}
}
\value{
It returns a vector with simulated values with length equal to the number of rows 
in the original data
}
\description{
Simulate values from an object of class gls. Unequal variances, 
as modeled using the \sQuote{weights} option are supported, and there is experimental
code for dealing with the \sQuote{correlation} structure. This generates just one simulation
from these type of models. To generate multiple simulations use \code{\link{simulate_lme}}
}
\details{
This function is based on \code{\link[nlme]{predict.gls}} function

It uses function \code{\link[MASS]{mvrnorm}} to generate new values for the coefficients
of the model using the Variance-Covariance matrix \code{\link{vcov}}. This variance-covariance matrix 
refers to the one for the parameters \sQuote{beta}, not the one for the residuals.
}
\examples{
\donttest{
require(nlme)
data(Orange)

fit.gls <- gls(circumference ~ age, data = Orange, 
               weights = varPower())

## Visualize covariance matrix
fit.gls.vc <- var_cov(fit.gls)
image(log(fit.gls.vc[,ncol(fit.gls.vc):1]))

sim <- simulate_gls(fit.gls)
}
}
\seealso{
\code{\link[nlme]{predict.gls}} \code{\link{simulate_lme}}
}
