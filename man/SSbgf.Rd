% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSbgf.R
\name{SSbgf}
\alias{SSbgf}
\alias{bgf}
\alias{bgf2}
\title{self start for Beta Growth Function}
\usage{
bgf(time, w.max, t.e, t.m)

SSbgf(time, w.max, t.e, t.m)

bgf2(time, w.max, w.b, t.e, t.m, t.b)
}
\arguments{
\item{time}{input vector (x) which is normally \sQuote{time}, the smallest value should be close to zero.}

\item{w.max}{value of weight or mass at its peak}

\item{t.e}{time at which the weight or mass reaches its peak.}

\item{t.m}{time at which half of the maximum weight or mass has been reached.}

\item{w.b}{weight or biomass at initial time}

\item{t.b}{initial time offset}
}
\value{
bgf: vector of the same length as x (time) using the beta growth function

bgf2: a numeric vector of the same length as x (time) containing parameter estimates for equation specified
}
\description{
Self starter for Beta Growth function with parameters w.max, t.m and t.e
}
\details{
For details see the publication by Yin et al. (2003) \dQuote{A Flexible Sigmoid Function of Determinate Growth}.

The form of the equation is: \deqn{w.max * (1 + (t.e - time)/(t.e - t.m)) * (time/t.e)^(t.e / (t.e - t.m))}.
Given this function weight is expected to decay and reach zero again at \eqn{2*t.e - t.m}.
}
\examples{
\donttest{
## See extended example in vignette 'nlraa-AgronJ-paper'
x <- seq(0, 17, by = 0.25)
y <- bgf(x, 5, 15, 7)
plot(x, y)
}
}
