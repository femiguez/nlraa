## Test for boot_nlme
require(nlme)
require(nlraa)
require(car)
require(boot)

run.boot.test <- FALSE

if(run.boot.test){
  
  set.seed(101)
  ## Simple test for barley and an object of class 'gnls'
  data(barley, package = "nlraa")

  ## Simplify the dataset to make the set up simpler
  barley2 <- subset(barley, year < 1974)

  fit.lp.gnls2 <- gnls(yield ~ SSlinp(NF, a, b, xs), data = barley2)

  intervals(fit.lp.gnls2)

  ## Compare this to the bootstrapping approach
  fit.lp.gnls2.bt <- boot_nlme(fit.lp.gnls2, R = 2000)

  summary(fit.lp.gnls2.bt)

  confint(fit.lp.gnls2.bt, type = "perc")

  ## hist(fit.lp.gnls2.bt, 1, ci = "perc")
  ## hist(fit.lp.gnls2.bt, 2, ci = "perc")
  ## hist(fit.lp.gnls2.bt, 3, ci = "perc")

  ## Testing the bootstrap function with an object which
  ## contains factors
  ## This does not work and I do not know how to handle the error
  ## [1] "approximate covariance matrix for parameter estimates not of full rank"
  set.seed(101)
  barley2$year.f <- as.factor(barley2$year)

  cfs <- coef(fit.lp.gnls2)

  fit.lp.gnls3 <- update(fit.lp.gnls2, 
                         params = list(a + b + xs ~ year.f),
                         start = c(cfs[1], 0, 0, 0, 
                                   cfs[2], 0, 0, 0,
                                   cfs[3], 0, 0, 0),
                         subset = 1:10)

  intervals(fit.lp.gnls3)

  fit.lp.gnls3.bt <- boot_nlme(fit.lp.gnls3, R = 3000)

  confint(fit.lp.gnls3.bt, type = "perc")
  
  ## hist(fit.lp.gnls3.bt, 1, ci = "perc")

  ## Testing the function for 'nlme'. It seems to work well because
  ## Error handling is somewhat different than for 'gnls'

  barley$year.f <- as.factor(barley$year)

  barleyG <- groupedData(yield ~ NF | year.f, data = barley)

  fitL.bar <- nlsList(yield ~ SSlinp(NF, a, b, xs), data = barleyG)

  fit.bar.nlme <- nlme(fitL.bar, random = pdDiag(a + b + xs ~ 1))

  ## Confidence intervals of the model fixed parameters
  intervals(fit.bar.nlme, which = "fixed")

  ## Bootstrap for the asymptote
  fna <- function(x) fixef(x)[1] + fixef(x)[2] * fixef(x)[3]

  fit.bar.nlme.bt <- boot_nlme(fit.bar.nlme, f = fna, R = 2000)

  confint(fit.bar.nlme.bt, type = "perc")

  hist(fit.bar.nlme.bt, 1, ci = "perc")

}