## Test for boot_nlme
require(nlme)
require(nlraa)
require(car)

run.boot.test <- FALSE

data(barley, package = "nlraa")
set.seed(101)
## Does Boot produce 'good' confidence intervals?

if(run.boot.test){
  
  fit.nls <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
  
  ## Profiled confidence intervals
  confint(fit.nls)

  ## Does bootstrap as implemented in car underestimates confidence intervals?  
  fit.nls.bt <- Boot(fit.nls)
  
  fit.gnls <- gnls(yield ~ SSlinp(NF, a, b, xs), data = barley)
  
  intervals(fit.gnls)
  
  system.time(fit.gnls.bt <- boot_nlme(fit.gnls, R = 999))
  
  fit.gnls.bt2 <- boot_nlme(fit.gnls, R = 999, psim = 0)
  
  confint(fit.gnls.bt)
  
  ## I do not know which one is "correct"
  ## Boot seems to have narrower confidence intervals
  ## It does not consider the uncertainty in the fitted values
}

if(run.boot.test){
  
  ## Simple test for barley and an object of class 'gnls'
  ## Simplify the dataset to make the set up simpler
  barley2 <- subset(barley, year < 1974)

  fit.lp.gnls2 <- gnls(yield ~ SSlinp(NF, a, b, xs), data = barley2)

  intervals(fit.lp.gnls2)

  ## Compare this to the bootstrapping approach
  fit.lp.gnls2.bt <- boot_nlme(fit.lp.gnls2, R = 2000)

  summary(fit.lp.gnls2.bt)

  confint(fit.lp.gnls2.bt, type = "perc")

  hist(fit.lp.gnls2.bt, 1, ci = "perc")
  hist(fit.lp.gnls2.bt, 2, ci = "perc")
  hist(fit.lp.gnls2.bt, 3, ci = "perc")

  ## Testing the bootstrap function with an object which
  ## contains factors
  
  barley2$year.f <- as.factor(barley2$year)

  cfs <- coef(fit.lp.gnls2)

  fit.lp.gnls3 <- update(fit.lp.gnls2, 
                         params = list(a + b + xs ~ year.f),
                         start = c(cfs[1], 0, 0, 0, 
                                   cfs[2], 0, 0, 0,
                                   cfs[3], 0, 0, 0))

  intervals(fit.lp.gnls3)

  fit.lp.gnls3.bt <- boot_nlme(fit.lp.gnls3, R = 3000)

  summary(fit.lp.gnls3.bt)
  
  confint(fit.lp.gnls3.bt, type = "perc")
  
  hist(fit.lp.gnls3.bt, 1, ci = "perc")
  hist(fit.lp.gnls3.bt, 2, ci = "perc")
  hist(fit.lp.gnls3.bt, 3, ci = "perc")

  ## Testing the function for 'nlme'. 
  barley$year.f <- as.factor(barley$year)

  barleyG <- groupedData(yield ~ NF | year.f, data = barley)

  fitL.bar <- nlsList(yield ~ SSlinp(NF, a, b, xs), data = barleyG)

  fit.bar.nlme <- nlme(fitL.bar, random = pdDiag(a + b + xs ~ 1))

  plot(augPred(fit.bar.nlme, level = 0:1))
  ## Confidence intervals of the model fixed parameters
  intervals(fit.bar.nlme, which = "fixed")

  ## Bootstrap for the asymptote
  fna <- function(x) fixef(x)[1] + fixef(x)[2] * fixef(x)[3]

  fit.bar.nlme.bt <- boot_nlme(fit.bar.nlme, f = fna, R = 2000)

  confint(fit.bar.nlme.bt, type = "perc")

  ## This is good because the estimate on the original data is
  ## 348.9912, which is in the middle of the confidence interval
  hist(fit.bar.nlme.bt, 1, ci = "perc")

}

## Trying to fit the model using different approaches
## This time with the SSasymp function

if(run.boot.test){
  
  # library(lme4)
  # library(rstanarm)
  
  fmL1 <- nlsList(yield ~ SSasymp(NF, Asym, R0, lrc), data = barleyG)
  
  fm1 <- nlme(fmL1)
  
  ## fm2 <- nlmer(yield ~ SSasymp(NF, Asym, R0, lrc) ~ (Asym | year.f) + (R0 | year.f), data = barley,
  ##             start = getInitial(yield ~ SSasymp(NF, Asym, R0, lrc), data = barley))
  ## This fails. Is nlmer usable?
  
  system.time(fm1.bt <- boot_nlme(fm1, R = 999, parallel = "multicore", ncpus = 2))
  
  summary(fm1.bt)
  
  confint(fm1.bt)
  
  hist(fm1.bt)
  
  ## I feel validated! These intervals are very, very close to the ones
  ## obtained through normal approximation

  ## rstanarm also failed with this example. Chains did not converge
  ## fm3 <- stan_nlmer(yield ~ SSasymp(NF, Asym, R0, lrc) ~ Asym + R0 + lrc | year.f, 
  ##                  data = barley, cores = 2)  
  
}
