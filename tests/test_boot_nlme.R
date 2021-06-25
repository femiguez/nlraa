## Test for boot_nlme
require(nlme)
require(nlraa)
require(car)

## I only run this on my computer and occasionally
run.boot.test <- Sys.info()[["user"]] == "fernandomiguez" && FALSE

data(barley, package = "nlraa")
set.seed(101)
## Does Boot produce 'good' confidence intervals?

if(run.boot.test){
  
  fit.nls <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
  
  ## Profiled confidence intervals
  confint(fit.nls)

  ## Does bootstrap as implemented in car underestimates confidence intervals?  
  fit.nls.bt <- Boot(fit.nls)
  ## It is not possible to parallelize the previous code, why?
  
  fit.gnls <- gnls(yield ~ SSlinp(NF, a, b, xs), data = barley)
  
  intervals(fit.gnls)
  
  ## This takes: 9.3 seconds
  system.time(fit.gnls.bt <- boot_nlme(fit.gnls, R = 999, cores = 4))
  
  ## This takes also <10 seconds 
  fit.gnls.bt2 <- boot_nlme(fit.gnls, R = 999, psim = 0, cores = 4)
  
  ## Compare confidence intervals
  confint(fit.nls) ## nls profile
  confint(fit.nls.bt) ## Boot
  confint(fit.gnls.bt) ## boot_nlme psim = 1
  confint(fit.gnls.bt2) ## boot_nlme psim = 0
  ## These intervals do not agree exaclty because different assumptions 
  ## are made and I'm running bootstrap a relatively small number of times
}

if(run.boot.test){
  
  ## Simple test for barley and an object of class 'gnls'
  ## Simplify the dataset to make the set up simpler
  barley2 <- subset(barley, year < 1974)

  fit.lp.gnls2 <- gnls(yield ~ SSlinp(NF, a, b, xs), data = barley2)

  intervals(fit.lp.gnls2)

  ## Compare this to the bootstrapping approach
  ## This take about ~13 seconds
  system.time(fit.lp.gnls2.bt <- boot_nlme(fit.lp.gnls2, R = 2000, cores = 4))

  summary(fit.lp.gnls2.bt) ## Bias is low, which is good

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

  ## This takes 13 seconds (Mac), not bad. Windows (21 seconds)
  system.time(fit.lp.gnls3.bt <- boot_nlme(fit.lp.gnls3, R = 3000, cores = 4))

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

  ## This takes much longer... 215-237 seconds on my laptop
  system.time(fit.bar.nlme.bt <- boot_nlme(fit.bar.nlme, f = fna, R = 2000, cores = 4))

  confint(fit.bar.nlme.bt, type = "perc")

  ## This is good because the estimate on the original data is
  ## 348.9912, which is in the middle of the confidence interval
  hist(fit.bar.nlme.bt, 1, ci = "perc")

}

## Trying to fit the model using different approaches
## This time with the SSasymp function

if(run.boot.test){
  
  fmL1 <- nlsList(yield ~ SSasymp(NF, Asym, R0, lrc), data = barleyG)
  
  fm1 <- nlme(fmL1)
  
  fm2 <- update(fm1, random = pdDiag(Asym + R0 + lrc ~ 1))
  ## Arguably the second model is better
  
  ## fm3 <- nlmer(yield ~ SSasymp(NF, Asym, R0, lrc) ~ (Asym | year.f) + (R0 | year.f), data = barley,
  ##             start = getInitial(yield ~ SSasymp(NF, Asym, R0, lrc), data = barley))
  ## This does not fail now, but it did before. 
  
  ## This is very slow... 17 seconds for 10 attempts
  ## system.time(fm1.bt <- boot_nlme(fm1, R = 10, cores = 4))
  ## What about the second model?
  ## This one takes about 39 sec (Mac and Windows?)
  system.time(fm2.bt <- boot_nlme(fm2, R = 1000, cores = 4))
  
  summary(fm2.bt)
  
  confint(fm2.bt)
  
  hist(fm2.bt, 1)
  hist(fm2.bt, 2)
  hist(fm2.bt, 3)
  
  ## I feel validated! These intervals are very, very close to the ones
  ## obtained through normal approximation
  ## I wrote the previous statement a long time aga, I don't know
  ## if it is still true (2020-5-22)

  ## rstanarm failed with this example. Chains did not converge (2020-05-22)
  ## fm3 <- stan_nlmer(yield ~ SSasymp(NF, Asym, R0, lrc) ~ Asym + R0 + lrc | year.f, 
  ##                  data = barley, cores = 2)  
  
  ## What is the impact of psim = 2 vs. psim = 3 on bootstrapping?
  data(barley, package = "nlraa")
  barley$year.f <- as.factor(barley$year)
  barleyG <- groupedData(yield ~ NF | year.f, data = barley)
  fmL1 <- nlsList(yield ~ SSasymp(NF, Asym, R0, lrc), data = barleyG)
  fm1 <- nlme(fmL1)
  fm2 <- update(fm1, random = pdDiag(Asym + R0 + lrc ~ 1))
  fm3 <- update(fm1, random = pdDiag(Asym + R0 ~ 1))

  anova(fm2, fm3)
  ## Bootstrap the random effects
  vcov_est <- function(x) diag(var_cov(x, type = "random"))
  
  system.time(fm3.vc.bt2 <- boot_nlme(fm3, vcov_est, psim = 2, cores = 4))
  
  fm3.vc.bt2
  
  hist(fm3.vc.bt2, 1, ci = "perc")
  hist(fm3.vc.bt2, 2, ci = "perc")
  
  system.time(fm3.vc.bt3 <- boot_nlme(fm3, vcov_est, psim = 3, cores = 4))

  fm3.vc.bt3
  
  hist(fm3.vc.bt3, 1, ci = "perc")
  hist(fm3.vc.bt3, 2, ci = "perc")
  
  ## Does this show that this method results in biased estimates of the
  ## variance components? Bootstrapping has substantial bias
  ## Is it because I attempt to fix the variance (I use empirical = TRUE in MASS::rmvnorm)
  ## I could either set empirical to FALSE or simulate new variance components...
  
}
