## Running other tests to see if the code actually works

test.other.examples <- FALSE

if(test.other.examples){

  require(nlraa)
  require(ggplot2)
  require(segmented)
  require(minpack.lm)
  require(car)
  require(nlstools)
  
  data(plant)
  
  ## For this first case I need to use nlsLM instead
  fit.rwc <- nlsLM(y ~ SSblin(time, a, b, xs, c), data = plant, subset = group == "RWC")
  fit.rkw <- nls(y ~ SSblin(time, a, b, xs, c), data = plant, subset = group == "RKW")
  fit.rkv <- nls(y ~ SSblin(time, a, b, xs, c), data = plant, subset = group == "RKV")
  
  fit.lm <- lm(y ~ time, data = subset(plant, group == "RWC"))
  fit.seg <- segmented(fit.lm)
  ## When comparing the outputs they differ in the interpretation of the
  ## second coefficient for SSblin it is for (x - xs) and for 'segmented'
  ## it must be with 'x' as a predictor, but I'm not sure
  
  ## plot for RWC
  ggplot(data = subset(plant, group == "RWC"), aes(x = time, y = y)) + 
    geom_point() + geom_line(aes(y = fitted(fit.rwc)))
  ## plot for RKW
  ggplot(data = subset(plant, group == "RKW"), aes(x = time, y = y)) + 
    geom_point() + geom_line(aes(y = fitted(fit.rkw)))
  ## plot for RKV
  ggplot(data = subset(plant, group == "RKV"), aes(x = time, y = y)) + 
    geom_point() + geom_line(aes(y = fitted(fit.rkv)))
  
  ## Data 'down'
  ## data(down)
  ## fit <- nls(births ~ SSricker(age, a, b), data = down)
  
  ## ggplot(data = down, aes(age, births)) + 
    ## geom_point() + 
    ## geom_line(aes(y = fitted(fit)))
  
  ## Testing with the swpg data
  data(swpg)
  ## blinear
  fit1 <- nls(lfgr ~ SSblin(ftsw, a, b, xs, c), data = swpg)
  ggplot(data = swpg, aes(x = ftsw, y = lfgr)) + 
         geom_point() + 
         geom_line(aes(y = fitted(fit1)))
  ## linear plateau
  fit2 <- nls(lfgr ~ SSlinp(ftsw, a, b, xs), data = swpg)
  ggplot(data = swpg, aes(x = ftsw, y = lfgr)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit2)))
  
  ## Some tests with datasets from package nlme
  require(nlme)
  data(Relaxin)
  fit.nlis <- nlsList(cAMP ~ SSquadp(conc, a, b, c, xs), data = Relaxin)
  fit.nlme <- nlme(fit.nlis, random = pdDiag(b + xs ~ 1))
  
  data(Fatigue)
  fit.nlis <- nlsList(relLength ~ SSexpf(cycles, a, c), data = Fatigue)
  fit.nlme <- nlme(fit.nlis, random = pdDiag(a + c ~ 1))
  
  data(Cefamandole)
  fit.nlis <- nlsList(conc ~ SSricker(Time, a, b), data = Cefamandole)
  fit.nlme <- nlme(fit.nlis, random = pdDiag(a + b ~ 1))
  
  data(Soybean)
  ## Default model
  fit.nlis <- nlsList(weight ~ SSbgf(Time, w.max, t.e, t.m), data = Soybean)
  ## The previous results in 7 warnings
  fit.soy.1 <- nlme(fit.nlis, random = pdDiag(w.max ~ 1))
  fit.nlis <- nlsList(weight ~ SSbgrp(Time, w.max, lt.e, ldt), data = Soybean)
  ## This results in just one
  fit.soy.2 <- nlme(fit.nlis, random = pdDiag(w.max + lt.e ~ 1))
  ## Four parameter
  fit.nlis <- nlsList(weight ~ SSbgf4(Time, w.max, t.e, t.m, t.b), data = Soybean)
  fit.soy.3 <- nlme(fit.nlis, random = pdDiag(w.max  ~ 1))
  ## Four parameter - reparameterized
  fit.nlis <- nlsList(weight ~ SSbg4rp(Time, w.max, lt.e, ldtm, ldtb), data = Soybean)
  fit.soy.4 <- nlme(fit.nlis, random = pdDiag(w.max + lt.e ~ 1))
  anova(fit.soy.1, fit.soy.2, fit.soy.3, fit.soy.4)
  ## The second model is better, but also the reparameterized versions
  ## allow for including the lt.e random effect which improves the fit
  
  ## Testing the behavior of the ratio function
  
  ## Changing a
  xx <- 1:100
  y1 <- ratio(xx, 1, 1, 1, 1)
  y2 <- ratio(xx, 0.5, 1, 1, 1)
  dat <- data.frame(xx = xx, y1 = y1, y2 = y2)
  ggplot(data = dat) + 
    geom_point(aes(x = xx, y = y1), color = "red") +
    geom_point(aes(x = xx, y = y2), color = "blue")
  ## If I decrease a, y decreases
  ## It makes sense because a is in the numerator
  
  ## Changing b
  y1 <- ratio(xx, 1, 1, 1, 1)
  y2 <- ratio(xx, 1, 0.5, 1, 1)
  dat <- data.frame(xx = xx, y1 = y1, y2 = y2)
  ggplot(data = dat) + 
    geom_point(aes(x = xx, y = y1), color = "red") +
    geom_point(aes(x = xx, y = y2), color = "blue")
  ## If I decrease b, y increases
  ## It makes sense because b is in the denominator
  ## When c and d = 1, then a/b is the asymptote
  
  ## Changing c
  y1 <- ratio(xx, 1, 1, 1, 1)
  y2 <- ratio(xx, 1, 1, 0.5, 1)
  dat <- data.frame(xx = xx, y1 = y1, y2 = y2)
  ggplot(data = dat) + 
    geom_point(aes(x = xx, y = y1), color = "red") +
    geom_point(aes(x = xx, y = y2), color = "blue")
  ## This turns the equation into an exponential decay
  
  ## Changing d, 1
  y1 <- ratio(xx, 1, 1, 1.1, 1.3)
  y2 <- ratio(xx, 1, 1, 1.1, 1.31)
  dat <- data.frame(xx = xx, y1 = y1, y2 = y2)
  ggplot(data = dat) + 
    geom_line(aes(x = xx, y = y1), color = "red") +
    geom_line(aes(x = xx, y = y2), color = "blue")
  ## Changing d, 2
  y1 <- ratio(xx, 1, 0.5, 1.1, 1.3)
  y2 <- ratio(xx, 1, 0.5, 1.1, 1.31)
  dat <- data.frame(xx = xx, y1 = y1, y2 = y2)
  ggplot(data = dat) + 
    geom_line(aes(x = xx, y = y1), color = "red") +
    geom_line(aes(x = xx, y = y2), color = "blue")
  
  ## Testing SSratio
  data(Theoph)
  Theoph.10 <- subset(Theoph, Subject == 10)
  fit.10 <- nlsLM(conc ~ SSratio(Time, a, b, c, d), data = Theoph.10)
  dat <- data.frame(x = Theoph.10$Time, y = Theoph.10$conc)
  ggplot(data = dat, aes(x = x, y = y)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit.10)))
  
  ## Dataset ChickWeight is an interesting one to analyze
  data(ChickWeight)
  chick.11 <- subset(ChickWeight, Chick == 11)
  fit.11.b <- nls(weight ~ bgf2(Time, w.max, w.b = 43, t.e, t.m, t.b = 0), data = chick.11, start = list(w.max = 170, t.e = 17, t.m = 10))
  fit.11.l <- nls(weight ~ SSfpl(Time, A, B, xmid, scal), data = chick.11)
  anova(fit.11.b, fit.11.l) ## RSS is lower for bgf2
  AIC(fit.11.l, fit.11.b) ## AIC is lower for bgf2
  
  ggplot(data = chick.11, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit.11.b)))
  
  ggplot(data = chick.11, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit.11.l)))
  
  chick.43 <- subset(ChickWeight, Chick == 43)
  fit.43.b <- nls(weight ~ SSbgf(Time, w.max, t.e, t.m), data = chick.43)
  fit.43.brp <- nls(weight ~ SSbgrp(Time, w.max, t.e, t.m), data = chick.43, subset = Time > 0)
  fit.43.b4rp <- nls(weight ~ SSbg4rp(Time, w.max, lt.e, ldtm, ldtb), data = chick.43)
  
  ggplot(data = chick.43, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit.43.b)))
  
  ggplot(data = subset(chick.43, Time > 0), aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit.43.brp)))
  
  ggplot(data = chick.43, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit.43.b4rp)))
  
  ## Looking at dataset 'Indometh'
  data(Indometh)
  
  fit <- nlsLM(conc ~ SSratio(time, a, b, c, d), data = Indometh)
  
  ggplot(data = Indometh, aes(x = time, y = conc)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  fit.nlis <- nlsLMList(conc ~ SSratio(time, a, b, c, d), data = Indometh)  
  ## Can't fit a reasonable model using nlme...
  
  ## Looking at 'Orange' dataset
  data(Orange)
  
  fit <- nlsList(circumference ~ SSbg4rp(age, w.max, lt.e, ldtm , ldtb), data = Orange)
  
  fit.nlme2 <- nlme(fit, random = list(w.max  ~ 1))
  
  plot(augPred(fit.nlme2, level = 0:1))
  
  ## Looking at dataset 'pressure'
  data(pressure)
  
  fit <- nls(pressure ~ SSexpf(temperature, a, c), data = pressure)
  
  ggplot(data = pressure, aes(x = temperature, y = pressure)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  ## Looking at dataset 'Loblloly'
  data(Loblolly)
  
  fit <- nlsLM(height ~ SSratio(age, a, b, c, d), data = Loblolly)
  
  ggplot(data = Loblolly, aes(x = age, y = height)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  fit <- nlsLM(height ~ SSblin(age, a, b, xs, c), data = Loblolly)
  
  ggplot(data = Loblolly, aes(x = age, y = height)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  ## Should create some examples from pacakge 'nlstools'
  require(nlstools)
  
  data(O2K)
  fit <- nls(VO2 ~ SSpquad(t, a, xs, b, c), data = O2K)
  
  plotfit(fit)
  
  ggplot(data = O2K, aes(x = t, y = VO2)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  data(vmkm)
  
  fit <- nlsLM(v ~ SSblin(S, a, b, xs, c), data = vmkm)
  
  ggplot(data = vmkm, aes(x = S, y = v)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  ## Let's review the datasets in NISTnls
  library(NISTnls)
  
  ## Bennett5 - not interested
  data(Chwirut1)
  fit <- nls(y ~ SSexpfp(x, a, c, d), data = Chwirut1)
  ggplot(data = Chwirut1, aes(x, y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  
  data(Chwirut2)
  fit <- nls(y ~ SSexpfp(x, a, c, d), data = Chwirut2)
  ggplot(data = Chwirut2, aes(x, y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  ## DanielWood - not interested
  ## ENSO - not interested
  data(Eckerle4)
  fit <- nls(y ~ SSbell(x, ymax, a, b, xc), data = Eckerle4)
  ggplot(data = Eckerle4, aes(x, y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  
  ## Gauss1, Gauss2, Gauss3 - nlraa cannot handle this!
  ## should use splines or gams instead
  
  data(Hahn1)
  fit <- nlsLM(y ~ SSratio(x, a, b, c, d), data = Hahn1, control = list(maxiter = 1e3))
  ggplot(data = Hahn1, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  ## Almost perfect fit, but there are complains
  fit2 <- nlsLM(y ~ SSratio(x, a, b, c, d), data = Hahn1, start = coef(fit))
  ggplot(data = Hahn1, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit2)))
  
  ## Kirby2 - maybe, but there is not much of an error
  ## Lanczos1, 2, and 3 - not too interesting
  data(MGH09)
  fit <- nls(y ~ SSblin(x, a, b, xs, c), data = MGH09)
  ggplot(data = MGH09, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  
  ## MGH10 - not interested
  ## MGH17 - challenging, but do not have a good candidate
  ## Misra1a, b, c and d - not interested
  ## Nelson - not sure what this one is about
  ## Ratkowsky2 and 3 - look like good candidates for logistic or Gompertz
  data(Ratkowsky2)
  fit1 <- nls(y ~ SSlogis(x, Asym, xmid, scal), data = Ratkowsky2)
  ggplot(data = Ratkowsky2, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit1)))
  
  ## What about a 5 parameter logistic?
  fit2 <- nlsLM(y ~ SSlogis5(x, asym1, asym2, xmid, iscal, theta), data = Ratkowsky2, control = list(maxiter = 1e3))
  ggplot(data = Ratkowsky2, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit2)))
  ## Comparing the models
  anova(fit1, fit2)
  ## fit2.bt <- Boot(fit2) It takes too long to run
  ## hist(fit2.bt)
  ## These shows that the parameters are not well constrained under this model
  
  data("Ratkowsky3")
  fit1 <- nls(y ~ SSlogis(x, asym, xmid, scal), data = Ratkowsky3)
  ggplot(data = Ratkowsky3, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit1)))
  
  fit2 <- nlsLM(y ~ SSlogis5(x, asym1, asym2, xmid, iscal, theta), data = Ratkowsky3, control = list(maxiter = 1e3))
  fit3 <- nls(y ~ SSgompertz(x, a, b, c), data = Ratkowsky3)
  ggplot(data = Ratkowsky3, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit2)))
  ggplot(data = Ratkowsky3, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit3)))
  
  ## Five-parameter logistic fits a little bit better
  anova(fit1, fit2, fit3)
  ## fit2.bt <- Boot(fit2) takes too long to run
  ## hist(fit2.bt)
  
  data(Roszman1)
  fit <- nls(y ~ SSexpf(x, a, c), data = Roszman1)
  ggplot(data = Roszman1, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  
  ## Blinear is not actually a terrible fit
  fit <- nls(y ~ SSblin(x, a, b, xs, c), data = Roszman1)
  ggplot(data = Roszman1, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  
  data(Thurber)
  ## Several possible models
  fit <- nls(y ~ SSfpl(x, a, b, c, d), data = Thurber)
  ggplot(data = Thurber, aes(x = x, y = y)) + geom_point() + geom_line(aes(y = fitted(fit)))
  
  ## Another idea is to show and compare the car::Boot and nlstools:nlsBoot functions
  ## These could be very useful, especially when the confidence interval function 'confint'
  ## fails
  
  ## Datasets from 'nlstools'
  data(L.minor)
  ## Blinear
  fit1 <- nls(rate ~ SSblin(conc, a, b, xs1, c), data = L.minor)
  ## MicMen
  fit2 <- nls(rate ~ SSmicmen(conc, Vm, K), data = L.minor)
  ## rational with minpack.lm
  fit3 <- nlsLM(rate ~ SSratio(conc, a, b, c, d), data = L.minor)
  
  ggplot(data = L.minor, aes(conc, rate)) + geom_point() + 
    geom_line(aes(y = fitted(fit1), color = "blinear")) + 
    geom_line(aes(y = fitted(fit2), color = "MicMen")) +
    geom_line(aes(y = fitted(fit3), color = "rational"))
  
  ## Oxygen uptake
  data("O2K")
  fit <- nls(VO2 ~ SSpquad(t, a, xs, b, c), data = O2K)
  ggplot(data = O2K, aes(x = t, y = VO2)) + geom_point() + geom_line(aes(y = fitted(fit)))
  
  ## nlstool::vmkm
  data(vmkm)
  fit1 <- nls(v ~ SSmicmen(S, Vm, K), data = vmkm)
  ggplot(data = vmkm, aes(x = S, y = v)) + geom_point() + geom_line(aes(y = fitted(fit1)))
  
  fit2 <- nls(v ~ SSasymp(S, Asym, R0, lrc), data = vmkm)
  ggplot(data = vmkm, aes(x = S, y = v)) + geom_point() + geom_line(aes(y = fitted(fit2)))
  
  ## Testing function SShill3
  ## With Soybean
  fit <- nls(weight ~ SShill3(Time, Ka, n, a), data = Soybean)
  
  ggplot(data = Soybean, aes(Time, weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  ## With Indometh
  fit <- nls(conc ~ SShill3(time, Ka, n, a), data = Indometh)
  
  ggplot(data = Indometh, aes(x = time, y = conc)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  ## With Loblolly
  fit <- nls(height ~ SShill3(age, Ka, n, a), data = Loblolly)
  
  ggplot(data = Loblolly, aes(x = age, y = height)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit)))
  
  ## With Thurber
  Thurber$x2 <- Thurber$x + 6
  fit <- nlsLM(y ~ SShill3(x2, Ka, n, a), data = Thurber)

  ggplot(data = Thurber, aes(x = x2, y = y)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fit))) 
}
