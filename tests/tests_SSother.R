## Running other tests to see if the code actually works

test.other.examples <- FALSE

if(test.other.examples){

  require(nlraa)
  require(ggplot2)
  require(segmented)
  require(minpack.lm)
  
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
}
