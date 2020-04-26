## Testing simulate_gnls
require(nlme)
require(nlraa)
require(ggplot2)
require(car)

run.simulate.gnls.test <- FALSE

if(run.simulate.gnls.test){
  
  set.seed(101)
  
  data(Orange)
  
  fitg <- gnls(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange, weights = varPower())
  
  fitg.bt0 <- boot_nlme(fitg, fitted, psim = 0)
  
  lwr0.q <- apply(t(fitg.bt0$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr0.q <- apply(t(fitg.bt0$t), 1, quantile, probs = 0.95, na.rm = TRUE)
 
  fitg.bt1 <- boot_nlme(fitg, fitted, psim = 1)
  
  lwr1.q <- apply(t(fitg.bt1$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr1.q <- apply(t(fitg.bt1$t), 1, quantile, probs = 0.95, na.rm = TRUE)
  
  ggplot() + 
    geom_point(data = Orange, aes(x = age, y = circumference)) + 
    geom_line(data = Orange, aes(x = age, y = fitted(fitg))) + 
    geom_ribbon(aes(x = Orange$age, ymin = lwr0.q, ymax = upr0.q), fill = "blue", alpha = 0.2) + 
    geom_ribbon(aes(x = Orange$age, ymin = lwr1.q, ymax = upr1.q), fill = "purple", alpha = 0.2)
  
  ## What about the model coefficients?
  fitg.cfs.bt0 <- boot_nlme(fitg, psim = 0) 
  fitg.cfs.bt1 <- boot_nlme(fitg, psim = 1) 

  summary(fitg.cfs.bt0)
  summary(fitg.cfs.bt1)
  
  ## Somewhat wider using bootstrap
  ## Wider when simulating the parameters
  confint(fitg.cfs.bt0)
  confint(fitg.cfs.bt1)
  intervals(fitg)
  
  hist(fitg.cfs.bt0)
  hist(fitg.cfs.bt1)
  
  sims0 <- simulate_gnls(fitg, psim = 0)
  sims1 <- simulate_gnls(fitg, psim = 1)
  sims2 <- simulate_gnls(fitg, psim = 2)
  
  Orange2 <- Orange
  Orange2$sims0 <- sims0
  Orange2$sims1 <- sims1
  Orange2$sims2 <- sims2
  
  ggplot() + 
    geom_point(data = Orange2, aes(x = age, y = circumference)) + 
    geom_line(data = Orange2, aes(x = age, y = sims0)) + 
    geom_point(data = Orange2, aes(x = age, y = sims1), color = "red") +
    geom_line(data = Orange2, aes(x = age, y = sims1), color = "red") +
    geom_point(data = Orange2, aes(x = age, y = sims2), color = "green3")  
    
  ## What about a nlme?
  
  fitL <- nlsList(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)
  
  fit.nm <- nlme(fitL, random = pdDiag(Asym + xmid + scal ~ 1))  
  
  plot(augPred(fit.nm, level = 0:1))
  
  ## Let's look at intervals first
  fit.nm.bt <- boot_nlme(fit.nm)

  confint(fit.nm.bt)  
  intervals(fit.nm)
  
  ## Fitted values...mmm?
  prdf <- function(x) predict(x, level = 0, newdata = data.frame(age = 50:1600)) 
  
  fit.nm.bt.f0 <- boot_nlme(fit.nm, prdf, psim = 0)
 
  lwr0.q <- apply(t(fit.nm.bt.f0$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr0.q <- apply(t(fit.nm.bt.f0$t), 1, quantile, probs = 0.95, na.rm = TRUE)
 
  Orange$fttd <- predict(fit.nm, level = 0)
 
  fit.nm.bt.f1 <- boot_nlme(fit.nm, prdf, psim = 1)
  
  lwr1.q <- apply(t(fit.nm.bt.f1$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr1.q <- apply(t(fit.nm.bt.f1$t), 1, quantile, probs = 0.95, na.rm = TRUE)
  
  ggplot() + 
    geom_point(data = Orange, aes(x = age, y = circumference)) + 
    geom_line(data = Orange, aes(x = age, y = fttd)) + 
    geom_ribbon(aes(x = 50:1600, ymin = lwr0.q, ymax = upr0.q), fill = "blue", alpha = 0.2) + 
    geom_ribbon(aes(x = 50:1600, ymin = lwr1.q, ymax = upr1.q), fill = "purple", alpha = 0.2)
    
  ## What about GAMS
  ## Can I use GAMS for just a few data points?
  ## GAMs do not work for this example, we need more data
  ## library(quantreg)
  
  ## fit.rq <- nlrq(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)
  
  # ggplot(data = Orange, aes(x = age, y = circumference)) + 
  #   geom_point() + 
  #   geom_smooth(method = "gam") + 
  #   geom_line(aes(x = age, y = fit.gm.prd$fit)) + 
  #   geom_ribbon(aes(x = age, ymin = fit.gm.prd$fit - fit.gm.prd$se.fit, ymax = fit.gm.prd$fit + fit.gm.prd$se.fit))
  # 
  
}