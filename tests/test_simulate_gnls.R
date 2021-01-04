## Testing simulate_gnls
require(nlme)
require(nlraa)
require(ggplot2)
require(car)

run.simulate.gnls.test <- Sys.info()[["user"]] == "fernandomiguez" && FALSE

if(run.simulate.gnls.test){
  
  set.seed(101)
  
  data(Orange)
  
  fitg <- gnls(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange, weights = varPower())
  
  system.time(fitg.bt0 <- boot_nlme(fitg, fitted, psim = 0, cores = 4))
  
  lwr0.q <- apply(t(fitg.bt0$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr0.q <- apply(t(fitg.bt0$t), 1, quantile, probs = 0.95, na.rm = TRUE)
 
  system.time(fitg.bt1 <- boot_nlme(fitg, fitted, psim = 1, cores = 4))
  
  lwr1.q <- apply(t(fitg.bt1$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr1.q <- apply(t(fitg.bt1$t), 1, quantile, probs = 0.95, na.rm = TRUE)
  
  ggplot() + 
    geom_point(data = Orange, aes(x = age, y = circumference)) + 
    geom_line(data = Orange, aes(x = age, y = fitted(fitg))) + 
    geom_ribbon(aes(x = Orange$age, ymin = lwr0.q, ymax = upr0.q), fill = "blue", alpha = 0.2) + 
    geom_ribbon(aes(x = Orange$age, ymin = lwr1.q, ymax = upr1.q), fill = "purple", alpha = 0.2)
  
  ## What about the model coefficients?
  system.time(fitg.cfs.bt0 <- boot_nlme(fitg, psim = 0, cores = 4))
  system.time(fitg.cfs.bt1 <- boot_nlme(fitg, psim = 1, cores = 4)) 

  summary(fitg.cfs.bt0)
  summary(fitg.cfs.bt1)
  
  ## Somewhat wider using bootstrap
  ## Wider when simulating the parameters
  confint(fitg.cfs.bt0)
  confint(fitg.cfs.bt1)
  intervals(fitg)
  
  hist(fitg.cfs.bt0)
  hist(fitg.cfs.bt1)
  
  sims0 <- simulate_gnls(fitg, psim = 0) ## Fitted values
  sims1 <- simulate_gnls(fitg, psim = 1) ## One simulation from the mean or population-level
  sims2 <- simulate_gnls(fitg, psim = 2) ## One simulation from the individual-level
  
  Orange2 <- Orange
  Orange2$sims0 <- sims0
  Orange2$sims1 <- sims1
  Orange2$sims2 <- sims2
  
  ggplot() + 
    geom_point(data = Orange2, aes(x = age, y = circumference)) + 
    geom_line(data = Orange2, aes(x = age, y = sims0)) + 
    geom_point(data = Orange2, aes(x = age, y = sims1), color = "red") +
    geom_point(data = Orange2, aes(x = age, y = sims2), color = "green3")  
    
  ## What about a nlme?
  
  fitL <- nlsList(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)
  
  fit.nm <- nlme(fitL, random = pdDiag(Asym + xmid + scal ~ 1))  
  
  plot(augPred(fit.nm, level = 0:1))
  
  ## Let's look at intervals first
  system.time(fit.nm.bt <- boot_nlme(fit.nm, cores = 4))

  confint(fit.nm.bt)  
  intervals(fit.nm, which = "fixed")
  
  ## Fitted values
  prdf <- function(x) predict(x, level = 0, newdata = data.frame(age = 50:1600)) 
  
  system.time(fit.nm.bt.f0 <- boot_nlme(fit.nm, prdf, psim = 0, cores = 4))
 
  lwr0.q <- apply(t(fit.nm.bt.f0$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr0.q <- apply(t(fit.nm.bt.f0$t), 1, quantile, probs = 0.95, na.rm = TRUE)
 
  Orange$fttd <- predict(fit.nm, level = 0)
 
  system.time(fit.nm.bt.f1 <- boot_nlme(fit.nm, prdf, psim = 1, cores = 4))
  
  lwr1.q <- apply(t(fit.nm.bt.f1$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  upr1.q <- apply(t(fit.nm.bt.f1$t), 1, quantile, probs = 0.95, na.rm = TRUE)
  
  ggplot() + 
    geom_point(data = Orange, aes(x = age, y = circumference)) + 
    geom_line(data = Orange, aes(x = age, y = fttd)) + 
    geom_ribbon(aes(x = 50:1600, ymin = lwr0.q, ymax = upr0.q), fill = "blue", alpha = 0.2) + 
    geom_ribbon(aes(x = 50:1600, ymin = lwr1.q, ymax = upr1.q), fill = "purple", alpha = 0.2)
    
  ## Testing the feature of sampling correlated errors
  ## The way to really test this would be to look at the correlation
  ## structure in the new simulations. This is not a good dataset for this
  ## Need to check in Pinheiro and Bates for a good example
  data("ChickWeight")
  
  fitcw <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), data = ChickWeight, 
                correlation = corCAR1(form = ~ Time | Chick),
                weights = varPower())
  
  fitcw.vc <- var_cov(fitcw) 
  image(fitcw.vc[,ncol(fitcw.vc):1], xaxt = "n", yaxt = "n")
  fitcw.vc2 <- fitcw.vc[1:36,1:36]
  image(log(fitcw.vc2[,ncol(fitcw.vc2):1]), xaxt = "n", yaxt = "n")
  ## In this case, this is much slower, because of the huge matrix involved
  fitcw.sim <- simulate_nlme(fitcw, nsim = 100, psim = 2, value = "data.frame")
 
  fitcw.sim$Chick_ID <- with(fitcw.sim, paste0(Chick,"_",ii))
  
  ## It looks like this approach does not capture the specific
  ## Chick-level effect, but it does seem to produce
  ## correlated lines.
  ## Perhaps I cannot expect to reproduce the original data for this 'gnls'
  ## object, but rather, I should check if it creates correlated errors
  ## similar in structure to the ones in the original dataset
  ggplot(data = fitcw.sim) + 
    facet_wrap(~ Chick) + 
    geom_line(aes(x = Time, y = sim.y, color = Chick, group = Chick_ID)) + 
    geom_point(aes(x = Time, y = weight)) +
    theme(legend.position = "none")
}
