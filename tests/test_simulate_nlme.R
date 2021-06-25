## Testing simulate_nlme levels
require(nlme)
require(nlraa)
require(ggplot2)
require(car)

run.test.simulate.nlme <- Sys.info()[["user"]] == "fernandomiguez" && FALSE

if(run.test.simulate.nlme){

  data(Orange)

  fitL <- nlsList(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)

  fmm <- nlme(fitL, random = pdDiag(Asym + xmid + scal ~ 1))

  ## This first simulation is at the level of population (fixed, level = 0)
  ## It also means that psim = 0, so there is no sampling from mvrnorm
  sim00 <- simulate_nlme(fmm, nsim = 100, psim = 0, level = 0, value = "data.frame")
  ## Level one means that we simlate at the level of individual Tree (or BLUP in this case)
  ## But psim = 0, so this is the same as predict(x, level = 0)
  sim01 <- simulate_nlme(fmm, nsim = 100, psim = 0, level = 1, value = "data.frame")

  ## What if I simulate?
  sim10 <- simulate_nlme(fmm, nsim = 100, psim = 1, level = 0, value = "data.frame")
  sim11 <- simulate_nlme(fmm, nsim = 100, psim = 1, level = 1, value = "data.frame")

  ## Level 2?
  ## Not sure if it is correct, but it seems to work as intended
  ## The line below does not make sense because we are adding error
  ## at the population level
  ## sim20 <- simulate_nlme(fmm, nsim = 100, psim = 2, level = 0)
  sim21 <- simulate_nlme(fmm, nsim = 100, psim = 2, level = 1, value = "data.frame")

  ggplot(data = sim00, aes(x = age, y = circumference)) + 
    geom_point() + 
    geom_line(aes(y = sim.y, group = ii)) + 
    ggtitle("psim = 0, level = 0")

  ggplot(data = sim01) + 
    geom_point(aes(x = age, y = circumference, color = Tree)) + 
    geom_line(aes(y = sim.y, x = age, color = Tree)) + 
    theme(legend.position = "none") + 
    ggtitle("psim = 0, level = 1")
    
  ggplot(data = sim10) + 
    geom_line(aes(y = sim.y, x = age, group = ii), color = "gray", alpha = 0.4) + 
    geom_point(aes(x = age, y = circumference)) + 
    theme(legend.position = "none") + 
    ggtitle("psim = 1, level = 0")

  sim11$Tree_ID <- with(sim11, paste0(Tree,"_",ii))
  
  ggplot(data = sim11) + 
    facet_wrap(~ Tree) + 
    geom_line(aes(y = sim.y, x = age, color = Tree, group = Tree_ID), alpha = 0.4) + 
    geom_point(aes(x = age, y = circumference)) + 
    theme(legend.position = "none") + 
    ggtitle("psim = 1, level = 1")
  
  ## Test with Soybean example
  data(Soybean)
  
  fmsL <- nlsList(SSlogis, Soybean)
  fmsm <- nlme(fmsL)  

  fxf <- fixef(fmsm)
  fmsm2 <- update(fmsm, fixed = list(Asym + xmid + scal ~ Variety),
                  start = c(fxf[1], 0, fxf[2], 0, fxf[3], 0),
                  random = pdDiag(Asym + xmid + scal ~ 1))
  
  fmsm.s <- simulate_nlme(fmsm2, nsim = 100, value = "data.frame")

  fmsm.s$ID <- with(fmsm.s, paste0(ii,"_",Plot))
  
  ggplot(data = fmsm.s, aes(x = Time, y = weight)) + 
    facet_grid(Variety  ~ Year) + 
    geom_line(aes(x = Time, y = sim.y, group = ID, color = Variety), alpha = 0.5) + 
    geom_point() 

  ## Test with CO2 example
  data(CO2)
  
  fmL <- nlsList(SSasymp, CO2)
  fmcm <- nlme(fmL, random = pdDiag(Asym + R0 ~ 1))
  
  fxf <- fixef(fmcm)
  fmcm2 <- update(fmcm, fixed = list(Asym + R0 + lrc ~ Treatment),
                  start = c(fxf[1], 0, fxf[2], 0, fxf[3], 0))

  fmcm2.s <- simulate_nlme(fmcm2, nsim = 100, value = "data.frame")
 
  fmcm2.s$ID <- with(fmcm2.s, paste0(ii,"_",Plant))
  
  ggplot(data = fmcm2.s, aes(x = conc, y = uptake)) + 
    facet_grid(Type ~ Treatment) + 
    geom_line(aes(x = conc, y = sim.y, group = ID, color = Treatment), alpha = 0.5) + 
    geom_point() 
  
  ## What about bootstrap?
  system.time(fmcm2.bt <- boot_nlme(fmcm2, cores = 4)) ## 35 seconds on my laptop MacBook Pro
  
  hist(fmcm2.bt, 2) 

  ## predictions
  ## Incorporate the fixed effect of type
  fxf2 <- fixef(fmcm2)
  fmcm3 <- update(fmcm2, fixed = list(Asym + R0 + lrc ~ Treatment + Type),
                  start = c(fxf2[1:2], 0, fxf2[3:4], 0, fxf2[5:6], 0),
                  weights = varPower())  
  
  fxf3 <- fixef(fmcm3)
  fmcm4 <- update(fmcm3, fixed = list(Asym + R0 + lrc ~ Treatment + Type + Treatment:Type),
                  start = c(fxf3[1:3], 0, fxf3[4:6], 0, fxf3[7:9], 0))
  
  prd_fun <- function(x) predict(x, level = 0)
  
  fmcm4.bt2 <- boot_nlme(fmcm4, prd_fun, cores = 4)
  
  CO2$lwr <- apply(t(fmcm4.bt2$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  CO2$upr <- apply(t(fmcm4.bt2$t), 1, quantile, probs = 0.95, na.rm = TRUE)

  CO2$ID <- with(CO2, paste0(Treatment,"_",Type))
  
  ggplot(data = CO2, aes(x = conc, y = uptake)) + 
    facet_grid( Type ~ Treatment) + 
    geom_ribbon(aes(x = conc, ymin = lwr, ymax = upr, group = ID, fill = Treatment), alpha = 0.5) + 
    geom_point() 
  
  ## Trying the nested approach
  fmcm4 <- update(fmcm2, random = list(Type = pdDiag(Asym + R0 ~ 1),
                                       Plot = list(Asym + R0 ~ 1)),
                  groups = ~Type/Plant)
  
  prd0 <- predict(fmcm4, newdata = CO2, level = 1)

  fmcm4.s <- simulate_nlme(fmcm4, nsim = 100, value = "data.frame")
  ## This shows that while predict.nlme fails, simulate_nlme works 
  ## (or seems to work)
  
  fmcm4.s$ID <- with(fmcm4.s, paste0(ii,"_",Plant))
  
  ggplot(data = fmcm4.s, aes(x = conc, y = uptake)) + 
    facet_grid( Type ~ Treatment) + 
    geom_line(aes(x = conc, y = sim.y, group = ID, color = Treatment), alpha = 0.5) + 
    geom_point() 
  
  ## GAMS don't quite get it right
  library(mgcv)
  CO2$id <- as.factor(with(CO2, paste0(Type, "_", Treatment)))
  fit.gam <- gam(uptake ~ Type*Treatment + s(conc, k = 3, by = id), data = CO2)
  # 
  ggplot(data = CO2, aes(x = conc, y = uptake)) + 
     facet_grid( Type ~ Treatment) + 
     geom_line(aes(y = fitted(fit.gam))) + 
     geom_point()
  
  ## Compared to gnls model
  fit.lis <- nlsLMList(uptake ~ SSnrh(conc, asym, phi, theta, rd), data = CO2)
  
  plot(fit.lis)
  
  fmco2 <- nlme(fit.lis, random = pdDiag(asym + phi + rd ~ 1))
  
  fmco2.2 <- update(fmco2, fixed = asym + phi + theta + rd ~ Type,
                    start = c(fixef(fmco2)[1], 0, 
                              fixef(fmco2)[2], 0, 
                              fixef(fmco2)[3], 0,
                              fixef(fmco2)[4], 0))

  fmco2.3 <- update(fmco2, fixed = asym + phi + theta + rd ~ Type + Treatment,
                    start = c(fixef(fmco2.2)[1:2], 0, 
                              fixef(fmco2.2)[3:4], 0, 
                              fixef(fmco2.2)[5:6], 0,
                              fixef(fmco2.2)[7:8], 0))
  
  ## Turns out this model has too many parameters and they are not stable
  ## for this reason the predictions are all over the place
  fmco2.4 <- update(fmco2, fixed = asym + phi + theta + rd ~ Type + Treatment + Type:Treatment,
                    start = c(fixef(fmco2.3)[1:3], 0, 
                              fixef(fmco2.3)[4:6], 0, 
                              fixef(fmco2.3)[7:9], 0,
                              fixef(fmco2.3)[10:12], 0))
  
  prd <- predict_nlme(fmco2.4, interval = "conf", level = 0.9)
  CO2A <- cbind(CO2, prd)
  
  CO2AS <- subset(CO2A, Type == "Quebec")
  ggplot(data = CO2AS, aes(x = conc, y = uptake)) + 
    facet_grid( Type ~ Treatment) + 
    geom_line(aes(y = Estimate)) +
    geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = "purple", alpha = 0.3) + 
    geom_point()
  
  ## A simpler model is better
  fmco2.5 <- update(fmco2.4, fixed = list(asym + phi ~ Type + Treatment + Type:Treatment,
                                        theta + rd ~ 1),
                    start = c(fixef(fmco2.3)[1:3], 0, 
                              fixef(fmco2.3)[4:6], 0, 
                              fixef(fmco2.3)[7],
                              fixef(fmco2.3)[10]))
  
  prd <- predict_nlme(fmco2.5, interval = "conf", level = 0.9)
  CO2A <- cbind(CO2, prd)
  
  ggplot(data = CO2A, aes(x = conc, y = uptake)) + 
    facet_grid( Type ~ Treatment) + 
    geom_line(aes(y = Estimate)) +
    geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = "purple", alpha = 0.3) + 
    geom_point()
  
  fit.gnls <- gnls(uptake ~ SSnrh(conc, asym, phi, theta, rd), data = CO2,
                   params = list(asym + phi ~ Type*Treatment, theta + rd ~ 1),
                   start = c(35.35, 0, 0, 0, 0.1577, 0, 0, 0, 0.955, 2.35))

  ## This shows that a 'good' nls model is better than a 'GAM' but 'gams'
  ## are much, much more flexible
  IC_tab(fit.gnls, fit.gam)
  
  ## Testing the implementation of psim = 3 for nlme objects
  ## Picking up the previous Soybean model
  data(Soybean)
  
  fmsL <- nlsList(SSlogis, Soybean)
  fmsm <- nlme(fmsL)  
  
  fxf <- fixef(fmsm)
  fmsm2 <- update(fmsm, fixed = list(Asym + xmid + scal ~ Variety),
                  start = c(fxf[1], 0, fxf[2], 0, fxf[3], 0),
                  random = pdDiag(Asym + xmid + scal ~ 1))
  
  ## Ok, this seems to work
  sim2.psim2 <- simulate_nlme(fmsm2, nsim = 1, 
                              psim = 2, level = 1,
                              value = "data.frame")
  
  sim2.psim2$id <- with(sim2.psim2, paste(Plot, Variety, sep = "_"))
  
  ggplot(data = sim2.psim2, aes(x = Time, y = weight)) + 
    geom_line(aes(x = Time, y = sim.y, group = id), color = "gray", alpha = 0.4) + 
    geom_point() 
  
  ## psim = 3 is not specific for a given subject
  sim2.psim3 <- simulate_nlme(fmsm2, nsim = 10, 
                              psim = 3, level = 1,
                              value = "data.frame")

  sim2.psim3$id <- with(sim2.psim3, paste(Plot, Variety, ii, sep = "_"))
  
  ggplot(data = sim2.psim3, aes(x = Time, y = weight)) + 
    geom_line(aes(x = Time, y = sim.y, group = id), color = "gray", alpha = 0.3) + 
    geom_point() 
  
  ## Since psim = 3 is not specific to an individual subject
  ## we should not expect for the lines to vary around the observed data
  ## in each individual panel. If that is what you want then
  ## you need to use psim = 2
  ggplot(data = sim2.psim3, aes(x = Time, y = weight)) + 
    facet_wrap(~ Variety * Plot) + 
    geom_line(aes(x = Time, y = sim.y, group = ii), color = "gray", alpha = 0.3) + 
    geom_point() 
  
  fmsm3 <- update(fmsm2, random = list(Year = pdDiag(Asym + xmid + scal ~ 1),
                                       Plot = pdDiag(Asym + xmid + scal ~ 1)),
                  groups = ~Year/Plot)
  
  ## I won't check if this gives you the right answer at the moment, but it runs
  sim3.psim3 <- simulate_nlme(fmsm3, nsim = 10, psim = 3, value = "data.frame")
  
  sim3.psim3$id <- with(sim3.psim3, paste(Plot, Variety, Year, ii, sep = "_"))
  
  ggplot(data = sim3.psim3, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(x = Time, y = sim.y, group = id), alpha = 0.1)
  
  ggplot(data = sim3.psim3, aes(x = Time, y = weight)) + 
    facet_wrap(~ Year * Variety) + 
    geom_point() + 
    geom_line(aes(x = Time, y = sim.y, group = id), alpha = 0.1)

  varT <- var(Soybean$weight)  
  
  sim3.psim2.2 <- simulate_nlme(fmsm3, nsim = 2e3, psim = 2)
  sim3.psim3.2 <- simulate_nlme(fmsm3, nsim = 2e3, psim = 3)
  
  varT.psim2 <- apply(sim3.psim2.2, 2, var, na.rm = TRUE)  
  varT.psim3 <- apply(sim3.psim3.2, 2, var, na.rm = TRUE)  
  
  ggplot() + 
    geom_density(aes(x = varT.psim2, color = "psim = 2")) + 
    geom_density(aes(x = varT.psim3, color = "psim = 3")) + 
    geom_vline(xintercept = varT)
  ## Unexpectedly, this shows that the variance for psim = 3
  ## is lower than the variance for psim2, perhaps I need to do more simulations?
  ## Try a different dataset
  
  data(Orange)
  fitL <- nlsList(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)
  fmm <- nlme(fitL, random = pdDiag(Asym + xmid + scal ~ 1))
  
  varT <- var(Orange$circumference)
  fmm.sim.psim2 <- simulate_nlme(fmm, nsim = 1e3, psim = 2)  
  fmm.sim.psim3 <- simulate_nlme(fmm, nsim = 1e3, psim = 3)  
 
  varT.psim2 <- apply(fmm.sim.psim2, 2, var, na.rm = TRUE)  
  varT.psim3 <- apply(fmm.sim.psim3, 2, var, na.rm = TRUE)  
  
  ggplot() + 
    geom_density(aes(x = varT.psim2, color = "psim = 2")) + 
    geom_density(aes(x = varT.psim3, color = "psim = 3")) + 
    geom_vline(xintercept = varT) 
  
  ## Prediction Interval for the mean of a NEW Tree
  nprd <- simulate_nlme(fmm, nsim = 1e3, psim = 3, value = "data.frame")  
  
  ## First step: remove the observation-level variation
  nprdA1 <- aggregate(sim.y ~ ii + age + Tree, data = nprd, FUN = mean)
  
  nprdA <- aggregate(sim.y ~ age, data = nprdA1, FUN = median)
  nprdA$lwr <- aggregate(sim.y ~ age, data = nprdA1, FUN = quantile, probs = 0.05)$sim.y
  nprdA$upr <- aggregate(sim.y ~ age, data = nprdA1, FUN = quantile, probs = 0.95)$sim.y

  ggplot() + 
    geom_point(data = Orange, aes(x = age, y = circumference, color = Tree)) + 
    geom_line(data = nprdA, aes(x = age, y = sim.y)) + 
    geom_ribbon(data = nprdA, aes(x = age, ymin = lwr, ymax = upr), alpha = 0.3) + 
    ggtitle("90% Prediction Interval for a NEW Tree regression")
}