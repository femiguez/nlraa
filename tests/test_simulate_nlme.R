## Testing simulate_nlme levels
require(nlme)
require(nlraa)
require(ggplot2)
require(car)

run.test.simulate.nlme <- FALSE

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
    geom_line(aes(y = y.sim, group = ii)) + 
    ggtitle("psim = 0, level = 0")

  ggplot(data = sim01) + 
    geom_point(aes(x = age, y = circumference, color = Tree)) + 
    geom_line(aes(y = y.sim, x = age, color = Tree)) + 
    theme(legend.position = "none") + 
    ggtitle("psim = 0, level = 1")
    
  ggplot(data = sim10) + 
    geom_line(aes(y = y.sim, x = age, group = ii), color = "gray", alpha = 0.4) + 
    geom_point(aes(x = age, y = circumference)) + 
    theme(legend.position = "none") + 
    ggtitle("psim = 1, level = 0")

  sim11$Tree_ID <- with(sim11, paste0(Tree,"_",ii))
  
  ggplot(data = sim11) + 
    facet_wrap(~ Tree) + 
    geom_line(aes(y = y.sim, x = age, color = Tree, group = Tree_ID), alpha = 0.4) + 
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
    geom_line(aes(x = Time, y = y.sim, group = ID, color = Variety), alpha = 0.5) + 
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
    geom_line(aes(x = conc, y = y.sim, group = ID, color = Treatment), alpha = 0.5) + 
    geom_point() 
  
  ## What about bootstrap?
  system.time(fmcm2.bt <- boot_nlme(fmcm2, parallel = "multicore", ncups = 4)) ## 35 seconds on my laptop MacBook Pro
  
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
  
  fmcm4.bt2 <- boot_nlme(fmcm4, prd_fun, parallel = "multicore", ncpus = 4)
  
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
    geom_line(aes(x = conc, y = y.sim, group = ID, color = Treatment), alpha = 0.5) + 
    geom_point() 
  
  ## GAMS
  ## library(mgcv)
  ## This actually does not provide the same level of fit
  # fit <- gam(uptake ~ Type*Treatment + conc + s(conc, k = 4), data = CO2)
  # 
  # ggplot(data = CO2, aes(x = conc, y = uptake)) + 
  #   facet_grid( Type ~ Treatment) + 
  #   geom_line(aes(y = fitted(fit))) + 
  #   geom_point()
}