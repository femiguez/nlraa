## Testing simulate_lme
require(nlme)
require(nlraa)
require(ggplot2)
require(car)

run.test.simulate.lme <- Sys.info()[["user"]] == "fernandomiguez"

if(run.test.simulate.lme){
  ## This first section is just to see if it works
  data(Orange)
  
  fm1 <- lme(circumference ~ age, random = ~ 1 | Tree, data = Orange)
  
  Orange$sim0 <- simulate_lme(fm1, psim = 0)
  Orange$sim1 <- simulate_lme(fm1, psim = 1)
  Orange$sim2 <- simulate_lme(fm1, psim = 2)
  Orange$sim3 <- simulate_lme(fm1, psim = 3)
  
  ggplot(data = Orange, aes(x = age, y = circumference)) +
    facet_wrap( ~ Tree) + 
    geom_point() + 
    geom_line(aes(y = sim0, color = "psim = 0")) + 
    geom_line(aes(y = sim1, color = "psim = 1")) + 
    geom_point(aes(y = sim2, color = "psim = 2")) + 
    geom_point(aes(y = sim3, color = "psim = 3"))
  
  data(Soybean)
  
  fm2 <- lme(weight ~ Time, random = ~ 1 | Year/Plot, data = Soybean)
  
  Soybean$sim0 <- simulate_lme(fm2, psim = 0)
  Soybean$sim1 <- simulate_lme(fm2, psim = 1)
  Soybean$sim2 <- simulate_lme(fm2, psim = 2)
  Soybean$sim3 <- simulate_lme(fm2, psim = 3)
  
  ggplot(data = Soybean, aes(x = Time, y = weight)) +
    facet_wrap( ~ Plot) + 
    geom_point() + 
    geom_line(aes(y = sim0, color = "psim = 0")) + 
    geom_line(aes(y = sim1, color = "psim = 1")) + 
    geom_point(aes(y = sim2, color = "psim = 2")) + 
    geom_point(aes(y = sim3, color = "psim = 3"))
  
  fm3 <- lme(weight ~ Time, random = list(Year = ~ 1, Plot = ~1), data = Soybean)
  
  data("Orthodont")
  
  fm4 <- lme(distance ~ age, random = ~ age | Subject, data = Orthodont)

  Orthodont$sim0 <- simulate_lme(fm4, psim = 0)
  Orthodont$sim1 <- simulate_lme(fm4, psim = 1)
  Orthodont$sim2 <- simulate_lme(fm4, psim = 2)
  Orthodont$sim3 <- simulate_lme(fm4, psim = 3)

  ggplot(data = Orthodont, aes(x = age, y = distance)) +
    facet_wrap( ~ Subject) + 
    geom_point() + 
    geom_line(aes(y = sim0, color = "psim = 0")) + 
    geom_line(aes(y = sim1, color = "psim = 1")) + 
    geom_point(aes(y = sim2, color = "psim = 2")) + 
    geom_point(aes(y = sim3, color = "psim = 3"))
  
  ## Assessing how reasonable the results are
  fm1 <- lme(circumference ~ age, random = ~ 1 | Tree, data = Orange)
  
  ## Identical to predicted, level = 0 
  Orange$sim0 <- simulate_lme(fm1, psim = 0, level = 0)
  
  ggplot(data = Orange, aes(x = age, y = circumference)) +
    geom_point() + 
    geom_line(aes(y = sim0))
  
  ## psim = 1, level = 0 (simulation of fixed effects only)
  sim.dat0 <- simulate_lme(fm1, nsim = 100, psim = 1, level = 0, value = "data.frame")
  
  ggplot(data = sim.dat0, aes(x = age, y = circumference)) +
    geom_line(aes(y = sim.y, group = ii), alpha = 0.3) + 
    geom_point() 
  
  ## psim = 1, level = 1, simulations at the subject level
  sim.dat1 <- simulate_lme(fm1, nsim = 100, value = "data.frame")
  
  sim.dat1$Tree_ii <- paste0(sim.dat1$Tree,"_",sim.dat1$ii)
  
  ggplot(data = sim.dat1, aes(x = age, y = circumference)) +
    facet_wrap(~Tree) + 
    geom_line(aes(y = sim.y, group = Tree_ii), color = "gray", alpha = 0.3) + 
    geom_point() 
  
  ## psim = 2, level = 1, simulations at the observation level
  sim.dat2 <- simulate_lme(fm1, nsim = 100, psim = 2, value = "data.frame")
  
  ggplot(data = sim.dat1, aes(x = age, y = circumference)) +
    facet_wrap(~Tree) + 
    geom_point(aes(y = sim.y), color = "gray", alpha = 0.3) + 
    geom_point() + theme_dark()
  
  ## psim = 3, level = 1, simulations at the observation level plus drawing from random effects
  sim.dat3 <- simulate_lme(fm1, nsim = 100, psim = 3, value = "data.frame")
  
  ggplot(data = sim.dat3, aes(x = age, y = circumference)) +
    facet_wrap(~Tree) + 
    geom_point(aes(y = sim.y), color = "gray", alpha = 0.3) + 
    geom_point() 
    
  ## Comparing the behavior of getVarCov and var_cov
  ## These give the same answer
  round((fm1.gvc.re <- getVarCov(fm1, type = "random.effects")))
  round((fm1.vc.re <- var_cov(fm1, type = "random")))
  ## Conditional or residual
  ## Same answer except that var_cov returns the full 
  ## matrix and getVarCov just for the first individual
  (fm1.gvc.cnd <- getVarCov(fm1, type = "conditional"))
  round((fm1.vc.rsd <- var_cov(fm1, type = "residual"))[1:7,1:7])
  ## Visualize
  fm1.vc.rea <- var_cov(fm1, type = "random", aug = TRUE)
  image(fm1.vc.rea[,ncol(fm1.vc.rea):1], main = "random only (augmented)", xaxt = "n", yaxt = "n")
  image(fm1.vc.rsd[,ncol(fm1.vc.rsd):1], main = "residual or conditional", xaxt = "n", yaxt = "n")
  ## Marginal or all
  (fm1.gvc.mrg <- getVarCov(fm1, type = "marginal"))
  round((fm1.vc.all <- var_cov(fm1, type = "all"))[1:7,1:7])
  ## Visualize
  image(fm1.vc.all[,ncol(fm1.vc.all):1], main = "all (random + residual) or marginal", 
        xaxt = "n", yaxt = "n")
  
  ## Testing the 'newdata' argument
  Orange.newdata <- expand.grid(Tree = unique(Orange$Tree),
                                age = seq(100, 1600, 50))
  
  sim.orange.newdata <- simulate_lme(fm1, nsim = 50, newdata = Orange.newdata)

  Orange.newdata$prd <- apply(sim.orange.newdata, 1, quantile, probs = 0.5)
  Orange.newdata$lwr <- apply(sim.orange.newdata, 1, quantile, probs = 0.05)
  Orange.newdata$upr <- apply(sim.orange.newdata, 1, quantile, probs = 0.95)
  
  ggplot() + 
    geom_point(data = Orange, aes(x = age, y = circumference, color = Tree)) + 
    geom_line(data = Orange.newdata, aes(x = age, y = prd, color = Tree)) + 
    geom_ribbon(data = Orange.newdata, aes(x = age, ymin = lwr, ymax = upr, fill = Tree), alpha = 0.1)
}

if(run.test.simulate.lme){
  
  data(barley, package = "nlraa")
  
  fm0 <- lme(yield ~ NF + I(NF^2), random = ~ 1 | year, data = barley)

  prd1 <- predict_lme(fm0, interval = "conf")
  barleyA1 <- cbind(barley, prd1)
  
  ggplot(data = barleyA1, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "blue", alpha = 0.3)
  
  ## Considering the effect of random effects
  ## Note: there is a bug, psim = 3 is not compatible with level = 0
  prd_fun0 <- function(x) predict(x, level = 0)
  ## boot1 <- boot_lme(fm0, prd_fun0, cores = 3)
  boot1 <- boot_lme(fm0, prd_fun0)
  
  barleyA2 <- cbind(barley, summary_simulate(t(boot1$t)))

  ggplot(data = barleyA2, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "blue", alpha = 0.3)
  
  ## Now resample the random effects
  ## boot2 <- boot_lme(fm0, prd_fun0, psim = 3, cores = 3)
  boot2 <- boot_lme(fm0, prd_fun0, psim = 3)
  
  barleyA3 <- cbind(barley, summary_simulate(t(boot2$t)))
  
  ggplot() + 
    geom_point(data = barley, aes(x = NF, y = yield)) + 
    geom_line(data = barleyA2, aes(x = NF, y = Estimate)) + 
    geom_ribbon(data = barleyA2, aes(x = NF, ymin = Q2.5, ymax = Q97.5, fill = "psim = 2"), alpha = 0.3) + 
    geom_ribbon(data = barleyA3, aes(x = NF, ymin = Q2.5, ymax = Q97.5, fill = "psim = 3"), alpha = 0.3)
  
  ### Bootstrapping the covariance parameter
  rand.int <- function(x) sqrt(var_cov(x, type = "random"))
  # boot.cvp.psim1 <- boot_lme(fm0, rand.int, cores = 3)
  # boot.cvp.psim3 <- boot_lme(fm0, rand.int, cores = 3, psim = 3)
  boot.cvp.psim1 <- boot_lme(fm0, rand.int)
  boot.cvp.psim3 <- boot_lme(fm0, rand.int, psim = 3)

  ## Clearly, psim 3 is necessary when bootstrapping the random effects
  hist(boot.cvp.psim1, ci = "perc")
  hist(boot.cvp.psim3, ci = "perc")
  
  ## This second model is better
  fm1 <- lme(yield ~ NF + I(NF^2), random = ~ NF | year, data = barley)
    
  anova(fm0, fm1)
  IC_tab(fm0, fm1)
  
  fm2 <- lme(yield ~ NF + I(NF^2), random = ~ NF + I(NF^2) | year, data = barley)

  anova(fm0, fm1, fm2)
  IC_tab(fm0, fm1, fm2)
  
  prd3 <- predict_lme(fm2, interval = "conf")
  barleyA3 <- cbind(barley, prd3)
  
  ggplot(data = barleyA3, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "blue", alpha = 0.3)
}

if(run.test.simulate.lme){
  
  ## Use datasets to evaluate the variance: total, fixed, random, residual
  data(Orange)
  
  tot.var <- var(Orange$circumference)
  
  fit0 <- lme(circumference ~ age + I(age^2) + I(age^3), 
              random = ~ 1 | Tree, data = Orange)
  
  sim1 <- simulate_lme(fit0, nsim = 1e3, psim = 2)
  sim2 <- simulate_lme(fit0, nsim = 1e3, psim = 3)
  
  varT1 <- apply(sim1, 2, var)
  varT2 <- apply(sim2, 2, var)
  
  ## This somewhat justifies the approach.
  ## The difference is that the value of the particular trees
  ## is very different in psim = 2 vs. psim = 3
  ggplot() + 
    geom_density(aes(x = varT1, color = "psim = 2")) + 
    geom_density(aes(x = varT2, color = "psim = 3")) +
    geom_vline(xintercept = tot.var) 
  
  sim11 <- simulate_lme(fit0, nsim = 1, psim = 2, value = "data.frame")
  
  ggplot(data = sim11, aes(x = age, y = circumference)) + 
    geom_point() + 
    facet_wrap(~ Tree) + 
    geom_point(aes(y = sim.y, color = "simulated"))
  
  sim12 <- simulate_lme(fit0, nsim = 1, psim = 3, value = "data.frame")
  
  ggplot(data = sim12, aes(x = age, y = circumference)) + 
    geom_point() + 
    facet_wrap(~ Tree) + 
    geom_point(aes(y = sim.y, color = "simulated"))
  
  R2M(fit0)
  
  fit1L <- nlsList(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)
  fit1 <- nlme(fit1L, random = pdDiag(Asym + xmid + scal ~ 1))
  
  sim.nl.11 <- simulate_nlme(fit1, nsim = 1e3, psim = 2)
  
  varT.nl.11 <- apply(sim.nl.11, 2, var)
  
  ggplot() + 
    geom_density(aes(x = varT.nl.11, color = "nlme - psim = 2")) + 
    geom_density(aes(x = varT1, color = "lme - psim = 2")) + 
    geom_vline(xintercept = tot.var) 
  
  simPI <- simulate_lme(fit0, nsim = 1e3, psim = 3)
  
  simPIs <- summary_simulate(simPI)
  simPIA <- cbind(Orange, simPIs)
  
  ggplot(data = simPIA) + 
    geom_point(aes(x = age, y = circumference, color = Tree)) + 
    geom_line(aes(x = age, y = Estimate )) + 
    geom_ribbon(aes(x = age, ymin = Q2.5, ymax = Q97.5), alpha = 0.2)

  fm0 <- nls(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)  
  
  fm0.prd <- predict2_nls(fm0, interval = "pred")
  Orange.fm0.A <- cbind(Orange, fm0.prd)
  
  ggplot(data = Orange.fm0.A, aes(x = age, y = circumference)) + 
    geom_point(aes(color = Tree)) + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2)
  
  fmg <- gnls(circumference ~ SSlogis(age, Asym, xmid, scal), 
              weights = varPower(),
              data = Orange)     
  
  fmg.conf <- predict_nlme(fmg, interval = "conf")
  fmg.prd <- predict_nlme(fmg, interval = "pred")
  Orange.fmg.A <- cbind(Orange, fmg.conf, 
                        data.frame(pQ2.5 = fmg.prd[,3], pQ97.5 = fmg.prd[,4]))

  ggplot(data = Orange.fmg.A, aes(x = age, y = circumference)) + 
    geom_point(aes(color = Tree)) + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "red", alpha = 0.2) + 
    geom_ribbon(aes(ymin = pQ2.5, ymax = pQ97.5), fill = "blue", alpha = 0.2)
    
}
