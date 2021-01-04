## Testing simulate_lme
require(nlme)
require(nlraa)
require(ggplot2)

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