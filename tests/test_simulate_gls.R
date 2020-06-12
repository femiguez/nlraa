## test simulate_lme and simulate_gls
require(nlme)
require(nlraa)
require(ggplot2)

run.simulate.gls.test <- FALSE

if(run.simulate.gls.test){
  
  data(Orthodont)
  
  fm1 <- gls(distance ~ age, data = Orthodont,
             correlation = corCAR1(form = ~ age | Subject)) ## Sub by default
  
  ## Visualize correlation
  fm1.vc <- var_cov(fm1)
  fm1.vc2 <- fm1.vc[1:28, 1:28]
  image(log(fm1.vc2[,ncol(fm1.vc2):1]), 
        main = "Covariance (residual or marginal)",
        xaxt = "n", yaxt = "n")
  
  ## auto-correlation
  plot(ACF(fm1))
  
  sim1 <- simulate_lme(fm1, nsim = 2e2, value = "data.frame")

  sim1$Sub_ii <- paste0(sim1$Subject,"_",sim1$ii)
  
  ggplot(data = sim1, aes(age, distance )) + 
    geom_line(aes(y = sim.y, group = Sub_ii), color = "grey",  
              alpha = 0.3) + 
    geom_point() + ggtitle("Orthodont")
  
  ## Trying the ChickWeight dataset
  data(ChickWeight)
  
  fm2 <- lme(weight ~ Time + I(Time^2) + I(Time^3), 
             data = ChickWeight, weights = varPower(),
             random = ~ Time)
  
  ## Level = 0 is at the population-level           
  sim2 <- simulate_lme(fm2, nsim = 200, level = 0, value = "data.frame")
  
  sim2$Chick_ii <- paste0(sim2$Chick,"_",sim2$ii)
  
  ggplot(data = sim2, aes(Time, weight)) + 
    geom_line(aes(y = sim.y, group = Chick_ii), color = "grey",  
              alpha = 0.3) + 
    geom_point()
  
  ## Test with the Orange dataset
  data("Orange")
  
  fm3 <- lme(circumference ~ age, data = Orange,
             random = ~ 1 | Tree)
  
  ## Small data set, no issues
  plot(fm3)
  
  sim3 <- simulate_lme(fm3, nsim = 200, value = "data.frame")
  
  sim3$Tree_ii <- paste0(sim3$Tree,"_",sim3$ii)
  
  ggplot(data = sim3, aes(age, circumference)) + 
    geom_line(aes(y = sim.y, group = Tree_ii), color = "grey",  
              alpha = 0.3) + 
    geom_point()
  
}