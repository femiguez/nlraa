## Testing both the psim = 3 and summary_simulate
require(nlme)
require(nlraa)
require(ggplot2)

run.sum.sim <- FALSE

if(run.sum.sim){

  data(Orange)
  fmL <- nlsList(circumference ~ SSlogis(age, Asym, xmid, scal), data = Orange)   
  fm1 <- nlme(fmL, random = pdDiag(Asym + xmid + scal ~ 1))  
  
  prd0 <- predict(fm1, level = 1)

  OrangeA.0 <- cbind(Orange, prd = prd0)  
  
  ## These are the BLUPs
  ggplot(data = OrangeA.0, aes(x = age, y = circumference, color = Tree)) + 
    geom_line(aes(y = prd)) + ggtitle("These are the BLUPs")
  
  ## Predictions only sampling from residuals
  sim2 <- simulate_nlme(fm1, nsim = 500, psim = 2)
  
  ss2 <- summary_simulate(sim2, data = Orange, by = ~ Tree + age)
  
  ggplot(ss2, aes(x = age, y = Estimate, color = Tree)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, color = NULL, fill = Tree), alpha = 0.2)
 
  ## Predictions sampling from residuals and random effects
  sim3 <- simulate_nlme(fm1, nsim = 1e3, psim = 3)

  ss3 <- summary_simulate(sim3, data = Orange, by = ~ Tree + age)

  ggplot(ss3, aes(x = age, y = Estimate, color = Tree)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, color = NULL, fill = Tree), alpha = 0.2) + 
    ggtitle("Here trees have lost their uniqueness compared \n to the previous plot")
  
  ggplot(ss3, aes(x = age, y = Estimate, color = Tree)) + 
    facet_wrap(~Tree) + 
    geom_line() + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, color = NULL, fill = Tree), alpha = 0.2) + 
    ggtitle("Here trees have lost their uniqueness compared \n to the previous plot")
  
  ## Using lme
  fm3 <- lme(circumference ~ age + I(age^2) + I(age^3), random = ~ 1 | Tree, data = Orange)

  sim4 <- simulate_lme(fm3, nsim = 1e3, psim = 2)
 
  ss4 <- summary_simulate(sim4, data = Orange, by = ~ Tree + age)
 
  ggplot(ss4, aes(x = age, y = Estimate, color = Tree)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, color = NULL, fill = Tree), alpha = 0.2)

  ## Predictions sampling from residuals and random effects
  sim5 <- simulate_lme(fm3, nsim = 1e3, psim = 3)

  ss5 <- summary_simulate(sim5, data = Orange, by = ~ Tree + age)
  
  ggplot(ss5, aes(x = age, y = Estimate, color = Tree)) + 
    geom_line() + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, color = NULL, fill = Tree), alpha = 0.2) + 
    ggtitle("Here trees have lost their uniqueness compared \n to the previous plot")   
  
}