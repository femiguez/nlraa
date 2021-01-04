## Testing the levels of simulation
require(ggplot2)
require(nlraa)

run.simulte.lm.test <- Sys.info()[["user"]] == "fernandomiguez"

if(run.simulte.lm.test){
  
  data(Orange)
  
  ## Simple simulation
  fit1 <- lm(circumference ~ age, data = Orange)
  sims1 <- simulate_lm(fit1, nsim = 100, value = "data.frame")
   
  ggplot(data = sims1) + 
    geom_line(aes(x = age, y = sim.y, group = ii), 
               color = "gray", alpha = 0.5) + 
    geom_point(aes(x = age, y = circumference)) + 
    geom_smooth(aes(x = age, y = circumference), method = "lm", 
                linetype = 2, se = FALSE)
  
  prd0 <- cbind(Orange, as.data.frame(predict(fit1, interval = "confidence")))
  
  ggplot() + 
    geom_ribbon(data = prd0, aes(x = age, ymin = lwr, ymax = upr),
                fill = "purple", alpha = 0.5) + 
    geom_point(data = Orange, aes(age, circumference))
    
  ## Testing quadratic
  fit2 <- lm(circumference ~ age + I(age^2), data = Orange)
  sims2 <- simulate_lm(fit2, nsim = 100, value = "data.frame")
  
  ggplot(data = sims2) + 
    geom_line(aes(x = age, y = sim.y, group = ii), 
              color = "gray", alpha = 0.5) + 
    geom_point(aes(x = age, y = circumference)) 
  
  ## Testing psim == 2
  sims3 <- simulate_lm(fit1, psim = 2, nsim = 100, value = "data.frame")
  
  ggplot(data = sims3) + 
    geom_line(aes(x = age, y = sim.y, group = ii), 
              color = "gray", alpha = 0.5) + 
    geom_point(aes(x = age, y = circumference)) 
  
  ## Testing psim == 3
  ## This level is very similar to the 'prediction' below
  sims4 <- simulate_lm(fit1, psim = 3, nsim = 100, value = "data.frame")
  
  sims4a <- aggregate(sim.y ~ age, data = sims4, FUN = quantile,
                      probs = 0.05)
  sims4b <- aggregate(sim.y ~ age, data = sims4, FUN = quantile,
                      probs = 0.95)
  sims4c <- merge(sims4a, sims4b, by = "age")
  
  names(sims4c) <- c("age","lwr","upr")
  
  ggplot() + 
    geom_ribbon(data = sims4c,
                aes(x = age, ymin = lwr, ymax = upr), 
                fill = "blue2", alpha = 0.3) + 
    geom_point(data = Orange, aes(x = age, y = circumference)) 

  prd <- predict(fit1, interval = "prediction")  
  prdd <- cbind(Orange, prd)
    
  ggplot() + 
    geom_ribbon(data = prdd, aes(x = age, ymin = lwr, ymax = upr), 
              fill = "purple", alpha = 0.3) + 
    geom_point(data = Orange, aes(x = age, y = circumference)) 
  
  ## I fitted this model with brms
  ## library(brms)
  ## obr <- brm(circumference ~ age, data = Orange)
  ## This is the same as interval 'confidence'
  ## plot(conditional_effects(obr), points = T)
  ## This is the same as interval 'prediction'
  ## plot(conditional_effects(obr, method = "posterior_predict"), points = T)
 
  ## Trying wild residuals
  ## Testing psim == 2
  sims3 <- simulate_lm(fit1, psim = 2, nsim = 100, value = "data.frame", resid.type = "wild")
  
  sims3$Tree_ID <- with(sims3, paste0(Tree,"_",ii))
  
  ggplot(data = sims3) + 
    geom_line(aes(x = age, y = sim.y, group = Tree_ID), 
              color = "gray", alpha = 0.5) + 
    geom_point(aes(x = age, y = circumference)) 
  
  ## Testing scoping issues
  ## This version does not work
  # f1 <- function(){
  #   dat <- data.frame(x = rnorm(10), y = rnorm(10))
  #   fm00 <- lm(y ~ x, data = dat)
  #   ans <- simulate_lm(fm00)
  #   ans
  # }
  
  f1 <- function(){
    dat <- data.frame(x = rnorm(10), y = rnorm(10))
    fm00 <- lm(y ~ x, data = dat)
    ans <- simulate_lm(fm00, data = dat)
    ans
  }
  ## Well, this seems to work because 'eval' is not inside a function
  ## such as 'getData'???
  res1 <- f1()
  
}