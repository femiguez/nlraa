## Testing the predict_varfun
require(nlme)
require(nlraa)
require(ggplot2)

run.test.predict_varFunc <- FALSE

if(run.test.predict_varFunc){
  
  data(Orthodont)
  fit.gls <- gls(distance ~ age, data = Orthodont, weights = varIdent(form = ~ 1 | Sex))
  orth.newdata <- data.frame(age = 8:15, Sex = "Male")
  
  ggplot(Orthodont, aes(age, distance, color = Sex)) + 
    geom_point()
  
  sim1 <- simulate_gls(fit.gls)
  sim2 <- simulate_gls(fit.gls, psim = 2)

  sim1.nd <- simulate_gls(fit.gls, newdata = orth.newdata)
  sim2.nd <- simulate_gls(fit.gls, psim = 2, newdata = orth.newdata)  
 
  orth.newdata.A <- cbind(orth.newdata, prd = sim2.nd)
  
  ggplot() + 
    geom_point(data = Orthodont, aes(age, distance, color = Sex)) + 
    geom_point(data = orth.newdata.A, aes(age, prd))
  
  ## Soybean dataset
  data(Soybean)
  
  ggplot(Soybean, aes(Time, weight, color = Variety)) + 
    geom_point()
  
  fm1 <- gls(weight ~ Time + I(Time^2), Soybean,
              weights = varPower())
  new.soy <- data.frame(Time = 14:84)
  
  sim2.nd <- simulate_gls(fm1, psim = 2, newdata = new.soy)
  new.soy.A <- cbind(new.soy, prd = sim2.nd)
  
  ggplot(data = new.soy.A, aes(x = Time, y = prd)) + 
    geom_point()
  
  fm2 <- gls(prd ~ Time + I(Time^2), new.soy.A,
             weights = varPower())
  
  stds1 <- sigma(fm1)/varWeights(fm1$modelStruct$varStruct)
  stds2 <- sigma(fm2)/varWeights(fm2$modelStruct$varStruct)

  ggplot() + 
    geom_density(aes(x = stds1, color = "Soybean")) + 
    geom_density(aes(x = stds2, color = "new.soy")) 
 
  fm3 <- gls(weight ~ Time + I(Time^2), Soybean,
             weights = varExp())    
  
  sim2.nd <- simulate_gls(fm3, psim = 2, newdata = new.soy) 
  new.soy.A <- cbind(new.soy, prd = sim2.nd) 

  ggplot(data = new.soy.A, aes(x = Time, y = prd)) + 
    geom_point()
  
  ## More complex variance structure with groups
  fm4.e <- gls(weight ~ Time + I(Time^2), Soybean,
             weights = varExp(form = ~ Time | Variety))
  
  fm4.p <- gls(weight ~ Time + I(Time^2), Soybean,
             weights = varPower(form = ~ Time | Variety))
  
  new.soy <- expand.grid(Time = 14:84, Variety = c("F", "P"))
  system.time(prd.e <- predict_gls(fm4.e, interval = "pred", newdata = new.soy))
  system.time(prd.p <- predict_gls(fm4.p, interval = "pred", newdata = new.soy))
}

if(run.test.predict_varFunc){
  ## Now for gnls
 
  data(Soybean)
  # variance increases with a power of the absolute fitted values
  fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
              weights = varPower())
 
  new.soy <- expand.grid(Time = 14:84)
  sim2.nd <- simulate_gnls(fm1, psim = 2, newdata = new.soy)
  new.soy.A <- cbind(new.soy, prd = sim2.nd)

  ggplot(data = new.soy.A, aes(x = Time, y = prd)) + 
    geom_point()   
  
  new.soy <- expand.grid(Time = 14:84, Variety = c("F", "P"))
  ## Changing the covariate
  fm2 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
              weights = varPower(form = ~ Time | Variety))
  
  sim2.nd <- simulate_gnls(fm2, psim = 2, newdata = new.soy) 
  new.soy.A <- cbind(new.soy, prd = sim2.nd)
  
  ggplot(data = new.soy.A, aes(x = Time, y = prd, color = Variety)) + 
    geom_point()   
  
  ## This takes about 42 seconds
  system.time(prds <- predict_gnls(fm2, interval = "pred", newdata = new.soy))
  
  new.soy.P <- cbind(new.soy, prds)

  ggplot(data = new.soy.P, aes(x = Time, y = Estimate, color = Variety)) + 
    geom_point() + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = Variety, color = NULL), alpha = 0.3)
    
}