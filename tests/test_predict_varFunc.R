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
  fm4 <- gls(weight ~ Time + I(Time^2), Soybean,
             weights = varExp(form = ~ Time | Variety))    
}

if(run.test.predict_varFunc){
  ## Now for gnls
 
  data(Soybean)
  # variance increases with a power of the absolute fitted values
  fm1 <- gnls(weight ~ SSlogis(Time, Asym, xmid, scal), Soybean,
              weights = varPower())
 
  new.soy <- data.frame(Time = 14:84)
  sim2.nd <- simulate_gnls(fm1, psim = 2, newdata = new.soy)
  new.soy.A <- cbind(new.soy, prd = sim2.nd)

  ggplot(data = new.soy.A, aes(x = Time, y = prd)) + 
    geom_point()   
}