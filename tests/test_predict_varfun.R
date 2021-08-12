## Testing the predict_varfun
require(nlme)
require(nlraa)
require(ggplot2)

run.test.predict_varfun <- FALSE

if(run.test.predict_varfun){
  
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
  
}