## Example
require(ggplot2)
require(nlme)
require(nlraa)

## I do not want to run this on CRAN servers
run.predict.nlme <- Sys.info()[["user"]] == "fernandomiguez" && FALSE

if(run.predict.nlme){
  
  data(Orange)
  
  ## All models should be fitted using Maximum Likelihood
  fm.L <- nlme(circumference ~ SSlogis(age, Asym, xmid, scal), 
               random = pdDiag(Asym + xmid + scal ~ 1), 
               method = "ML", data = Orange)
  fm.G <- nlme(circumference ~ SSgompertz(age, Asym, b2, b3), 
               random = pdDiag(Asym + b2 + b3 ~ 1), 
               method = "ML", data = Orange)
  fm.F <- nlme(circumference ~ SSfpl(age, A, B, xmid, scal), 
               random = pdDiag(A + B + xmid + scal ~ 1), 
               method = "ML", data = Orange)
  fm.B <- nlme(circumference ~ SSbg4rp(age, w.max, lt.e, ldtm, ldtb), 
               random = pdDiag(w.max + lt.e + ldtm + ldtb ~ 1), 
               method = "ML", data = Orange)
  
  ## Print the table with weights
  IC_tab(fm.L, fm.G, fm.F, fm.B)
  
  ## Each model prediction is weighted according to their AIC values
  prd <- predict_nlme(fm.L, fm.G, fm.F, fm.B)
  
  ggplot(data = Orange, aes(x = age, y = circumference)) + 
    geom_point() + 
    geom_line(aes(y = predict(fm.L, level = 0), color = "Logistic")) +
    geom_line(aes(y = predict(fm.G, level = 0), color = "Gompertz")) +
    geom_line(aes(y = predict(fm.F, level = 0), color = "4P-Logistic")) +  
    geom_line(aes(y = predict(fm.B, level = 0), color = "Beta")) +
    geom_line(aes(y = prd, color = "Avg. Model"), size = 1.2)
  
  ## Testing confidence intervals
  prdc <- predict_nlme(fm.L, fm.G, fm.F, fm.B, interval = "conf")
  OrangeA <- cbind(Orange, prdc)
  
  ggplot(data = OrangeA, aes(x = age, y = circumference)) + 
    geom_point() + 
    geom_line(aes(y = predict(fm.L, level = 0), color = "Logistic")) +
    geom_line(aes(y = predict(fm.G, level = 0), color = "Gompertz")) +
    geom_line(aes(y = predict(fm.F, level = 0), color = "4P-Logistic")) +  
    geom_line(aes(y = predict(fm.B, level = 0), color = "Beta")) +
    geom_line(aes(y = prd, color = "Avg. Model"), size = 1.2) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3)
  
  ## Confidence at the level of individual or Tree in this case
  prdc1 <- predict_nlme(fm.L, fm.G, fm.F, fm.B, interval = "conf", plevel = 1, level = 0.9)
  OrangeAC <- cbind(Orange, prdc1)
  
  ggplot(data = OrangeAC, aes(x = age, y = circumference, color = Tree)) + 
    facet_wrap(~ Tree, ncol = 5) + 
    geom_point() + 
    geom_line(aes(y = Estimate), size = 1.2) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95, color = NULL, fill = Tree), alpha = 0.3) + 
    ggtitle("90% Confidence bands at the level of 'Tree'")
  
  ## Predictions only make sense at the level of individual or Tree in this case
  prdp <- predict_nlme(fm.L, fm.G, fm.F, fm.B, interval = "pred", plevel = 1, level = 0.9)
  OrangeAP <- cbind(Orange, prdp)
  
  ggplot(data = OrangeAP, aes(x = age, y = circumference, color = Tree)) + 
    facet_wrap(~ Tree, ncol = 5) + 
    geom_point() + 
    geom_line(aes(y = Estimate), size = 1.2) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95, color = NULL, fill = Tree), alpha = 0.3) + 
    ggtitle("90% Prediciton bands at the level of 'Tree'")
  
  ## Example using Soybean
  data(Soybean)
  
  ## All models should be fitted using Maximum Likelihood
  fms.L <- nlme(weight ~ SSlogis(Time, Asym, xmid, scal), 
                random = pdDiag(Asym + xmid ~ 1), 
                weights = varPower(),
                method = "ML", data = Soybean)
  fms.G <- nlme(weight ~ SSgompertz(Time, Asym, b2, b3), 
                random = pdDiag(Asym + b2 ~ 1), 
                weights = varPower(),
                method = "ML", data = Soybean)
  fms.B <- nlme(weight ~ SSbg4rp(Time, w.max, lt.e, ldtm, ldtb), 
                random = pdDiag(w.max + lt.e ~ 1), 
                weights = varPower(),
                method = "ML", data = Soybean) 
  fms.W <- nlme(weight ~ SSweibull(Time, Asym, Drop, lrc, pwr),
                 random = pdDiag(Asym + lrc ~ 1),
                 weights = varPower(),
                 method = "ML", data = Soybean)
  fms.E <- nlme(weight ~ SSexplin(Time, cm, rm, tb),
                random = pdDiag(cm + rm ~ 1),
                weights = varPower(),
                method = "ML", data = Soybean)
  fms.T <- nlme(weight ~ SStrlin(Time, a, b, xs1, c, xs2, d),
                random = pdDiag(b + xs1 + c + d ~ 1),
                weights = varPower(),
                method = "ML", data = Soybean)
  
  IC_tab(fms.L, fms.G, fms.B, fms.W, fms.E, fms.T)
  
  ## The beta growth function is the best one by far
  ## No need to perform averaging
  prdsc <- predict_nlme(fms.B, interval = "conf", level = 0.9)
  SoybeanAC <- cbind(Soybean, prdsc)
  
  ggplot(data = SoybeanAC, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = Estimate), size = 1.2) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = "purple", alpha = 0.3) + 
    ggtitle("90% Confidence bands ")

  prdsp <- predict_nlme(fms.B, interval = "pred", level = 0.9, plevel = 1)
  SoybeanAP <- cbind(Soybean, prdsp)
  
  ggplot(data = SoybeanAP, aes(x = Time, y = weight)) + 
    facet_wrap(~ Plot) + 
    geom_point() + 
    geom_line(aes(y = Estimate), size = 1.2) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = "purple", alpha = 0.3) + 
    ggtitle("90% Prediction bands for each plot ")
  
  ## Example with ChickWeight
  data(ChickWeight)
  
  fmc.L <- nlme(weight ~ SSlogis(Time, Asym, xmid, scal), 
                random = pdDiag(Asym + xmid ~ 1), 
                method = "ML", data = ChickWeight)
  
  fmc.G <- nlme(weight ~ SSgompertz(Time, Asym, b2, b3), 
                random = pdDiag(Asym + b2 ~ 1), 
                method = "ML", data = ChickWeight)
  
  fmc.F <- nlme(weight ~ SSfpl(Time, A, B, xmid, scal), 
                random = pdDiag(A + B + xmid ~ 1), 
                method = "ML", data = ChickWeight)

  fmc.B <- nlme(weight ~ SSbg4rp(Time, w.max, lt.e, ldtm, ldtb), 
                random = pdDiag(w.max + lt.e ~ 1), 
                method = "ML", data = ChickWeight)
  
  ## FPL FTW  
  IC_tab(fmc.L, fmc.G, fmc.F, fmc.B)

  prdsch <- predict_nlme(fmc.F, interval = "conf", level = 0.90)
  ChickWeightAC <- cbind(ChickWeight, prdsch)  
  
  ggplot(data = ChickWeightAC, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = Estimate), size = 1.2) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = "purple", alpha = 0.3) + 
    ggtitle("90% Confidence bands")

  prdschp <- predict_nlme(fmc.F, interval = "pred", level = 0.99, plevel = 1)
  ChickWeightAP <- cbind(ChickWeight, prdschp)  
  
  ggplot(data = ChickWeightAP, aes(x = Time, y = weight)) + 
    facet_wrap(~Chick) + 
    geom_point() + 
    geom_line(aes(y = Estimate), size = 1.2) + 
    geom_ribbon(aes(ymin = Q0.5, ymax = Q99.5), fill = "purple", alpha = 0.3) + 
    ggtitle("99% Prediction bands")
  
}
