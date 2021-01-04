# Testing simulate_gam
require(ggplot2)
require(nlraa)
require(nlme)
require(mgcv)

if(Sys.info()[["user"]] == "fernandomiguez"){

  y <- c(12, 14, 33, 50, 67, 74, 123, 141, 165, 204, 253, 246, 240)
  t <- 1:13

  dat <- data.frame(y = y, t = t)

  ggplot(data = dat, aes(x = t, y = y)) + geom_point()

  ## From page 132 form GAM book by Simon Wood
  m1 <- gam(y ~ t + I(t^2), data = dat, family = poisson)

  ggplot(data = dat, aes(x = t, y = y)) + 
    geom_point() + 
    geom_line(aes(y = fitted(m1)))

  m1.sim <- simulate_gam(m1, nsim = 1e3)

  m1.sims <- summary_simulate(m1.sim)
  datA <- cbind(dat, m1.sims)

  ggplot(data = datA, aes(x = t, y = y)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) + 
    ggtitle("Simulate + summary_simulate method")

  ## Compare to native method from mgcv package
  m1.prd <- predict(m1, se.fit = TRUE, type = "response")

  m1.prdd <- data.frame(dat, prd = m1.prd$fit, 
                        lwr = m1.prd$fit - 1.96 * m1.prd$se.fit, 
                        upr = m1.prd$fit + 1.96 * m1.prd$se.fit)

  ggplot(data = m1.prdd, aes(x = t, y = y)) + 
    geom_point() + 
    geom_line(aes(y = prd)) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "purple", alpha = 0.3) + 
    ggtitle("Built-in predict.gam method")

  ## Prediction. The default method uses 'simulate.lm' under the hood
  m1.simP <- simulate_gam(m1, nsim = 1e3, psim = 2)

  m1.simPs <- summary_simulate(m1.simP)
  datAP <- cbind(dat, m1.simPs)

  ggplot(data = datAP, aes(x = t, y = y)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3)

  ## This method includes the uncertainty in the coefficients and the residuals
  m1.simP2 <- simulate_gam(m1, nsim = 1e3, psim = 2, resid.type = "resample")

  m1.simP2s <- summary_simulate(m1.simP2)
  datAP2 <- cbind(dat, m1.simP2s)

  ggplot(data = datAP2, aes(x = t, y = y)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3)

  ## Wild residuals 
  m1.simPW <- simulate_gam(m1, nsim = 1e3, psim = 2, resid.type = "wild")

  m1.simPWs <- summary_simulate(m1.simPW)
  datAPW <- cbind(dat, m1.simPWs)

  ggplot(data = datAPW, aes(x = t, y = y)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3)
  
  ## Dataset Soybean
  data(Soybean)
  
  fms.G <- gam(weight ~ Time + s(Time), data = Soybean)
  fms.C <- lm(weight ~ Time + I(Time^2) + I(Time^3), data = Soybean)
  fms.B <- nls(weight ~ SSbgrp(Time, w.max, lt.e, ldt), data = Soybean)

  ## AIC strongly favors GAM and BIC strongly favors beta
  IC_tab(fms.G, fms.C, fms.B, criteria = "AIC")
  IC_tab(fms.G, fms.C, fms.B, criteria = "BIC")
  
  ggplot(data = Soybean, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fms.C), color = "Cubic")) + 
    geom_line(aes(y = fitted(fms.G), color = "GAM")) + 
    geom_line(aes(y = fitted(fms.B), color = "Beta"))
  
  ## Prediction based on GAM
  prd <- predict_gam(fms.G, interval = "confidence")
  
  SoybeanAG <- cbind(Soybean, prd)
  
  ggplot(data = SoybeanAG, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fms.G), color = "GAM")) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3)

  ## Prediction based on beta
  prdc <- predict_nls(fms.B, interval = "confidence")
  prdp <- predict_nls(fms.B, interval = "prediction")
  colnames(prdp) <- paste0("p", c("Estimate", "Est.Error", "Q2.5", "Q97.5"))
  
  SoybeanAB <- cbind(Soybean, prdc, prdp)
  
  ggplot(data = SoybeanAB, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fms.B), color = "beta")) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) + 
    geom_ribbon(aes(ymin = pQ2.5, ymax = pQ97.5), fill = "purple", alpha = 0.2)
  
  ## Clearly a gnls model would be better
  fms.Bg <- gnls(weight ~ SSbgrp(Time, w.max, lt.e, ldt), data = Soybean,
                 weights = varPower())

  IC_tab(fms.G, fms.C, fms.B, fms.Bg, criteria = "AIC")    
  
  prdc.bg <- predict_nlme(fms.Bg, interval = "confidence")
  prdp.bg <- predict_nlme(fms.Bg, interval = "prediction")  
  colnames(prdp.bg) <- paste0("p", c("Estimate", "Est.Error", "Q2.5", "Q97.5"))

  SoybeanBG <- cbind(Soybean, prdc.bg, prdp.bg)
  
  ggplot(data = SoybeanBG, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fms.Bg), color = "mean")) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, color = "confidence"), fill = "purple", alpha = 0.3) + 
    geom_ribbon(aes(ymin = pQ2.5, ymax = pQ97.5, color = "prediction"), fill = "purple", alpha = 0.2) + 
    ggtitle("Beta growth with increasing variance is the best model")
  
  ## Looking at data barley
  data(barley)
  fgb0 <- gam(yield ~ s(NF, k = 5), data = barley)
  
  fgb0p <- predict_gam(fgb0, interval = "conf")
  barleyA <- cbind(barley, fgb0p)

  prd <- predict(fgb0, se = TRUE)
  barleyG <- barley
  barleyG$fit <- prd$fit
  barleyG$lwr <- prd$fit - 1.96 * prd$se.fit
  barleyG$upr <- prd$fit + 1.96 * prd$se.fit
  
  barleyAG <- merge(barleyA, barleyG)
  
  ## It is hard to see, but they are identical
  ggplot(data = barleyAG, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = fit)) + 
    geom_line(aes(y = Estimate), linetype = 2) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.2) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "yellow", alpha = 0.1)
  
  ## What if we include a random effect of year?
  barley$year.f <- as.factor(barley$year)
  
  fgb1 <- gam(yield ~ s(NF, k = 5) + s(year.f, bs = "re"), data = barley)    

  fgb1p <- predict_gam(fgb1, interval = "conf")
  barleyA <- cbind(barley, fgb1p)

  prd <- predict(fgb1, se = TRUE)
  barleyG <- barley
  barleyG$fit <- prd$fit
  barleyG$lwr <- prd$fit - 1.96 * prd$se.fit
  barleyG$upr <- prd$fit + 1.96 * prd$se.fit
  
  barleyAG <- merge(barleyA, barleyG)
  
  ## It is hard to see, but they are identical
  ggplot(data = barleyAG, aes(x = NF, y = yield)) + 
    facet_wrap(~year.f) + 
    geom_point() + 
    geom_line(aes(y = fit)) + 
    geom_line(aes(y = Estimate), linetype = 2) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "blue", alpha = 0.3) + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), fill = "white", alpha = 0.3)
  
  f1 <- function(){
    dat <- data.frame(x = rnorm(10), y = rnorm(10))
    fm00 <- mgcv::gam(y ~ x, data = dat)
    ans <- simulate_gam(fm00)
    ans
  }
  ## For GAM objects it works regardless because I use 'predict.gam' under the hood and method = "lpmatrix"
  res1 <- f1()

}

