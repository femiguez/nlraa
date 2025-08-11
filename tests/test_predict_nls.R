## Example of a bad model vs. uncertainty vs. model averaging
require(nlraa)
packageVersion("nlraa")
require(car)
require(ggplot2)

run.predict.nls <- Sys.info()[["user"]] == "fernandomiguez" && FALSE

if(run.predict.nls){
 
  data(barley, package = "nlraa")
  
  ggplot(data = barley, aes(x = NF, y = yield)) + 
    geom_point() + 
    xlab("NF (g/m2)") + ylab("Yield (g/m2)") + 
    ggtitle("Barley yield response to N fertilizer")
  
  ## This is not a 'good' model but we'll go with it
  fm.LP <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
  
  sim.LP <- simulate_nls(fm.LP, nsim = 1e3)
  
  ## Does predict work for a single model?
  prd.LP <- predict_nls(fm.LP)
  prd.LP.ci <- predict_nls(fm.LP, interval = "confidence")
  prd.LP.pi <- predict_nls(fm.LP, interval = "prediction")
  
  ggplot(data = barley, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm.LP))) + 
    geom_vline(xintercept = coef(fm.LP)[3]) + 
    xlab("NF (g/m2)") + ylab("Yield (g/m2)") + 
    ggtitle("Linear-plateau fit with break-point")
  
  barleyA <- cbind(barley, summary_simulate(sim.LP, probs = c(0.05, 0.95)))
  
  fm.LP.bt <- boot_nls(fm.LP) ## Bootstrap
  fm.LP.bt.ci <- confint(fm.LP.bt) ## Bootstrap CI
  fm.LP.ci <- confint(fm.LP) ## Profiled CI
  
  ggplot(data = barleyA, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm.LP))) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), alpha = 0.3, fill = "purple") + 
    geom_vline(xintercept = fm.LP.bt$t0[3]) +
    geom_errorbarh(aes(y = 100, xmin = fm.LP.ci[3,1], xmax = fm.LP.ci[3,2], 
                       color = "profiled"), color = "blue") + 
    geom_errorbarh(aes(y = 50, xmin = fm.LP.bt.ci[3,1], xmax = fm.LP.bt.ci[3,2], 
                       color = "bootstrap"), color = "purple") + 
    geom_text(aes(x = 13, y = 100, label = "profiled"), color = "blue") +  
    geom_text(aes(x = 13, y = 50, label = "bootstrap"), color = "purple") +
    xlab("NF (g/m2)") + ylab("Yield (g/m2)") + 
    ggtitle("90% uncertainty bands and intervals for the break-point")
  
  ## What if we fit several models?
  fm.L <- lm(yield ~ NF, data = barley)
  fm.Q <- lm(yield ~ NF + I(NF^2), data = barley)
  fm.A <- nls(yield ~ SSasymp(NF, Asym, R0, lrc), data = barley)
  fm.BL <- nls(yield ~ SSblin(NF, a, b, xs, c), data = barley)
  
  print(IC_tab(fm.L, fm.Q, fm.A, fm.LP, fm.BL), digits = 2)
  
  ggplot(data = barley, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm.L), color = "Linear")) +
    geom_line(aes(y = fitted(fm.Q), color = "Quadratic")) +
    geom_line(aes(y = fitted(fm.A), color = "Asymptotic")) +  
    geom_line(aes(y = fitted(fm.LP), color = "Linear-plateau")) + 
    geom_line(aes(y = fitted(fm.BL), color = "Bi-linear")) + 
    xlab("NF (g/m2)") + ylab("Yield (g/m2)") + 
    ggtitle("Different model fits")
  
  ## Each model prediction is weighted using the AIC values
  prd <- predict_nls(fm.L, fm.Q, fm.A, fm.LP, fm.BL)
  prdc <- predict_nls(fm.L, fm.Q, fm.A, fm.LP, fm.BL, interval = "confidence")
  prdp <- predict_nls(fm.L, fm.Q, fm.A, fm.LP, fm.BL, interval = "prediction")
  
  ggplot(data = barley, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm.L), color = "Linear")) +
    geom_line(aes(y = fitted(fm.Q), color = "Quadratic")) +
    geom_line(aes(y = fitted(fm.A), color = "Asymptotic")) +  
    geom_line(aes(y = fitted(fm.LP), color = "Linear-plateau")) + 
    geom_line(aes(y = fitted(fm.BL), color = "Bi-linear")) + 
    geom_line(aes(y = prd, color = "Avg. Model"), size = 1.2, color = "black") + 
    xlab("NF (g/m2)") + ylab("Yield (g/m2)") + 
    ggtitle("Different model fits and average model weighted by AIC")
  
  ggplot(data = barley, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm.L), color = "Linear")) +
    geom_line(aes(y = fitted(fm.Q), color = "Quadratic")) +
    geom_line(aes(y = fitted(fm.A), color = "Asymptotic")) +  
    geom_line(aes(y = fitted(fm.LP), color = "Linear-plateau")) + 
    geom_line(aes(y = fitted(fm.BL), color = "Bi-linear")) + 
    geom_line(aes(y = prd, color = "Avg. Model"), size = 1.2, color = "black") + 
    geom_ribbon(aes(ymin = prdc[,3], ymax = prdc[,4]), 
                fill = "purple", alpha = 0.3) + 
    geom_ribbon(aes(ymin = prdp[,3], ymax = prdp[,4]), 
                fill = "purple", alpha = 0.1) + 
    xlab("NF (g/m2)") + ylab("Yield (g/m2)") + 
    ggtitle("Model fits, 90% uncertainty bands for confidence and prediction")

  ## Do GAMs work?
  require(mgcv)
  
  fm.L <- lm(yield ~ NF, data = barley)
  fm.Q <- lm(yield ~ NF + I(NF^2), data = barley)
  fm.C <- lm(yield ~ NF + I(NF^2) + I(NF^3), data = barley)
  fm.A <- nls(yield ~ SSasymp(NF, Asym, R0, lrc), data = barley)
  fm.LP <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley)
  fm.G <- gam(yield ~ NF + s(NF, k = 3), data = barley)

  fm.Gs <- simulate_lm(fm.G, nsim = 1e3)
  fm.Gss <- summary_simulate(fm.Gs, probs = c(0.05, 0.95))
  barleyAS <- cbind(barley, fm.Gss)
  
  ## The default predict method for GAMs does not produce intervals
  ## But we can generate them
  fm.Gp <- predict(fm.G, se.fit = TRUE)
  qnt <- qt(0.05, 72)
  fm.Gpd <- data.frame(prd = fm.Gp$fit, 
                       lwr = fm.Gp$fit + qnt * fm.Gp$se.fit,
                       upr = fm.Gp$fit - qnt * fm.Gp$se.fit)
  ## These intervals are almost exactly the same as the ones
  ## obtained through simulation
  
  print(IC_tab(fm.L, fm.Q, fm.C, fm.A, fm.LP, fm.G), digits = 2)
  
  fm.prd <- predict_nls(fm.L, fm.Q, fm.C, fm.A, fm.LP, fm.G)
  
  ggplot(data = barleyAS, aes(x = NF, y = yield)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm.G), color = "gam")) + 
    geom_line(aes(y = fitted(fm.C), color = "cubic")) + 
    geom_line(aes(y = Estimate, color = "simulate_lm")) + 
    geom_line(aes(y = fm.prd, color = "Avg. Model")) + 
    geom_ribbon(aes(ymin = Q5, ymax = Q95), fill = "purple", alpha = 0.3) +
    ggtitle("90% bands based on simulation")

}

### Testing predict2_nls and also using newdata with a function which is not an SS ----

if(run.predict.nls){
 
  require(ggplot2)
  require(nlme)
  data(Soybean)
  
  SoyF <- subset(Soybean, Variety == "F" & Year == 1988)
  fm1 <- nls(weight ~ SSlogis(Time, Asym, xmid, scal), data = SoyF)
  ## The SSlogis also supplies analytical derivatives
  ## therefore the predict function returns the gradient too
  prd1 <- predict(fm1, newdata = SoyF)
  
  ## Gradient
  head(attr(prd1, "gradient"))
  ## Prediction method using gradient
  prds <- predict2_nls(fm1, interval = "conf")
  SoyFA <- cbind(SoyF, prds)
  ggplot(data = SoyFA, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) +
    ggtitle("95% Confidence Bands")
  
  ### Without using a SS function
  getInitial(weight ~ SSlogis(Time, Asym, xmid, scal), data = SoyF)
  
  fm11 <- nls(weight ~ Asym / (1 + exp((xmid - Time)/scal)), 
              start = c(Asym = 21, xmid = 45, scal = 10),
              data = SoyF)
  
  SoyF2 <- subset(Soybean, Variety == "F" & Year == 1989)

  #### Using Monte Carlo method  
  prds2 <- predict_nls(fm11, interval = "conf", newdata = SoyF2)
  SoyFA2 <- cbind(SoyF2, prds2)
  ggplot(data = SoyFA2, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) +
    ggtitle("Newdata: 95% Confidence Bands")
  
  #### Using Delta method
  prds3 <- predict2_nls(fm11, interval = "conf", newdata = SoyF2)
  SoyFA3 <- cbind(SoyF2, prds3)
  ggplot(data = SoyFA3, aes(x = Time, y = weight)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) +
    ggtitle("Newdata: 95% Confidence Bands")

  #### Monte Carlo vs. Delta Method
  ggplot() +
    geom_point(aes(x = SoyFA2$Estimate, y = SoyFA3$Estimate)) + 
    xlab("Monte Carlo method") + 
    ylab("Delta method") + 
    geom_abline(intercept = 0, slope = 1) + 
    ggtitle("Estimates are identical as they should be")

  ggplot() +
    geom_point(aes(x = SoyFA2$Q2.5, y = SoyFA3$Q2.5)) + 
    xlab("Monte Carlo method") + 
    ylab("Delta method") + 
    geom_abline(intercept = 0, slope = 1) + 
    ggtitle("Lower bounds are similar")
  
  ggplot() +
    geom_point(aes(x = SoyFA2$Q97.5, y = SoyFA3$Q97.5)) + 
    xlab("Monte Carlo method") + 
    ylab("Delta method") + 
    geom_abline(intercept = 0, slope = 1) + 
    ggtitle("Upper bounds are similar")
  
  ggplot() + 
    geom_point(aes(x = SoyFA2$Time, y = SoyFA2$Q2.5, color = "Monte Carlo")) + 
    geom_point(aes(x = SoyFA3$Time, y = SoyFA3$Q2.5, color = "Delta method")) + 
    geom_point(aes(x = SoyFA2$Time, y = SoyFA2$Q97.5, color = "Monte Carlo")) + 
    geom_point(aes(x = SoyFA3$Time, y = SoyFA3$Q97.5, color = "Delta method")) + 
    geom_line(aes(x = SoyFA2$Time, y = SoyFA2$Q2.5, color = "Monte Carlo")) + 
    geom_line(aes(x = SoyFA3$Time, y = SoyFA3$Q2.5, color = "Delta method")) + 
    geom_line(aes(x = SoyFA2$Time, y = SoyFA2$Q97.5, color = "Monte Carlo")) + 
    geom_line(aes(x = SoyFA3$Time, y = SoyFA3$Q97.5, color = "Delta method")) + 
    xlab("Time") + ylab("weight")
  
  #### It appears that Monte Carlo is narrower so it might be underestimating
  #### the uncertainty
  #### Another example
  data(Orange)
  head(Orange)
  
  Orange1 <- subset(Orange, Tree == 1)  
  
  fm111 <- nls(circumference ~ Asym / (1 + exp((xmid - age)/scal)),
               start = c(Asym = 145, xmid = 922, scal = 200),
               data = Orange1)
  
  prds111 <- predict2_nls(fm111, interval = "conf")
  Orange1A <- cbind(Orange1, prds111)
  
  ggplot(data = Orange1A, aes(x = age, y = circumference)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm111))) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 1/3) + 
    ggtitle("Using the Delta method but no newdata")
  
  
  prds112 <- predict_nls(fm111, interval = "conf")
  Orange1A2 <- cbind(Orange1, prds112)
  
  ggplot(data = Orange1A2, aes(x = age, y = circumference)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm111))) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 1/3) + 
    ggtitle("Using the Monte Carlo method but no newdata")
  
  fgm <- gam(circumference ~ s(age, k = 7), data = Orange1)
  
  prds113 <- predict_gam(fgm, interval = "conf")
  Orange1A3 <- cbind(Orange1, prds113)
  
  fm <- lm(circumference ~ age + I(age^2), data = Orange1)
  
  prds114 <- predict_nls(fm, interval = "conf")
  Orange1A4 <- cbind(Orange1, prds114)
  
  ggplot() + 
    geom_point(aes(x = Orange1A2$age, y = Orange1A2$Q2.5, color = "Monte Carlo")) + 
    geom_point(aes(x = Orange1A$age, y = Orange1A$Q2.5, color = "Delta method")) + 
    geom_point(aes(x = Orange1A3$age, y = Orange1A3$Q2.5, color = "GAM")) + 
    geom_point(aes(x = Orange1A4$age, y = Orange1A4$Q2.5, color = "LM")) + 
    geom_point(aes(x = Orange1A2$age, y = Orange1A2$Q97.5, color = "Monte Carlo")) + 
    geom_point(aes(x = Orange1A$age, y = Orange1A$Q97.5, color = "Delta method")) + 
    geom_point(aes(x = Orange1A3$age, y = Orange1A3$Q97.5, color = "GAM")) + 
    geom_point(aes(x = Orange1A4$age, y = Orange1A4$Q97.5, color = "LM")) + 
    geom_line(aes(x = Orange1A2$age, y = Orange1A2$Q2.5, color = "Monte Carlo")) + 
    geom_line(aes(x = Orange1A$age, y = Orange1A$Q2.5, color = "Delta method")) + 
    geom_line(aes(x = Orange1A3$age, y = Orange1A3$Q2.5, color = "GAM")) + 
    geom_line(aes(x = Orange1A4$age, y = Orange1A4$Q2.5, color = "LM")) + 
    geom_line(aes(x = Orange1A2$age, y = Orange1A2$Q97.5, color = "Monte Carlo")) + 
    geom_line(aes(x = Orange1A$age, y = Orange1A$Q97.5, color = "Delta method")) + 
    geom_line(aes(x = Orange1A3$age, y = Orange1A3$Q97.5, color = "GAM")) + 
    geom_line(aes(x = Orange1A4$age, y = Orange1A4$Q97.5, color = "LM")) + 
    xlab("age") + ylab("circumference")
  
  ### With New data
  Orange2 <- subset(Orange, Tree == 2)
  
  prds121 <- predict2_nls(fm111, interval = "conf", newdata = Orange2)
  Orange2A <- cbind(Orange2, prds121)
    
  ggplot(data = Orange2A, aes(x = age, y = circumference)) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 1/3) + 
    ggtitle("Using the Delta method with newdata")
  
  prds122 <- predict_nls(fm111, interval = "conf")
  Orange2A2 <- cbind(Orange1, prds122)
  
  ggplot(data = Orange2A2, aes(x = age, y = circumference)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fm111))) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 1/3) + 
    ggtitle("Using the Monte Carlo method with newdata")
  
  prds123 <- predict_gam(fgm, interval = "conf", newdata = Orange2)
  Orange2A3 <- cbind(Orange2, prds123)
  
  prds124 <- predict_nls(fm, interval = "conf")
  Orange2A4 <- cbind(Orange2, prds124)
  
  ggplot() + 
    geom_point(aes(x = Orange2A2$age, y = Orange2A2$Q2.5, color = "Monte Carlo")) + 
    geom_point(aes(x = Orange2A$age, y = Orange2A$Q2.5, color = "Delta method")) + 
    geom_point(aes(x = Orange2A3$age, y = Orange2A3$Q2.5, color = "GAM")) + 
    geom_point(aes(x = Orange2A4$age, y = Orange2A4$Q2.5, color = "LM")) + 
    geom_point(aes(x = Orange2A2$age, y = Orange2A2$Q97.5, color = "Monte Carlo")) + 
    geom_point(aes(x = Orange2A$age, y = Orange2A$Q97.5, color = "Delta method")) + 
    geom_point(aes(x = Orange2A3$age, y = Orange2A3$Q97.5, color = "GAM")) + 
    geom_point(aes(x = Orange2A4$age, y = Orange2A4$Q97.5, color = "LM")) + 
    geom_line(aes(x = Orange2A2$age, y = Orange2A2$Q2.5, color = "Monte Carlo")) + 
    geom_line(aes(x = Orange2A$age, y = Orange2A$Q2.5, color = "Delta method")) + 
    geom_line(aes(x = Orange2A3$age, y = Orange2A3$Q2.5, color = "GAM")) + 
    geom_line(aes(x = Orange2A4$age, y = Orange2A4$Q2.5, color = "LM")) + 
    geom_line(aes(x = Orange2A2$age, y = Orange2A2$Q97.5, color = "Monte Carlo")) + 
    geom_line(aes(x = Orange2A$age, y = Orange2A$Q97.5, color = "Delta method")) + 
    geom_line(aes(x = Orange2A3$age, y = Orange2A3$Q97.5, color = "GAM")) + 
    geom_line(aes(x = Orange2A4$age, y = Orange2A4$Q97.5, color = "LM")) + 
    xlab("age") + ylab("circumference")
  
}
