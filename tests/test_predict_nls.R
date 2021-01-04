## Example of a bad model vs. uncertainty vs. model averaging
require(nlraa)
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

