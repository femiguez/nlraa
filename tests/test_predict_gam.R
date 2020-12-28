## Testing predict_gam
require(nlraa)
library(mgcv)
require(minpack.lm)
require(ggplot2)

if(Sys.info()[["user"]] == "fernandomiguez"){
 
  data("maizeleafext") 
  
  ggplot(data = maizeleafext, aes(x = temp, y = rate)) + 
    geom_point()
  
  fmg <- gam(rate ~ temp + s(temp, k = 9), data = maizeleafext)
  fmb <- nlsLM(rate ~ SSbeta5(temp, mu, tb, a, tc, b), data = maizeleafext)
  ##fmr <- nlsLM(rate ~ SSratio(temp, a, b, c, d), data = maizeleafext)
  
  IC_tab(fmg, fmb)
  
  fmgp <- predict_gam(fmg, interval = "conf")
  fmbp <- predict_gam(fmb, interval = "conf")
  
  maizeleafextAG <- cbind(maizeleafext, fmgp)
  maizeleafextAB <- cbind(maizeleafext, fmbp)
  
  ggplot() + 
    geom_point(data = maizeleafextAG, aes(x = temp, y = rate)) + 
    geom_line(data = maizeleafextAG, aes(x = temp, y = Estimate, color = "GAM")) + 
    geom_line(data = maizeleafextAB, aes(x = temp, y = Estimate, color = "Beta5")) + 
    geom_ribbon(data = maizeleafextAG, 
                aes(x = temp, ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3) + 
    geom_ribbon(data = maizeleafextAB, 
                aes(x = temp, ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3)
  
  data(swpg)
  
  ggplot(data = swpg, aes(x = ftsw, y = lfgr)) + 
    geom_point()

  fsw.G <- gam(lfgr ~ ftsw + s(ftsw), data = swpg)    
  fsw.LP <- nls(lfgr ~ SSlinp(ftsw, a, b, xs), data = swpg)
  fsw.A <- nls(lfgr ~ SSasymp(ftsw, Asym, R0, lrc), data = swpg)
  fsw.C <- lm(lfgr ~ poly(ftsw, 3), data = swpg)
  
  IC_tab(fsw.G, fsw.LP, fsw.A, fsw.C, criteria = "BIC")

  prd <- predict_nls(fsw.G, fsw.LP, fsw.A, fsw.C, criteria = "BIC", interval = "conf")
  
  swpgA <- cbind(swpg, prd)
  
  ggplot(data = swpgA, aes(x = ftsw, y = lfgr)) + 
    geom_point() + 
    geom_line(aes(y = fitted(fsw.G), color = "GAM")) + 
    geom_line(aes(y = fitted(fsw.LP), color = "LP")) + 
    geom_line(aes(y = fitted(fsw.A), color = "A")) +
    geom_line(aes(y = fitted(fsw.C), color = "C")) +
    geom_line(aes(y = Estimate, color = "Avg. model"), size = 1.2) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.3)
  
  data(sm)

  fsm.G <- gam(Yield ~ DOY*Crop*Input + s(DOY), data = sm)
  
  fsm.Gp <- predict_gam(fsm.G, interval = "conf")
  smAG <- cbind(sm, fsm.Gp)
  
  ggplot(data = smAG, aes(x = DOY, y = Yield, color = Crop)) + 
    facet_wrap(~Input) + 
    geom_point() + 
    geom_line(aes(y = Estimate)) + 
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = Crop, color = NULL), alpha = 0.3)
  
}