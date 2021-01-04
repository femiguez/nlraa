## For this test I'm planning to use the lfmc dataset
## The goals is to test both simulate and bootstrap
require(nlme)
require(nlraa)
require(ggplot2)

run.sim.boot.lfmc <- Sys.info()[["user"]] == "fernandomiguez" && FALSE

if(run.sim.boot.lfmc){
  
  data(lfmc, package = "nlraa")
  
  lfmc <- droplevels(subset(lfmc, leaf.type != "Grass E"))

  lfmcG <- groupedData(lfmc ~ time | group, data = lfmc)
  
  fitL <- nlsList(lfmc ~ SSdlf(time, asym, a2, xmid, scal), data = lfmcG)  

  plot(fitL)
  plot(intervals(fitL))

  fm <- nlme(fitL, random = pdDiag(asym + a2 + xmid + scal ~ 1))
  
  fm1 <- update(fm, random = pdDiag(asym + xmid ~ 1))
  ## I'll ignore site for this example
  fxf <- fixef(fm1)
  
  fm2 <- update(fm1, fixed = list(asym + a2 + xmid + scal ~ leaf.type),
                weights = varPower(),
                start = c(fxf[1], rep(0, 2), fxf[2], rep(0, 2), fxf[3], rep(0, 2), fxf[4], rep(0, 2)))
  
  anova(fm2)
  
  ## Simulate from this model 
  ## psim = 1, level = 0
  ## This is at the level of the mean function for each level of factor 'leaf.type'
  sim10 <- simulate_nlme(fm2, nsim = 100, level = 0, value = "data.frame") 

  ggplot(data = sim10, aes(x = time, y = lfmc)) + 
    facet_wrap(~ leaf.type) + 
    geom_line(aes(x = time, y = sim.y, group = ii, color = leaf.type)) + 
    geom_point()

  ## This is at the level of each 'group' 
  sim11 <- simulate_nlme(fm2, nsim = 100, level = 1, value = "data.frame") 
  
  sim11$ID <- with(sim11, paste0(ii,"_",group))

  ## This is interesting because it shows why we need bootstrapping, right?                   
  ggplot(data = sim11, aes(x = time, y = lfmc)) + 
    facet_wrap(~ leaf.type) + 
    geom_line(aes(x = time, y = sim.y, group = ID, color = leaf.type)) + 
    geom_point()

  ## This is at the level of observation 
  ## We can't go deeper than level 1 (fm2$dims$Q)
  sim21 <- simulate_nlme(fm2, nsim = 100, psim = 2, level = 1, value = "data.frame") 
  
  ggplot(data = sim21, aes(x = time, y = lfmc)) + 
    facet_wrap(~ leaf.type) + 
    geom_point(aes(x = time, y = sim.y, color = leaf.type)) + 
    geom_point()
  
  ## Bootstrapping section
  prd_fun <- function(x) predict(x, level = 0)
  ## This takes about 111 seconds
  system.time(bprd <- boot_nlme(fm2, prd_fun, cores = 4))
  
  lfmcG$lwr.q <- apply(t(bprd$t), 1, quantile, probs = 0.05, na.rm = TRUE)
  lfmcG$upr.q <- apply(t(bprd$t), 1, quantile, probs = 0.95, na.rm = TRUE)

  lfmcG$prd <- predict(fm2, level = 0)
  
  ggplot(data = lfmcG, aes(x = time, y = lfmc)) + 
    facet_wrap(~ leaf.type) + 
    geom_ribbon(aes(x = time, ymin = lwr.q, ymax = upr.q, fill = leaf.type), alpha = 0.5) + 
    geom_point() + 
    geom_line(aes(y = prd)) + 
    ggtitle("Bootstrapped fits at the fixed level of species")
  
}