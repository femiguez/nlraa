## test simulate_lme and simulate_gls
require(nlme)
require(nlraa)
require(ggplot2)

run.simulate.gls.test <- Sys.info()[["user"]] == "fernandomiguez"

if(run.simulate.gls.test){
  
  data(Orthodont)
  
  fm1 <- gls(distance ~ age, data = Orthodont,
             correlation = corCAR1(form = ~ age | Subject)) ## Sub by default
  
  ## Visualize correlation
  fm1.vc <- var_cov(fm1)
  fm1.vc2 <- fm1.vc[1:28, 1:28]
  image(log(fm1.vc2[,ncol(fm1.vc2):1]), 
        main = "Covariance (residual or marginal)",
        xaxt = "n", yaxt = "n")
  
  ## auto-correlation
  plot(ACF(fm1))
  
  sim1 <- simulate_lme(fm1, nsim = 2e2, value = "data.frame")

  sim1$Sub_ii <- paste0(sim1$Subject,"_",sim1$ii)
  
  ggplot(data = sim1, aes(age, distance )) + 
    geom_line(aes(y = sim.y, group = Sub_ii), color = "grey",  
              alpha = 0.3) + 
    geom_point() + ggtitle("Orthodont")
  
  ## Trying the ChickWeight dataset
  data(ChickWeight)
  
  fm2 <- lme(weight ~ Time + I(Time^2) + I(Time^3), 
             data = ChickWeight, weights = varPower(),
             random = ~ Time)
  
  ## Level = 0 is at the population-level           
  sim2 <- simulate_lme(fm2, nsim = 200, level = 0, value = "data.frame")
  
  sim2$Chick_ii <- paste0(sim2$Chick,"_",sim2$ii)
  
  ggplot(data = sim2, aes(Time, weight)) + 
    geom_line(aes(y = sim.y, group = Chick_ii), color = "grey",  
              alpha = 0.3) + 
    geom_point()
  
  ## Test with the Orange dataset
  data("Orange")
  
  fm3 <- lme(circumference ~ age, data = Orange,
             random = ~ 1 | Tree)
  
  ## Small data set, no issues
  plot(fm3)
  
  sim3 <- simulate_lme(fm3, nsim = 200, value = "data.frame")
  
  sim3$Tree_ii <- paste0(sim3$Tree,"_",sim3$ii)
  
  ggplot(data = sim3, aes(age, circumference)) + 
    geom_line(aes(y = sim.y, group = Tree_ii), color = "grey",  
              alpha = 0.3) + 
    geom_point()
  
  ## Testing scoping issues
  ## Script from bug report
  simulation_summary_fixed <- function(duration_per_period = 6L) {

    treatment <- as.factor(c(0,1))
    time <- 1:20
    cluster <- as.factor(c("a", "b", "c"))
    SWdesign <- expand.grid(treatment = treatment, time = time, cluster = cluster)
    SWdesign$id <- with(SWdesign, paste0(cluster, "_", treatment))
    SWdesign$outcome <- with(SWdesign, as.numeric(treatment) + 
                               (1 + as.numeric(treatment)) * time + 
                               scale(as.numeric(cluster)) * 0.5 + 
                               rnorm(nrow(SWdesign)))
    
    cor_init <- corAR1(form= ~ time|id )
    cor_init <- Initialize(cor_init, data = SWdesign )
    
    SWdesign$outcome <- rnorm( nrow(SWdesign) )
    
    local.gls <- gls( outcome ~ treatment*time, correlation=cor_init, data=SWdesign, method="ML") 
    local.lme <- lme( outcome ~ treatment*time, random = ~ 1 | cluster, 
                      data=SWdesign, method="ML") 

    ans1 <- simulate_lme(local.lme, psim=2, data = SWdesign)
    ans2 <- simulate_gls(local.gls, psim=2, data = SWdesign)
  }
  
  ## If this runs without errors it means the above works
  simulation_summary_fixed()

}