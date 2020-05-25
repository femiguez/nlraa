## Test for the var_cov function
require(nlme)
require(nlraa)
## Tip: the image function can be used to visualize these matrices
## image(x[,ncol(x):1])

run.test.var.cov <- FALSE

if(run.test.var.cov){

  data(ChickWeight)

  fm0 <- lm(weight ~ Time, data = ChickWeight)
  vm0 <- var_cov(fm0)

  ## First model with no modeling of the Variance-Covariance
  fit0 <- gls(weight ~ Time, data = ChickWeight)
  v0 <- var_cov(fit0)
  ## getVarCov cannot handle this case

  ## Only modeling the diagonal (weights)
  fit1 <- gls(weight ~ Time, data = ChickWeight, weights = varPower())
  v1 <- var_cov(fit1)
  ## getVarCov cannot handle this case

  ## Only the correlation structure is defined and there are no groups
  fit2 <- gls(weight ~ Time, data = ChickWeight, correlation = corAR1())
  v2 <- var_cov(fit2)
  ## getVarCov cannot handle this case

  ## Only the correlation structure is defined and there are groups present
  fit3 <- gls(weight ~ Time, data = ChickWeight, correlation = corCAR1(form = ~ Time | Chick))
  v3 <- var_cov(fit3)
  ## getVarCov CAN handle this case, but it returns only the one for the first subject

  ## There are both weights and correlations
  fit4 <- gls(weight ~ Time, data = ChickWeight, 
               weights = varPower(),
               correlation = corCAR1(form = ~ Time | Chick))
  v4 <- var_cov(fit4)
  ## getVarCov cannot handle this case

  ## Note: I get a halving factor error when I try to use SSlogis on the
  ## ChickWeight data
  ## What about fitting a gnls?
  data(CO2)
  fitn0 <- gnls(uptake ~ SSasymp(conc, Asym, lrc, c0), data = CO2)
  vn0 <- var_cov(fitn0)
  ## getVarCov cannot handle this case

  fitn1 <- gnls(uptake ~ SSasymp(conc, Asym, lrc, c0), data = CO2,
                weights = varPower())
  vn1 <- var_cov(fitn1)
  ## getVarCov cannot handle this case

  fitn2 <- gnls(uptake ~ SSasymp(conc, Asym, lrc, c0), data = CO2,
                correlation = corAR1())
  vn2 <- var_cov(fitn2)
  ## getVarCov cannot handle this case

  fitn3 <- gnls(uptake ~ SSasymp(conc, Asym, lrc, c0), data = CO2,
                correlation = corCAR1(form = ~ conc | Plant))
  vn3 <- var_cov(fitn3)
  ## getVarCov CAN handle this case, but it returns the matrix for just the
  ## first subject by default

  fitn4 <- gnls(uptake ~ SSasymp(conc, Asym, lrc, c0), data = CO2,
                weights = varPower(),
                correlation = corCAR1(form = ~ conc | Plant))
  vn4 <- var_cov(fitn4)
  ## getVarCov returns an empty object in this case

  ## Testing an nlme object
  ## getVarCov cannot handle any of these cases
  fitnm0 <- nlsList(uptake ~ SSasymp(conc, Asym, lrc, c0), data = CO2)

  fitnm1 <- nlme(fitnm0, random = pdDiag(Asym + lrc + c0 ~ 1))
  vnm1 <- var_cov(fitnm1)

  fitnm2 <- update(fitnm1, weights = varPower())
  vnm2 <- var_cov(fitnm2)

  fitnm3 <- update(fitnm2, correlation = corCAR1(form = ~ conc))
  vnm3 <- var_cov(fitnm3)

  ## Test for whether my var_cov code agrees with getVarCov
  ChickWeight <- ChickWeight[order(ChickWeight$Chick),]
  fitcw.lme <- gls(weight ~ Time, data = ChickWeight, 
                   correlation = corAR1(form = ~ Time | Chick))

  ## This is per individual and it fails if I try to get both
  gt.vc.1 <- getVarCov(fitcw.lme, 1, type = "marginal")
  gt.vc.2 <- getVarCov(fitcw.lme, 2, type = "marginal")
  ## The same as
  vc <- var_cov(fitcw.lme)
  vc[1:9,1:9]
  
  vc2 <- vc[1:50,1:50]
  ## This looks nice
  image(vc2[,ncol(vc2):1])    
 
  fitcw.lme2 <- gls(weight ~ Time, data = ChickWeight, 
                    weights = varPower(),
                    correlation = corAR1(form = ~ Time | Chick))

  ## This is per individual and it fails if I try to get both
  ## This code returns errors or warnings
  ## gt2.vc.1 <- getVarCov(fitcw.lme2, 1, type = "marginal")
  ## gt2.vc.2 <- getVarCov(fitcw.lme2, 2, type = "marginal")
    
  vc2.2 <- var_cov(fitcw.lme2)
  vc2.2[1:9,1:9]
  
  vc2.3 <- vc2.2[1:50, 1:50]
  image(vc2.3[,ncol(vc2.3):1])       
}