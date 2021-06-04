## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 6)
library(ggplot2)
library(nlraa)

## ----apropos------------------------------------------------------------------
apropos("^SS")

## ----sm-----------------------------------------------------------------------
## Sorghum and Maize dataset
data(sm)
ggplot(data = sm, aes(x = DOY, y = Yield, color = Crop)) + 
  geom_point() + 
  facet_wrap(~ Input)

## ----lfmc---------------------------------------------------------------------
## Live fuel moisture content
data(lfmc)
ggplot(data = lfmc, aes(x = time, y = lfmc, color = leaf.type)) +
  geom_point() + 
  ylab("Live fuel moisture content (%)")

## ----swpg---------------------------------------------------------------------
## Soil water and plant growth
data(swpg)
ggplot(data = swpg, aes(x = ftsw, y = lfgr)) +
  geom_point() + 
  xlab("Fraction Transpirable Soil Water") + 
  ylab("Relative Leaf Growth")

## ----barley-------------------------------------------------------------------
## Response of barley to nitrogen fertilizer
## There is a barley dataset also in package 'lattice'
data(barley, package = "nlraa")
ggplot(data = barley, aes(x = NF, y = yield, color = as.factor(year))) +
  geom_point() +
  xlab("Nitrogen fertilizer (g/m^2)") +
  ylab("Grain (g/m^2)")

## ----maizeleafext-------------------------------------------------------------
## Response of barley to nitrogen fertilizer
## There is a barley dataset also in package 'lattice'
data(maizeleafext, package = "nlraa")
ggplot(data = maizeleafext, aes(x = temp, y = rate)) +
  geom_point() + geom_line() + 
  xlab("Temperature (C)") +
  ylab("Leaf Extension Rate (relative)")

## ---- eval = FALSE------------------------------------------------------------
#  ## Error in nls(y ~ SSratio(x, a, b, c, d), data = dat) :
#  ##  step factor 0.000488281 reduced below 'minFactor' of 0.000976562

## ---- eval = FALSE------------------------------------------------------------
#  ## Error in qr.default(.swts * gr) :
#  ##  NA/NaN/Inf in foreign function call (arg 1)

## ----barleyG------------------------------------------------------------------
library(nlme)
data(barley, package = "nlraa")
barley$yearf <- as.factor(barley$year)
barleyG <- groupedData(yield ~ NF | yearf, data = barley)

## ----barleyG-mixed------------------------------------------------------------
## Fit the nonlinear model for each year
fit.nlis <- nlsList(yield ~ SSasymp(NF, Asym, R0, lrc), data = barleyG)
## Use this to fit a nonlinear mixed model
fit.nlme <- nlme(fit.nlis)
## Investigate residuals
plot(fit.nlme)
## Look at predictions
plot(augPred(fit.nlme, level = 0:1))
## Compute confidence intervals
intervals(fit.nlme)
## A simpler model is possible...

