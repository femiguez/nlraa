## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 6)
library(lattice)
library(nlme)
library(ggplot2)
library(nlraa)

## ----strsm--------------------------------------------------------------------
data(sm)
str(sm)
head(sm)

## ----sm-ggplot, echo = FALSE--------------------------------------------------
ggplot(data = sm, aes(y = Yield, x = DOY)) +
   facet_grid(. ~ Input) +
   geom_point(aes(fill=Crop, shape=Crop), size=2) +
   scale_shape_manual(values=c(24,21,1)) +
   scale_fill_manual(values = c("grey","black","black")) +
   scale_x_continuous("Day of the Year") +
   scale_y_continuous("Dry biomass (Mg/ha)") +
   theme_bw()

## ----create-eu----------------------------------------------------------------
sm$eu <- with(sm, factor(Block):factor(Input):factor(Crop))
sm2 <- subset(sm, DOY != 141)

## ----grouped-data-------------------------------------------------------------
smG <- groupedData(Yield ~ DOY | eu, data = sm2)

## ----nls-list-sm--------------------------------------------------------------
fit.lis <- nlsList(Yield ~ SSbgf(DOY, w.max, t.e, t.m), data = smG)
## But this works better
## Added 2020/1/2
fit.lis.rp <- nlsList(Yield ~ SSbgrp(DOY, w.max, lt.e, ldt), data = smG) 

## ----nls-list-plot, echo = FALSE----------------------------------------------
print(plot(fit.lis))
print(plot(intervals(fit.lis)))

## ----relax-control------------------------------------------------------------
fit.me <- nlme(fit.lis, control = list(maxIter = 100, msMaxIter = 300, pnlsMaxIter = 20))

## ----plot-resid-nlme, echo = FALSE--------------------------------------------
print(plot(fit.me))
print(plot(augPred(fit.me, level = 0:1)))

## ----bgf2---------------------------------------------------------------------
fit.lis2 <- nlsList(Yield ~ bgf2(DOY, w.max, w.b = 0, t.e, t.m, t.b = 141),
                    data = smG,
                    start = c(w.max = 30, t.e=280, t.m=240))

## ----plot-bgf2, echo = FALSE--------------------------------------------------
print(plot(fit.lis2))

## ----nlme-update--------------------------------------------------------------
fit.me2 <- nlme(fit.lis2)
## Error message, but the next model is the one we care about
fit2.me2 <- update(fit.me2, random = pdDiag(w.max + t.e + t.m ~ 1))
anova(fit.me2, fit2.me2)
## The second model is simpler and it seems to be marginally better than 
## the orginial, but we need to keep in mind that the simpler model
## converges much more easily

## ----nlme-update-two----------------------------------------------------------
fe <- fixef(fit2.me2) ## Some starting values with visual help
fit3.me2 <- update(fit2.me2, fixed = list(w.max + t.e + t.m ~ Crop),
                  start = c(fe[1], -10, 20, fe[2], -40, 0, fe[3], -40, 0))
## We next include the Input
fe2 <- fixef(fit3.me2)
fit4.me2 <- update(fit3.me2, fixed = list(w.max + t.e + t.m
                               ~ Crop + Input),
                  start = c(fe2[1:3], 0, fe2[4:6], 0, fe2[7:9], 0))
## and the interaction
fe3 <- fixef(fit4.me2)
fit5.me2 <- update(fit4.me2,
                   fixed = list(w.max + t.e + t.m
                     ~ Crop + Input + Crop:Input),
                  start = c(fe3[1:4], 0, 0,
                            fe3[5:8], 0, 0,
                            fe3[9:12], 0, 0))

## ----fit5-plot, echo = FALSE--------------------------------------------------
print(plot(fit5.me2))

## ----fit6-and-fit7------------------------------------------------------------
fit6.me2 <- update(fit5.me2,
                   weights = varPower(form = ~ fitted(.) | Crop))

fit7.me2 <- update(fit6.me2, weights = varPower(form = ~ fitted(.)))

anova(fit6.me2, fit7.me2)

## ----fit6.me2-----------------------------------------------------------------
fit6.me2

## ----gnls---------------------------------------------------------------------
## Random effects are almost zero
fit8.me2 <- gnls(Yield ~ bgf2(DOY, w.max, t.e, t.m, w.b=0, t.b=141),
                 data = smG,
                 params = list(w.max + t.e + t.m ~ Crop + Input
                                                   + Crop:Input),
                 weights = varPower(form = ~ fitted(.) | Crop),
                 start = fixef(fit7.me2))
anova(fit6.me2, fit8.me2)

## ----anova-fit8---------------------------------------------------------------
anova(fit8.me2)

## ----plot-fit8----------------------------------------------------------------
print(plot(fit8.me2))

## ----fit8-fitted--------------------------------------------------------------
smG$prds <- fitted(fit8.me2)

doys <- 168:303
ndat <- expand.grid(DOY=doys, Crop= unique(smG$Crop), Input=c(1,2))
ndat$preds <- predict(fit8.me2, newdata = ndat)

## Here I'm just removing prediction for maize that go beyond
## day of the year 270
ndat2 <- ndat
ndat2[ndat2$Crop == "M" & ndat2$DOY > 270,"preds"] <- NA
ndat2 <- na.omit(ndat2)

## ----fit8-fitted-plot, echo = FALSE-------------------------------------------
 ggplot(data = smG, aes(y = Yield, x = DOY)) +
  facet_grid(. ~ Input) +
   geom_point(aes(fill=Crop, shape=Crop), size=2) +
   geom_line(aes(x = DOY, y = preds, linetype = Crop), data=ndat2) +
   scale_shape_manual(values=c(24,21,1)) +
   scale_fill_manual(values = c("grey","black","black")) +
   scale_x_continuous("Day of the Year") +
   scale_y_continuous("Dry biomass (Mg/ha)") +
   theme_bw()

