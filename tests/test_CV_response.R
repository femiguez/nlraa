## Fitting nlme
library(nlme)
set.seed(101)

dat <- expand.grid(time = 1:10, rep = 1:4, trt = letters[1:6])
dat$eu <- with(dat, paste0(trt,"_",rep))

dat$y <- c(replicate(24, SSlogis(1:10, 10, 5, 1.5) + rnorm(10, sd = 0.25)))

datG <- groupedData(y ~ time | eu, data = dat)

plot(datG)

fitL <- nlsList(SSlogis, datG)

fnm1 <- nlme(fitL, random = pdDiag(Asym + xmid + scal ~ 1))

fxf <- fixef(fnm1)

fnm2 <- update(fnm1, fixed = Asym + xmid + scal ~ trt, 
               start = c(fxf[1], rep(0, 5), 
                         fxf[2], rep(0, 5),
                         fxf[3], rep(0, 5)))

anova(fnm2)

plot(fnm2)

