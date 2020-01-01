## Testing all SS functions in nlraa package
require(nlraa)

## 1. SSbgf
## This function won't be tested here as it is used extensively in the 
## vignette nlraa::nlraa-AgronJ-paper
## However, in the future I will create a version reparameterized in 
## terms of unconstrained parameters, because the condition t.m < t.e 
## is not guranteed

## 2. SSbgf4 
## It would also be beneficial to reparameterize this function
## for routine work in terms of unconstrained parameters
data(sm)
#' ## Let's just pick one crop
sm2 <- subset(sm, Crop == "M")
## For this particular problem it is easier to 'fix' t.b and w.b
fit <- nls(Yield ~ bgf2(DOY, w.max, w.b = 0, t.e, t.m, t.b = 141), 
            data = sm2, start = list(w.max = 16, t.e= 240, t.m = 200))

## 3. SSdlf 
## Extended example in vignette 'nlraa-Oddi-LFMC'

## 4. SSricker
set.seed(123)
x <- 1:30
y <- 30 * x * exp(-0.3 * x) + rnorm(30, 0, 0.25)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSricker(x, a, b), data = dat)
confint(fit)

## 5. SSprofd
## I'm not a huge fan of this function as I think the dlf is more stable
set.seed(1234)
x <- 1:10
y <- profd(x, 0.3, 0.05, 0.5, 4) + rnorm(10, 0, 0.01)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSprofd(x, a, b, c, d), data = dat)

## 6. SSlinp
x <- 1:30
y <- linp(x, 0, 1, 20) + rnorm(30, 0, 0.5)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSlinp(x, a, b, xs), data = dat)
confint(fit)

## 7. SSplin
set.seed(123)
x <- 1:30
y <- plin(x, 10, 20, 1) + rnorm(30, 0, 0.5)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSplin(x, a, xs, b), data = dat)
confint(fit)

## 8. SSquadp
set.seed(123)
x <- 1:30
y <- quadp(x, 5, 1.7, -0.04, 20) + rnorm(30, 0, 0.6)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSquadp(x, a, b, c, xs), data = dat, algorithm = "port")
## It's strange but confint will return NAs unless level is 0.5
confint(fit, level = 0.5)

## 9. SSpquad
set.seed(12345)
x <- 1:40
y <- pquad(x, 5, 20, 1.7, -0.04) + rnorm(40, 0, 0.6)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSpquad(x, a, xs, b, c), data = dat)
confint(fit)

## 10. SSblin
set.seed(123)
x <- 1:30
y <- blin(x, 0, 0.75, 15, 1.75) + rnorm(30, 0, 0.5)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSblin(x, a, b, xs, c), data = dat)
## Confidence intervals do not work for this example...?

## 11. SSexpf
set.seed(1234)
x <- 1:15
y <- expf(x, 10, -0.3) + rnorm(15, 0, 0.2)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSexpf(x, a, c), data = dat)
confint(fit)

## 12. SSexpfp
set.seed(12345)
x <- 1:30
y <- expfp(x, 10, 0.1, 15) + rnorm(30, 0, 1.5)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSexpfp(x, a, c, xs), data = dat)
confint(fit)

## 13. SSpexpf
set.seed(1234)
x <- 1:30
y <- pexpf(x, 20, 15, -0.2) + rnorm(30, 0, 1)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSpexpf(x, a, xs, c), data = dat)

## 14. SSbgrp
x <- 1:30
y <- bgrp(x, 20, log(25), log(5)) + rnorm(30, 0, 1)
dat <- data.frame(x = x, y = y)
fit <- nls(y ~ SSbgrp(x, w.max, lt.e, ldt), data = dat)
exp(confint(fit)[2:3,])
