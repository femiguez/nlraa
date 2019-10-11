## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 7, fig.height = 6)
library(ggplot2)
library(nlraa)

## ----apropos-------------------------------------------------------------
apropos("^SS")

## ----sm------------------------------------------------------------------
## Sorghum and Maize dataset
data(sm)
ggplot(data = sm, aes(x = DOY, y = Yield, color = Crop)) + 
  geom_point() + 
  facet_wrap(~ Input)

## ----lfmc----------------------------------------------------------------
## Live fuel moisture content
data(lfmc)
ggplot(data = lfmc, aes(x = time, y = lfmc, color = leaf.type)) +
  geom_point() + 
  ylab("Live fuel moisture content (%)")

## ----swpg----------------------------------------------------------------
## Soil water and plant growth
data(swpg)
ggplot(data = swpg, aes(x = ftsw, y = lfgr)) +
  geom_point() + 
  xlab("Fraction Transpirable Soil Water") + 
  ylab("Relative Leaf Growth")

