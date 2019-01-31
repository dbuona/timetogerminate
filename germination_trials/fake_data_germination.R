###fake data for germination mdodels
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(drc)
library(tidyverse)
##use lognoral to simulate the data
set.seed(613613)


time<- seq(0, 25, by=3)

germ<-round(rlnorm(9, meanlog = 0, sdlog = 1))
cum.germ<-cumsum(germ)
goo<-data.frame(cbind(time,germ,cum.germ))
goo$unit<-"A"

mod<-drm(cum.germ~time,data=goo,fct = LL.3(), type ="continuous")
summary(mod)
plot(mod)
