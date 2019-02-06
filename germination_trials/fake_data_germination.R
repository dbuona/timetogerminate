###fake data for germination mdodels
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(drc)
library(tidyverse)
library(MCMCglmm)
library(shinystan)
library(rstan)
options(mc.cores = parallel::detectCores())

setwd("~/Documents/git/timetogerminate/germination_trials")

##use lognoral to simulate the data
set.seed(613613)


###This make data folowing the log logistic function
germ<-function(t,d,b,t50){
y<- d/(1+((t/t50)^b))
return(data.frame(time=t, y=y)) 
  
}

germ(seq(0,24,by=3),rtnorm(1,12,1,lower=8,upper=20),rnorm(1,-5,0.1),rnorm(1,15,1))


##3 petridishes of the same treatment, there is probably a loop or apply function for this

A<-germ(seq(0,24,by=3),rtnorm(1,12,1,lower=8,upper=20),rnorm(1,-5,0.5),rnorm(1,15,2))
A$dish<-"A"
B<-germ(seq(0,24,by=3),rtnorm(1,12,1,lower=8,upper=20),rnorm(1,-5,0.5),rnorm(1,15,2))
B$dish<-"B"
C<-germ(seq(0,24,by=3),rtnorm(1,12,1,lower=8,upper=20),rnorm(1,-5,0.5),rnorm(1,15,2))
C$dish<-"C"

#make your data
d<-rbind(A,B,C)
d$y<-round(d$y)
  
mod<-drm(y~time,fct=LL.3(),data=d,type="continuous")
summary(mod)
plot(mod,ylim=c(0,20),xlim=c(0,24),log="",pch=16,type="all")


###now try it in stan

data.list <- with(d, 
                  list(Y=y, 
                    t = time, 
                    N = nrow(d)
                  ))

germ.mod = stan('stan/fakeseedmodel.stan', data = data.list,
                  iter = 2500, warmup=1500)

mod.sum<-summary(germ.mod)$summary
mod.sum[c("d", "beta", "t50", "sigma"),]


G<-germ(seq(0,24,by=3),rtnorm(1,18,1,lower=8,upper=20),rnorm(1,-6,0.5),rnorm(1,10,1))
G$dish<-"G"
H<-germ(seq(0,24,by=3),rtnorm(1,18,1,lower=8,upper=20),rnorm(1,-6,0.5),rnorm(1,10,1))
H$dish<-"H"
I<-germ(seq(0,24,by=3),rtnorm(1,18,1,lower=8,upper=20),rnorm(1,-6,0.5),rnorm(1,10,1))
I$dish<-"I"
dd<-rbind(G,H,I)
dd$chill<-1
dd$y<-round(dd$y)

d$chill<-0

dat<-rbind(d,dd)

mod2<-drm(y~time,factor(chill),fct=LL.3(),data=dat,type="continuous")
summary(mod2)
plot(mod2,ylim=c(0,20),xlim=c(0,24),log="",pch=16,type="all")



data.list2 <- with(dat, 
                  list(Y=y, 
                       t = time,
                       chill=chill,
                       N = nrow(dat)
                  )
)

germ.mod2 = stan('stan/fakeseed_wchill.stan', data = data.list2,
                iter = 3000, warmup=2000) #this model is fishy
mod.sum2<-summary(germ.mod2)$summary
mod.sum2[c("d", "beta", "t50", "sigma",
           "b_chill_beta","b_chill_t50","b_chill_d"),]
