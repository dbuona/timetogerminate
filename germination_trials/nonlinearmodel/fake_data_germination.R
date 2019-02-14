###fake data for germination mdodels
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(drc)
library(tidyverse)
library(MCMCglmm)
library(shinystan)
library(rstan)
library(rstanarm)
library(bayesplot)
options(mc.cores = parallel::detectCores())

setwd("~/Documents/git/timetogerminate/germination_trials")

##use lognoral to simulate the data
set.seed(613613)



###This make data folowing the log logistic function
germ<-function(t,d,b,t50){
y<- d/(1+((t/t50)^b))
return(data.frame(time=t, y=y))
  } 


##no chill, cool inc

A<-c(rep("1",9),c(rep("2",9),rep("3",9),rep("4",9),rep("5",9),rep("6",9),rep("7",9),rep("8",9),rep("9",9),rep("10",9))) ###generate replicate ids/time tratment

###it seems apply functions cant take paramenters from the envirnment

##same parements in all functions
time<-seq(0,24,by=3)
d.low<-0 
d.up<-20
#code: paremeter.chill(0-90).inc(0,1).

d.0.0.mu<-4
d.0.0.sd<-1
t50.0.0.mu<-20
t50.0.0.sd<-3
beta.0.0.mu<--3
beta.0.0.sd<-0.2

d.0.1.mu<-5
d.0.1.sd<-1
t50.0.1.mu<-18
t50.0.1.sd<-3
beta.0.1.mu<--3.5
beta.0.1.sd<-0.2


d.14.0.mu<-8
d.14.0.sd<-1
t50.14.0.mu<-15
t50.14.0.sd<-3
beta.14.0.mu<--4
beta.14.0.sd<-0.2

d.14.1.mu<-9
d.14.1.sd<-1
t50.14.1.mu<-13
t50.14.1.sd<-3
beta.14.1.mu<--4.5
beta.14.1.sd<-0.2

d.28.0.mu<-12
d.28.0.sd<-1
t50.28.0.mu<-11
t50.28.0.sd<-3
beta.28.0.mu<--5
beta.28.0.sd<-0.2

d.28.1.mu<-13
d.28.1.sd<-1
t50.28.1.mu<-10
t50.28.1.sd<-3
beta.28.1.mu<--5.5
beta.28.1.sd<-0.2

d.42.0.mu<-16
d.42.0.sd<-1
t50.42.0.mu<-8
t50.42.0.sd<-2
beta.42.0.mu<--6
beta.42.0.sd<-0.2

d.42.1.mu<-17
d.42.1.sd<-1
t50.42.1.mu<-7
t50.42.1.sd<-2
beta.42.1.mu<--6.5
beta.42.1.sd<-0.2

c0.0<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.0.0.mu,d.0.0.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.0.0.mu,t50.0.0.sd))^rnorm(1,beta.0.0.mu,beta.0.0.sd)))
  return(data.frame(time=time, y=y))
})
c0.0<-do.call(rbind.data.frame, c0.0)
c0.0$Chill<-0
c0.0$inc<-0
c0.0$ID<-A

c0.1<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.0.1.mu,d.0.1.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.0.1.mu,t50.0.1.sd))^rnorm(1,beta.0.1.mu,beta.0.1.sd)))
  return(data.frame(time=time, y=y))
})
c0.1<-do.call(rbind.data.frame, c0.1)
c0.1$Chill<-0
c0.1$inc<-1
c0.1$ID<-A

c14.0<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.14.0.mu,d.14.0.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.14.0.mu,t50.14.0.sd))^rnorm(1,beta.14.0.mu,beta.14.0.sd)))
  return(data.frame(time=time, y=y))
})
c14.0<-do.call(rbind.data.frame, c14.0)
c14.0$Chill<-14
c14.0$inc<-0
c14.0$ID<-A

c14.1<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.14.1.mu,d.14.1.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.14.1.mu,t50.14.1.sd))^rnorm(1,beta.14.1.mu,beta.14.1.sd)))
  return(data.frame(time=time, y=y))
})
c14.1<-do.call(rbind.data.frame, c14.1)
c14.1$Chill<-14
c14.1$inc<-1
c14.1$ID<-A



c28.0<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.28.0.mu,d.28.0.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.28.0.mu,t50.28.0.sd))^rnorm(1,beta.28.0.mu,beta.28.0.sd)))
  return(data.frame(time=time, y=y))
})
c28.0<-do.call(rbind.data.frame, c28.0)
c28.0$Chill<-28
c28.0$inc<-0
c28.0$ID<-A


c28.1<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.28.1.mu,d.28.1.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.28.1.mu,t50.28.1.sd))^rnorm(1,beta.28.1.mu,beta.28.1.sd)))
  return(data.frame(time=time, y=y))
})
c28.1<-do.call(rbind.data.frame, c28.1)
c28.1$Chill<-28
c28.1$inc<-1
c28.1$ID<-A

c42.0<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.42.0.mu,d.42.0.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.42.0.mu,t50.42.0.sd))^rnorm(1,beta.42.0.mu,beta.42.0.sd)))
  return(data.frame(time=time, y=y))
})
c42.0<-do.call(rbind.data.frame, c42.0)
c42.0$Chill<-42
c42.0$inc<-0
c42.0$ID<-A


c42.1<-lapply(A,function(d,t,b,t50){
  y<- rtnorm(1,d.42.1.mu,d.42.1.sd,lower=d.low,upper=d.up)/(1+((time/rnorm(1,t50.42.1.mu,t50.42.1.sd))^rnorm(1,beta.42.1.mu,beta.42.1.sd)))
  return(data.frame(time=time, y=y))
})
c42.1<-do.call(rbind.data.frame, c42.1)
c42.1$Chill<-42
c42.1$inc<-1
c42.1$ID<-A

test.dat<-rbind(c0.0,c0.1,c14.0,c14.1,c28.0,c28.1,c42.0,c42.1)

test.dat$y<-round(test.dat$y)
  
mod3<-drm(y~time,factor(Chill):factor(inc),fct=LL.3(),data=test.dat,type="continuous")
summary(mod3)
plot(mod3,ylim=c(0,20),xlim=c(0,24),log="",pch=16)


#make sure inc is a factor
#test.dat$inc<-as.factor(test.dat$inc)

data.list3 <- with(test.dat, 
                   list(Y=y, 
                        t = time,
                        chill=Chill,
                        inc=inc,
                        N = nrow(test.dat)
                   )
)

germ.mod3 = stan('stan/fakeseed_wchillandinc.stan', data = data.list3,
                 iter = 3000, warmup=2000) 

launch_shinystan(germ.mod3)
mod.sum3<-summary(germ.mod3)$summary
mod.sum3[c("a_beta","a_t50","a_d", "b_chill_beta","b_chill_t50","b_chill_d","b_inc_beta","b_inc_t50","b_inc_d","sigma"),]

cL<-7
cH<-14
iL<-20
iH<-25

### beta
9.0780587+cL* 0.3670658+iL*0.1044341
9.0780587+cL* 0.3670658+iH*0.1044341
9.0780587+cH* 0.3670658+iL*0.1044341 #underestimiating
9.0780587+cH* 0.3670658+iH*0.1044341 ##underestimating

germ.mod4 = stan('stan/fakeseed_winters.stan', data = data.list3,
                 iter = 3000, warmup=2000) ### dosent converge, but maybe cause data list is too small
#Chain 4: Rejecting initial value:
#  Chain 4:   Error evaluating the log probability at the initial value.
#Chain 4: Exception: normal_lpdf: Location parameter[7] is nan, but must be finite!  (in 'modelf34a78075689_fakeseed_winters' at line 100)
summary(germ.mod4)

