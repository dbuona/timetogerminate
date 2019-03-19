##Began by Dan and Lizzie late Feb 2019
###updated most recently by Dan on March 1 2019
##Purpose is to simulat germiantion data a kind to Dan's trial

##next step, a thoughtful way to add error between replicates
#1. I could add error on each parementer, I think using rnorm
# or just add error to the final why value
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")

library(rstan)
library(tidyr)
library(drc)
library(dplyr)
library(shinystan)
library(extraDistr)

######Chilling only
time<-seq(0,24,by=3) #time of each trial
treat<-c(0,1) # level of chilling, continuous data
sigma_y <- 0.1

t50.a<-20 #intercept of t50
t50.b<--5 # slope of t50 with chilling

beta.a<--5 #intercept of beta (shape paramenter)
beta.b<--2 # slope of beta with chilling

d.a<-5 # intercept of d (maximum germination)
d.b<-10 #slope of d with chilling

repz<-seq(1,50,by=1) ## number of replicates




####Force and chill
time<-seq(0,24,by=3) #time of each trial
treat<-c(0,1) # level of chilling, continuous data
force<-c(0,1)# 2 levels of forching, low/high

##parameters
t50.a<-20 #intercept of t50
t50.b<--5 # slope of t50 with chilling
t50.f<--1 #slope of t50 with forcing
sigma_y <- 0.1

beta.a<--5 #intercept of beta (shape paramenter)
beta.b<--2 # slope of beta with chilling
beta.f<--1.5 # slope of beta with forcing


d.a<-5 # intercept of d (maximum germination)
d.b<-10 #slope of d with chilling


repz<-seq(1,50,by=1) ## number of replicates

df<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),ID=numeric())  ##3why is this breaking?

for (i in c(1:length(treat))){
  y <- c()
    for(k in c(1:length(repz))){ 
      y<-(d.b*treat[i]+d.a)/(1+((time/(t50.b*treat[i]+t50.a))^(beta.b*treat[i]+beta.a)))
      dfhere <- data.frame(time=time, y=rtnorm(length(y),y,sigma_y,a=0,b=Inf),chilltreat=rep(treat[i], length(y)),ID=rep(repz[k],length(y))) ## make a data frame for each level, this over rights so
      
      df <- rbind(df, dfhere) ## rbind it here for safty
    }
  }

ploty<-ggplot(df,aes(time,y))+geom_point(aes(color=as.factor(chilltreat)))
ploty
ploty+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat))) ### plot point with mean lines

data.list<-with(df,
                list(Y=y,
                     t=time,
                     chill=chilltreat,
                     N=nrow(df)
                )
)

germ.mod.chill = stan('stan/fakeseed_chillonly.stan', data = data.list,
                                iter = 3000, warmup=2000) 

bin.sum<-summary(germ.mod.chill)$summary
bin.sum[c("a_beta","a_t50","a_d","b_chill_beta","b_chill_t50","b_chill_d","sigma"),] ###returns the right paraments but seems to stuggle


stop("not an error")


##This is the base loop 
df<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),forcetreat=numeric(),ID=numeric())  ##3why is this breaking?

  for (i in c(1:length(treat))){
   y <- c()
  for (l in c(1:length(force))){
    y<-c()
    for(k in c(1:length(repz))){ 
      y<-(d.b*treat[i]+d.f*force[l]+d.a)/(1+((time/(t50.b*treat[i]+t50.f*force[l]+t50.a))^(beta.b*treat[i]+beta.f*force[l]+beta.a)))
  dfhere <- data.frame(time=time, y=rtnorm(length(y),y,sigma_y,a=0,b=Inf),chilltreat=rep(treat[i], length(y)),forcetreat=rep(force[l],length(y)),ID=rep(repz[k],length(y))) ## make a data frame for each level, this over rights so

   df <- rbind(df, dfhere) ## rbind it here for safty
}
}
}

##df$y<-round(df$y) #maybe I don't want this if I am trying to get the model to return the error
 ##plot the data
ploty<-ggplot(df,aes(time,y))+geom_point(aes(color=as.factor(chilltreat),shape=as.factor(forcetreat)))
ploty
ploty+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat),linetype=as.factor(forcetreat))) ### plot point with mean lines




##model is good when when chill is 0 and focing is 0,1

data.list<-with(df,
                list(Y=y,
                     t=time,
                     warm=forcetreat,
                     chill=chilltreat,
                     N=nrow(df)
                     )
                )

germ.mod.warmchill.noint = stan('stan/fakeseed_forcechill_noint.stan', data = data.list,
                         iter = 2000, warmup=1500) 

bin.sum<-summary(germ.mod.warmchill.noint)$summary
bin.sum[c("a_beta","a_t50","a_d","b_warm_beta","b_warm_t50","b_warm_d","b_chill_beta","b_chill_t50","b_chill_d","sigma"),] ###whoohoo it works! but does it work with logitsitc?

## this model(treatments binary) returns the proper parementser when Y~ normal(yhat,sigma), but the model with chilling as continuous won't run well, rejects initial values and divergent transitions
## binary is still okay, but parameters drift off when Y~ lognormal(yhat,sigma)





### now try it with continuous values for chillin
treat2<-c(0,1,2,4,5,6,7,8,9,10) # level of chilling, continuous data

t50.b2<--1
beta.b2<--0.75
d.b2<-1.2

df2<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),forcetreat=numeric(),ID=numeric())  ##3why is this breaking?

for (i in c(1:length(treat2))){
  y <- c()
  for (l in c(1:length(force))){
    y<-c()
    for(k in c(1:length(repz))){ 
      y2<-(d.b2*treat2[i]+d.f*force[l]+d.a)/(1+((time/(t50.b2*treat2[i]+t50.f*force[l]+t50.a))^(beta.b2*treat2[i]+beta.f*force[l]+beta.a)))
      dfhere2 <- data.frame(time=time, y=rtnorm(length(y2),y2,sigma_y,a=0,b=Inf),chilltreat=rep(treat2[i], length(y2)),forcetreat=rep(force[l],length(y2)),ID=rep(repz[k],length(y2))) ## make a data frame for each level, this over rights so
      
      df2 <- rbind(df2, dfhere2) ## rbind it here for safty
    }
  }
}

ploty2<-ggplot(df2,aes(time,y))+geom_point(aes(color=as.factor(chilltreat),shape=as.factor(forcetreat)))
ploty2
ploty2+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat),linetype=as.factor(forcetreat))) ### plot point with mean lines

chillonly<-filter(df2,forcetreat==1)

data.list.cont<-with(chillonly,
                list(Y=y,
                     t=time,
                     chill=chilltreat,
                     N=nrow(chillonly)
                )
)
germ.mod.warmchill.noint.cont = stan('stan/fakeseed_chillonly.stan', data = data.list.cont,
                                iter = 2000, warmup=1200) 
cont.sum<-summary(germ.mod.warmchill.noint.cont)$summary #585 divergent transitions :(

cont.sum[c("a_beta","a_t50","a_d","b_warm_beta","b_warm_t50","b_warm_d","b_chill_beta","b_chill_t50","b_chill_d","sigma"),]
save.image("fake_germ_models")
#this happens 
#Chain 4: Rejecting initial value:
#  Chain 4:   Error evaluating the log probability at the initial value.
#Chain 4: Exception: normal_lpdf: Location parameter[452] is nan, but must be finite!  (in 'modelfd5bc445cb0_fakeseed_forcechill_noint' at line 56)
##forum say; (c) that's a sign that your constraints don't match your model.  if a value of the parameters satisfies the declared constraints, it should have positive density (finite log density). 
