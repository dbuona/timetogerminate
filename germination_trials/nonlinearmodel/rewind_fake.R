##Began by Dan and Lizzie late Feb 2019
###updated most recently by Dan on March 1 2019
##Purpose is to simulat germiantion data a kind to Dan's trial

##next step, a thoughtful way to add error between replicates
#1. I could add error on each parementer, I think using rnorm
# or just add error to the final why value
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

setwd("~/Documents/git/timetogerminate/germination_trials")

library(rstan)
library(tidyr)
library(drc)
library(dplyr)
library(shinystan)
time<-seq(0,24,by=3) #time of each trial
treat<-c(0,14,28,35,42,49,56,63,77,91) # level of chilling, continuous data
force<-c(0,1)# 2 levels of forching, low/high

##parameters
t50.a<-20 #intercept of t50
t50.b<--0.09 # slope of t50 with chilling
t50.f<--1 #slope of t50 with forcing


beta.a<--4 #intercept of beta (shape paramenter)
beta.b<--0.05 # slope of beta with chilling
beta.f<--1.5 # slope of beta with forcing

d.a<-0 # intercept of d (maximum germination)
d.b<-0.19 #slope of d with chilling
d.f<-2 # slope of d with forcing

repz<-seq(1,10,by=1) ## number of replicates


# Trying running just this:  ##DB, rthis doesn't run, but I think it was just trouble shooting
#l <- 1
#  y<-(d.b*treat[i]+d.a)/(1+((time/(t50.a)^beta.b)))
#  dfhere <- data.frame(time=time, y=y,treat=rep(treat[i], length(y)))
#  df <- rbind(df, dfhere)
# This tells you that in each loop of l you get ALL the time runs ... so this means the way your time loop works it does the above NINE TIMES
# 9 * 9 * 2 treats = 162
# So, now we know the time loop is the problem.... and actually it make me think that you just don't need it ...


##This is the base loop 
df<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),forcetreat=numeric(),ID=numeric())  ##3why is this breaking?

  for (i in c(1:length(treat))){
   y <- c()
  for (l in c(1:length(force))){
    y<-c()
    for(k in c(1:length(repz))){ 
      d.error<-rnorm(1,0,0.5)
      t50.error<-rnorm(1,0,0.5)
      beta.error<-rnorm(1,0,0.5)
      y<-(d.b*treat[i]+d.f*force[l]+d.a+d.error)/(1+((time/(t50.b*treat[i]+t50.f*force[l]+t50.a+t50.error))^(beta.b*treat[i]+beta.f*force[l]+beta.a+beta.error)))
  dfhere <- data.frame(time=time, y=y,chilltreat=rep(treat[i], length(y)),forcetreat=rep(force[l],length(y)),ID=rep(repz[k],length(y))) ## make a data frame for each level, this over rights so
  df <- rbind(df, dfhere) ## rbind it here for safty
}
}
}

##df$y<-round(df$y) #maybe I don't want this if I am trying to get the model to return the error
 ##plot the data
ploty<-ggplot(df,aes(time,y))+geom_point(aes(color=as.factor(chilltreat),shape=as.factor(forcetreat)))
ploty
ploty+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat),linetype=as.factor(forcetreat))) ### plot point with mean lines

chillonly<-filter(df,forcetreat==0)### try to model just chilling for now.

###try it in a stan model
data.list.chill <- with(chillonly, 
                   list(Y=y, 
                        t = time,
                        chill=chilltreat,
                        N = nrow(chillonly)
                   )
)


germ.mod.chillonly = stan('stan/fakeseed_wchill_emw.stan', data = data.list.chill,
                 iter = 3000, warmup=2000) 
summary(germ.mod.chillonly) # 700+ divergent transitions. 
launch_shinystan(germ.mod.chillonly)
