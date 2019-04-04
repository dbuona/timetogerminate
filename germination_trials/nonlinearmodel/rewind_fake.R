##Began by Dan and Lizzie late Feb 2019
###updated most recently by Dan on April 4 2019
##Purpose is to simulat germiantion data a kind to Dan's trial

rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

loadmodels <- FALSE

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")

if(loadmodels){
load("fake_germ_models") # emw: What is this, do I need it? DB...no its only so I can load models but none currently work.
}

library(rstan)
library(tidyr)
library(drc)
library(dplyr)
library(shinystan)
library(extraDistr)

######Chilling only
time<-seq(0,24,by=1) #time of each trial
treat<-c(0,1) # level of chilling, continuous data
sigma_y <- 0.01 ## small signma

#Start with data where only t50 changes
t50.a<-15 #intercept of t50
t50.b<--5 # slope of t50 with chilling

beta.a<-6 #intercept of beta (shape paramenter)
d.a<-0.9 # intercept of d (maximum germination %)

repz<-seq(1,50,by=1) ## number of replicates

### generate fake data
df<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),ID=numeric())  ##generate fake data

for (i in c(1:length(treat))){
  y <- c()
    for(k in c(1:length(repz))){ 
      #y<-(d.b*treat[i]+d.a)/(1+((time/(t50.b*treat[i]+t50.a))^-(beta.b*treat[i]+beta.a)))
      y<-d.a/(1+((time/(t50.b*treat[i]+t50.a))^-(beta.a)))
      dfhere <- data.frame(time=time, y=rtnorm(length(y),y,sigma_y,a=0,b=Inf),chilltreat=rep(treat[i], length(y)),ID=rep(repz[k],length(y)))
      
      df <- rbind(df, dfhere) ## rbind it here for safty
    }
  }

###plot the trials
ploty<-ggplot(df,aes(time,y))+geom_point(aes(color=as.factor(chilltreat)))
ploty
ploty+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat))) ### plot with mean lines

df.adj<-df
df.adj$time<-ifelse(df.adj$time==0,0.0001,df.adj$time) ###the log logistic can't handle time values of zero. I am not sure if this is the right way to handle this



### PART 1: NO TREATMENTS ###########
notreat<-filter(df.adj,chilltreat==1) ## do this for models with out any predictors
ggplot(notreat,aes(time,y))+geom_point() #plot single treatment

data.list.notreat<-with(notreat, # prepare for stan
                list(Y=y,
                     t=time,
                     N=nrow(notreat)
                     ))
                
mod1= stan('stan/fakeseedmodel.stan', data = data.list.notreat,
                      iter = 3000, warmup=2000) #### stan model

mod1.sum<-summary(mod1)$summary
mod1.sum[c("beta","t50","d","sigma"),] ### these parameters are correct

##compare with drc
mod1.drc<-drm(y~time,fct=LL.3(),data=notreat,type="continuous")
summary(mod1.drc)
# mattches
######################M#########

###Part II chilling (0,1) alters  t50 and beta but not d#################################
                                                          
data.list<-with(df.adj,
                list(Y=y,
                     t=time,
                     chill=chilltreat,
                     N=nrow(df.adj)
                )
)

mod3 = stan('stan/fakeseedgoodchill.stan', data = data.list,  ###13 divergent transitions
                      iter = 3000, warmup=2200) 
mod3.sum<-summary(mod3)$summary
mod3.sum[c("a_t50","d","beta","b_t50","sigma"),] ###but very wrong parameter estimate for a_t50 and b_t50 and sigma.

mod3.drc<-drm(y~time,chilltreat,fct=LL.3(),data=df,type="continuous") ### check it against drc
summary(mod3.drc) ### this work

##########


