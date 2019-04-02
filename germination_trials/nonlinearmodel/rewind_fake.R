##Began by Dan and Lizzie late Feb 2019
###updated most recently by Dan on March 1 2019
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
sigma_y <- 0.01 ## small signmal

t50.a<-15 #intercept of t50
t50.b<--5 # slope of t50 with chilling

beta.a<-6 #intercept of beta (shape paramenter)
beta.b<-2 # slope of beta with chilling

d.a<-5 # intercept of d (maximum germination)
d.b<-10 #slope of d with chilling

repz<-seq(1,50,by=1) ## number of replicates

###for starters assme a 100% germination scenario, no (D parameter)
df<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),ID=numeric())  ##generate fake data

for (i in c(1:length(treat))){
  y <- c()
    for(k in c(1:length(repz))){ 
      #y<-(d.b*treat[i]+d.a)/(1+((time/(t50.b*treat[i]+t50.a))^-(beta.b*treat[i]+beta.a)))
      y<-1/(1+((time/(t50.b*treat[i]+t50.a))^-(beta.b*treat[i]+beta.a)))
      dfhere <- data.frame(time=time, y=rtnorm(length(y),y,sigma_y,a=0,b=1),chilltreat=rep(treat[i], length(y)),ID=rep(repz[k],length(y))) ## make a data frame for each level, this over rights so
      
      df <- rbind(df, dfhere) ## rbind it here for safty
    }
  }

df$uniqueID<-paste(df$chilltreat,df$ID,sep="-")
ploty<-ggplot(df,aes(time,y))+geom_point(aes(color=as.factor(chilltreat)))
ploty
ploty+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat))) ### plot point with mean lines these datas look good


### PART 1: NO TREATMENTS ###########
notreat<-filter(df,chilltreat==1) ## do this for models with out any predictors)
ggplot(notreat,aes(time,y))+geom_point()
data.list.notreat<-with(notreat,
                list(Y=y,
                     t=time,
                     N=nrow(notreat)
                     ))
                
###below are the same model coded in different ways. mod 1 ia good and returns proper parameners, mod 2 is baaad. 
mod1= stan('stan/fakeseedmodel.stan', data = data.list.notreat,
                      iter = 3000, warmup=2000) ###good mdel
mod1.sum<-summary(mod1)$summary

mod1.sum[c("beta","t50","d","sigma"),] 

#mod2= stan('stan/altfakeseed.model.stan', data = data.list.notreat,
 #          iter = 3000, warmup=2000) ##1930 divergent transitions

##compare with drc

mod1.drc<-drm(y~time,fct=LL.3(),data=notreat,type="continuous")
summary(mod1.drc)
#
######################M#########
###Part II chilling (0,1) alters  t50 and beta but not d#################################
df.adj<-df
df.adj$time<-ifelse(df.adj$time==0,0.0001,df.adj$time) ### I think the log logistic distrubtion can't handle when time=0 because it creata log of (0) situation
                                                        ###this alwows germ mo    
data.list<-with(df.adj,
                list(Y=y,
                     t=time,
                     chill=chilltreat,
                     N=nrow(df.adj)
                )
)

###below are the same model coded in different ways. 
mod3 = stan('stan/fakeseedgoodchill.stan', data = data.list,                        ##4/1/19 
                      iter = 3000, warmup=2000, control = list(max_treedepth = 15)) #16 divergent transitions after warmup
                                    #3176  transitions after warmup that exceeded the maximum treedepth
                                    # bad Rhats, but I broke it again worse trying to imrpove it. Doesn't run in 24 hours

mod3.sum<-summary(mod3)$summary
mod3.sum[c("a_beta","a_t50","d","b_beta","b_t50","sigma"),] 

mod3.drc<-drm(y~time,chilltreat,fct=LL.3(),data=df,type="continuous")
summary(mod3.drc)

mod4 = stan('stan/fakeseed_chillonly.stan', data = data.list,
                                iter = 3000, warmup=2000,control = list(max_treedepth = 15)) ###  250 divergent transitions
                                                                                              ##better rhats
mod4.sum<-summary(mod4)$summary
mod4.sum[c("a_beta","a_t50","b_chill_beta","b_chill_t50","sigma"),]                                                          #
launch_shinystan(mod4)

# emw: this model returns 3864 div trans for me and the rhat values are awful -- it is not converging at all --- are you sure it is running for you?
#DB agreeed! 

bin.sum<-summary(germ.mod.chill)$summary
bin.sum[c("a_beta","a_t50","a_d","b_chill_beta","b_chill_t50","b_chill_d","sigma"),] 







stop("not an error") # emw ... assuming I don't run anything after this? Many errors below. 


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
                         iter = 6000, warmup=4000) 

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
