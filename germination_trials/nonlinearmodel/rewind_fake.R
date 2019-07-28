##Began by Dan and Lizzie late Feb 2019
###updated most recently by Dan on April 4 2019
##Purpose is to simulat germiantion data a kind to Dan's trial

rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

loadmodels <- FALSE
runpart1 <- FALSE

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

##The code below genrates fake data with 10 levels of chilling and 2 levels of forcing, with no interactions between treatments.
time <-seq(0,24,by=1) #time of each trial
chilltreat <- c(0,1,2,3,4,5,6,7,8,9,10) # level of chilling, continuous data
forcetreat<-c(0,1) #levels of forgince
sigma_y <- 0.01 ## small signma


t50.a <- 15 #intercept of t50
t50.cb <- -1 # slope of t50 with chilling
t50.fb<--.5 # slope of t50 with forcing

beta.a <- 4 #intercept of beta (shape paramenter)
beta.cb<-.5 #slope of beta chilling
beta.fb<-1 #slope of beta forcing

d.a <- 0.4 # intercept of d (maximum germination %)
d.cb<-0.03 # slope of d chilling
d.fb<-0.3 # slope of d forcing

#Interactions
d.cxf<--0.02
t50.cxf<-.5
b.cxf<--.2

repz<-seq(1,50,by=1) ## number of replicates



###make a fake data list with above paramenters
df2<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),forcetreat=numeric(),ID=numeric())  ##generate fake data

for (j in c(1:length(forcetreat))){
  y <- c()
for (i in c(1:length(chilltreat))){
  y <- c()
  for(k in c(1:length(repz))){ 
    y<-(d.fb*forcetreat[j]+d.cb*chilltreat[i]+d.cxf*(forcetreat[j]*chilltreat[i])+d.a)/(1+((time/(t50.fb*forcetreat[j]+t50.cb*chilltreat[i]+t50.cxf*(forcetreat[j]*chilltreat[i])+t50.a))^-(beta.fb*forcetreat[j]+beta.cb*chilltreat[i]+b.cxf*(forcetreat[j]*chilltreat[i])+beta.a)))
    #y<-(d.b*treat[i]+d.a)/(1+((time/(t50.b*treat[i]+t50.a))^-(beta.a)))
    dfhere2 <- data.frame(time=time, y=rtnorm(length(y),y,sigma_y,a=0,b=Inf),chilltreat=rep(chilltreat[i], length(y)),forcetreat=rep(forcetreat[j], length(y)),ID=rep(repz[k],length(y)))
    
    df2 <- rbind(df2, dfhere2) ## rbind it here for safty
  }
}
}


ploty2<-ggplot(df2,aes(time,y))+geom_point(aes(color=as.factor(chilltreat),shape=as.factor(forcetreat))) #plot fake data
ploty2+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat),linetype=as.factor(forcetreat))) # plot fake data with average lines

df2$y<-ifelse(df2$y>=1,1,df2$y) # correct any data what was more than 100% germiantion to 100%

ploty2<-ggplot(df2,aes(time,y))+geom_point(aes(color=as.factor(chilltreat),shape=as.factor(forcetreat))) #plot fake data
ploty2+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat),linetype=as.factor(forcetreat))) # plot fake data with average lines

df.adj2<-df2
df.adj2$time<-ifelse(df.adj2$time==0,0.0001,df.adj2$time) ### change time=0 to .0001 because models struggle to fit zero values


data.list2<-with(df.adj2,
                list(Y=y,
                     t=time,
                     chill=chilltreat,
                     force=forcetreat,
                     N=nrow(df.adj2) # datalist
                )
)

mod4.alt = stan('stan/fakeseedgoodchill_alt2param.stan', data = data.list2, 
                iter = 6000, warmup=5000, chain=1) # try model on one chaian
summary(mod4.alt)$summary[c("a_t50","a_d","a_beta","bf_beta","bf_t50","bf_d","bc_beta","bc_t50","bc_d","sigma"),]

mod4.alt.mega = stan('stan/fakeseedgoodchill_alt2param.stan', data = data.list2, 
                     iter = 10000, warmup=9000, chain=4) # try full model

summary(mod4.alt.mega)$summary[c("a_t50","a_d","a_beta","bf_beta","bf_t50","bf_d","bc_beta","bc_t50","bc_d","sigma"),]
launch_shinystan(mod4.alt.mega) # model fits

#interaction
modintxn = stan('stan/fakeseedgoodchill_winters.stan', data = data.list2, 
                iter = 6000, warmup=5000, chain=1)
summary(modintxn)$summary[c("a_d","bc_d","bf_d","a_beta","bc_beta","bf_beta","a_t50","bc_t50","bf_t50","inter_d","inter_beta","inter_t50","sigma"),]

Y_mean <- rstan::extract(modintxn, pars=c("Y_mean"))
Y_mean_cred <- apply(Y_mean$Y_mean, 2, quantile, c(0.05, 0.95))
Y_mean_mean <- apply(Y_mean$Y_mean, 2, mean)

Y_pred <- rstan::extract(modintxn, "Y_pred")
Y_pred_cred <- apply(Y_pred$Y_pred, 2, quantile, c(0.05, 0.95))
Y_pred_mean <- apply(Y_pred$Y_pred, 2, mean)
df.adj2$y ~ df.adj2$time
plot(c(0,20),c(0,1), xlab="time", ylab="germination", 
     ylim=c(0, 1),xlim=c(0,20), main="Non-linear Growth Curve")
lines(df.adj2$time, Y_mean_mean)
points(df.adj2$time, Y_pred_mean, pch=19, col=4)
lines(df.adj2$time, Y_mean_cred[1,], col=4)
lines(df.adj2$time, Y_mean_cred[2,], col=4)
lines(df.adj2$time, Y_pred_cred[1,], col=2)
lines(df.adj2$time, Y_pred_cred[2,], col=2)
legend(x="bottomright", bty="n", lwd=2, lty=c(NA, NA, 1, 1,1),
       legend=c("observation", "prediction", "mean prediction",
                "90% mean cred. interval", "90% pred. cred. interval"),
       col=c(1,1,1,4,2),  pch=c(1, 19, NA, NA, NA))


modintxn.mega = stan('stan/fakeseedgoodchill_winters.stan', data = data.list2, 
                iter = 12000, warmup=11000, chain=4)
summary(modintxn.mega)$summary[c("a_d","bc_d","bf_d","a_beta","bc_beta","bf_beta","a_t50","bc_t50","bf_t50","inter_d","inter_beta","inter_t50","sigma"),] #horrible Rhats



save.image("fake_germ_models") 
stop("Not an error, just everything below is older code that was used to build current edition")
#######################################Older editiions below
### generate fake data with only chilling and t50 beta
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


if(runpart1){
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
  # matches
}

###### Priors #########
hist(rnorm(1000,0,3))
hist(rbeta(1000,2,2))
###Part II chilling (0,1) alters  t50 and beta but not d#################################

data.list<-with(df.adj,
                list(Y=y,
                     t=time,
                     chill=chilltreat,
                     N=nrow(df.adj)
                )
)

mod3 = stan('stan/fakeseedgoodchill.stan', data = data.list, 
            iter = 6000, warmup=5000, chain=1)

mod3.alt = stan('stan/fakeseedgoodchill_alt.stan', data = data.list, 
                iter = 6000, warmup=5000, chain=1) # this runs, but I did get divergent transitions sometimes when increasing chain number unless I had high iter
summary(mod3.alt)$summary[c("a_t50","d","beta","b_t50","sigma"),]

mod3.alt.mega = stan('stan/fakeseedgoodchill_alt.stan', data = data.list, 
                     iter = 10000, warmup=9000, chain=4)
summary(mod3.alt.mega)$summary[c("a_t50","d","beta","b_t50","sigma"),]

mod3.sum<-summary(mod3)$summary
mod3.sum[c("a_t50","d","beta","b_t50","sigma"),] ### but very wrong parameter estimate for a_t50 and b_t50 and sigma.
launch_shinystan(mod3)


mod3.drc<-drm(y~time,chilltreat,fct=LL.3(),data=df,type="continuous") ### check it against drc
summary(mod3.drc) ### this work

##########
