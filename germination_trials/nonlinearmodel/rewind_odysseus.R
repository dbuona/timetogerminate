rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(rstan)
#library(tidyr)
#library(drc)
#library(dplyr)
#library(shinystan)
#library(extraDistr)
### 9/12/19 new plan: for real data. run a chilling model on subset of data at each forcing level, so we just need to simulate data for chilling
time <-seq(0,24,by=1) #time of each trial
chilltreat <- c(0,1,2,3,4,5,6,7,8,9,10)
sigma_y <- 0.01 ## small signma




t50.a <- 20 #intercept of t50
t50.cb <- -1.5 # slope of t50 with chilling

beta.a <- 4 #intercept of beta (shape paramenter)
beta.cb<-.5 #slope of beta chilling

d.a <- 0.2 # intercept of d (maximum germination %)
d.cb<-0.07 # slope of d chilling

repz<-seq(1,50,by=1)

df2<-data.frame(time=numeric(), y=numeric(),chilltreat=numeric(),ID=numeric())  ##generate fake data


for (i in c(1:length(chilltreat))){
  y <- c()
  for(k in c(1:length(repz))){ 
    y<-(d.cb*chilltreat[i]+d.a)/(1+((time/(t50.cb*chilltreat[i]+t50.a))^-(beta.cb*chilltreat[i]+beta.a)))
    #y<-(d.b*treat[i]+d.a)/(1+((time/(t50.b*treat[i]+t50.a))^-(beta.a)))
    dfhere2 <- data.frame(time=time, y=rtnorm(length(y),y,sigma_y,a=0,b=Inf),chilltreat=rep(chilltreat[i], length(y)),ID=rep(repz[k],length(y)))
    
    df2 <- rbind(df2, dfhere2) ## rbind it here for safty
  }
}


#ploty2<-ggplot(df2,aes(time,y))+geom_point(aes(color=as.factor(chilltreat))) #plot fake data
#ploty2+geom_line(stat = "summary", fun.y = mean, aes(color=as.factor(chilltreat))) # plot fake data with average lines


df.adj2<-df2
df.adj2$time<-ifelse(df.adj2$time==0,0.0001,df.adj2$time) ### change time=0 to .0001 because models struggle to fit zero values


data.list2<-with(df.adj2,
                 list(Y=y,
                      t=time,
                      chill=chilltreat,
                      N=nrow(df.adj2) # datalist
                 ))



modchill.mega = stan('stan/fakeseedgoodchill_alt.stan', data = data.list2, 
                     iter = 10000, warmup=9000, chain=4)

save(modchill.mega, file="/n/home04/dbuonaiuto/odysseus/modchill.mega.rds")