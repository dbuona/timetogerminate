###a script to amke my data survival formated
###need to add a way to elminate seeds that rotted.
rm(list=ls()) 
options(stringsAsFactors = FALSE)


library(ggplot2)
library(dplyr)
library("arm")
library(brms)
library(ggstance)
setwd("~/Documents/git/timetogerminate/germination_trials/")
load("bernoullis")

d<-read.csv("survival_analysis/surival_dat_nointerval.csv")
d<- filter(d,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea","Carex grayi"))

d$chillweeks<-d$chill_time/7 # make chilling weeks instead of days

d$force<-NA # make forcing numeric
d<- within(d, force[INC=="L"]<-0)
d<- within(d, force[INC=="H"]<-1)

mod1<-brm(germinated~chillweeks*force+(chillweeks*force|Taxa),data=d, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)

beta<-as.data.frame(coef(mod1))
hi<-1
lo<-0
###for [pollination syndrome]
beta<-as.data.frame(coef(mod1))

new.data<-data.frame(Taxa=rep(unique(d$Taxa),each=20),chillweeks=rep(unique(d$chillweeks),20),force=rep(c(0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1),20))
goot<-fitted(mod1,newdata = new.data,summary = TRUE)
prediction<-cbind(new.data,goot)                     
ggplot(prediction,aes(chillweeks,Estimate))+geom_point()+stat_smooth(method=loess,aes(color=Taxa,fill=Taxa))+facet_wrap(~as.factor(force))+ggthemes::theme_base(base_size = 11)

unique(d$chillweeks)


d.rest<-filter(d,chillweeks %in% c(0,6,13))
modrest<-brm(germinated~chillweeks*force+(chillweeks*force|Taxa),data=d.rest, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)

betas2<-as.data.frame(coef(modrest))
coef(mod1)
save.image("bernoullis")

beta$Taxa.Estimate.chillweeks
betas2$Taxa.Estimate.chillweeks

new.data<-data.frame(Taxa=rep(unique(d$Taxa),5),chillweeks=rep(c(0,2,6,9,13),each=10), force=rep(0,50))
threetreat<-fitted(modrest,newdata = new.data,summary = TRUE)
threetreat<-cbind(new.data,threetreat)
fulltreat<-fitted(mod1,newdata = new.data,summary = TRUE)
fulltreat<-cbind(new.data,fulltreat)

threetreat$design<-"classify"
fulltreat$design<-"regression"

comps<-rbind(fulltreat,threetreat)

ggplot(comps,aes(chillweeks,Estimate))+geom_point(aes(color=design),size=2)+geom_errorbar(aes(ymin=Q2.5,ymax=Q97.5,color=design),width=0)+facet_wrap(~Taxa)+scale_y_continuous(breaks=(c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1)))
?geom_errorbar()
