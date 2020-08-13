rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")

library(rstan)
library(tidyr)
library(dplyr)
library(brms)
library(ggplot2)
realdat<-read.csv("input/daily_dat_nointerval.csv")

##clean data
realdat$germ_perc<-NA
realdat$germ_perc<-realdat$germ_num/realdat$tot_seed
realdat$germ_perc<-ifelse(realdat$germ_perc>1,1,realdat$germ_perc)

##make chilling numeric
realdat$chill_time<-NA
realdat<- within(realdat, chill_time[COLD=="0" ]<-0)
realdat<- within(realdat, chill_time[COLD=="A" ]<-14)
realdat<- within(realdat, chill_time[COLD=="B" ]<-28)
realdat<- within(realdat, chill_time[COLD=="C" ]<-35)
realdat<- within(realdat, chill_time[COLD=="D" ]<-42)
realdat<- within(realdat, chill_time[COLD=="E" ]<-49)
realdat<- within(realdat, chill_time[COLD=="f" ]<-56)
realdat<- within(realdat, chill_time[COLD=="G" ]<-63)
realdat<- within(realdat, chill_time[COLD=="H" ]<-77)
realdat<- within(realdat, chill_time[COLD=="i" ]<-91)
realdat$chillweeks<-realdat$chill_time/7 # make chilling weeks instead of days

realdat$force<-NA # make forcing numeric
realdat<- within(realdat, force[INC=="L"]<-0)
realdat<- within(realdat, force[INC=="H"]<-5)
realdat<- filter(realdat,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea","Carex grayi","Silene stellata"))
realdat$DAY<-ifelse(realdat$DAY==0,0.0001,realdat$DAY)
###select plates where there is at least 25% final germination
fgp.dat<-filter(realdat,DAY==25)
fgp.dat<-filter(fgp.dat,germ_perc>=.5)

plates<-unique(fgp.dat$plate_num)
df<-data.frame(plate_num=numeric(),Taxa=character(),force=numeric(),chillweeks=numeric(),T25=numeric())

for (p in seq_along(plates)){
  dataoneplate <- subset(realdat, plate_num==plates[p])
  t50 <- sum(dataoneplate$DAY*dataoneplate$germ.daily)/20
  dfhere<-data.frame(platenum=plates[p],Taxa=dataoneplate$Taxa,force=dataoneplate$force, chillweeks=dataoneplate$chillweeks,T50=t50)
  df <- rbind(df, dfhere)
}


df<-df %>% distinct()
sumt50<- df %>% group_by(Taxa,force,chillweeks) %>% summarize(meant50=mean(T50),sdt50=sd(T50)) 

maxminmeans<- sumt50 %>% group_by(Taxa) %>% summarise(mint50=min(meant50),maxt50=max(meant50))
class(df$chillweeks)

ggplot(df,aes(chillweeks,T50))+geom_point(aes(color=Taxa,shape=as.factor(force)))+
  geom_smooth(method=lm,aes(color=Taxa,fill=Taxa),alpha=0.2)+facet_wrap(~as.factor(force))+ggthemes::theme_base()

ggplot(df,aes(chillweeks,T50))+stat_summary(aes(color=Taxa,shape=as.factor(force)))+
  facet_wrap(~as.factor(force))+ggthemes::theme_base()
##center predictors



                     
mod.T50<-brm(T50~chillweeks*force+(1+chillweeks*force|Taxa),data=df)## maybe I dont want to pool because some dont germinate at low chilling but others do
coef(mod.T50)
new.data<-data.frame(Taxa=rep(unique(df$Taxa),each=6),
                     force=rep(c(0,0,0,5,5,5),9),
                     chillweeks=rep(c(4,7,13),9))

prediction<-predict(mod.T50,newdata=new.data,probs = c(.25,.75))
predy<-cbind(new.data,prediction)
predy<-filter(predy,Taxa!="Anemone virginana")
ggplot(predy,aes(as.factor(chillweeks),Estimate))+geom_point(aes(color=Taxa),position=position_dodge(width = .2,))+facet_wrap(~as.factor(force))+geom_errorbar(aes(ymin=Q25,ymax=Q75,color=Taxa,width=0),position=position_dodge(width = .2))+ggthemes::theme_base()


forest<-filter(predy, Taxa %in% c("Cryptotaenia canadensis","Eurbia diviricata","Hesperis matronalis","Polygonum virginiatum"))
field<-filter(predy, !Taxa %in% c("Cryptotaenia canadensis","Eurbia diviricata","Hesperis matronalis","Polygonum virginiatum"))

ggplot(forest,aes(as.factor(chillweeks),Estimate))+geom_point(aes(color=Taxa),position=position_dodge(width = .2,))+facet_wrap(~as.factor(force))+geom_errorbar(aes(ymin=Q25,ymax=Q75,color=Taxa,width=0),position=position_dodge(width = .2))+ggthemes::theme_base()

ggplot(field,aes(as.factor(chillweeks),Estimate))+geom_point(aes(color=Taxa),position=position_dodge(width = .2,))+facet_wrap(~as.factor(force))+geom_errorbar(aes(ymin=Q25,ymax=Q75,color=Taxa,width=0),position=position_dodge(width = .2))+ggthemes::theme_base()

ggplot()+geom_point(data=sum.crpy,aes(chillweeks,meant50,shape=as.factor(force)),color="blue")+
  geom_point(data=cryp.pred,aes(chillweeks,Estimate,shape=as.factor(force)))
