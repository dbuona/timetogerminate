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
library(bayesplot)
library(lme4)
library(tibble)
library(segmented)
library(RColorBrewer)
library(brms)
library(grid)
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

realdat$DAY<-ifelse(realdat$DAY==0,0.0001,realdat$DAY) #elimiate 0 values for log logistic dist

realdatshorty<- filter(realdat,Taxa %in% c("Cryptotaenia canadensis","Polygonum virginiatum","Eurbia diviricata")) #These 3 sp show a strong chill responmse

fgp.dat<-filter(realdatshorty,DAY==25) ### this makes a dataset of only final germination percentate

#full.fgp.mod<-brms::brm(germ_perc~force*chillweeks+(force*chillweeks|Taxa),data=fgp.dat,iter=4000,warmup=3000)


ggplot(fgp.dat,aes(chillweeks,germ_perc))+
  geom_smooth(method="lm",level=0.9,aes(color=as.factor(force),fill = as.factor(force)))+
  facet_wrap(~Taxa)+scale_fill_manual(values=c("royalblue","firebrick1"))+
  scale_color_manual(values=c("royalblue","firebrick1"))+geom_hline(aes(yintercept=0),color="black")+geom_hline(aes(yintercept=1),color="black")+
  ylim(-.3,1.3)+theme_linedraw()+geom_point(aes(color=as.factor(force)),size=0.5)+geom_vline(aes(xintercept=8),color="gray",size=2)


fgp.dat.cool<-filter(fgp.dat,force==0)
ggplot(fgp.dat.cool,aes(chillweeks,germ_perc))+stat_summary(aes(color=Taxa))+geom_vline(aes(xintercept=7))

fgp.dat.cool.before<-filter(fgp.dat.cool,chillweeks<=7)
fgp.dat.cool.after<-filter(fgp.dat.cool,chillweeks>=7)

fgp.dat.cool$thresh<-ifelse(fgp.dat.cool$chillweeks>7.9,"after","before")
ggplot(data=fgp.dat.cool,aes(x=chillweeks,y=germ_perc))+
  geom_smooth(method="lm",aes(color=Taxa,linetype=thresh,fill=Taxa),alpha=0.1,size=.3)+scale_linetype_manual(values=c("solid","solid"))+
  stat_summary(aes(color=Taxa))+geom_vline(aes(xintercept=7))

a<-ggplot()+
  geom_smooth(data=fgp.dat.cool,method="lm",aes(x=chillweeks,y=germ_perc,color=Taxa,fill=Taxa),alpha=0.1,size=.3)+geom_hline(yintercept=1)+
  geom_hline(yintercept=0)+stat_summary(data=fgp.dat.cool,aes(x=chillweeks,y=germ_perc,color=Taxa))+geom_vline(aes(xintercept=7),linetype="dotted")+
  ylim(0,1.2)
  
b<-ggplot()+
  geom_smooth(data=fgp.dat.cool.before,method="lm",aes(x=chillweeks,y=germ_perc,color=Taxa,fill=Taxa),alpha=0.1,size=.3)+
  geom_smooth(data=fgp.dat.cool.after,method="lm",aes(x=chillweeks,y=germ_perc,color=Taxa,fill=Taxa),alpha=0.1,size=.3)+
  stat_summary(data=fgp.dat.cool,aes(x=chillweeks,y=germ_perc,color=Taxa))+geom_vline(aes(xintercept=7),linetype="dotted")+geom_hline(yintercept=1)+
  geom_hline(yintercept=0)+ylim(0,1.2)

jpeg("figures/germ_perc_plot.jpeg",width = 7,height = 5,units = "in",res=150)
ggpubr::ggarrange(a,b,common.legend = TRUE)
dev.off()
fgp.dat.crypto<-filter(fgp.dat.cool,Taxa=="Cryptotaenia canadensis")
fgp.dat.poly<-filter(fgp.dat.cool,Taxa=="Polygonum virginiatum")
fgp.dat.eury<-filter(fgp.dat.cool,Taxa=="Eurbia diviricata")
cryptoafter<-filter(fgp.dat.crypto,chillweeks>=7)
cryptobefore<-filter(fgp.dat.crypto,chillweeks<7)

polyafter<-filter(fgp.dat.poly,chillweeks>=7)
polybefore<-filter(fgp.dat.poly,chillweeks<7)
euryafter<-filter(fgp.dat.eury,chillweeks>=7)
eurybefore<-filter(fgp.dat.eury,chillweeks<7)


summary(lm(germ_perc~chillweeks,data=fgp.dat.crypto))
summary(lm(germ_perc~chillweeks,data=cryptobefore))
summary(lm(germ_perc~chillweeks,data=cryptoafter))

summary(lm(germ_perc~chillweeks,data=fgp.dat.poly))
summary(lm(germ_perc~chillweeks,data=polybefore))
summary(lm(germ_perc~chillweeks,data=polyafter))

summary(lm(germ_perc~chillweeks,data=fgp.dat.eury))
summary(lm(germ_perc~chillweeks,data=eurybefore))
summary(lm(germ_perc~chillweeks,data=euryafter))



crypto.full.after<-filter(realdat,Taxa=="Cryptotaenia canadensis")
crypto.full.after<-filter(crypto.full.after,chillweeks>=7)
poly.full.after<-filter(realdat,Taxa=="Polygonum virginiatum")
poly.full.after<-filter(poly.full.after,chillweeks>=7)

shorty.after<-filter(realdatshorty,force==0)
shorty.after<-filter(shorty.after,chillweeks>=7)
plates<-unique(shorty.after$plate_num)
df<-data.frame(plate_num=numeric(),chillweeks=numeric(),Taxa=character(), MGT=numeric())

for (p in seq_along(plates)){
dataoneplate <- subset(shorty.after, plate_num==plates[p])
MGT<-sum(dataoneplate$DAY*dataoneplate$germ.daily)/20
dfhere <- data.frame(plate_num=plates[p],chillweeks=dataoneplate$chillweeks, Taxa=dataoneplate$Taxa,MGT=MGT)
   df <- rbind(df, dfhere) ## rbind it here for safty
}


plates2<-unique(poly.full.after$plate_num)
df2<-data.frame(plate_num=numeric(),chillweeks=numeric(),force=numeric(), MGT=numeric())

for (p in seq_along(plates2)){
  dataoneplate <- subset(poly.full.after, plate_num==plates2[p])
  MGT[p]<-sum(dataoneplate$DAY*dataoneplate$germ.daily)/20
  dfhere2 <- data.frame(plate_num=plates2[p],chillweeks=dataoneplate$chillweeks, force=dataoneplate$force,MGT=MGT[p])
  df2 <- rbind(df2, dfhere2) ## rbind it here for safty
}


df<-df %>% distinct()
df2<-df2 %>% distinct()

df$adjchill<-df$chillweeks-7
df2$adjchill<-df2$chillweeks-7

d<-ggplot()+ geom_smooth(data=fgp.dat.cool.after,method="lm",aes(x=chillweeks,y=germ_perc,color=Taxa,fill=Taxa),alpha=0.1,size=.3)+
  stat_summary(data=fgp.dat.cool.after,aes(x=chillweeks,y=germ_perc,color=Taxa),size=.1)+geom_vline(aes(xintercept=7),linetype="dotted")+geom_hline(yintercept=1)+
  geom_hline(yintercept=0)+ylim(0,1.1)+theme(legend.position="none")
c<-ggplot()+geom_smooth(data=df,method="lm",aes(x=chillweeks,y=MGT,color=Taxa,fill=Taxa),alpha=0.1,size=.7)+geom_point(data=df,aes(x=chillweeks,y=MGT,color=Taxa))

vp <- viewport(width = 0.3, height = 0.4, x = 0.35, y = .99,just=c("left","top"))
jpeg("figures/MGT_plot.jpeg",width = 7,height = 5,units = "in",res=150)
c
print(d, vp = vp )
dev.off()
summary(lm(MGT~adjchill*force,data=df))
summary(lm(MGT~adjchill*force,data=df2))
###Identify breakpoints in linear relationship for forest species
unique(fgp.dat$Taxa)
Cc<-filter(fgp.dat,Taxa=="Cryptotaenia canadensis")
Pv<-filter(fgp.dat,Taxa=="Polygonum virginiatum")
Ed<-filter(fgp.dat,Taxa=="Eurbia diviricata")
Hm<-filter(fgp.dat,Taxa=="Hesperis matronalis")
Ss<-filter(fgp.dat,Taxa=="Silene stellata")
Av<-filter(fgp.dat,Taxa=="Anemone virginana")
Td<-filter(fgp.dat,Taxa=="Thalictrum dioicum")

CC.mod.lm<-lm(germ_perc~chillweeks*force,Cc)
Pv.mod.lm<-lm(germ_perc~chillweeks*force,Pv)
Ed.mod.lm<-lm(germ_perc~chillweeks*force,Ed)
Hm.mod.lm<-lm(germ_perc~chillweeks*force,Hm)
Ss.mod.lm<-lm(germ_perc~chillweeks*force,Ss)
Av.mod.lm<-lm(germ_perc~chillweeks*force,Av)
Td.mod.lm<-lm(germ_perc~chillweeks*force,Td)

CC.break.lm<-segmented(CC.mod.lm,seg.Z=~chillweeks,npsi=1,seed=613)
PV.break.lm<-segmented(Pv.mod.lm,seg.Z=~chillweeks,npsi=2,seed=613)
ED.break.lm<-segmented(Ed.mod.lm,seg.Z=~chillweeks,npsi=1,seed=613)
HM.break.lm<-segmented(Hm.mod.lm,seg.Z=~chillweeks,npsi=1,seed=613)
SS.break.lm<-segmented(Ss.mod.lm,seg.Z=~chillweeks,npsi=1,seed=613)
AV.break.lm<-segmented(Av.mod.lm,seg.Z=~chillweeks,npsi=1,seed=613)
TD.break.lm<-segmented(Td.mod.lm,seg.Z=~chillweeks,npsi=1,seed=613)
brewer.pal(n=7, "Dark2")
summary(CC.break.lm)[9] #adjr2 from 0.36-0.89
summary(PV.break.lm)[9]
summary(ED.break.lm)[9]
summary(HM.break.lm)[9]

summary(AV.break.lm)[9]
summary(TD.break.lm)[9]

#uuper break point
summary(CC.break.lm) #7.552
summary(PV.break.lm) #7.2
summary(ED.break.lm) #7
summary(HM.break.lm) #2.36
summary(AV.break.lm) #6.5
summary(TD.break.lm) # 7.27

###most break around 7.5

###Run a model to show that chill doesn't generally affect final germ percentage after 7.5 weeks 
realdat.aftr<-filter(realdat, chillweeks>7) ## run this model to confirm chillweeks dont affect final germ perc after 7 weeks
ger.perc.upper<-brms::brm(germ_perc~chillweeks*force+(1|Taxa),data=realdat.aftr)
summary(ger.perc.upper)

####check individual species for same as above
Cc.maxy<-dplyr::filter(Cc,chillweeks>7.5)
range(Cc.maxy$germ_perc)
summary(lm(germ_perc~chillweeks*force,Cc.maxy))
##
Pv.maxy<-dplyr::filter(Pv,chillweeks>7.2)
summary(lm(germ_perc~chillweeks*force,Pv.maxy))##no effect of chilll

Ed.maxy<-dplyr::filter(Ed,chillweeks>7)
summary(lm(germ_perc~chillweeks*force,Ed.maxy))

### Now investigate germination time
d<-read.csv("survival_analysis/surival_dat_nointerval.csv")
d$chill_time<-NA
d<- within(d, chill_time[COLD=="0" ]<-0)
d<- within(d, chill_time[COLD=="A" ]<-14)
d<- within(d, chill_time[COLD=="B" ]<-28)
d<- within(d, chill_time[COLD=="C" ]<-35)
d<- within(d, chill_time[COLD=="D" ]<-42)
d<- within(d, chill_time[COLD=="E" ]<-49)
d<- within(d, chill_time[COLD=="f" ]<-56)
d<- within(d, chill_time[COLD=="G" ]<-63)
d<- within(d, chill_time[COLD=="H" ]<-77)
d<- within(d, chill_time[COLD=="i" ]<-91)
d$chillweeks<-d$chill_time/7 # make chilling weeks instead of days

d$force<-NA # make forcing numeric
d<- within(d, force[INC=="L"]<-0)
d<- within(d, force[INC=="H"]<-1)
d$censored<-ifelse(d$germinated==1,0,1)
d$censored<-ifelse(d$DAY==0.001 & d$germinated==1,-1,d$censored)

d.aftr<-filter(d,chillweeks>7.5) ## ronly after 7.5 weeks
d.aftr.1<-filter(d.aftr,germinated==1) #remove spcies that didnt germinate
ger.time.upper<-brms::brm(DAY~chillweeks*force+(1|Taxa),data=d.aftr.1) ### effectively a mgt model
summary(ger.time.upper)

priorz<-brms::get_prior(DAY | cens(censored)~chillweeks*force+(1|Taxa),data=d.aftr,family= lognormal, inits = "0")
ger.survial.upper<-brms::brm(DAY | cens(censored)~chillweeks*force+(1|Taxa), 
    data=d.aftr, family =lognormal, inits = "0" ,prior=priorz,iter=6000,warmup = 5000, chains=4) ##
summary(ger.survial.upper)
Cc.time<-filter(d,Taxa=="Cryptotaenia canadensis")
Cc.time<-filter(Cc.time, chillweeks>7.5)
Cc.time<-filter(Cc.time,germinated==1)
summary(lm(DAY~chillweeks*force,data=Cc.time))## yes

Pv.time<-filter(d,Taxa=="Polygonum virginiatum")
Pv.time<-filter(Pv.time, chillweeks>7.2)
Pv.time<-filter(Pv.time,germinated==1)
summary(lm(DAY~chillweeks*force,data=Pv.time))## yes

Ed.time<-filter(d,Taxa=="Eurbia diviricata")
Ed.time<-filter(Ed.time, chillweeks>7)
Ed.time<-filter(Ed.time,germinated==1)
summary(lm(DAY~chillweeks*force,data=Ed.time))## yes

Hm.time<-filter(d,Taxa=="Hesperis matronalis")
Hm.time<-filter(Hm.time, chillweeks>2.3)
Hm.time<-filter(Hm.time,germinated==1)
summary(lm(DAY~chillweeks*force,data=Hm.time))## yes

Td.time<-filter(d,Taxa=="Thalictrum dioicum")
Td.time<-filter(Td.time, chillweeks>7.2)
Td.time<-filter(Td.time,germinated==1)
summary(lm(DAY~chillweeks*force,data=Td.time))


range(Pv.maxy$germ_perc)
summary(lm(germ_perc~chillweeks*force,Pv.maxy))


plot.segmented(CC.break.lm,conf.level =.9,shade=FALSE,col="#D95F02",ylim=c(0,1.2),lwd=4,ylab="Effect of stratification on germination %",xlab="Weeks of stratification")
plot.segmented(PV.break.lm,conf.level =.9,add=TRUE,col="#66A61E",lwd=4)
plot.segmented(ED.break.lm,conf.level =.9,add=TRUE,col="#7570B3",lwd=4)
plot.segmented(HM.break.lm,conf.level =.9,add=TRUE,col="#E7298A",lwd=4)
plot.segmented(SS.break.lm,conf.level =.9,add=TRUE,col="#E6AB02",lwd=4)
plot.segmented(AV.break.lm,conf.level =.9,add=TRUE,col="#1B9E77",lwd=4)
plot.segmented(TD.break.lm,conf.level =.9,add=TRUE,col="#A6761D",lwd=4)


###nonlinerar

CC.mod.glm<-glm(germ_perc~chillweeks*force,family=binomial,Cc)
Pv.mod.glm<-glm(germ_perc~chillweeks*force,family=binomial,Pv)
Ed.mod.glm<-glm(germ_perc~chillweeks*force,family=binomial,Ed)
Hm.mod.glm<-glm(germ_perc~chillweeks*force,family=binomial,Hm)
Ss.mod.glm<-glm(germ_perc~chillweeks*force,family=binomial,Ss)
Av.mod.glm<-glm(germ_perc~chillweeks*force,family=binomial,Av)
Td.mod.glm<-glm(germ_perc~chillweeks*force,family=binomial,Td)

CC.break.glm<-segmented(CC.mod.glm,seg.Z=~chillweeks,npsi=2,seed=613)
PV.break.glm<-segmented(Pv.mod.glm,seg.Z=~chillweeks,npsi=2,seed=613)
ED.break.glm<-segmented(Ed.mod.glm,seg.Z=~chillweeks,npsi=2,seed=613)
HM.break.glm<-segmented(Hm.mod.glm,seg.Z=~chillweeks,npsi=2,seed=613)
SS.break.glm<-segmented(Ss.mod.glm,seg.Z=~chillweeks,npsi=2,seed=613)
AV.break.glm<-segmented(Av.mod.glm,seg.Z=~chillweeks,npsi=2,seed=613)
TD.break.glm<-segmented(Td.mod.glm,seg.Z=~chillweeks,npsi=2,seed=613)
brewer.pal(n=7, "Dark2")
#uuper break point
summary(CC.break.glm)#7.03
summary(PV.break.glm) #7.095
summary(ED.break.glm)#7
summary(HM.break.glm) #9.3
summary(AV.break.glm) #6.5
summary(TD.break.glm) # 7.27


plot.segmented(CC.break.glm,conf.level =.9,shade=FALSE,col="#D95F02",ylim=c(-6,7),lwd=4)
plot.segmented(PV.break.glm,conf.level =.9,add=TRUE,col="#66A61E",lwd=4)
plot.segmented(ED.break.glm,conf.level =.9,add=TRUE,col="#7570B3",lwd=4)
plot.segmented(HM.break.glm,conf.level =.9,add=TRUE,col="#E7298A",lwd=4)
plot.segmented(SS.break.glm,conf.level =.9,add=TRUE,col="#E6AB02",lwd=4)
plot.segmented(AV.break.glm,conf.level =.9,add=TRUE,col="#1B9E77",lwd=4)
plot.segmented(TD.break.glm,conf.level =.9,add=TRUE,col="#A6761D",lwd=4)
?plot.segmented()
