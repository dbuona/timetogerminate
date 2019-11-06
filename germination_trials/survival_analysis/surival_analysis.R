###a script to amke my data survival formated
###need to add a way to elminate seeds that rotted.
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(survival)
library(ggplot2)
library(dplyr)
library("ggfortify")
library("coxme")
library("icenReg")
library("emmeans")
library(drc)
library(lme4)
library(brms)
library(ggstance)
setwd("~/Documents/git/timetogerminate/germination_trials/input")

d<-read.csv("..//survival_analysis/surival_dat_nointerval.csv")
d$DAY<-ifelse(d$DAY==0,0.00001,d$DAY)
realdat<-d
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
realdat<- within(realdat, force[INC=="H"]<-1)
###brms is backwords from other survival
realdat$censored<-ifelse(realdat$germinated==1,0,1)
realdat<- filter(realdat,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea"))
unique(realdat$DAY)
realdat$censored<-ifelse(realdat$DAY==0.00001 & realdat$germinated==1,-1,realdat$censored)



priorz<-get_prior(DAY | cens(censored)~chillweeks*force+(chillweeks*force|Taxa),data=realdat,family= lognormal, inits = "0")

fit2 <- brm(DAY | cens(censored)~chillweeks*force+(chillweeks*force|Taxa), 
            data=realdat, family =   lognormal, inits = "0" ,prior=priorz,iter=8000,warmup = 7000, chains=4) ##

summary(fit2)
ranef(fit2)
WAIC(fit2)

###predict and plot
Field<-c("Oenethera biennis","Silene vulgaris","Asclepias syriaca")
Forest<-c("Anemone virginana","Cryptotaenia canadensis","Eurbia diviricata","Hesperis matronalis","Polygonum virginiatum","Silene stellata","Thalictrum dioicum")
new.data <- data.frame(Taxa = factor(c(rep(unique(realdat$Taxa),10))),
                                chillweeks = c(rep(unique(realdat$chillweeks),each=11)),
                                force = c(rep(c(0,1),each=110)))
daty<-(predict(fit2, newdata = new.data))
daty<-cbind(daty,new.data)
Forestsp<-dplyr::filter(daty,Taxa %in% c(Forest))
Forestsp<-dplyr::filter(Forestsp,chillweeks %in% c(0,4,6,9,13))

pd=position_dodgev(height=20)
jpeg("AFT.forest.jpeg",width = 10, height = 6, units = 'in', res=300)
ggplot(Forestsp,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
xlim(1,120)+ylim(-40,40)+geom_point(position=pd,aes(color=Taxa),size=4)+
  geom_vline(aes(xintercept=Estimate, color= Taxa))+
facet_grid(as.factor(force)~as.factor(chillweeks))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")
dev.off()

jpeg("AFT.forest.zoom.jpeg",width = 10, height = 6, units = 'in', res=300)
ggplot(Forestsp,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
  xlim(0,30)+ylim(-40,40)+geom_point(position=pd,aes(color=Taxa),size=4)+
  geom_vline(aes(xintercept=Estimate, color= Taxa))+
  facet_grid(as.factor(force)~as.factor(chillweeks))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")
dev.off()


jpeg("AFT.forestsp.change.jpeg",width = 12, height = 6, units = 'in', res=300)
ggplot(Forestsp,aes(x=Estimate,y=0, color=as.factor(chillweeks)))+
  xlim(-10,90)+ylim(-40,40)+geom_point(aes(color=as.factor(chillweeks)),size=4)+
  geom_vline(aes(xintercept=Estimate, color= as.factor(chillweeks)))+
  facet_grid(as.factor(force)~Taxa)+theme_linedraw()+
  scale_colour_brewer( type = "seq", palette = 16, direction = 2)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")
dev.off()

ggplot(Forestsp,aes(x=chillweeks,y=Estimate, color=Taxa))+
  geom_smooth(method="lm",level=.5)+ylim(-100,1800)+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  facet_wrap(~as.factor(force))+theme_linedraw()

Fieldtsp<-dplyr::filter(daty,Taxa %in% c(Field))
Fieldtsp<-dplyr::filter(Fieldtsp,chillweeks %in% c(0,4,6,9,13))
ggplot(Fieldtsp,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
  xlim(-10,120)+ylim(-40,80)+geom_point(aes(color=Taxa),size=4)+
  geom_vline(aes(xintercept=Estimate, color= Taxa))+
  facet_grid(as.factor(force)~as.factor(chillweeks))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = 1, direction = 1)


priorzPv<-get_prior(DAY | cens(censored)~chillweeks*force,data=Pv,family= lognormal, inits = "0")

fitPv <- brm(DAY | cens(censored)~chillweeks*force, 
            data=Pv, family =   lognormal, inits = "0" ,prior=priorz,iter=6000,warmup = 5000, chains=2) ##
summary(fitPv)


###nonlinear mixed effect model: 
#1ignore those that didn't germinate

head(brms::kidney)
fit1 <- brm(time | cens(censored) ~ age + sex + disease, 
            data = brms::kidney, family = weibull, inits = "0")
summary(fit1) 
plot(fit1)
head(Cryp.C)
