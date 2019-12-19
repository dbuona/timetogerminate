###a script to amke my data survival formated
###need to add a way to elminate seeds that rotted.
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(survival)
library(ggplot2)
library(dplyr)
library("ggfortify")
library("icenReg")
library("emmeans")
library(drc)
library(lme4)
library(brms)
library(ggstance)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
load("survmodel")

d<-read.csv("..//survival_analysis/surival_dat_nointerval.csv")
d$DAY<-ifelse(d$DAY==0,0.00001,d$DAY)


d$chillweeks<-d$chill_time/7 # make chilling weeks instead of days

d$force<-NA # make forcing numeric
d<- within(d, force[INC=="L"]<-0)
d<- within(d, force[INC=="H"]<-5)
###brms is backwords from other survival
d$censored<-ifelse(d$germinated==1,0,1)
d<- filter(d,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea"))

d$censored<-ifelse(d$DAY==0.00001 & d$germinated==1,-1,d$censored)


##try weibull

priorz.wei<-get_prior(DAY | cens(censored)~chillweeks+force+(chillweeks+force|Taxa),data=d,family= weibull())

fit.wei.all<- brm(DAY | cens(censored)~chillweeks+force+(chillweeks+force|Taxa), data=d, family =   weibull(),inits=0 ,prior=priorz.wei,iter=4000,warmup = 3000, chains=4) 




summary(fit.wei.all)
coef(fit.wei.all,probs =c(0.25,.75))

brmsfamily("weibull")

exp(3.009708) #74 days
exp(0.11582958)#1.12
exp(4.092365) ##59
exp(4.544279) ## 94

exp(4.307620-0.24874228*9) 7.9

pred.weeks<-c(4,8,12,16)
pred.force<-c(0,5)
unique(d$Taxa)
new.data <- data.frame(Taxa=c(rep(c("Cryptotaenia canadensis","Polygonum virginiatum","Eurbia diviricata","Hesperis matronalis"),each=4,2)),
  chillweeks = c(rep(pred.weeks,8)),
                       force = c(rep(pred.force,each=16)))

daty.wei.all<-predict(fit.wei.all,probs =c(0.25,.75),newdata=new.data)### something is wrong with error
daty.wei<-cbind(daty.wei.all,new.data)
dev.new()

pd<-position_dodge(width=0.4)
ggplot(daty.wei,aes(x=chillweeks,y=Estimate))+geom_point(aes(color=Taxa),position=pd)+
  ggplot2::geom_errorbar(aes(ymax=Q75,ymin=Q25,color=Taxa),width=0,position=pd)+
  facet_wrap(~as.factor(force))+theme_bw()+ylim(0,25)
 

  geom_vline(aes(xintercept=Estimate))+
  facet_grid(as.factor(force)~as.factor(chillweeks))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")

dev.off()


###predict and plot
Field<-c("Oenethera biennis","Silene vulgaris","Asclepias syriaca")
Forest<-c("Anemone virginana","Carex grayi","Cryptotaenia canadensis","Eurbia diviricata","Hesperis matronalis","Polygonum virginiatum","Silene stellata","Thalictrum dioicum")

daty.wei<-predict(fit.wei,probs =c(0.25,.75),newdata=new.data)### something is wrong with error
daty.wei<-cbind(daty.wei,new.data)

Forestsp<-dplyr::filter(daty,Taxa %in% c(Forest))
Forestsp.wei<-dplyr::filter(daty.wei,Taxa %in% c(Forest))

library("timelineS")

pd=position_dodgev(height=0)
jpeg("AFT.forest.jpeg",width = 10, height = 6, units = 'in', res=300)
ggplot(Forestsp,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
xlim(-20,120)+ylim(-5,20)+geom_point(position=pd,aes(color=Taxa),size=4)+
  geom_vline(aes(xintercept=Estimate, color= Taxa))+
facet_grid(as.factor(force)~as.factor(chillweeks))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")
dev.off()

ggplot(Forestsp.wei,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
  xlim(-0,35)+ylim(-5,20)+geom_point(position=pd,aes(color=Taxa),size=4)+
  geom_vline(aes(xintercept=Estimate, color= Taxa))+
  facet_grid(as.factor(force)~as.factor(chillweeks))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")

ggplot(Forestsp.wei,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
  xlim(-0,35)+ylim(-5,20)+geom_point(position=pd,aes(color=Taxa),size=4)+
  geom_vline(aes(xintercept=Estimate, color= Taxa))+geom_vline(aes(xintercept=Q2.5,color=Taxa),linetype="dotted")+
  facet_grid(as.factor(force)~as.factor(chillweeks))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")



jpeg("AFT.forest.zoom.jpeg",width = 10, height = 6, units = 'in', res=300)
ggplot(Forestsp,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
  xlim(0,35)+ylim(-10,20)+geom_point(position=pd,aes(color=Taxa),size=4)+
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




###not sure how to interpret these models on what scale but we need a cure funciton anyway

save.image("survmodel")
