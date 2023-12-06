######## May 5 emprical plots for priorty effect project

rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()
library(drc)
library(dplyr)
if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")


library(ggplot2)
realdat<-read.csv("input/daily_dat_nointerval.csv")
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
realdat<- dplyr::filter(realdat,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea", "Carex grayi"))
realdat$germ_perc<-NA
realdat$germ_perc<-realdat$germ_num/realdat$tot_seed
realdat$germ_perc<-ifelse(realdat$germ_perc>1,1,realdat$germ_perc)
range(realdat$germ_perc)


unique(realdat$Taxa)
realdat$taxa<-NA
realdat<- within(realdat, taxa[Taxa=="Oenethera biennis" ]<-"Oenothera \nbiennis")
realdat<- within(realdat, taxa[Taxa=="Cryptotaenia canadensis" ]<-"Cryptotaenia \ncanadensis")
realdat<- within(realdat, taxa[Taxa=="Hesperis matronalis" ]<-"Hesperis \nmatronalis")
realdat<- within(realdat, taxa[Taxa=="Polygonum virginiatum" ]<-"Persicaria \nvirginiana")
realdat<- within(realdat, taxa[Taxa=="Asclepias syriaca" ]<-"Asclepias \nsyriaca")
realdat<- within(realdat, taxa[Taxa=="Silene vulgaris" ]<-"Silene \nvulgaris")
realdat<- within(realdat, taxa[Taxa=="Silene stellata" ]<-"Silene \nstellata")
realdat<- within(realdat, taxa[Taxa=="Eurbia diviricata" ]<-"Eurybia \ndivaricata")
realdat<- within(realdat, taxa[Taxa=="Thalictrum dioicum" ]<-"Thalictrum \ndioicum")
realdat<- within(realdat, taxa[Taxa=="Anemone virginana" ]<-"Anemone \nvirginiana")

realdat$incubation<-ifelse(realdat$INC=="L","20/10","25/15")

jitter <- position_jitter(width = 0.1, height = 0.01)

jpeg("figures/priority_effect_model/germ_courses.jpeg",height=10,width=11,unit='in',res=300)
ggplot(realdat,aes(DAY,germ_perc))+
  geom_smooth(aes(group=incubation),span=.85,size=.3,se=FALSE)+
  stat_summary(aes(shape=incubation),position=jitter,size=0.2)+
  facet_grid(taxa~chillweeks)+ggthemes::theme_few(base_size = 11)+ scale_x_continuous(breaks=c(0,25))+scale_y_continuous(breaks=c(0,.5,1))+
  ylab("germination %")+xlab("day of experiment")+ theme(strip.text.y = element_text(face = "italic"))

dev.off()


unique(realdat$Taxa)
drc1<-dplyr::filter(realdat, Taxa %in% c("Eurbia diviricata","Polygonum virginiatum"))        

#drc2<-filter(realdat, !Taxa %in% c("Polygonum virginiatum", "Anemone virginana","Cryptotaenia canadensis",
#                                  "Eurbia diviricata","Thalictrum dioicum","Silene stellata" ))

unique(drc1$chillweeks)
drc2<-dplyr::filter(drc1, chillweeks %in% c(4))    
drc2<-dplyr::filter(drc2, INC %in% c("L")) 

drc3<-dplyr::filter(drc1, chillweeks %in% c(4))    
drc3<-dplyr::filter(drc3, INC %in% c("H")) 

drc4<-dplyr::filter(drc1, chillweeks %in% c(8))    
drc4<-dplyr::filter(drc4, INC %in% c("H")) 

drc5<-dplyr::filter(drc1, chillweeks %in% c(8))    
drc5<-dplyr::filter(drc5, INC %in% c("L")) 

drc6<-dplyr::filter(drc1, chillweeks %in% c(13))    
drc6<-dplyr::filter(drc6, INC %in% c("H")) 

drc7<-dplyr::filter(drc1, chillweeks %in% c(13))    
drc7<-dplyr::filter(drc7, INC %in% c("L")) 






sds<-drm(germ_perc~DAY,factor(taxa), data=drc2,fct = LL.3(fixed = c(NA,NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
#sds2<-drm(germ_perc~DAY,factor(taxa), data=drc3,fct = LL.3(fixed = c(NA,NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
#sds3<-drm(germ_perc~DAY,factor(taxa), data=drc4,fct = LL.3(fixed = c(NA,NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
sds4<-drm(germ_perc~DAY,factor(taxa), data=drc5,fct = LL.3(fixed = c(NA,NA, NA), names = c("b", "gmax", "e50")), type ="continuous")


sdsa<-drm(germ_perc~DAY,factor(taxa), data=drc2,fct = LL.3(fixed = c(NA,.8, NA), names = c("b", "gmax", "e50")), type ="continuous")

sds4a<-drm(germ_perc~DAY,factor(taxa), data=drc5,fct = LL.3(fixed = c(NA,.8, NA), names = c("b", "gmax", "e50")), type ="continuous")


par(mfrow=c(2,2))
par(mar = c(4.1, 4.1, .5, .5))
plot(sds,log="",xlim=c(0,120),cex=.3,cex.axis = .8,ylim=c(0,1),xt=c(20,40,60,80,100,120),ylab="germination %",xlab="time (days)",legend=FALSE)
#plot(sds2,log="",xlim=c(0,50),ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50))
plot(sds4,log="",xlim=c(0,120),cex=.3,ylim=c(0,1),xt=c(20,40,60,80,100,120),ylab="germination %",xlab="time (days)",legend=FALSE)


plot(sdsa,log="",xlim=c(0,120),cex=.3,cex.axis = .8,ylim=c(0,1),xt=c(20,40,60,80,100,120),ylab="germination %",xlab="time (days)",legend=FALSE)
#plot(sds2,log="",xlim=c(0,50),ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50))
plot(sds4a,log="",xlim=c(0,120),cex=.3,ylim=c(0,1),xt=c(20,40,60,80,100,120),ylab="germination %",xlab="time (days)",legend=FALSE)

coef(sds)
coef(sdsa)

coef(sds4a)
coef(sds4)
dev.off()



sds5<-drm(germ_perc~DAY,factor(taxa), data=drc6,fct = LL.3(fixed = c(NA,NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
sds6<-drm(germ_perc~DAY,factor(taxa), data=drc7,fct = LL.3(fixed = c(NA,NA, NA), names = c("b", "gmax", "e50")), type ="continuous")

jpeg("figures/priority_effect_model/examp_Ed_andPv.jpeg",height=5,width=5,unit='in',res=200)
par(mfrow=c(2,1))
par(mar = c(4.1, 4.1, .5, .5))
plot(sds,log="",xlim=c(0,60),cex=.3,cex.axis = .8,ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50,55,60),ylab="germination %",xlab="time (days)",legend=FALSE)
#plot(sds2,log="",xlim=c(0,50),ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50))
plot(sds4,log="",xlim=c(0,60),cex=.3,ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50,55,60),ylab="germination %",xlab="time (days)",legend=FALSE)

dev.off()



plot(sds3,log="",xlim=c(0,50),ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50),legend=FALSE)
plot(sds6,log="",xlim=c(0,50),ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50),legend=FALSE)
plot(sds5,log="",xlim=c(0,50),ylim=c(0,1),xt=c(5,10,15,20,25,30,35,40,45,50),legend=FALSE)





library(tidyverse)
bb<-as.data.frame(coef(b))
goo<-summary(b)
goo<-as.data.frame(goo[[3]])
goo<-rownames_to_column(goo,"descriptor")                  
#?separate()
goo<-separate(goo,col = descriptor,sep=":",c("param","chill_weeks","Taxa"))
goo$chill_weeks<-as.numeric(goo$chill_weeks)
goo2<-filter(goo,param=="gmax")
goo3<-filter(goo,param=="e50")
 ?resp_mi()

ggplot(goo2,aes(chill_weeks,Estimate))+geom_point(aes(color=Taxa))+geom_smooth(method="lm",aes(color=Taxa))
ggplot(goo3,aes(chill_weeks,Estimate))+geom_point(aes(color=Taxa))+geom_smooth(method="lm",aes(color=Taxa))




library(brms)
##mi model
colnames(goo3)[5]<-"se"
colnames(goo2)[5]<-"se"


gmax.mod<-brm(
  bf(Estimate|se(se,sigma=TRUE)~chill_weeks*Taxa,
     phi ~chill_weeks*Taxa),
  data = goo2,
  family = Beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr")



get_prior(Estimate|se(se,sigma=TRUE)~chill_weeks*Taxa,data=goo2)


gmax.mod<-brm(Estimate|se(se,sigma=TRUE)~chill_weeks*Taxa,data=goo2,control=list(adapt_delta=.99),warmup=3000,iter=4000)

median(goo2$se)

newbers<-data.frame(Taxa=rep(unique(goo3$Taxa),2),se=(.37),chill_weeks=c(4,4,4,4,10,10,10,10))
newbers2<-data.frame(Taxa=rep(unique(goo3$Taxa),2),se=(.03),chill_weeks=c(4,4,4,4,10,10,10,10))
ePred<-predict(e50.mod,newdata = newbers,probs = c(.25,.75,.05,.95))
gPred<-predict(gmax.mod,newdata = newbers2,probs = c(.25,.75,.05,.95))

plot<-cbind(newbers,ePred)
plot2<-cbind(newbers2,gPred)

p1<-ggplot(plot,aes(chill_weeks,Estimate))+geom_point(aes(color=Taxa),size=3.5)+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=Taxa),width=0)+
  geom_smooth(method="lm",aes(color=Taxa),size=0.2)+scale_x_continuous(name="weeks of chilling",breaks=c(4,10))+
  scale_y_continuous(name = "t50",breaks=c(0,2,4,6,8,10,12,14,16,18))+ggthemes::theme_few()+geom_hline(yintercept=0)+scale_color_viridis_d(option="B",begin = .2,end=.8)

p2<-ggplot(plot2,aes(chill_weeks,Estimate))+geom_point(aes(color=Taxa),size=3)+geom_errorbar(aes(ymin=`Q25`,ymax=`Q75`,color=Taxa),width=0)+
  geom_smooth(method="lm",aes(color=Taxa),size=0.2)+scale_x_continuous(name="weeks of chilling",breaks=c(4,10))+
  scale_y_continuous(name = "Gmax")+ggthemes::theme_few()+geom_hline(yintercept=0)+scale_color_viridis_d(option="B",begin = .2,end=.8)

ggpubr::ggarrange(p2,p1,common.legend=TRUE,legend="right")


library(tidybayes)

