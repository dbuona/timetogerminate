######## May 5 emprical plots for priorty effect project

rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()


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

jitter <- position_jitter(width = 0.1, height = 0.01)
ggplot(realdat,aes(DAY,germ_perc))+
  stat_summary(aes(group=chillweeks,color=chillweeks),position=jitter,size=0.1)+
  #geom_smooth(aes(group=chillweeks,color=chillweeks),span=.85,size=.5,se=FALSE)+
  facet_grid(INC~Taxa)+ggthemes::theme_few(base_size = 9)+ scale_color_distiller(,type = "seq",
                                                                           direction = 1,
                                                                       palette = "Greys")

unique(realdat$Taxa)
drc1<-filter(realdat, Taxa %in% c("Polygonum virginiatum", "Anemone virginana","Cryptotaenia canadensis",
                                  "Eurbia diviricata"))        

#drc2<-filter(realdat, !Taxa %in% c("Polygonum virginiatum", "Anemone virginana","Cryptotaenia canadensis",
#                                  "Eurbia diviricata","Thalictrum dioicum","Silene stellata" ))

unique(drc1$chillweeks)
drc1<-filter(drc1, chillweeks %in% c(4,5,6,7,8,9,11,13))    

drc1<-filter(drc1, INC %in% c("L")) 



b<-drm(germ_perc~DAY,factor(chillweeks):factor(Taxa), data=drc1,fct = LL.3(fixed = c(NA,NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
plot(b)

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

e50.mod<-brm(Estimate|se(se,sigma=TRUE)~chill_weeks*Taxa,data=goo3,control=list(adapt_delta=.99),warmup=3000,iter=4000)


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

