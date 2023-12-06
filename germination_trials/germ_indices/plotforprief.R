rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")

library("brms")
library(tidyverse)
library(lme4)
library(broom)
library(ggstance)
library(RColorBrewer)
#load("GImodels")
d.dat<-read.csv("daily_dat_nointerval.csv",header=TRUE)

goober<-d.dat

goober$chill_time<-NA
goober<- within(goober, chill_time[COLD=="0" ]<-0)
goober<- within(goober, chill_time[COLD=="A" ]<-14)
goober<- within(goober, chill_time[COLD=="B" ]<-28)
goober<- within(goober, chill_time[COLD=="C" ]<-35)
goober<- within(goober, chill_time[COLD=="D" ]<-42)
goober<- within(goober, chill_time[COLD=="E" ]<-49)
goober<- within(goober, chill_time[COLD=="f" ]<-56)
goober<- within(goober, chill_time[COLD=="G" ]<-63)
goober<- within(goober, chill_time[COLD=="H" ]<-77)
goober<- within(goober, chill_time[COLD=="i" ]<-91)
goober$chillweeks<-goober$chill_time/7

goober$warmT<-NA
goober<- within(goober, warmT[INC=="H" ]<-25)
goober<- within(goober, warmT[INC=="L" ]<-20)

goober$germ_perc<-NA
goober$germ_perc<-goober$germ_num/goober$tot_seed
goober$germ_perc<-ifelse(goober$germ_perc>1,1,goober$germ_perc)
goober$incubation<-ifelse(goober$INC=="L","cool","warm")

germ.perc<-goober %>% filter(DAY==25)
#germ.perc$germ_perc<-ifelse(germ.perc$germ_perc==1,.9999999999,germ.perc$germ_perc)
fgp.dat<-dplyr::select(germ.perc,plate_num,germ_num)
colnames(fgp.dat)[2]<-"finalgerm"
goober<-left_join(goober,fgp.dat)
plates<-unique(germ.perc$plate_num)

df<-data.frame()

for (p in seq_along(plates)){
  dataoneplate <- dplyr::filter(goober, plate_num==plates[p])
  MGT <- sum(dataoneplate$DAY*dataoneplate$germ.daily)/unique(dataoneplate$finalgerm)
  dfhere<-data.frame(plate_num=plates[p],Taxa=unique(dataoneplate$Taxa),INC=unique(dataoneplate$incubation), chillweeks=unique(dataoneplate$chillweeks),MGT=MGT)
  df <- rbind(df, dfhere)
}






#germ.perc$MGT<-ifelse(germ.perc$MGT=="Inf",NA,germ.perc$MGT)
#germ.perc$T50<-ifelse(germ.perc$T50=="Inf",NA,germ.perc$MGT)
germ.perc<-left_join(germ.perc,df,by="plate_num")
round(cor(germ.perc$germ_perc,germ.perc$MGT,use = "pairwise.complete.obs",method = "pearson"),2)

goodsps<-filter(germ.perc,!Taxa.x%in%c("Impatiens capensis","Carex grisea","Phlox cuspidata","Carex grayi"))#,"Oenethera biennis","Silene vulgaris","Asclepias syriaca","Hesperis matronalis","Anemone virginana"))
round(cor(goodsps$germ_perc,goodsps$MGT,use = "pairwise.complete.obs",method = "pearson"),2)

a<-ggplot(goodsps,aes(chillweeks.x,germ_perc))+stat_summary(aes(shape=incubation))+facet_wrap(~Taxa.x,nrow=1)+geom_smooth(se = FALSE,span=1.1)+ggthemes::theme_few()
b<-ggplot(goodsps,aes(chillweeks.x,MGT))+stat_summary(aes(shape=incubation))+
  facet_wrap(~Taxa.x,nrow=2)+geom_smooth(se = FALSE,span=1.9,aes(color=incubation))+
  ggthemes::theme_few()+xlab("weeks of chilling exposure")+theme(legend.position = "top")+theme(legend.title = element_blank())

a<-ggplot(goodsps,aes(chillweeks.x,MGT))+stat_summary(aes(color=Taxa.x))+geom_smooth(se = FALSE,span=1.1,aes(color=Taxa.x))+ggthemes::theme_few()+scale_color_viridis_d(option="turbo")
b<-ggplot(goodsps,aes(chillweeks.x,germ_perc))+stat_summary(aes(color=Taxa.x))+
  geom_smooth(se = FALSE,span=0.9,aes(color=Taxa.x))+
  ggthemes::theme_few()+scale_color_viridis_d(option="turbo")



c<-ggplot(goodsps,aes(MGT,germ_perc))+geom_point(aes(color=chillweeks.x))+facet_wrap(~Taxa.x,ncol=2)+geom_smooth(method="lm")+scale_x_discrete()


comp<-filter(goodsps,Taxa.x %in% c("Cryptotaenia canadensis","Polygonum virginiatum"))

comp<-filter(comp, chillweeks.x %in% c(5,13))
cc<-ggplot(comp,aes(chillweeks.x,MGT))+stat_summary(aes(color=Taxa.x))+facet_grid(incubation~chillweeks.x,scales="free_x")+
  scale_color_viridis_d()+ggthemes::theme_few()+scale_x_discrete()+scale_y_continuous(n.breaks = 7)+xlab("")+
  theme(legend.position = "top")+theme(legend.title = element_blank())


goodsps2<-filter(germ.perc,Taxa.x%in%c("Cryptotaenia canadensis","Thalictrum dioicum","Polygonum virginiatum","Eurbia diviricata","Silene stellata"))#,"Oenethera biennis","Silene vulgaris","Asclepias syriaca","Hesperis matronalis","Anemone virginana"))
ggplot(goodsps2,aes(chillweeks.x,MGT))+stat_summary(aes(color=Taxa.x))+facet_grid(incubation~chillweeks.x,scales="free_x")+
  scale_color_viridis_d()+ggthemes::theme_few()+scale_x_discrete()+scale_y_continuous(n.breaks = 7)+xlab("")+
  theme(legend.position = "top")+theme(legend.title = element_blank())


jpeg("..//figures/empricalplots2.jpeg",width=11, height=7,unit='in',res=200)
ggpubr::ggarrange(b,cc,ncol=2,widths=c(.7,.3),common.legend = FALSE,labels = c("a)","b)"))
  
dev.off()
fgp<-brm(
  bf(germ_perc ~ chillweeks*incubation+(chillweeks*incubation|Taxa),
     phi ~chillweeks*incubation+(chillweeks*incubation|Taxa),
zi ~1),
  data = germ.perc,
  family = zero_inflated_beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,backend = "cmdstanr")

coef(fgp)
