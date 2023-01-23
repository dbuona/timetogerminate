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
germ.perc$final_perc<-germ.perc$germ_num/germ.perc$tot_seed

#invasive.perc<-filter(germ.perc,Taxa %in% c("Hesperis matronalis","Cryptotaenia canadensis"))
#FGP<-invasive.perc %>% dplyr::group_by(Taxa,chillweeks,warmT) %>% summarise(meanFGP=mean(final_perc,na.rm=TRUE),sdFGP=sd(final_perc,na.rm=TRUE))




goodsp<-filter(goober, Taxa %in% c("Hesperis matronalis","Cryptotaenia canadensis"))

plot1<-ggplot(goodsp,aes(DAY,germ_perc))+stat_summary(aes(color=Taxa,shape=Taxa),size=.5)+ stat_smooth(aes(color=Taxa,fill=Taxa),size=0.3,se=FALSE)+
  facet_grid(chillweeks~incubation)+
  ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(begin=0,end=.5)+scale_fill_viridis_d(begin=0,end=.5)+xlab("Day of trial")+ylab("Germination percentatge")+
   theme(legend.text = element_text(face = "italic"),legend.position = "top")

germ.perc<-filter(germ.perc, Taxa %in% c("Hesperis matronalis","Cryptotaenia canadensis"))
germ.perc$final_perc<-ifelse(germ.perc$final_perc>1,1,germ.perc$final_perc)
plot1a<-ggplot(germ.perc,aes(chillweeks,final_perc*100))+geom_point(aes(color=Taxa,shape=Taxa),size=0.9)+stat_smooth(aes(color=Taxa,fill=Taxa),size=0.5,method="lm",se=FALSE)+facet_wrap(~incubation)+
  ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(begin=0,end=.5)+scale_fill_viridis_d(begin=0,end=.5)+xlab("")+ylab("Final germination percentatge")+
  theme(legend.text = element_text(face = "italic"))+theme(legend.position = "none")+scale_y_continuous(limits=c(0,100))


####calculate GI germination index Kader 2005
#d<-d %>% group_by(plate_num) %>%  mutate(weight = row_number())
#d<-d %>% group_by(plate_num)  %>%  mutate(weight = rev(weight))

###or is it better for the actualy
goodsp2<-goodsp %>% group_by(plate_num)  %>%  mutate(weight = rev(DAY))

goodsp2$weighting<-goodsp2$weight*goodsp2$germ.daily

GIs<-goodsp2 %>% group_by(Taxa,warmT,chillweeks,COLD,incubation,plate_num) %>% summarise(GI=sum(weighting))                                          
GI2<-GIs %>%group_by(Taxa,warmT,chillweeks,COLD,incubation) %>% summarise(mean=mean(GI))
GI3<-GIs %>%group_by(Taxa,warmT,chillweeks,COLD,incubation) %>% summarise(sd=sd(GI))

GIforplot<-left_join(GI2,GI3)
plot2<-ggplot(GIs,aes(chillweeks,GI))+geom_point(aes(color=Taxa),size=0.9)+geom_smooth(method="lm",aes(color=Taxa,fill=Taxa,shape=Taxa),size=0.5,se=FALSE)+facet_wrap(~incubation)+
  ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(begin=0,end=.5)+scale_fill_viridis_d(begin=0,end=.5)+xlab("")+ylab("Germination index")+theme(legend.position = "none")

ggplot(fgp.dat,aes(chillweeks,final_perc))+geom_point(aes(color=Taxa))+geom_smooth(method="lm",aes(color=Taxa,fill=Taxa),size=0.3)+facet_wrap(~INC)+
  ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(begin=0,end=.5)+scale_fill_viridis_d(begin=0,end=.5)+xlab("Weeks of cold stratification")+ylab("Germination index")+theme(legend.position = "none")


fgp.dat<-filter(goodsp,DAY==25)
#fgp.dat<-filter(fgp.dat,germ_perc>=.5)
fgp.dat<-dplyr::select(fgp.dat,plate_num,germ_num)
colnames(fgp.dat)[2]<-"finalgerm"
goodsp<-left_join(goodsp,fgp.dat)

plates<-unique(fgp.dat$plate_num)
df<-data.frame(plate_num=numeric(),Taxa=character(),INC=character(),chillweeks=numeric(),T50=numeric())

for (p in seq_along(plates)){
  dataoneplate <- dplyr::filter(goodsp, plate_num==plates[p])
  t50 <- sum(dataoneplate$DAY*dataoneplate$germ.daily)/dataoneplate$finalgerm
  dfhere<-data.frame(plate_num=plates[p],Taxa=unique(dataoneplate$Taxa),INC=unique(dataoneplate$incubation), chillweeks=unique(dataoneplate$chillweeks),T50=t50)
  df <- rbind(df, dfhere)
}

plot3<-ggplot(df,aes(chillweeks,T50))+geom_point(aes(color=Taxa,shape=Taxa),size=0.9)+stat_smooth(method="lm",aes(color=Taxa),size=0.5,se=FALSE)+facet_wrap(~INC)+
  ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(begin=0,end=.5)+scale_fill_viridis_d(begin=0,end=.5)+xlab("Weeks of cold stratification")+ylab("Mean germination time")+
  theme(legend.position = "none")


####MGT
MGT<- df %>% group_by(Taxa,chillweeks,INC)%>%summarize(MGT=mean(T50,na.rm=TRUE),sd=sd(T50,na.rm=TRUE))
MGT$MGT<-round(MGT$MGT,digits = 2)
MGT$sd<-round(MGT$sd,digits = 1)

MGT$MGTsd<-paste(MGT$MGT,MGT$sd,sep = " (")
MGT<-dplyr::select(MGT,Taxa,chillweeks,INC,MGTsd)

MGT<-spread(MGT,Taxa,MGTsd)



time<-goodsp %>% dplyr::group_by(plate_num,Taxa) %>% dplyr::summarize(MGT=mean(DAY))

germ.perc<-filter(germ.perc,Taxa %in% c("Hesperis matronalis","Cryptotaenia canadensis"))
germ.perc$final_perc<-ifelse(germ.perc$final_perc>1.0,1.0,germ.perc$final_perc)

FGP<-dplyr::select(germ.perc,Taxa,chillweeks,INC,final_perc)
FGP<- FGP %>% group_by(Taxa,chillweeks,INC)%>%summarize(FGP=mean(final_perc,na.rm=TRUE),sd=sd(final_perc,na.rm=TRUE))
FGP$FGP<-round(FGP$FGP,digits = 2)
FGP$sd<-round(FGP$sd,digits = 1)

FGP$FGPsd<-paste(FGP$FGP,FGP$sd,sep = " (")
FGP<-dplyr::select(FGP,Taxa,chillweeks,INC,FGPsd)

FGP<-spread(FGP,Taxa,FGPsd)


table.data<-cbind(FGP,MGT)
table.data<-table.data[c(1,2,3,4,7,8)]
write.csv(table.data,"descriptive.csv")
xtable(table.data)

ggplot(germ.perc,aes(chillweeks,final_perc))+geom_point(aes(color=Taxa))+stat_smooth(method="glm",aes(color=Taxa))+facet_wrap(~INC)+
  ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(begin=0,end=.5)+scale_fill_viridis_d(begin=0,end=.5)+xlab("Weeks of cold stratification")+ylab("FGP")+
  theme(legend.position = "none")+ylim(0,1)



 

plot4<-ggpubr::ggarrange(plot1a,plot3,plot2,ncol=3,nrow=1,labels = c("b)","c)","d)"))

jpeg("..//figures/crp_hesp2.jpeg",width=7,height=10,units='in',res=400)
plot1
dev.off()
jpeg("..//figures/crp_hesp1.jpeg",width=10,height=10,units='in',res=400)
ggpubr::ggarrange(plot1,plot4,nrow=2,ncol=1,common.legend = TRUE,heights=c(.9,.4),labels=c("a)"))
dev.off()
###center predictors
GIs$inc.cent<-GIs$warmT-mean(GIs$warmT)
GIs$chill.cent<-GIs$chill_time-mean(GIs$chill_time)

#z-score predictors
GIs$inc.z<-(GIs$warmT-mean(GIs$warmT))/sd(GIs$warmT)
GIs$chill.z<-(GIs$chill_time-mean(GIs$chill_time))/sd(GIs$chill_time)

#
mod<-brm(GI~warmT*chill_time+(warmT*chill_time|Taxa),data=GIs, iter=3000)

mod.z<-brm(GI~inc.z*chill.z+(inc.z*chill.z|Taxa),data=GIs,
  iter=3000)

#prior_summary(mod.z)

summary(mod,prob=.95)
summary(mod.z)
mod.z
coef(mod.z)
coef(mod)
#####BRROM

pd=position_dodgev(height=0.4)
ggplot(sults,aes(estimate,term))+geom_point(aes(color=level),position=pd)+geom_errorbarh(aes(color=level,xmin=lower,xmax=upper),position=pd,width=.1)


extract_coefs<-function(x){
gather(rownames_to_column(dplyr::select(as.data.frame(coef(x)),contains("Estimate")),"Taxa"),term,Estimate,2:5)
  }

extract_CIlow<-function(x){
  gather(rownames_to_column(dplyr::select(as.data.frame(coef(x)),contains("Q2.5")),"Taxa"),term,lower,2:5)
}

extract_CIhigh<-function(x){
  gather(rownames_to_column(dplyr::select(as.data.frame(coef(x)),contains("Q9")),"Taxa"),term,upper,2:5)
}


d<-extract_coefs(mod)
l<-extract_CIlow(mod)
h<-extract_CIhigh(mod)

d<-d %>% separate(term,c("xtra","term"),sep="Estimate.")
l<-l %>% separate(term,c("xtra","term"),sep="Q2.5.")
h<-h %>% separate(term,c("xtra","term"),sep="Q97.5.")

d<-d%>% dplyr::select(-xtra)
l<-l%>% dplyr::select(-xtra)
h<-h %>%dplyr::select(-xtra)

dat<-left_join(h,l)
dat<-left_join(dat,d)

daty<-filter(dat,term!="Intercept")

pd=position_dodgev(height=0.6)
ggplot(daty,aes(Estimate,term))+geom_point(aes(color=Taxa),position=pd)+geom_errorbarh(aes(color=Taxa,xmin=lower,xmax=upper),position=pd,width=.1)+geom_vline(xintercept = 0)+theme_bw()
ggplot(dat,aes(Estimate,term))+geom_point

dat2<-dat %>% dplyr::select(-upper,-lower)
  
dat2<-spread(dat2,term,Estimate)
dat2$start<-dat2$Intercept+(20*dat2$warmT)
dat2$end<-dat2$start+(dat2$chill_time*91)

ggplot(dat2,aes(0,start))+geom_segment(aes(x=0,xend=91,y=start, yend=end,col=Taxa,linetype=Taxa),size=1.5)+scale_linetype_manual(values = c(rep("solid", 7), rep("dotdash", 7)))+scale_color_manual(values = c(brewer.pal(7, "Set1"), brewer.pal(7 ,"Set1")))+theme_linedraw()+ggtitle("Projected GI values at 20 C")

### ploting
unique(GIs$Taxa)
GI2s<-filter(GIs,Taxa %in% c("Anemone virginana"     ,  "Asclepias syriaca" ,"Cryptotaenia canadensis","Eurbia diviricata","Hesperis matronalis","Oenethera biennis","Polygonum virginiatum" ,"Silene vulgaris","Silene stellata"))

GIwet<-filter(GIforplot,Taxa %in% c("Cryptotaenia canadensis","Hesperis matronalis","Polygonum virginiatum","Eurbia diviricata","Anemone virginana"))
GIdry<-filter(GIforplot,Taxa %in% c("Silene vulgaris","Asclepias syriaca"))

ggplot(GIwet,aes(reorder(Taxa,-mean),mean))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+geom_errorbar(aes(x=Taxa, ymin=mean-sd, ymax=mean+sd),width=0.4)+facet_grid(INC~COLD)+theme(axis.text.x = element_blank())
ggplot(GIdry,aes(reorder(Taxa,-mean),mean))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+geom_errorbar(aes(x=Taxa, ymin=mean-sd, ymax=mean+sd),width=0.4)+facet_grid(INC~COLD)+theme(axis.text.x = element_blank())


ggplot(germ.perc,aes(reorder(Taxa,-final_perc),final_perc))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+facet_grid(INC~COLD)+theme(axis.text.x = element_blank())
####germination percentage
### modeling for moldy seeds
save.image("GImodels")
