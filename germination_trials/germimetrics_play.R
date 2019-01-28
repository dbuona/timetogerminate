rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")

library("brms")
library(tidyverse)
library(lme4)
library(broom)
library(ggstance)
load("GImodels")
d<-read.csv("daily_dat_nointerval.csv",header=TRUE)
goober<-d

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

goober$warmT<-NA
goober<- within(goober, warmT[INC=="H" ]<-25)
goober<- within(goober, warmT[INC=="L" ]<-20)

d<-goober

####calculate GI germination index Kader 2005
#d<-d %>% group_by(plate_num) %>%  mutate(weight = row_number())
#d<-d %>% group_by(plate_num)  %>%  mutate(weight = rev(weight))

###or is it better for the actualy
d<-d %>% group_by(plate_num)  %>%  mutate(weight = rev(DAY))

d$weighting<-d$weight*d$germ.daily
GIs<-d %>% group_by(Taxa,warmT,chill_time,COLD,INC,plate_num) %>% summarise(GI=sum(weighting))                                          
GI2<-GIs %>%group_by(Taxa,warmT,chill_time,COLD,INC) %>% summarise(mean=mean(GI))
GI3<-GIs %>%group_by(Taxa,warmT,chill_time,COLD,INC) %>% summarise(sd=sd(GI))

GIforplot<-left_join(GI2,GI3)

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
save.image("GImodels")
#prior_summary(mod.z)

summary(mod,prob=.95)
summary(mod.z)
mod.z
coef(mod.z)
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


d<-extract_coefs(mod.z)
l<-extract_CIlow(mod.z)
h<-extract_CIhigh(mod.z)

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
dat2$start<-dat2$Intercept+(20*dat2$inc.z)


### ploting
unique(GIs$Taxa)
GI2s<-filter(GIs,Taxa %in% c("Anemone virginana"     ,  "Asclepias syriaca" ,"Cryptotaenia canadensis","Eurbia diviricata","Hesperis matronalis","Oenethera biennis","Polygonum virginiatum" ,"Silene vulgaris","Silene stellata"))

GIwet<-filter(GIforplot,Taxa %in% c("Cryptotaenia canadensis","Hesperis matronalis","Polygonum virginiatum","Eurbia diviricata","Anemone virginana"))
GIdry<-filter(GIforplot,Taxa %in% c("Silene vulgaris","Asclepias syriaca"))

ggplot(GIwet,aes(reorder(Taxa,-mean),mean))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+geom_errorbar(aes(x=Taxa, ymin=mean-sd, ymax=mean+sd),width=0.4)+facet_grid(INC~COLD)+theme(axis.text.x = element_blank())
ggplot(GIdry,aes(reorder(Taxa,-mean),mean))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+geom_errorbar(aes(x=Taxa, ymin=mean-sd, ymax=mean+sd),width=0.4)+facet_grid(INC~COLD)+theme(axis.text.x = element_blank())

####germination percentage
### modeling for moldy seeds

