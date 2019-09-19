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
realdat<- within(realdat, force[INC=="H"]<-1)

realdat$DAY<-ifelse(realdat$DAY==0,0.0001,realdat$DAY) #elimiate 0 values for log logistic dist

realdatshorty<- filter(realdat,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea"))
#jpeg("figures/raw_germy.jpeg",res=200,width=1500,height=800)
#ploty<-ggplot(realdatshorty,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) #plot fake data
#ploty+geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#dev.off()

candies<-filter(realdatshorty, Taxa %in% c("Hesperis matronalis", "Asclepias syriaca"))
#jpeg("figures/woodspecies.jpeg",res=200,width=1500,height=800)
candies<-filter(candies,chillweeks==5)

ggplot(candies,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
drm(germ_perc~DAY,factor(INC):factor(Taxa), data=candies,fct = LL.3(), type ="continuous")

candies2<-filter(realdatshorty, Taxa %in% c("Polygonum virginiatum", "Eurbia diviricata"))
ggplot(candies2,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
candies2<-filter(candies2,chillweeks %in% c(6,9))
drm(germ_perc~DAY,factor(INC):factor(Taxa), data=candies2,fct = LL.3(), type ="continuous")


candies3<-filter(realdatshorty, Taxa %in% c("Cryptotaenia canadensis", "Polygonum virginiatum"))
ggplot(candies3,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
candies3<-filter(candies3,chillweeks %in% c(8))
#candies3<-filter(candies3,force==0)
drm(germ_perc~DAY,factor(force):factor(Taxa), data=candies3,fct = LL.3(), type ="continuous")

candies4<-filter(realdatshorty, Taxa %in% c("Asclepias syriaca","Silene vulgaris"))
ggplot(candies4,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
candies4<-filter(candies4,chillweeks %in% c(6))
drm(germ_perc~DAY,factor(INC):factor(Taxa), data=candies4,fct = LL.3(), type ="continuous")
#dev.off()
#"Hesperis matronalis", "Asclepias syriaca","Oenethera biennis","Eurbia diviricata","Hesperis matronalis","Cryptotaenia canadensis","Eurbia diviricata"
specieslist<-sort(unique(realdat$Taxa))
X<-split(realdat, with(realdat, realdat$Taxa), drop = TRUE)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:20]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)

crypto.cold<-filter(`Cryptotaenia canadensis`,INC=="L")

data.cryto.cold<-with(crypto.cold,
               list(Y=germ_perc,
                    t=DAY,
                    chill=chillweeks,
                    N=nrow(crypto.cold)
               )
)


mod.crypto.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.cryto.cold, 
                    iter = 9000, warmup=7000 , chain=4) ## run model
summary(mod.crypto.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

crypto.warm<-filter(`Cryptotaenia canadensis`,INC=="H")

data.cryto.warm<-with(crypto.warm,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(crypto.warm)
                      )
)

mod.crypto.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.cryto.warm, 
                      iter = 4000, warmup=3000 , chain=4) ## run model
summary(mod.crypto.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]
