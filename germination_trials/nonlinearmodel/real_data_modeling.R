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

#candies<-filter(realdatshorty, Taxa %in% c("Hesperis matronalis", "Asclepias syriaca"))
#jpeg("figures/woodspecies.jpeg",res=200,width=1500,height=800)
#candies<-filter(candies,chillweeks==5)

#ggplot(candies,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#drm(germ_perc~DAY,factor(INC):factor(Taxa), data=candies,fct = LL.3(), type ="continuous")

#candies2<-filter(realdatshorty, Taxa %in% c("Polygonum virginiatum", "Eurbia diviricata"))
#ggplot(candies2,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#candies2<-filter(candies2,chillweeks %in% c(6,9))
#candies2<-filter(candies2,force==1)
#drm(germ_perc~DAY,factor(chillweeks):factor(Taxa), data=candies2,fct = LL.3(), type ="continuous")


#candies3<-filter(realdatshorty, Taxa %in% c("Hesperis matronalis", "Cryptotaenia canadensis"))
#ggplot(candies3,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#candies3<-filter(candies3,chillweeks %in% c(5,8))
#candies3<-filter(candies3,force==0)
#drm(germ_perc~DAY,factor(chillweeks):factor(Taxa), data=candies3,fct = LL.3(c(NA,.8,NA)), type ="continuous")

#candies4<-filter(realdatshorty, Taxa %in% c(,"))
#ggplot(candies4,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#candies4<-filter(candies4,chillweeks %in% c(6))
#drm(germ_perc~DAY,factor(INC):factor(Taxa), data=candies4,fct = LL.3(), type ="continuous")
#dev.off()
#"Hesperis matronalis", "Asclepias syriaca","Oenethera biennis","Eurbia diviricata","Hesperis matronalis","Cryptotaenia canadensis","Eurbia diviricata"

specieslist<-sort(unique(realdat$Taxa))
X<-split(realdat, with(realdat, realdat$Taxa), drop = TRUE)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:20]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)


##Anemone
anemo.cold<-filter(`Anemone virginana`,INC=="L")

data.anemo.cold<-with(anemo.cold,
               list(Y=germ_perc,
                    t=DAY,
                    chill=chillweeks,
                    N=nrow(anemo.cold)
               )
)

mod.anemo.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.anemo.cold, 
                    iter = 9000, warmup=8000 , chain=4) ## 40 divergent transition
summary(mod.anemo.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

anemo.warm<-filter(`Anemone virginana`,INC=="H")

data.anemo.warm<-with(anemo.warm,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(anemo.warm)
                      )
)

mod.anemo.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.anemo.warm, 
                      iter = 9000, warmup=8000 , chain=4) ## 3 divergent transitions
summary(mod.anemo.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]


##Asclepias
asclep.cold<-filter(`Asclepias syriaca`,INC=="L")

data.asclep.cold<-with(asclep.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(asclep.cold)
                      )
)
mod.asclep.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.asclep.cold, 
                     iter = 9000, warmup=8000 , chain=4) ## 0 divergent transition
summary(mod.asclep.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

asclep.warm<-filter(`Asclepias syriaca`,INC=="H")

data.asclep.warm<-with(asclep.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(asclep.warm)
                       )
)
mod.asclep.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.asclep.warm, 
                      iter = 9000, warmup=8000 , chain=4) ## 0 divergent transition, but high bad r hats and low neffs
summary(mod.asclep.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

##Cryptotaenia canadensis
crypto.cold<-filter(`Cryptotaenia canadensis`,INC=="L")

data.crypto.cold<-with(crypto.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(crypto.cold)
                      )
)


mod.crypto.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.crypto.cold, 
                     iter = 9000, warmup=8000 , chain=4) ## 0 divergent transition, but bad rhats and neff
summary(mod.crypto.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

crypto.warm<-filter(`Cryptotaenia canadensis`,INC=="H")

data.crypto.warm<-with(crypto.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(crypto.warm)
                       )
)


mod.crypto.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.crypto.warm, 
                      iter = 9000, warmup=8000 , chain=4) ## 0 divergent transition, good rhats
summary(mod.crypto.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

##Eurybia
eury.cold<-filter(`Eurbia diviricata`,INC=="L")

data.eury.cold<-with(eury.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(eury.cold)
                      )
)

mod.eury.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.eury.cold, 
                     iter = 9000, warmup=8000 , chain=4) ## 224 divergent transition
summary(mod.eury.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

eury.warm<-filter(`Eurbia diviricata`,INC=="H")

data.eury.warm<-with(eury.warm,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(eury.warm)
                     )
)

mod.eury.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.eury.warm, 
                    iter = 9000, warmup=8000 , chain=4) ## 1173 divergent transition
summary(mod.eury.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

##Hesperis
hesper.cold<-filter(`Hesperis matronalis`,INC=="L")

data.hesper.cold<-with(hesper.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(hesper.cold)
                      )
)

mod.hesper.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.hesper.cold, 
                     iter = 9000, warmup=8000 , chain=4) ## 603 divergent transition
summary(mod.hesper.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

hesper.warm<-filter(`Hesperis matronalis`,INC=="H")

data.hesper.warm<-with(hesper.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(hesper.warm)
                       )
)

mod.hesper.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.hesper.warm, 
                      iter = 9000, warmup=8000 , chain=4) ## 839 divergent transition
summary(mod.hesper.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

##Oenethera
oene.cold<-filter(`Oenethera biennis`,INC=="L")

data.oene.cold<-with(oene.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(oene.cold)
                      )
)

mod.oene.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.oene.cold, 
                     iter = 9000, warmup=8000 , chain=4) ##  divergent transition
summary(mod.oene.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

oene.warm<-filter(`Oenethera biennis`,INC=="H")

data.oene.warm<-with(oene.warm,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(oene.warm)
                     )
)

mod.oene.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.oene.warm, 
                    iter = 9000, warmup=8000 , chain=4) ## 0 divergent transition
summary(mod.oene.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

##Polygonum
poly.cold<-filter(`Polygonum virginiatum`,INC=="L")

data.poly.cold<-with(poly.cold,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(poly.cold)
                     )
)

mod.poly.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.poly.cold, 
                    iter = 9000, warmup=8000 , chain=4) ## 3 divergent transition
summary(mod.poly.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

poly.warm<-filter(`Polygonum virginiatum`,INC=="H")

data.poly.warm<-with(poly.warm,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(poly.warm)
                     )
)

mod.poly.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.poly.warm, 
                    iter = 9000, warmup=8000 , chain=4) ##  divergent transition
summary(mod.poly.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

##silene stellata
stella.cold<-filter(`Silene stellata`,INC=="L")

data.stella.cold<-with(stella.cold,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(stella.cold)
                     )
)

mod.stella.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.stella.cold, 
                    iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.stella.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

stella.warm<-filter(`Silene stellata`,INC=="H")

data.stella.warm<-with(stella.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(stella.warm)
                       )
)

mod.stella.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.stella.warm, 
                      iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.stella.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

##silene vulgaris
vulga.cold<-filter(`Silene vulgaris`,INC=="L")

data.vulga.cold<-with(vulga.cold,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(vulga.cold)
                       )
)

mod.vulga.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.vulga.cold, 
                      iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.vulga.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

vulga.warm<-filter(`Silene vulgaris`,INC=="H")

data.vulga.warm<-with(vulga.warm,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(vulga.warm)
                      )
)

mod.vulga.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.vulga.warm, 
                     iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.vulga.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]
save.image(file="realgermers.Rdata")
