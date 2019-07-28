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

specieslist<-sort(unique(realdat$Taxa))
X<-split(realdat, with(realdat, realdat$Taxa), drop = TRUE)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:20]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)

data.silene.vulgaris<-with(`Silene vulgaris`,
               list(Y=germ_perc,
                    t=DAY,
                    chill=chillweeks,
                    force=force,
                    N=nrow(`Silene vulgaris`)
               )
)

mod.silene.vulgaris= stan('stan/fakeseedgoodchill_winters.stan', data = data.silene.vulgaris, 
                    iter = 10000, warmup=9000 , chain=4,control=list(adapt_delta=0.95)) ## run model
summary(mod.silene.vulgaris)$summary[c("a_d","bc_d","bf_d","a_beta","bc_beta","bf_beta","a_t50","bc_t50","bf_t50","inter_d","inter_beta","inter_t50","sigma"),]


save.image("real_germ_models") 



