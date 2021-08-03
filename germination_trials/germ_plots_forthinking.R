## Started 2 August 20221 ##
## By Lizzie (for now) ##

# Copied some code from germ_indices/fgpcode.R 
# The column "germ_perc_daily" is the daily germination percentages and "germ_perc"has the cumulative percentages so you could subset to the last day of the trial (column "DAY" == 25, not to be confused with column "day") to get the final germination percentages which I did for a subset of species in line 52.

rm(list=ls()) 
options(stringsAsFactors = FALSE)

if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/others/dan/timetogerminate/germination_trials") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("~/Documents/git/timetogerminate/germination_trials/")

library(tidyr)
library(plyr)
library(dplyr)
library(ggplot2)
library(grid)

realdat<-read.csv("input/daily_dat_nointerval.csv")

##clean data
realdat$germ_perc <- NA
realdat$germ_perc <- realdat$germ_num/realdat$tot_seed
realdat$germ_perc_daily <- realdat$germ.daily/realdat$tot_seed
realdat$germ_perc <- ifelse(realdat$germ_perc>1,1,realdat$germ_perc)

##make chilling numeric
realdat$chill_time<-NA
realdat <- within(realdat, chill_time[COLD=="0" ]<-0)
realdat <- within(realdat, chill_time[COLD=="A" ]<-14)
realdat <- within(realdat, chill_time[COLD=="B" ]<-28)
realdat <- within(realdat, chill_time[COLD=="C" ]<-35)
realdat <- within(realdat, chill_time[COLD=="D" ]<-42)
realdat <- within(realdat, chill_time[COLD=="E" ]<-49)
realdat <- within(realdat, chill_time[COLD=="f" ]<-56)
realdat <- within(realdat, chill_time[COLD=="G" ]<-63)
realdat <- within(realdat, chill_time[COLD=="H" ]<-77)
realdat <- within(realdat, chill_time[COLD=="i" ]<-91)
realdat$chillweeks <-realdat$chill_time/7 # make chilling weeks instead of days

realdat$force <-NA # make forcing numeric
realdat <- within(realdat, force[INC=="L"]<-0)
realdat <- within(realdat, force[INC=="H"]<-5)

realdatfin <- subset(realdat, DAY==25) # get the final germination fraction

# summarize ...
dsum <-
      ddply(realdat, c("DAY", "Taxa", "chillweeks", "force"), summarise,
      mean = mean(germ_perc_daily),
      sd = sd(germ_perc_daily),
      sem = sd(germ_perc_daily)/sqrt(length(germ_perc_daily)))

dsumfin <-
      ddply(realdatfin, c("Taxa", "chillweeks", "force"), summarise,
      mean = mean(germ_perc),
      sd = sd(germ_perc),
      sem = sd(germ_perc)/sqrt(length(germ_perc)))


# plotting
ggplot(realdat, aes(x=germ_perc_daily, fill=Taxa, color=Taxa)) +
    geom_density(alpha=0.2) +
    facet_wrap(.~chillweeks)

ggplot(realdat, aes(x=DAY, y=germ_perc_daily, group=Taxa, color=Taxa)) +
    geom_line() +
    facet_wrap(.~force*chillweeks)

ggplot(dsum, aes(x=DAY, y=mean, group=Taxa, color=Taxa)) +
    geom_line() +
    facet_wrap(.~force*chillweeks) +
    xlab("Days after forcing") + ylab("Daily percent germinated (mean)")

dsumforce0 <- subset(dsum, force==0)
dsumforce5 <- subset(dsum, force==5)

ggplot(dsumforce0, aes(x=DAY, y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_point() +
    geom_ribbon(aes(ymin=(dsumforce0$mean-dsumforce0$sem), ymax=(dsumforce0$mean+dsumforce0$sem)), alpha=0.1) + 
    facet_wrap(.~Taxa*chillweeks)  +
    xlab("Days after forcing (at 0)") + ylab("Daily percent germinated (mean +/- SE)")

ggplot(dsumforce5, aes(x=DAY, y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_point() +
    geom_ribbon(aes(ymin=(dsumforce5$mean-dsumforce5$sem), ymax=(dsumforce5$mean+dsumforce5$sem)), alpha=0.1) + 
    facet_wrap(.~Taxa*chillweeks)  +
    xlab("Days after forcing (at 5)") + ylab("Daily percent germinated (mean +/- SE)")

ggplot(dsumfin, aes(x=as.numeric(chillweeks), y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_line() +
    geom_ribbon(aes(ymin=(dsumfin$mean-dsumfin$sem), ymax=(dsumfin$mean+dsumfin$sem)), alpha=0.1) + 
    facet_wrap(.~Taxa*force) +
    xlab("Weeks of chilling)") + ylab("Final total percent germinated (mean +/- SE)")
