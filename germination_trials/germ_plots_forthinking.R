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

maxdaily <- realdat %>%
  select(Taxa, chillweeks, force, DAY, germ_perc_daily)%>%
  group_by(Taxa,chillweeks,force)%>% filter(germ_perc_daily== max(germ_perc_daily)) %>%
  group_by(Taxa,chillweeks,force,germ_perc_daily) %>% summarise(mean_day=mean(DAY),sd_day=sd(DAY))
maxdaily$maxgermin20 <- maxdaily$germ_perc_daily*20

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

totalgermbyplate <-
      ddply(realdat, c("Taxa", "chillweeks", "force", "plate_num"), summarise,
      sum = sum(germ_perc_daily))

totalgerm <-
      ddply(totalgermbyplate, c("Taxa", "chillweeks", "force"), summarise,
      mean = mean(sum),
      sd = sd(sum),
      sem = sd(sum)/sqrt(length(sum)))


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

## BIG plots (next two)
# free-form, but you can read species names
ggplot(dsumforce0, aes(x=DAY, y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_point() +
    geom_ribbon(aes(ymin=(dsumforce0$mean-dsumforce0$sem), ymax=(dsumforce0$mean+dsumforce0$sem)), alpha=0.1) + 
    facet_wrap(.~Taxa*chillweeks)  +
    xlab("Days after forcing (at 0)") + ylab("Daily percent germinated (mean +/- SE)")

ggplot(dsumforce5, aes(x=DAY, y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_point() +
    geom_ribbon(aes(ymin=(dsumforce5$mean-dsumforce5$sem), ymax=(dsumforce5$mean+dsumforce5$sem)), alpha=0.1) + 
    facet_wrap(.~Taxa*chillweeks) +
    xlab("Days after forcing (at 5)") + ylab("Daily percent germinated (mean +/- SE)")

# rows are each species
ggplot(dsumforce5, aes(x=DAY, y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_point() +
    geom_ribbon(aes(ymin=(dsumforce5$mean-dsumforce5$sem), ymax=(dsumforce5$mean+dsumforce5$sem)), alpha=0.1) + 
    facet_grid(Taxa~chillweeks)  +
    xlab("Days after forcing (at 5)") + ylab("Daily percent germinated (mean +/- SE)")

ggplot(dsumforce0, aes(x=DAY, y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_point() +
    geom_ribbon(aes(ymin=(dsumforce0$mean-dsumforce0$sem), ymax=(dsumforce0$mean+dsumforce0$sem)), alpha=0.1) + 
    facet_grid(Taxa~chillweeks)  +
    xlab("Days after forcing (at 5)") + ylab("Daily percent germinated (mean +/- SE)")

ggplot(dsum, aes(x=DAY, y=mean, group=as.factor(force), color=as.factor(force))) +
    geom_line() +
    facet_grid(Taxa~chillweeks)  +
    xlab("Days after forcing") + ylab("Daily percent germinated (mean +/- SE)")
## End BIG plots

ggplot(dsumfin, aes(x=as.numeric(chillweeks), y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_line() +
    geom_ribbon(aes(ymin=(dsumfin$mean-dsumfin$sem), ymax=(dsumfin$mean+dsumfin$sem)), alpha=0.1) + 
    facet_wrap(.~Taxa*force) +
    xlab("Weeks of chilling)") + ylab("Final total percent germinated (mean +/- SE)")

ggplot(totalgerm, aes(x=as.numeric(chillweeks), y=mean, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_point() +
    geom_ribbon(aes(ymin=(totalgerm$mean-totalgerm$sem), ymax=(totalgerm$mean+totalgerm$sem)), alpha=0.1) + 
    facet_wrap(.~Taxa*force) +
    xlab("Weeks of chilling") + ylab("Total germinated (mean +/- SE)")


# max germination, first plot has just day of max germ, second one also plots that value
ggplot(maxdaily, aes(x=as.numeric(chillweeks), y=mean_day, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_line() +
    facet_wrap(.~Taxa*force) +
    xlab("Weeks of chilling") + ylab("Day of max germination)")

ggplot(maxdaily, aes(x=as.numeric(chillweeks), y=mean_day, group=Taxa, fill=Taxa, color=Taxa)) +
    geom_line() +
    geom_line(aes(y=maxgermin20), linetype = "dashed") +
    facet_wrap(.~Taxa*force) +
    xlab("Weeks of chilling)") +
    scale_y_continuous(name = "Day of max germination (solid line)", sec.axis = sec_axis(~.*1, name="Max germination (20X, dashed line)"))
