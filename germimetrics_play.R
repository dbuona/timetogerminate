rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")

library("brms")

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

####calculate GI germination index
d<-d %>% group_by(plate_num) %>%  mutate(weight = row_number())
d<-d %>% group_by(plate_num)  %>%  mutate(weight = rev(weight))


d$weighting<-d$weight*d$germ.daily
GIs<-d %>% group_by(Taxa,warmT,chill_time,COLD,INC) %>% summarise(GI=sum(weighting))                                          

mod<-brm(GI~warmT+chill_time+(warmT+chill_time|Taxa),data=GIs)
summary(mod)
coef(mod)



unique(GIs$Taxa)
GI2s<-filter(GIs,Taxa %in% c("Anemone virginana"     ,  "Asclepias syriaca" ,"Cryptotaenia canadensis","Eurbia diviricata","Hesperis matronalis","Oenethera biennis","Polygonum virginiatum" ,"Silene vulgaris","Silene stellata"))

GIwet<-filter(GI2s,Taxa %in% c("Cryptotaenia canadensis","Hesperis matronalis","Polygonum virginiatum","Eurbia diviricata","Anemone virginana"))
GIdry<-filter(GI2s,Taxa %in% c("Silene vulgaris","Asclepias syriaca"))

ggplot(GIwet,aes(reorder(Taxa,-GI),GI))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+facet_grid(INC~COLD)+theme(axis.text.x = element_blank())
ggplot(GIdry,aes(reorder(Taxa,-GI),GI))+geom_bar(stat = "identity",position="dodge",aes(fill=Taxa))+facet_grid(INC~COLD)+theme(axis.text.x = element_blank())


####germination percentage
dt<-read.csv("cumulative_data.csv")
final<-filter(dt,DAY==25)
final$percent<-final$germ_num/final$tot_seed
