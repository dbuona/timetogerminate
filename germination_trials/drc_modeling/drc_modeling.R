####drc models Dec 17, 2018
##drc treats right centered data is if they will oneday gertminatrion, this is not biologically fair due to dormancy cycling
#but if I dont include these data, the model has no away of accounting for final germination percentage hmm
#also models with with left censoring struggle to converged, although I think its actually having trouble with models with treeatments that have all zero, +++actually its both

##current idea, could ignore left censor and let the model estimate the upper bound, thats probably accuraten for most species except Carex.

rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
library(tidyverse)
library(lubridate)
library("Hmisc")
library(drc)
library(brms)
load("drcmods")
d<-read.csv("germ_data_forDRC.csv",header= TRUE)

d<-filter(d, !Taxa %in% c("Phlox cuspidata","Impatiens capensis", "Carex grisea"))
specieslist<-sort(unique(d$Taxa))
d$DAY<-ifelse(d$DAY==-Inf,NA,d$DAY)###I think this might remove all left censoring
#d$DAY<-ifelse(d$DAY==-Inf,0,d$DAY)# wont converge if left censoring is included

###remove infinity rows with no left of right censoring
nonstartcens<-filter(d,END==0 & germination==0)
noendcense<-filter(d,END==Inf & germination==0)         
d<-anti_join(d, nonstartcens)
d<-anti_join(d,noendcense)

###model without right censoring
#d.nocen<-filter(d,END!=Inf)
#average ungerminated
#dormy<-dplyr::filter(d,END==Inf)
#mean(dormy$germination) #7.1 so on average germination percent was ~60

### make each species its own data framed 
transit<-filter(d, DAY!=25)
transit<-transit %>% group_by(plate_num) %>% summarise(mean=mean(germination))
transit<-filter(transit,mean>0)
germed<-unique(transit$plate_num)

d<-d %>%filter(plate_num %in% c(germed))
X<-split(d, with(d, d$Taxa), drop = TRUE)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:14]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)


mody<-function(x) {
  the_fit <- drm(germination~DAY+END, data=x ,fct = LL.3(), type ="event",upperl=c(NA,1,NA))
  setNames(data.frame(t(coef(the_fit))), c("b", "d","e"))
"predict"(the_fit)
  }

mody.aher<- function(x){
  the_fit <- drm(germination~DAY+END, data=x ,fct = LL.3(), type ="event")
  setNames(data.frame(t(coef(the_fit))), c("b", "d","e"))
}

###Think about changing things in drmc


As<-`Asclepias syriaca` %>% group_by(Taxa,COLD,INC) %>%do(mody(.))
Ed<-`Eurbia diviricata`%>% group_by(Taxa,COLD,INC) %>%do(mody(.))
Sv<-`Silene vulgaris`%>% group_by(Taxa,COLD,INC) %>%do(mody(.))
Av<-`Anemone virginana`%>% group_by(Taxa,COLD,INC) %>%do(mody(.))
Pv<-`Polygonum virginiatum` %>% group_by(Taxa,COLD,INC) %>%do(mody(.))
Td<-`Thalictrum dioicum`%>% group_by(Taxa,COLD,INC) %>%do(mody.aher(.))
OB<-`Oenethera biennis`%>% group_by(Taxa,COLD,INC) %>%do(mody(.))###weird d
CC<-C.canadensis %>% group_by(Taxa,COLD,INC) %>%do(mody(.))
Ss<-S.stellata%>%group_by(Taxa,COLD,INC) %>%do(mody.aher(.))




goober<-rbind(As,Ed,Sv,Av,Pv,Td,CC,Ss)

goober$chill_time<-NA
goober<- within(goober, chill_time[COLD=="O" ]<-0)
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

goober$warm.cent<-goober$warmT-mean(goober$warmT)
goober$chill.cent<-goober$chill_time-mean(goober$chill_time)
  
full.dat<-goober

jpeg("..//rawplot1.jpeg")
ggplot(full.dat,aes(chill_time,e))+geom_point(aes(color=Taxa,shape=INC))+geom_smooth(method="loess",se=FALSE,aes(color=Taxa,linetype=INC))+ylim(0,27)+theme_light()
dev.off()
jpeg("..//rawplot2.jpeg")
ggplot(full.dat,aes(chill_time,e))+geom_point(aes(color=Taxa,shape=INC))+geom_smooth(method="glm",se=FALSE,aes(color=Taxa,linetype=INC))+ylim(0,27)+theme_light()
dev.off()
jpeg("..//rawplot3.jpeg")
ggplot(full.dat,aes(chill_time,d))+geom_point(aes(color=Taxa,shape=INC))+geom_smooth(method="loess",se=FALSE,aes(color=Taxa,linetype=INC))+theme_light()
dev.off()
jpeg("..//rawplot4.jpeg")
ggplot(full.dat,aes(chill_time,d))+geom_point(aes(color=Taxa,shape=INC))+geom_smooth(method="glm",se=FALSE,aes(color=Taxa,linetype=INC))+theme_light()
dev.off()



mods<-brms::brm(e~chill.cent*warm.cent+(chill.cent*warm.cent|Taxa),data=full.dat)
coef(mods)
mods2<-brms::brm(e~chill_time*warmT+(chill_time*warmT|Taxa),data=full.dat)
coef(mods2)
#Hm<- ###not working
lmer()

###not running: both carex, impatiens and phlox 

##Scratch to find out what is bad
C.L<-filter(`Cryptotaenia canadensis`,INC!="H")
C.H<-filter(`Cryptotaenia canadensis`,INC=="H")
C.H<-filter(C.H,COLD!="A")
C.canadensis<-rbind(C.H,C.L)
##########
S.L<-filter(`Silene stellata`,INC!="H")
S.H<-filter(`Silene stellata`,INC=="H")
S.H<-filter(S.H,COLD!="G")
S.stellata<-rbind(S.L,S.H)
################



###can't figure out why hesperis wont converge###########3
H.H<-filter(`Hesperis matronalis`,INC=="H")
H.H<-filter(H.H, COLD=="i") #BDEGHi
drm(germination~DAY+END,COLD, data=H.H,fct = LL.2(), type ="event")


H.L<-filter(`Hesperis matronalis`,INC!="H")
H.L<-filter(H.L, COLD=="B") 
drm(germination~DAY+END,COLD, data=S.H,fct = LL.3(), type ="event",upperl=c(NA,1,0))


save.image("drcmods")






