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

d<-read.csv("germ_data_forDRC.csv",header= TRUE)

unique(d$COLD)
d<-filter(d, !Taxa %in% c("Phlox cuspidata","Impatiens capensis", "Carex grisea"))
specieslist<-sort(unique(d$Taxa))

###remove infinity rows with no left of right censoring
nonstartcens<-filter(d,END==0 & germination==0)
noendcense<-filter(d,END==Inf & germination==0)         
d<-anti_join(d, nonstartcens)
d<-anti_join(d,noendcense)

####model without right censoring
#d.nocen<-filter(d,DAY!=-Inf)
#d.nocen<-filter(d,END!=Inf)
#average ungerminated
#dormy<-dplyr::filter(d,END==Inf)
#mean(dormy$germination) #7.1 so on average germination percent was ~60

### make each species its own data frame
X<-split(d, with(d, d$Taxa), drop = TRUE)


Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:14]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)

####Asclepias first
As.global<-drm(germination~DAY+END, data=`Asclepias syriaca`,fct = LL.3(c(NA,NA,NA)), type ="event")
summary(As.global)

As.global<-drm(germination~DAY+END,factor(INC):factor(COLD), data=`Asclepias syriaca`,fct = LL.3(c(NA,.95,NA)), type ="event")
summary(As.global)
ED(As.global,c(50),"delta")

As.global.est<-drm(germination~DAY+END,factor(INC):factor(COLD), data=`Asclepias syriaca`,fct = LL.3(c(NA,NA,NA)), type ="event")
summary(As.global)
ED(As.global,c(50),"delta")

AsO.mod<-drm(germination~DAY+END,INC, data=AsO,fct = LL.3(c(NA,NA,NA)), type ="event")
ED(AsO.mod,c(50),"delta")

summary(AsO.mod)

##anemone
Av.global<-drm(germination~DAY+END, data=`Anemone virginana`,fct = LL.3(), type ="event")
summary(Av.global)

AvO<-dplyr::filter(`Anemone virginana`,COLD=="i")

AvO.mod<-drm(germination~DAY+END,INC, data=AvO,fct = LL.3(c(NA,.99,NA)), type ="event")
ED(AvO.mod,c(50),"delta")

##oenethera
`Oenethera biennis`<-filter(`Oenethera biennis`,DAY!=-Inf) ###remove the one row that is left censored
On.global<-drm(germination~DAY+END, data=`Oenethera biennis`,fct = LL.3(), type ="event")
summary(On.global)

OnO<-dplyr::filter(`Oenethera biennis`,COLD=="i")

OnO.mod<-drm(germination~DAY+END,factor(INC):factor(COLD), data=`Oenethera biennis`,fct = LL.3(c(NA,.90,NA)), type ="event")
ED(OnO.mod,c(50),"delta")

##crypto
`Cryptotaenia canadensis`<-filter(`Cryptotaenia canadensis`,DAY!=-Inf) ###remove the one row that is left censored
CR.global<-drm(germination~DAY+END, INC,data=`Cryptotaenia canadensis`,fct = LL.3(c(NA,NA,NA)), type ="event")
summary(CR.global)
ED(CR.global,c(50),"delta")
plot(CR.global)
segments(4.5769,0,4.5769,0.2,col="red")
segments(5.385048,0,5.385048,0.38,col="blue")

CR.global.foxed<-drm(germination~DAY+END, INC,data=`Cryptotaenia canadensis`,fct = LL.3(c(NA,1,NA)), type ="event")
summary(CR.global.foxed)
plot(CR.global.foxed)
segments(4.59,0,4.59,0.2,col="red")
segments(5.385048,0,5.385048,0.38,col="blue")




CR.H<-dplyr::filter(`Cryptotaenia canadensis`,INC=="H")

CR.H.mod<-drm(germination~DAY+END,COLD, data=CR.H,fct = LL.3(c(NA,NA,NA)), type ="event")
ED(CR.H.mod,c(50),"delta")
summary(CR.H.mod)
CR.L<-dplyr::filter(`Cryptotaenia canadensis`,INC=="L")
CR.L.mod<-drm(germination~DAY+END,COLD, data=CR.L,fct = LL.3(c(NA,.73,NA)), type ="event")
ED(CR.L.mod,c(50),"delta")



####old attempt to model each treatment together
Modelit<-function(x,y,z){y<-drm(germination~DAY+END, factor(INC):factor(COLD), data = x, fct = LL.2(), type ="event")
z<-as.data.frame(ED(y,c(50)))}


##nothing converges
Av<-Modelit(`Anemone virginana`)
As<-Modelit(`Asclepias syriaca`)


Cc<-Modelit(`Cryptotaenia canadensis`)##no converge
Ed<-Modelit(`Eurbia diviricata`)###
Hm<-Modelit(`Hesperis matronalis`)##
OB<-Modelit(`Oenethera biennis`)## convergese without censoring
Pv<-Modelit(`Polygonum virginiatum`)###
Ss<-Modelit(`Silene stellata`) ##convergest with out censoring, but estimates seem pretty small
Sv<-Modelit(`Silene vulgaris`) ##converges without censoring
Td<-Modelit(`Thalictrum dioicum`)###

