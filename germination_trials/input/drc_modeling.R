####drc models Dec 17, 2018
##drc treats right centered data is if they will oneday gertminatrion, this is not biologically fair due to dormancy cycling
#but if I dont include these data, the model has no away of accounting for final germination percentage hmm
#also models with with left censoring struggle to converged, although I think its actually having trouble with models with treeatments that have all zero, +++actually its both


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

 
goop<-dplyr::filter(`Cryptotaenia canadensis`, COLD=="G" & INC=="L")


goo<-drm(germination~DAY+END,, data = goop,fct = LL.3(), type ="event") ##problem if this model estimates 
goo50<-drm(germination~DAY+END,factor(INC):factor(COLD), data = goop,fct = LL.3(c(NA,.5,NA)), type ="event") 
goo100<-drm(germination~DAY+END, data = goop,fct = LL.2(), type ="event") ##problem if this model estimates 

ED(goo,c(50),"delta")
summary(goo)
ED(goo50,c(50),"delta")
ED(goo100,c(50),"delta")

summary(goo)
plot(goo100,xlim=c(0,100), ylim=c(0,1),col="darkgreen",pch=20)
axis(1, at=seq(0, 25, by=5), labels = TRUE)
plot(goo,add=TRUE,xlim=c(0,100), ylim=c(0,1),col="red",pch=20)
plot(goo50,add=TRUE,xlim=c(0,100), ylim=c(0,1),col="blue",pch=20)



ED(goo,(50),"delta")
drm(germination~DAY+END,factor(INC):factor(COLD), data = `Av.f`, fct = LL.4(), type ="event", start = c(coef(goo)[1:18],rep(-Inf,18), coef(goo)[19:36]))

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

