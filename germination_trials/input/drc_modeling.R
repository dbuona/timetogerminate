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
d$DAY<-ifelse(d$DAY==-Inf,0,d$DAY)

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

### make each species its own data frame

X<-split(d, with(d, d$Taxa), drop = TRUE)


Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:14]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)

####CR
CC.L.G<-dplyr::filter(`Cryptotaenia canadensis`,INC=="L",COLD=="G")
CC.L.G<-filter(CC.L.G,END!=0)
hmm<-drm(germination~DAY+END, data=CC.L.G,fct = LL.3(c(NA,NA,NA)), type ="event",)

?drm()
As<-filter(`Asclepias syriaca`,INC=="H")
As2<-filter(`Asclepias syriaca`,INC=="L")
As.global<-drm(germination~DAY+END,factor(COLD):factor(INC), data=`Asclepias syriaca` ,fct = LL.3(c(NA,NA,NA)), type ="event",upperl =c(NA,1,NA))
As.global2<-drm(germination~DAY+END,factor(COLD), data=As2 ,fct = LL.3(c(NA,NA,NA)), type ="event",upperl =c(NA,1,NA))


summary(As.global)
coef(As.global)
As.global$dataList
par(mfrow=c(1,1))
plot(As.global, xlab = "Time", ylab = "Proportion germinated", 
     xlim=c(0, 25), ylim=c(0, 1), log="", lwd=1, cex=1.2, col="black",legendPos ="" )



legend("bottomright", inset=c(-0.2,0), legend=c("A","B"), pch=c(1,3), title="Group")

plot(As.global2,add=FALSE, xlab = "Time", ylab = "Proportion germinated", 
      xlim=c(0, 25), ylim=c(0, 1), log="", lwd=1, cex=1.2,col="blue")
legend("topright", inset=c(-0.2,0), legend=c("A","B"), pch=c(1,3), title="Group")



ED(As.global,c(50),"delta")

Cry<-filter(`Cryptotaenia canadensis`,END!="0")
Cry<-filter(Cry,COLD!="O")
Cry<-filter(Cry,COLD!="A")

Cc.global<-drm(germination~DAY+END,factor(INC):factor(COLD), data=Cry,fct = LL.3(c(NA,.90,NA)), type ="event")
ED(Cc.global,c(50),"delta")

Av.global<-drm(germination~DAY+END,factor(INC):factor(COLD), data=`Anemone virginana`,fct = LL.3(c(NA,NA,NA)), type ="event")
summary(Av.global)
ED(Av.global,c(50),"delta")

Sil.v<-filter(`Silene vulgaris`,END!=0)
Sv.global<-drm(germination~DAY+END,factor(INC):factor(COLD), data=Sil.v,fct = LL.3(c(NA,.90,NA)), type ="event")
ED(Sv.global,c(50),"delta")


Pol<-filter(`Polygonum virginiatum`,END!=0)
Pol<-filter(Pol,COLD != "O")
Pol<-filter(Pol,COLD != "A")
Pv.global<-drm(germination~DAY+END,factor(INC):factor(COLD), data=Pol,fct = LL.3(c(NA,.90,NA)), type ="event")
ED(Pv.global,c(50),"delta")

AvO.mod<-drm(germination~DAY+END,factor(INC):factor(COLD), data=`Anemone virginana`,fct = LL.3(c(NA,.99,NA)), type ="event")
ED(AvO.mod,c(50),"delta")

##oenethera
 ###remove the one row that is left censored
On.global<-drm(germination~DAY+END, data=`Oenethera biennis`,fct = LL.3(), type ="event")
summary(On.global)

OnO<-dplyr::filter(`Oenethera biennis`,COLD=="i")

OnO.mod<-drm(germination~DAY+END,factor(INC):factor(COLD), data=`Oenethera biennis`,fct = LL.3(c(NA,.90,NA)), type ="event")
ED(OnO.mod,c(50),"delta")

##crypto
`Cryptotaenia canadensis`<-filter(`Cryptotaenia canadensis`,DAY!=0) ###remove the one row that is left censored
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

