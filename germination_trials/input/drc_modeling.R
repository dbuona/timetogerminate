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
noendcense<-filter(d,END==Inf,germination==0)         
d<-anti_join(d, nonstartcens)
d<-anti_join(d,noendcense)

####model without censoring
d.nocen<-filter(d,END!=Inf)
d.nocen<-filter(d.nocen,DAY!=-Inf)

### make each species its own data frame
X<-split(d, with(d, d$Taxa), drop = TRUE)


Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:14]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)

#=============================drc example
germLL.2 <- drm(germinated ~ start + end, species:factor(temp), 
 data = germination[c(1:23, 25:61, 63:192), ], fct = LL.2(), type = "event")
plot(germLL.2, ylim=c(0, 1.5), legendPos=c(2.5,1.5))  # plotting the fitted curves and the data
 summary(germLL.2)
 
 
 germLL.3 <- drm(germinated~start+end, species:factor(temp), 
                 data = germination[c(1:23, 25:61, 63:192), ], fct = LL.3(), type = "event",
                  start = c(coef(germLL.2)[1:13], rep(0.7,13), coef(germLL.2)[14:26]), 
                  upper = c(rep(Inf, 13), rep(1, 13), rep(Inf, 13)))
#+++++++++++++++++++++++END Drc

Poly.L<-dplyr::filter(`Cryptotaenia canadensis`,INC=="L")
PolyL.G.<-dplyr::filter(Poly.L,COLD=="H")

goo<-drm(germination~DAY+END, data = PolyL.G., fct = LL.2(), type ="event")
summary(goo)
plot(goo,xlim=c(1,25), ylim=c(0,1.2),legendPos=c(2.5,1))
ED(goo,(90))
drm(germination~DAY+END,factor(INC):factor(COLD), data = `Av.f`, fct = LL.3(), type ="event", start = c(coef(goo)[1:18],rep(-Inf,18), coef(goo)[19:36]))

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

