
rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
library(tidyverse)
library(lubridate)
library("Hmisc")
library(drc)

d<-read.csv("germ_data_forDRC.csv",header= TRUE)

d<-filter(d, !Taxa %in% c("Phlox cuspidata","Impatiens capensis", "Carex grisea"))
specieslist<-sort(unique(d$Taxa))

#<- within(d, germ[datasetID=="caffara11b" & response.time==""]
          
        

d<-anti_join(d, noncens)

nonstartsense<-filter(d,end==0 & germination!=0)


### make each species its own data frame
X<-split(d, with(d, d$Taxa  ), drop = TRUE)


Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:14]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)





drm(germination~DAY+END, factor(INC):factor(COLD), data = `Anemone virginana`, fct = LL.2(), type ="event")



Modelit<-function(x,y,z){y<-drm(germination~start+end, factor(INC), data = x, fct = LL.2(), type ="event")
z<-as.data.frame(ED(y,c(50)))}


##nothing converges
Av<-Modelit(`Anemone virginana`)
As<-Modelit(`Asclepias syriaca`)
Cc<-Modelit(`Cryptotaenia canadensis`)
Ed<-Modelit(`Eurbia diviricata`)
Hm<-Modelit(`Hesperis matronalis`)
OB<-Modelit(`Oenethera biennis`)
Pv<-Modelit(`Polygonum virginiatum`)
Ss<-Modelit(`Silene stellata`)
Sv<-Modelit(`Silene vulgaris`)
Td<-Modelit(`Thalictrum dioicum`)

