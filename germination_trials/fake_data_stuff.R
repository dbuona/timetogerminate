rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

install.packages("drcSeedGerm")
library(drc)
### right and interval censoring only=goo

###no left or right censoring
start<-c(0,1,2,5,8,10,13)
goo3<-as.data.frame(start)
goo3$end<-c(1,2,5,8,10,13,15)
goo3$count<-c(0,1,1,4,6,0,0) 
#########

##trying to add left censoring, goo2
start<-c(NA,0,1,2,5,8,10,13,15)
goo2<-as.data.frame(start)
goo2$end<-c(0,1,2,5,8,10,13,15,Inf)
goo2$count<-c(10,0,1,1,4,6,0,0,3) 

start<-c(0,1,2,5,8,10,13,15)
goo<-as.data.frame(start)
goo$end<-c(1,2,5,8,10,13,15,Inf)
goo$count<-c(0,1,1,4,6,0,0,3) 

a<-drm(count~start+end,data=goo2,fct=LL.4(c(NA,NA,.9,NA),names=c("Slope","Lower Limit","Upper Limit", "ED50")),type="event")
a2<-drm(count~start+end,data=goo2,fct=LL.4(c(NA,NA,NA,NA),names=c("Slope","Lower Limit","Upper Limit", "ED50")),type="event")
a3<-drm(count~start+end,data=goo3,fct=LL.4(c(NA,NA,NA,NA),names=c("Slope","Lower Limit","Upper Limit", "ED50")),type="event")

b<-drm(count~start+end,data=goo,fct=LL.4(c(NA,NA,.9,NA),names=c("Slope","Lower Limit","Upper Limit", "ED50")),type="event")
b2<-drm(count~start+end,data=goo,fct=LL.4(c(NA,NA,NA,NA),names=c("Slope","Lower Limit","Upper Limit", "ED50")),type="event")

summary(a)
summary(b)
summary(a2)
summary(b2)
plot(a,xlim=c(0,50), ylim=c(-.5,1),col="darkgreen",pch=20)
plot(b,xlim=c(0,50), ylim=c(-.5,1),add=TRUE,col="red",pch=20)

###other optio 