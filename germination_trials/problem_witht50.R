library(drc)
rm(list=ls()) 
options(stringsAsFactors = FALSE)

start<-c(0,2,4,6,8,10)
end<-c(2,4,6,8,10,Inf)
germination<-c(1,2,4,2,1,4)
dater<-data.frame(start,end,germination)

mod<-drm(germination~start+end, data = dater,fct = LL.2(), type ="event",robust=TRUE) 
summary(mod)
mod2<-drm(germination~start+end, data = dater,fct = LL.3(), type ="event") 
summary(mod2)
mod3<-drm(germination~start+end, data = dater,fct = LL.3(c(NA,.2,NA)), type ="event") 
summary(mod3)

plot(mod, xlim=c(0,50),ylim=c(0,1),col="red") ### if there is no upper parament it estimate complete germination and t50 to that estimatn
segments(6.434,0,6.434,0.5,col="red")
segments(1,0.5,10,0.5,col="red")


plot(mod2,xlim=c(0,50),ylim=c(0,1),col="blue",add=TRUE)
segments(5.48948,0,5.48948,0.5,col="blue")
segments(1,.85/2,10,.85/2,col="blue")


plot(mod3,xlim=c(0,50),ylim=c(0,1),col="green",add=TRUE)
segments(4.75864,0,4.75864,0.2,col="green")
segments(1,0.1,10,0.1,col="green")

germination2<-c(0,0,1,2,2,7)
dater2<-data.frame(start,end,germination2)


mod0<-drm(germination2~start+end, data = dater2,fct = LL.2(), type ="event") 
summary(mod0)
plot(mod0,xlim=c(0,50),ylim=c(0,1))

mod20<-drm(germination2~start+end, data = dater2,fct = LL.3(), type ="event") 
summary(mod20)

plot(mod0,xlim=c(0,50),ylim=c(0,1),col="red")
segments(10.7224,0,10.7224,0.5,col="red")
segments(5,0.5,20,0.5,col="red")
plot(mod20,xlim=c(0,50),ylim=c(0,1),col="blue", add=TRUE)
segments(3,.52/2,15,.52/2,col="blue")
segments(8.06,0,8.06,.3,col="blue")     
###NOW TRY SOMETHING WITH LEFT CENSORING
start<-c(NA,0,2,4,6,8,10) ### NA just omits the line of data
end<-c(0,2,4,6,8,10,Inf)
germination<-c(2,1,2,4,2,1,4)
lefty<-data.frame(start,end,germination)
mod<-drm(germination~start+end, data = lefty,fct = LL.3(), type ="event",robust=TRUE) 
mod2<-drm(germination~start+end, data = lefty,fct = LL.3(), type ="event")
plot(mod,xlim=c(0,50),ylim=c(0,1),col="red")
plot(mod2,add=TRUE,xlim=c(0,50),ylim=c(0,1),col="blue")
ED(mod,c(5,50))
ED(mod2,c(5,50))
