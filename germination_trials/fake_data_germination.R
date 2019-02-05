###fake data for germination mdodels
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

library(drc)
library(tidyverse)
library(MCMCglmm)
##use lognoral to simulate the data
set.seed(613613)


###This make data folowing the log logistic function
germ<-function(t,d,b,t50){
y<- d/(1+((t/t50)^b))
return(data.frame(time=t, y=y)) 
  
}

germ(seq(0,24,by=3),rtnorm(1,.6,0.1,lower=0.4,upper=1),rnorm(1,-5,0.1),rnorm(1,15,1))


##3 petridishes of the same treatment, there is probably a loop or apply function for this

A<-germ(seq(0,24,by=3),rtnorm(1,.6,0.1,lower=0.4,upper=1),rnorm(1,-10,0.1),rnorm(1,15,1))
A$dish<-"A"
B<-germ(seq(0,24,by=3),rtnorm(1,.6,0.1,lower=0.4,upper=1),rnorm(1,-10,0.1),rnorm(1,15,1))
B$dish<-"B"
C<-germ(seq(0,24,by=3),rtnorm(1,0.6,0.1,lower=0.4,upper=1),rnorm(1,-10,0.1),rnorm(1,15,1))
C$dish<-"C"

#make your data
d<-rbind(A,B,C)
  
mod<-drm(y~time,fct=LL.3(),data=d,type="continuous")
lines(d$time,predict(mod),lty=2,col="red",lwd=3)
summary(mod)
plot(mod,ylim=c(0,1),xlim=c(0,24),log="",pch=16,type="all")

###now now to change the chilling

?pch
