rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

time<-seq(0,24,by=3)
treat<-c(0,1)

##parameters
t50.a<-15
t50.b<--6
beta.b<--4

d.a<-6
d.b<-10




for (i in c(1:length(treat))){
for (l in time){
  y<-(d.b*treat[i]+d.a)/(1+((time/(t50.a)^beta.b)))
dfhere<-data.frame(time=time, y=y,treat=treat)
}
}

dfhere$y<-round(dfhere$y)

