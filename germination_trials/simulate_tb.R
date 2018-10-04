###practice calculating Tb

###these values simulate ascepias in  A
temp<-c(20,20,20,15,15,15)
t25<-c(4,4,5,8,6,7)
t50<-c(4,4,6,10,10,11)
t75<-c(6,8,10,11,13,25)



recip50<-c(1/t50)
recip75<-c(1/t75)
recip25<-c(1/t25)

p.1<-plot(temp,,xlim=c(0,20),ylim=c(0,0.3))+points(temp,recip25, col="blue")+points(temp,recip75, col="green")+points(temp,recip50,col="red")+abline(lm(recip50~temp),col="red")+abline(lm(recip25~temp),col="blue")+abline(lm(recip75~temp),col="green")

a<-lm(recip25~temp)
xcept25<--(a$coef[1])/a$coef[2]

b<-lm(recip50~temp)
xcept50<--(b$coef[1])/b$coef[2]

c<-lm(recip75~temp)
xcept75<--(c$coef[1])/c$coef[2]

mean(c(xcept25,xcept50,xcept75)) #8.22

####### with more chill
temp<-c(20,20,20,15,15,15)
t25<-c(4,5,6,8,8,7)
t50<-c(6,7,6,9,9,11)
t75<-c(7,9,9,11,11,15)

recip501<-c(1/t50)
recip751<-c(1/t75)
recip251<-c(1/t25)

p.2<-plot(temp,recip50,xlim=c(0,20),ylim=c(0,0.28))+points(temp,recip25, col="blue")+points(temp,recip75, col="green")+points(temp,recip50,col="red")+abline(lm(recip50~temp),col="red")+abline(lm(recip25~temp),col="blue")+abline(lm(recip75~temp),col="green")+points(temp,recip251, col="darkblue")+points(temp,recip751, col="darkgreen")+points(temp,recip501,col="brown")+abline(lm(recip501~temp),col="red",lty="dashed")+abline(lm(recip251~temp),col="blue",lty="dashed")+abline(lm(recip751~temp),col="green",lty="dashed")

a<-lm(recip251~temp)
xcept25<--(a$coef[1])/a$coef[2]

b<-lm(recip501~temp)
xcept50<--(b$coef[1])/b$coef[2]

c<-lm(recip751~temp)
xcept75<--(c$coef[1])/c$coef[2]

mean(c(xcept25,xcept50,xcept75)) #4.222

Tb<-c(8.22, 5.32)
chilling<-c(14,28)
lm(Tb~chilling)




######


