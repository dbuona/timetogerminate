###practice calculating Tb

###these values simulate ascepias in  A
temp<-c(20,20,20,15,15,15)
t50<-c(4,4,6,10,10,10)
t75<-c(6,8,10,11,13,25)
t25<-c(4,4,4,5,6,7)


recip50<-c(1/t50)
recip75<-c(1/t75)
recip25<-c(1/t25)

plot(temp,,xlim=c(0,20),ylim=c(0,0.3))+points(temp,recip25, col="blue")+points(temp,recip75, col="green")+points(temp,recip50,col="red")+abline(lm(recip50~temp),col="red")+abline(lm(recip25~temp),col="blue")+abline(lm(recip75~temp),col="green")

a<-lm(recip25~temp)
xcept25<--(a$coef[1])/a$coef[2]

b<-lm(recip50~temp)
xcept50<--(b$coef[1])/b$coef[2]

c<-lm(recip75~temp)
xcept75<--(c$coef[1])/c$coef[2]

mean(c(xcept25,xcept50,xcept75)) #8.22

####### with more chill
temp<-c(20,20,20,15,15,15)
t25<-c(4,6,6,8,8,6)
t50<-c(6,7,7,9,9,9)
t75<-c(7,9,9,11,11,11)

recip50<-c(1/t50)
recip75<-c(1/t75)
recip25<-c(1/t25)

plot(temp,,xlim=c(-15,20),ylim=c(0,0.3))+points(temp,recip25, col="blue")+points(temp,recip75, col="green")+points(temp,recip50,col="red")+abline(lm(recip50~temp),col="red")+abline(lm(recip25~temp),col="blue")+abline(lm(recip75~temp),col="green")

a<-lm(recip25~temp)
xcept25<--(a$coef[1])/a$coef[2]

b<-lm(recip50~temp)
xcept50<--(b$coef[1])/b$coef[2]

c<-lm(recip75~temp)
xcept75<--(c$coef[1])/c$coef[2]

mean(c(xcept25,xcept50,xcept75)) #1.244

Tb<-c(8.22,1.244)
chilling<-c(14,28)
lm(Tb~chilling)

#(Intercept)     chilling  
#15.1960      -0.4983 


######


