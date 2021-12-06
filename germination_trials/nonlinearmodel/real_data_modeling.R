rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")

library(rstan)
library(tidyr)
library(drc)
library(dplyr)
library(shinystan)
library(extraDistr)
library(bayesplot)
library(ggplot2)
realdat<-read.csv("input/daily_dat_nointerval.csv")

##clean data
realdat$germ_perc<-NA
realdat$germ_perc<-realdat$germ_num/realdat$tot_seed
realdat$germ_perc<-ifelse(realdat$germ_perc>1,1,realdat$germ_perc)

##make chilling numeric
realdat$chill_time<-NA
realdat<- within(realdat, chill_time[COLD=="0" ]<-0)
realdat<- within(realdat, chill_time[COLD=="A" ]<-14)
realdat<- within(realdat, chill_time[COLD=="B" ]<-28)
realdat<- within(realdat, chill_time[COLD=="C" ]<-35)
realdat<- within(realdat, chill_time[COLD=="D" ]<-42)
realdat<- within(realdat, chill_time[COLD=="E" ]<-49)
realdat<- within(realdat, chill_time[COLD=="f" ]<-56)
realdat<- within(realdat, chill_time[COLD=="G" ]<-63)
realdat<- within(realdat, chill_time[COLD=="H" ]<-77)
realdat<- within(realdat, chill_time[COLD=="i" ]<-91)
realdat$chillweeks<-realdat$chill_time/7 # make chilling weeks instead of days

realdat$force<-NA # make forcing numeric
realdat<- within(realdat, force[INC=="L"]<-0)
realdat<- within(realdat, force[INC=="H"]<-1)

write.csv(realdat,"germ_perc_data.csv", row.names = FALSE)


realdat$DAY<-ifelse(realdat$DAY==0,0.0001,realdat$DAY) #elimiate 0 values for log logistic dist

realdatshorty<- filter(realdat,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea"))
#jpeg("figures/raw_germy.jpeg",res=200,width=1500,height=800)
#ploty<-ggplot(realdatshorty,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) #plot fake data
#ploty+geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#dev.off()

#candies<-filter(realdatshorty, Taxa %in% c("Hesperis matronalis", "Asclepias syriaca"))
#jpeg("figures/woodspecies.jpeg",res=200,width=1500,height=800)
#candies<-filter(candies,chillweeks==5)

#ggplot(candies,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#drm(germ_perc~DAY,factor(INC):factor(Taxa), data=candies,fct = LL.3(), type ="continuous")

#candies2<-filter(realdatshorty, Taxa %in% c("Polygonum virginiatum", "Cryptotaenia canadensis"))
#ggplot(candies2,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#candies2<-filter(candies2,chillweeks %in% c(6,9,13))
#candies2<-filter(candies2,force==1)
#drm(germ_perc~DAY,factor(chillweeks):factor(Taxa), data=candies2,fct = LL.3(), type ="continuous")


candies3<-filter(realdatshorty, Taxa %in% c("Hesperis matronalis","Cryptotaenia canadensis"))
ggplot(candies3,aes(DAY,germ_perc,color=Taxa))+geom_point(aes(color=Taxa,alpha=as.factor(chillweeks)),size=0.2,shape=1)+geom_line(stat = "summary",aes(alpha=as.factor(chillweeks),size=as.factor(chillweeks)))+
  facet_wrap(~force)+theme_minimal(base_size = 6)+scale_size_manual(values=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,2))

ggplot(candies3,aes(DAY,germ_perc,color=Taxa))+geom_point(aes(color=Taxa,alpha=as.factor(chillweeks)),size=0.2,shape=1)+stat_summary(geom="line", fun="mean",aes(alpha=as.factor(chillweeks),size=as.factor(chillweeks)))+
  facet_wrap(~force)+ggthemes::theme_few(base_size = 6)+scale_size_manual(values=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,2))+scale_colour_viridis_d()

#+geom_point(aes(color=Taxa,alpha=chillweeks,group=chillweeks),size=0.2,shape=1)

ggplot(candies3,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+stat_smooth(aes(color=Taxa,fill=Taxa))+facet_grid(chillweeks~INC)+
ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(option="turbo")+scale_fill_viridis_d(option="turbo")+xlab("Day of experiment")+ylab("Germination percentatge")  

ggplot(candies3,aes(DAY,germ_perc))+geom_point(aes(color=Taxa,alpha=as.factor(chillweeks)),size=0.2,shape=1)+stat_smooth(aes(color=Taxa,alpha=as.factor(chillweeks),size=as.factor(chillweeks)),se=FALSE)+facet_wrap(~INC)+
  ggthemes::theme_few()+scale_color_viridis_d(option="turbo")+scale_fill_viridis_d(option="turbo")+scale_size_manual(values=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,2))  

candies3<-filter(candies3,chillweeks %in% c(5,9,13))
candies3<-filter(candies3,force==0)
#candies3<-filter(candies3,chillweeks==13)
mody<-drm(germ_perc~DAY,factor(chillweeks),data=candies3,fct = LL.3(c(NA,NA,NA)), type ="continuous")
unique(`Cryptotaenia canadensis`$chillweeks)
plot(mody,xlim=c(0.1,25))


candies4<-filter(realdatshorty, Taxa %in% c("Hesperis matronalis"))
#ggplot(candies3,aes(DAY,germ_perc,color=Taxa))+geom_point(aes(color=Taxa,alpha=as.factor(chillweeks)),size=0.2,shape=1)+geom_line(stat = "summary",aes(alpha=as.factor(chillweeks),size=as.factor(chillweeks)))+
# facet_wrap(~force)+theme_minimal(base_size = 6)+scale_size_manual(values=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,2))

#ggplot(candies3,aes(DAY,germ_perc,color=Taxa))+geom_point(aes(color=Taxa,alpha=as.factor(chillweeks)),size=0.2,shape=1)+stat_summary(geom="line", fun="mean",aes(alpha=as.factor(chillweeks),size=as.factor(chillweeks)))+
# facet_wrap(~force)+ggthemes::theme_few(base_size = 6)+scale_size_manual(values=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,2))+scale_colour_viridis_d()

#+geom_point(aes(color=Taxa,alpha=chillweeks,group=chillweeks),size=0.2,shape=1)

#ggplot(candies3,aes(DAY,germ_perc))+geom_point(aes(color=Taxa,alpha=c),size=0.2,shape=1)+stat_smooth(aes(color=Taxa,fill=Taxa))+facet_grid(chillweeks~INC)+
#ggthemes::theme_few()+scale_color_viridis_d(option="turbo")+scale_fill_viridis_d(option="turbo")  

#ggplot(candies3,aes(DAY,germ_perc))+geom_point(aes(color=Taxa,alpha=as.factor(chillweeks)),size=0.2,shape=1)+stat_smooth(aes(color=Taxa,alpha=as.factor(chillweeks),size=as.factor(chillweeks)),se=FALSE)+facet_wrap(~INC)+
# ggthemes::theme_few()+scale_color_viridis_d(option="turbo")+scale_fill_viridis_d(option="turbo")+scale_size_manual(values=c(.1,.2,.3,.4,.5,.6,.7,.8,.9,2))  

candies4<-filter(candies4,chillweeks %in% c(5,6,7))
candies4<-filter(candies4,force==0)
mody2<-drm(germ_perc~DAY,factor(chillweeks),data=candies4,fct = LL.3(c(NA,NA,NA)), type ="continuous")
unique(`Cryptotaenia canadensis`$chillweeks)
plot(mody2,col="blue",add=TRUE)



#candies4<-filter(realdatshorty, Taxa %in% c(,"))
#ggplot(candies4,aes(DAY,germ_perc))+geom_point(aes(color=Taxa),size=0.2,shape=1)+facet_grid(force~chillweeks) +geom_line(stat = "summary", fun.y = mean, aes(color=Taxa),size=1.2)+theme_minimal(base_size = 6)
#candies4<-filter(candies4,chillweeks %in% c(6))
#drm(germ_perc~DAY,factor(INC):factor(Taxa), data=candies4,fct = LL.3(), type ="continuous")
#dev.off()
#"Hesperis matronalis", "Asclepias syriaca","Oenethera biennis","Eurbia diviricata","Hesperis matronalis","Cryptotaenia canadensis","Eurbia diviricata"

specieslist<-sort(unique(realdat$Taxa))
X<-split(realdat, with(realdat, realdat$Taxa), drop = TRUE)
Y <- lapply(seq_along(X), function(x) as.data.frame(X[[x]])[, 1:20]) 
names(Y) <-(c(specieslist))
list2env(Y, envir = .GlobalEnv)




##Anemone
anemo.cold<-filter(`Anemone virginana`,INC=="L")
modanemo.cold<-drm(germ_perc~DAY,factor(chillweeks), data=anemo.cold,fct = LL.3(fixed = c(-60, NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
plot(modanemo.cold,xlim=c(10,25))

asclep.cold<-filter(`Asclepias syriaca`,INC=="H")
modasclep.cold<-drm(germ_perc~DAY,factor(chillweeks), data=asclep.cold,fct = LL.3(fixed = c(-60, NA, NA), names = c("b", "gmax", "e50")), type ="continuous")
plot(modasclep.cold,xlim=c(0,25))
goober<-as.data.frame(coef(modasclep.cold))


Sil<-drm(germ_perc~DAY,factor(chillweeks), data=`Silene vulgaris`,fct = LL.3(fixed = c(-60, NA, NA), names = c("b", "gmax", "e50")), type ="continuous")

data.anemo.cold<-with(anemo.cold,
               list(Y=germ_perc,
                    t=DAY,
                    chill=chillweeks,
                    N=nrow(anemo.cold)
               )
)





mod.anemo.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.anemo.cold, 
                    iter = 10000, warmup=9000 , chain=4) ## 40 divergent transition
summary(mod.anemo.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

anemo.warm<-filter(`Anemone virginana`,INC=="H")

data.anemo.warm<-with(anemo.warm,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(anemo.warm)
                      )
)

mod.anemo.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.anemo.warm, 
                      iter = 10000, warmup=9000 , chain=4) ## 3 divergent transitions
summary(mod.anemo.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]


####ppcheck
jpeg("nonlinearmodel/pp_checks/anemone.jpeg",width = 5, height = 6, units = 'in', res = 300)
y.anemo.warm<-data.anemo.warm$Y
y_pred.anemo.warm <- rstan::extract(mod.anemo.warm, "Y_pred")
y.anemo.cool<-data.anemo.cold$Y
y_pred.anemo.cool <- rstan::extract(mod.anemo.cold, "Y_pred")

par(mfrow=c(2,2))
hist(y.anemo.warm, breaks=10, xlab="real data germination response", main="Anemone warm")
hist(y_pred.anemo.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.anemo.cool, breaks=10, xlab="real data germination response", main="Anemone cool")
hist(y_pred.anemo.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")
dev.off()
##Asclepias
asclep.cold<-filter(`Asclepias syriaca`,INC=="L")

data.asclep.cold<-with(asclep.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(asclep.cold)
                      )
)
mod.asclep.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.asclep.cold, 
                     iter = 10000, warmup=9000 , chain=4) ## 2 divergent transition
summary(mod.asclep.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

asclep.warm<-filter(`Asclepias syriaca`,INC=="H")

data.asclep.warm<-with(asclep.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(asclep.warm)
                       )
)
mod.asclep.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.asclep.warm, 
                      iter = 10000, warmup=9000 , chain=4) ## 0 divergent transition, but high bad r hats and low neffs
summary(mod.asclep.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]


####ppcheck
y.asclep.warm<-data.asclep.warm$Y
y_pred.asclep.warm <- rstan::extract(mod.asclep.warm, "Y_pred")
y.asclep.cool<-data.asclep.cold$Y
y_pred.asclep.cool <- rstan::extract(mod.asclep.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/asclepias.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(asclep.warm$germ_perc, breaks=10, xlab="real data germination response", main="Asclepias warm")
hist(y_pred.asclep.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(asclep.cold$germ_perc, breaks=10, xlab="real data germination response", main="Asclepias cool")
hist(y_pred.asclep.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

dev.off()

##Cryptotaenia canadensis
crypto.cold<-filter(`Cryptotaenia canadensis`,INC=="L")

data.crypto.cold<-with(crypto.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(crypto.cold)
                      )
)


mod.crypto.cold= stan('stan/2param.stan', data = data.crypto.cold, 
                     iter = 6000, warmup=5000 , chain=2) 

## 0 divergent transition, but bad rhats and neff
summary(mod.crypto.cold)$summary[c("a_beta","b_beta","a_t50","b_t50","sigma"),]



crypto.warm<-filter(`Cryptotaenia canadensis`,INC=="H")

data.crypto.warm<-with(crypto.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(crypto.warm)
                       )
)


mod.crypto.warm= stan('stan/2param.stan', data = data.crypto.warm, 
                      iter = 9000, warmup=8000 , chain=2) ## 0 divergent transition, good rhats

summary(mod.crypto.warm)$summary[c("a_beta","b_beta","a_t50","b_t50","sigma"),]

y.crypto.warm<-data.crypto.warm$Y
y_pred.crypto.warm <- rstan::extract(mod.crypto.warm, "Y_pred")
y.crypto.cool<-data.crypto.cold$Y
y_pred.crypto.cool <- rstan::extract(mod.crypto.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/crypto.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.crypto.warm, breaks=10, xlab="real data germination response", main="Cryptotaenia warm")
hist(y_pred.crypto.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.crypto.cool, breaks=10, xlab="real data germination response", main="Crypto cool")
hist(y_pred.crypto.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

dev.off()



##Eurybia
eury.cold<-filter(`Eurbia diviricata`,INC=="L")

data.eury.cold<-with(eury.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(eury.cold)
                      )
)

mod.eury.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.eury.cold, 
                     iter = 9000, warmup=8000 , chain=4) ## 224 divergent transition
summary(mod.eury.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

eury.warm<-filter(`Eurbia diviricata`,INC=="H")

data.eury.warm<-with(eury.warm,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(eury.warm)
                     )
)

mod.eury.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.eury.warm, 
                    iter = 9000, warmup=8000 , chain=4) ## 1173 divergent transition
summary(mod.eury.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]



y.eury.warm<-data.eury.warm$Y
y_pred.eury.warm <- rstan::extract(mod.eury.warm, "Y_pred")
y.eury.cool<-data.eury.cold$Y
y_pred.eury.cool <- rstan::extract(mod.eury.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/eury.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.eury.warm, breaks=10, xlab="real data germination response", main="Eurybia warm")
hist(y_pred.eury.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.eury.cool, breaks=10, xlab="real data germination response", main="Eurybia cool")
hist(y_pred.eury.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

dev.off()

##Hesperis
hesper.cold<-filter(`Hesperis matronalis`,INC=="L")

data.hesper.cold<-with(hesper.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(hesper.cold)
                      )
)

mod.hesper.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.hesper.cold, 
                     iter = 9000, warmup=8000 , chain=4) ## 603 divergent transition
summary(mod.hesper.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

hesper.warm<-filter(`Hesperis matronalis`,INC=="H")

data.hesper.warm<-with(hesper.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(hesper.warm)
                       )
)

mod.hesper.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.hesper.warm, 
                      iter = 9000, warmup=8000 , chain=4) ## 839 divergent transition
summary(mod.hesper.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]


y.hesper.warm<-data.hesper.warm$Y
y_pred.hesper.warm <- rstan::extract(mod.hesper.warm, "Y_pred")
y.hesper.cool<-data.hesper.cold$Y
y_pred.hesper.cool <- rstan::extract(mod.hesper.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/hesper.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.hesper.warm, breaks=10, xlab="real data germination response", main="Hesperis warm")
hist(y_pred.hesper.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.hesper.cool, breaks=10, xlab="real data germination response", main="Hesperis cool")
hist(y_pred.hesper.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

dev.off()


##Oenethera
oene.cold<-filter(`Oenethera biennis`,INC=="L")

data.oene.cold<-with(oene.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(oene.cold)
                      )
)

mod.oene.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.oene.cold, 
                     iter = 9000, warmup=8000 , chain=4) ##  divergent transition
summary(mod.oene.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

oene.warm<-filter(`Oenethera biennis`,INC=="H")

data.oene.warm<-with(oene.warm,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(oene.warm)
                     )
)

mod.oene.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.oene.warm, 
                    iter = 9000, warmup=8000 , chain=4) ## 0 divergent transition
summary(mod.oene.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

y.oene.warm<-data.oene.warm$Y
y_pred.oene.warm <- rstan::extract(mod.oene.warm, "Y_pred")
y.oene.cool<-data.oene.cold$Y
y_pred.oene.cool <- rstan::extract(mod.oene.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/oene.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.oene.warm, breaks=10, xlab="real data germination response", main="oenethera warm")
hist(y_pred.oene.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.oene.cool, breaks=10, xlab="real data germination response", main="oenethera cool")
hist(y_pred.oene.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

dev.off()



##Polygonum
poly.cold<-filter(`Polygonum virginiatum`,INC=="L")

data.poly.cold<-with(poly.cold,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(poly.cold)
                     )
)

mod.poly.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.poly.cold, 
                    iter = 9000, warmup=8000 , chain=4) ## 3 divergent transition
summary(mod.poly.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

poly.warm<-filter(`Polygonum virginiatum`,INC=="H")

data.poly.warm<-with(poly.warm,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(poly.warm)
                     )
)

mod.poly.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.poly.warm, 
                    iter = 9000, warmup=8000 , chain=4) ##  divergent transition
summary(mod.poly.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

y.poly.warm<-data.poly.warm$Y
y_pred.poly.warm <- rstan::extract(mod.poly.warm, "Y_pred")
y.poly.cool<-data.poly.cold$Y
y_pred.poly.cool <- rstan::extract(mod.poly.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/poly.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.poly.warm, breaks=10, xlab="real data germination response", main="polygonum warm")
hist(y_pred.poly.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.poly.cool, breaks=10, xlab="real data germination response", main="polygonum cool")
hist(y_pred.poly.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

dev.off()

##silene stellata
stella.cold<-filter(`Silene stellata`,INC=="L")

data.stella.cold<-with(stella.cold,
                     list(Y=germ_perc,
                          t=DAY,
                          chill=chillweeks,
                          N=nrow(stella.cold)
                     )
)

mod.stella.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.stella.cold, 
                    iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.stella.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

stella.warm<-filter(`Silene stellata`,INC=="H")

data.stella.warm<-with(stella.warm,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(stella.warm)
                       )
)

mod.stella.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.stella.warm, 
                      iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.stella.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]


y.stella.warm<-data.stella.warm$Y
y_pred.stella.warm <- rstan::extract(mod.stella.warm, "Y_pred")
y.stella.cool<-data.stella.cold$Y
y_pred.stella.cool <- rstan::extract(mod.stella.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/stella.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.stella.warm, breaks=10, xlab="real data germination response", main="S. stellata warm")
hist(y_pred.stella.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.stella.cool, breaks=10, xlab="real data germination response", main="S.stellata cool")
hist(y_pred.stella.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")
dev.off()

##silene vulgaris
vulga.cold<-filter(`Silene vulgaris`,INC=="L")

data.vulga.cold<-with(vulga.cold,
                       list(Y=germ_perc,
                            t=DAY,
                            chill=chillweeks,
                            N=nrow(vulga.cold)
                       )
)

mod.vulga.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.vulga.cold, 
                      iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.vulga.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

vulga.warm<-filter(`Silene vulgaris`,INC=="H")

data.vulga.warm<-with(vulga.warm,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(vulga.warm)
                      )
)

mod.vulga.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.vulga.warm, 
                     iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.vulga.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

y.vulga.warm<-data.vulga.warm$Y
y_pred.vulga.warm <- rstan::extract(mod.vulga.warm, "Y_pred")
y.vulga.cool<-data.vulga.cold$Y
y_pred.vulga.cool <- rstan::extract(mod.vulga.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/vulga.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.vulga.warm, breaks=10, xlab="real data germination response", main="S. vulgaris warm")
hist(y_pred.vulga.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.vulga.cool, breaks=10, xlab="real data germination response", main="S.vulgaris cool")
hist(y_pred.vulga.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")
dev.off()


# Thalictrum
thali.cold<-filter(`Thalictrum dioicum`,INC=="L")

data.thali.cold<-with(thali.cold,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(thali.cold)
                      )
)

mod.thali.cold= stan('stan/fakeseedgoodchill_alt.stan', data = data.thali.cold, 
                     iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.thali.cold)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

thali.warm<-filter(`Thalictrum dioicum`,INC=="H")

data.thali.warm<-with(thali.warm,
                      list(Y=germ_perc,
                           t=DAY,
                           chill=chillweeks,
                           N=nrow(thali.warm)
                      )
)

mod.thali.warm= stan('stan/fakeseedgoodchill_alt.stan', data = data.thali.warm, 
                     iter = 9000, warmup=8000 , chain=4) ## 
summary(mod.thali.warm)$summary[c("a_d","b_d","a_beta","b_beta","a_t50","b_t50","sigma"),]

####ppcheck
y.thali.warm<-data.thali.warm$Y
y_pred.thali.warm <- rstan::extract(mod.thali.warm, "Y_pred")
y.thali.cool<-data.thali.cold$Y
y_pred.thali.cool <- rstan::extract(mod.thali.cold, "Y_pred")
jpeg("nonlinearmodel/pp_checks/thali.jpeg",width = 5, height = 6, units = 'in', res = 300)
par(mfrow=c(2,2))
hist(y.thali.warm, breaks=10, xlab="real data germination response", main="Thalictrum warm")
hist(y_pred.thali.warm[[1]][1,], breaks=10, xlab="PPC germ perc", main="")

hist(y.thali.cool, breaks=10, xlab="real data germination response", main="Thalictrum cool")
hist(y_pred.thali.cool[[1]][1,], breaks=10, xlab="PPC germ perc", main="")
dev.off()
  
save.image(file="realgermers.Rdata")
