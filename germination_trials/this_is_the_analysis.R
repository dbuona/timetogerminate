
### prediction: Cliamte change will change the competitive dynamics  between species through
##temprotal effects

rm(list=ls()) 
options(stringsAsFactors = FALSE)


library(ggplot2)
library(dplyr)
library("arm")
library(brms)
library(ggstance)

if(length(grep("Lizzie", getwd())>0)) {
setwd("~/Documents/git/projects/others/dan/timetogerminate/germination_trials") } else if (length(grep("boomer", getwd()))>0) {setwd("boom/boom")
}  else setwd("~/Documents/git/timetogerminate/germination_trials/")


load("goingforitgerm")
##likelihood of germination== propogule pressure
d<-read.csv("survival_analysis/surival_dat_nointerval.csv")
unique(d$Taxa)


d$chillweeks<-d$chill_time/7 # make chilling weeks instead of days

d$force<-NA # make forcing numeric
d<- within(d, force[INC=="L"]<-0)
d<- within(d, force[INC=="H"]<-1)

###prepfor survival analysis
d$censored<-ifelse(d$germinated==1,0,1)
d$DAY<-ifelse(d$DAY==0,0.00001,d$DAY)
d$censored<-ifelse(d$DAY==0.00001 & d$germinated==1,-1,d$censored)


Pv<- filter(d,Taxa %in% c("Polygonum virginiatum"))
Cc<- filter(d,Taxa %in% c("Cryptotaenia canadensis"))
Ed<-filter(d,Taxa %in% c("Eurbia diviricata"))
Td<-filter(d,Taxa=="Thalictrum dioicum")
Av<-filter(d,Taxa=="Anemone virginana")
Ss<-filter(d,Taxa=="Silene stellata")
mod1<-brm(germinated~chill_time*force,data=Pv, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)
mod2<-brm(germinated~chill_time*force,data=Cc, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)
mod3<-brm(germinated~chill_time*force,data=Ed, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)
mod4<-brm(germinated~chill_time*force,data=Td, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)
mod5<-brm(germinated~chill_time*force,data=Av, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)
mod6<-brm(germinated~chill_time*force,data=Ss, family =bernoulli(link="logit"),chains=4, control=list(adapt_delta=0.99),iter=4000,warmup = 3000)

new.data<-data.frame(chill_time=rep(c(79.2, 47.73),2), force=rep(c(0,1),each=2))


goot<-fitted(mod1,newdata = new.data,summary = TRUE)
goot2<-fitted(mod2,newdata = new.data,summary = TRUE)
goot3<-fitted(mod3,newdata = new.data,summary = TRUE)
goot4<-fitted(mod4,newdata = new.data,summary = TRUE)
goot5<-fitted(mod5,newdata = new.data,summary = TRUE)
goot6<-fitted(mod6,newdata = new.data,summary = TRUE)

prediction<-cbind(new.data,goot) 
prediction$Taxa<-"P.virginiatum"

prediction2<-cbind(new.data,goot2) 
prediction2$Taxa<-"C.canadensis"

prediction3<-cbind(new.data,goot3) 
prediction3$Taxa<-"E.divaricata"

prediction4<-cbind(new.data,goot4) 
prediction4$Taxa<-"T.dioicum"

prediction5<-cbind(new.data,goot5) 
prediction5$Taxa<-"A.virginiana"

prediction6<-cbind(new.data,goot6) 
prediction6$Taxa<-"S.stellata"


predictionall<-rbind(prediction, prediction2,prediction3,prediction4)

png("figures/germpercs.png",width = 5,height = 5,units = "in",res = 300)
ggplot(predictionall,aes(as.factor(force),Estimate))+
  geom_bar(stat="identity",position="dodge",aes(fill=Taxa))+
  ggthemes::theme_base()+facet_wrap(~chill_time)+scale_fill_brewer( type = "qual", palette = "Dark2", direction = 1)

dev.off()
new.data2<-


ggplot(predictionall,aes(scenario,Estimate))+
  stat_summary(aes(color=Taxa))+
  ggthemes::theme_base()
###survivial analysis


fit.wei.Pv<- brm(DAY | cens(censored)~chill_time+force, data=Pv, family =   weibull(),inits=0 ,iter=4000,warmup = 3000, chains=4) 
fit.wei.Cc<- brm(DAY | cens(censored)~chill_time+force, data=Cc, family =   weibull(),inits=0 ,iter=4000,warmup = 3000, chains=4) 
fit.wei.Ed<- brm(DAY | cens(censored)~chill_time+force, data=Ed, family =   weibull(),inits=0 ,iter=4000,warmup = 3000, chains=4) 
dev.off()
daty.pv<-fitted(fit.wei.Pv,probs =c(0.25,.75),newdata=new.data)
daty.pv<-cbind(daty.pv,new.data)
daty.pv$Taxa<-"P.virginiatum"

daty.cc<-fitted(fit.wei.Cc,probs =c(0.25,.75),newdata=new.data)
daty.cc<-cbind(daty.cc,new.data)
daty.cc$Taxa<-"C.canadensis"


daty.ed<-fitted(fit.wei.Ed,probs =c(0.25,.75),newdata=new.data)
daty.ed<-cbind(daty.ed,new.data)
daty.ed$Taxa<-"E.divaricata"

daty<-rbind(daty.cc,daty.ed,daty.pv)
pd<-position_dodge(width=0.4)


png("figures/germtime.png",width =7,height = 5,units = "in",res = 300)
ggplot(daty,aes(x=Estimate,y=0, color=Taxa,label=Taxa))+
  #xlim(-20,120)+
  ylim(-5,5)+#geom_point(position=pd,aes(color=Taxa),size=4)+
  geom_vline(aes(xintercept=Estimate, color= Taxa),size=1.4)+
  geom_errorbarh(aes(xmin=Q25,xmax=Q75,color=Taxa),height=1,alpha=0.8)+
  facet_grid(as.factor(force)~as.factor(chill_time))+theme_linedraw()+
  scale_colour_brewer( type = "qual", palette = "Dark2", direction = 1)+
  scale_fill_brewer( type = "qual", palette = "Dark2", direction = 1)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank())+
  labs(x="Model Estimated Days to 50% Germination")+scale_x_continuous(breaks=c(0,5,10,20,30,40,50,75,100,125))
dev.off()
?geom_rect()

#####predict number od year germ is above 90%, betwenn 50 and 90 and less that 50 
new.data.2<-data.frame(chill_time=c(rnorm(100,79.2,26.63),rnorm(100,47.73,25.7)),force=rep(c(0,1),each=100),scenario=rep(c("ambient","warming"),each=100))
gooter<-fitted(mod1,newdata = new.data.2,summary = TRUE)
gooter2<-fitted(mod2,newdata = new.data.2,summary = TRUE)
gooter3<-fitted(mod3,newdata = new.data.2,summary = TRUE)
gooter4<-fitted(mod4,newdata = new.data.2,summary = TRUE)
gooter5<-fitted(mod5,newdata = new.data.2,summary = TRUE)
gooter6<-fitted(mod6,newdata = new.data.2,summary = TRUE)

Cprediction<-cbind(new.data.2,gooter) 
Cprediction$Taxa<-"P.virginiatum"


Cprediction2<-cbind(new.data.2,gooter2) 
Cprediction2$Taxa<-"C.canadensis"

Cprediction3<-cbind(new.data.2,gooter3) 
Cprediction3$Taxa<-"E.divaricata"

Cprediction4<-cbind(new.data.2,gooter4) 
Cprediction4$Taxa<-"T.dioicum"

Cprediction5<-cbind(new.data.2,gooter5) 
Cprediction5$Taxa<-"A.virginiana"

Cprediction6<-cbind(new.data.2,gooter6) 
Cprediction6$Taxa<-"S.stellata"

Cpredictionall<-rbind(Cprediction, Cprediction2,Cprediction3)

Cpredictionall$germination<-NA
Cpredictionall$germination[which(Cpredictionall$Estimate>=.80)]<-">80%"
Cpredictionall$germination[which(Cpredictionall$Estimate<.80 & Cpredictionall$Estimate>=.50 )]<-"50-80%"
Cpredictionall$germination[which(Cpredictionall$Estimate<.50 & Cpredictionall$Estimate>=.20 )]<-"20-50%"
Cpredictionall$germination[which(Cpredictionall$Estimate<.20)]<-"<20%"

#write.csv(Cpredictionall,"dogoober.csv") ### do some excell stuff



png("figures/germpercs_wvar.png",width = 5,height = 5,units = "in",res = 300)
ggplot(Cpredictionall,aes(scenario,Estimate))+
  geom_bar(stat="identity",position="dodge",aes(fill=Taxa))+
  ggthemes::theme_base()+scale_fill_brewer( type = "qual", palette = "Dark2", direction = 1)

dev.off()

png("figures/germfreq.png",width = 8,height = 5,units = "in",res = 300)
Cpredictionall %>%

    
  arrange(germination) %>%
  mutate(germination = factor(germination, levels=c(">80%","50-80%","20-50%","<20%"))) %>%
ggplot(aes(germination,fill=Taxa))+geom_histogram(stat="count",position="dodge")+facet_wrap(~scenario)+
  ggthemes::theme_base()+xlab("Germination")+scale_fill_brewer( type = "qual", palette = "Dark2", direction = 1)

dev.off()
save.image("goingforitgerm")


PVse<-read.csv("dogoober.csv")
PVse$ignore<-100*PVse$Estimate
PVse$year<-rep(1:100,6)
colnames(PVse)
ggplot(PVse,aes(year,germinated))+geom_line(aes(color=Taxa))+facet_wrap(~scenario)+geom_hline(yintercept = 100)+  scale_color_brewer( type = "qual", palette = "Dark2", direction = 1)+ggthemes::theme_calc()


ggplot(PVse,aes(year,seedbank))+geom_line(aes(color=Taxa))+facet_wrap(~scenario)+geom_hline(yintercept = 100)+ggthemes::theme_calc()


##simular seed bank dynamics
