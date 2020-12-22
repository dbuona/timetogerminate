###a script to amke my data survival formated
###need to add a way to elminate seeds that rotted.
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(survival)
library(ggplot2)
library(dplyr)
library("ggfortify")
library("icenReg")
library("emmeans")
library(drc)
library(lme4)
library(brms)
library(ggstance)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
load("survmodel.Rda")

d<-read.csv("..//survival_analysis/surival_dat_nointerval.csv")
d$DAY<-ifelse(d$DAY==0,0.00001,d$DAY)


d$chillweeks<-d$chill_time/7 # make chilling weeks instead of days

d$force<-NA # make forcing numeric
d<- within(d, force[INC=="L"]<-0)
d<- within(d, force[INC=="H"]<-5)
###brms is backwords from other survival
d$censored<-ifelse(d$germinated==1,0,1)
d<- filter(d,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea"))

d$censored<-ifelse(d$DAY==0.00001 & d$germinated==1,-1,d$censored)


##try weibull

priorz.wei<-get_prior(DAY | cens(censored)~chillweeks+force+(chillweeks+force|Taxa),data=d,family= weibull())

fit.wei.all<- brm(DAY | cens(censored)~chillweeks+force+(chillweeks+force|Taxa), data=d, family =   weibull(),inits=0 ,prior=priorz.wei,iter=4000,warmup = 3000, chains=4) 




summary(fit.wei.all)
coef(fit.wei.all,probs =c(0.55,0.25,.75,.945))

brmsfamily("weibull")

###########################
###### make plot 1############
pred.weeks<-c(1:16)
pred.force<-c(0,5)
unique(d$Taxa)
new.data <- data.frame(Taxa=c(rep(c(unique(d$Taxa)),each=16,2)),
  chillweeks = c(rep(pred.weeks,22)), force = c(rep(pred.force,each=352)))



daty.wei.all<-fitted(fit.wei.all,probs =c(0.055,0.25,.75,.945),newdata=new.data)### something is wrong with error
daty.wei<-cbind(daty.wei.all,new.data)
daty.wei$incubation<-ifelse(daty.wei$force==0,"20/10","25/15")



pd<-position_dodge(width=0.1)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
a<-ggplot()+
  geom_line(data=daty.wei,aes(x=chillweeks,y=Estimate,color=Taxa),size=1)+

 geom_line(data=daty.wei,aes(x=chillweeks,y=Q25,color=Taxa),size=.5,alpha=0.5)+
  geom_line(data=daty.wei,aes(x=chillweeks,y=Q75,color=Taxa),size=.5,alpha=0.5)+
  geom_line(data=daty.wei,aes(x=chillweeks,y=Q5.5,color=Taxa),size=.5,alpha=0.5,linetype="dashed")+
  geom_line(data=daty.wei,aes(x=chillweeks,y=Q94.5,color=Taxa),size=.5,alpha=0.5,linetype="dashed")+##geom_point(aes(color=Taxa),position=pd)+
  #ggplot2::geom_errorbar(aes(ymax=Q75,ymin=Q25,color=Taxa),width=0,position=pd)+
  facet_wrap(~incubation,nrow = 2)+ggthemes::theme_base(base_size = 11)+scale_color_manual(values = mycolors)+ylim(0,60)+
  labs(y="Model Estimated Days to 50% Germination",x="Weeks of cold stratification")+theme(legend.text = element_text(face = "italic"))
 

###############################
###### make plot 2############
##########################
pred.weeks2<-c(8,12)
pred.force2<-c(5,0)
unique(d$Taxa)
new.data2 <- data.frame(Taxa=c(rep(c(unique(d$Taxa)), 2)),
                       chillweeks = c(rep(pred.weeks2,each=11)), force = c(rep(pred.force2,each=11)))

daty.wei.all2<-fitted(fit.wei.all,probs =c(0.055,0.25,.75,.945),newdata=new.data2)### something is wrong with error
daty.wei2<-cbind(daty.wei.all2,new.data2)

daty.wei2$scenario<-ifelse(daty.wei2$chillweeks==12,"average","warming")
b<-ggplot(daty.wei2,aes(scenario,Estimate))+geom_point(aes(color=Taxa),position=pd)+
  ggplot2::geom_errorbar(aes(ymax=Q75,ymin=Q25,color=Taxa),width=0,position=pd)+
  ggplot2::geom_errorbar(aes(ymax=Q5.5,ymin=Q94.5,color=Taxa),width=0,linetype="dashed",position=pd)+
  ylim(0,42)+
  ggthemes::theme_base(base_size = 11)+scale_color_manual(values = mycolors)+
  geom_line(aes(x=as.factor(scenario),y=Estimate,group=Taxa,color=Taxa),linetype="dotted")+
  labs(y="Model Estimated Days to 50% Germination",x="scenario")

jpeg("..//figures/AFTplots.jpeg",height=8,width=10, units="in",res=300)
ggpubr::ggarrange(a,b,common.legend = TRUE,nrow=1,widths = c(.7,.3),legend="bottom",labels = c("a)","b)"))
dev.off()

###########table###############
Taxa=unique(daty.wei2$Taxa)
HP<-data.frame(Taxa=unique(daty.wei2$Taxa),HP_warm=c(
daty.wei2[1,1]-daty.wei2[1,1],
daty.wei2[1,1]-daty.wei2[2,1],
daty.wei2[1,1]-daty.wei2[3,1],
daty.wei2[1,1]-daty.wei2[4,1],
daty.wei2[1,1]-daty.wei2[5,1],
daty.wei2[1,1]-daty.wei2[6,1],
daty.wei2[1,1]-daty.wei2[7,1],
daty.wei2[1,1]-daty.wei2[8,1],
daty.wei2[1,1]-daty.wei2[9,1],
daty.wei2[1,1]-daty.wei2[10,1],
daty.wei2[1,1]-daty.wei2[11,1]))

HP$SVwarm<-c(daty.wei2[2,1]-daty.wei2[1,1],
daty.wei2[2,1]-daty.wei2[2,1],
daty.wei2[2,1]-daty.wei2[3,1],
daty.wei2[2,1]-daty.wei2[4,1],
daty.wei2[2,1]-daty.wei2[5,1],
daty.wei2[2,1]-daty.wei2[6,1],
daty.wei2[2,1]-daty.wei2[7,1],
daty.wei2[2,1]-daty.wei2[8,1],
daty.wei2[2,1]-daty.wei2[9,1],
daty.wei2[2,1]-daty.wei2[10,1],
daty.wei2[2,1]-daty.wei2[11,1])

HP$OBwarm<-c(daty.wei2[3,1]-daty.wei2[1,1],
             daty.wei2[3,1]-daty.wei2[2,1],
             daty.wei2[3,1]-daty.wei2[3,1],
             daty.wei2[3,1]-daty.wei2[4,1],
             daty.wei2[3,1]-daty.wei2[5,1],
             daty.wei2[3,1]-daty.wei2[6,1],
             daty.wei2[3,1]-daty.wei2[7,1],
             daty.wei2[3,1]-daty.wei2[8,1],
             daty.wei2[3,1]-daty.wei2[9,1],
             daty.wei2[3,1]-daty.wei2[10,1],
             daty.wei2[3,1]-daty.wei2[11,1])

HP$ASwarm<-c(daty.wei2[4,1]-daty.wei2[1,1],
             daty.wei2[4,1]-daty.wei2[2,1],
             daty.wei2[4,1]-daty.wei2[3,1],
             daty.wei2[4,1]-daty.wei2[4,1],
             daty.wei2[4,1]-daty.wei2[5,1],
             daty.wei2[4,1]-daty.wei2[6,1],
             daty.wei2[4,1]-daty.wei2[7,1],
             daty.wei2[4,1]-daty.wei2[8,1],
             daty.wei2[4,1]-daty.wei2[9,1],
             daty.wei2[4,1]-daty.wei2[10,1],
             daty.wei2[4,1]-daty.wei2[11,1])

HP$SSwarm<-c(daty.wei2[5,1]-daty.wei2[1,1],
             daty.wei2[5,1]-daty.wei2[2,1],
             daty.wei2[5,1]-daty.wei2[3,1],
             daty.wei2[5,1]-daty.wei2[4,1],
             daty.wei2[5,1]-daty.wei2[5,1],
             daty.wei2[5,1]-daty.wei2[6,1],
             daty.wei2[5,1]-daty.wei2[7,1],
             daty.wei2[5,1]-daty.wei2[8,1],
             daty.wei2[5,1]-daty.wei2[9,1],
             daty.wei2[5,1]-daty.wei2[10,1],
             daty.wei2[5,1]-daty.wei2[11,1])

HP$EDwarm<-c(daty.wei2[6,1]-daty.wei2[1,1],
             daty.wei2[6,1]-daty.wei2[2,1],
             daty.wei2[6,1]-daty.wei2[3,1],
             daty.wei2[6,1]-daty.wei2[4,1],
             daty.wei2[6,1]-daty.wei2[5,1],
             daty.wei2[6,1]-daty.wei2[6,1],
             daty.wei2[6,1]-daty.wei2[7,1],
             daty.wei2[6,1]-daty.wei2[8,1],
             daty.wei2[6,1]-daty.wei2[9,1],
             daty.wei2[6,1]-daty.wei2[10,1],
             daty.wei2[6,1]-daty.wei2[11,1])
  
HP$AVwarm<-c(daty.wei2[7,1]-daty.wei2[1,1],
             daty.wei2[7,1]-daty.wei2[2,1],
             daty.wei2[7,1]-daty.wei2[3,1],
             daty.wei2[7,1]-daty.wei2[4,1],
             daty.wei2[7,1]-daty.wei2[5,1],
             daty.wei2[7,1]-daty.wei2[6,1],
             daty.wei2[7,1]-daty.wei2[7,1],
             daty.wei2[7,1]-daty.wei2[8,1],
             daty.wei2[7,1]-daty.wei2[9,1],
             daty.wei2[7,1]-daty.wei2[10,1],
             daty.wei2[7,1]-daty.wei2[11,1])


HP$CCwarm<-c(daty.wei2[8,1]-daty.wei2[1,1],
             daty.wei2[8,1]-daty.wei2[2,1],
             daty.wei2[8,1]-daty.wei2[3,1],
             daty.wei2[8,1]-daty.wei2[4,1],
             daty.wei2[8,1]-daty.wei2[5,1],
             daty.wei2[8,1]-daty.wei2[6,1],
             daty.wei2[8,1]-daty.wei2[7,1],
             daty.wei2[8,1]-daty.wei2[8,1],
             daty.wei2[8,1]-daty.wei2[9,1],
             daty.wei2[8,1]-daty.wei2[10,1],
             daty.wei2[8,1]-daty.wei2[11,1])

HP$CGwarm<-c(daty.wei2[9,1]-daty.wei2[1,1],
             daty.wei2[9,1]-daty.wei2[2,1],
             daty.wei2[9,1]-daty.wei2[3,1],
             daty.wei2[9,1]-daty.wei2[4,1],
             daty.wei2[9,1]-daty.wei2[5,1],
             daty.wei2[9,1]-daty.wei2[6,1],
             daty.wei2[9,1]-daty.wei2[7,1],
             daty.wei2[9,1]-daty.wei2[8,1],
             daty.wei2[9,1]-daty.wei2[9,1],
             daty.wei2[9,1]-daty.wei2[10,1],
             daty.wei2[9,1]-daty.wei2[11,1])

HP$TDwarm<-c(daty.wei2[10,1]-daty.wei2[1,1],
             daty.wei2[10,1]-daty.wei2[2,1],
             daty.wei2[10,1]-daty.wei2[3,1],
             daty.wei2[10,1]-daty.wei2[4,1],
             daty.wei2[10,1]-daty.wei2[5,1],
             daty.wei2[10,1]-daty.wei2[6,1],
             daty.wei2[10,1]-daty.wei2[7,1],
             daty.wei2[10,1]-daty.wei2[8,1],
             daty.wei2[10,1]-daty.wei2[9,1],
             daty.wei2[10,1]-daty.wei2[10,1],
             daty.wei2[10,1]-daty.wei2[11,1])

HP$PVwarm<-c(daty.wei2[11,1]-daty.wei2[1,1],
             daty.wei2[11,1]-daty.wei2[2,1],
             daty.wei2[11,1]-daty.wei2[3,1],
             daty.wei2[11,1]-daty.wei2[4,1],
             daty.wei2[11,1]-daty.wei2[5,1],
             daty.wei2[11,1]-daty.wei2[6,1],
             daty.wei2[11,1]-daty.wei2[7,1],
             daty.wei2[11,1]-daty.wei2[8,1],
             daty.wei2[11,1]-daty.wei2[9,1],
             daty.wei2[11,1]-daty.wei2[10,1],
             daty.wei2[11,1]-daty.wei2[11,1])


#######
HP2<-data.frame(Taxa=unique(daty.wei2$Taxa),HP_warm=c(
  daty.wei2[12,1]-daty.wei2[12,1],
  daty.wei2[12,1]-daty.wei2[13,1],
  daty.wei2[12,1]-daty.wei2[14,1],
  daty.wei2[12,1]-daty.wei2[15,1],
  daty.wei2[12,1]-daty.wei2[16,1],
  daty.wei2[12,1]-daty.wei2[17,1],
  daty.wei2[12,1]-daty.wei2[18,1],
  daty.wei2[12,1]-daty.wei2[19,1],
  daty.wei2[12,1]-daty.wei2[20,1],
  daty.wei2[12,1]-daty.wei2[21,1],
  daty.wei2[12,1]-daty.wei2[22,1]))

HP2$SVwarm<-c(daty.wei2[13,1]-daty.wei2[12,1],
             daty.wei2[13,1]-daty.wei2[13,1],
             daty.wei2[13,1]-daty.wei2[14,1],
             daty.wei2[13,1]-daty.wei2[15,1],
             daty.wei2[13,1]-daty.wei2[16,1],
             daty.wei2[13,1]-daty.wei2[17,1],
             daty.wei2[13,1]-daty.wei2[18,1],
             daty.wei2[13,1]-daty.wei2[19,1],
             daty.wei2[13,1]-daty.wei2[20,1],
             daty.wei2[13,1]-daty.wei2[21,1],
             daty.wei2[13,1]-daty.wei2[22,1])

HP2$OBwarm<-c(daty.wei2[14,1]-daty.wei2[12,1],
             daty.wei2[14,1]-daty.wei2[13,1],
             daty.wei2[14,1]-daty.wei2[14,1],
             daty.wei2[14,1]-daty.wei2[15,1],
             daty.wei2[14,1]-daty.wei2[16,1],
             daty.wei2[14,1]-daty.wei2[17,1],
             daty.wei2[14,1]-daty.wei2[18,1],
             daty.wei2[14,1]-daty.wei2[19,1],
             daty.wei2[14,1]-daty.wei2[20,1],
             daty.wei2[14,1]-daty.wei2[21,1],
             daty.wei2[14,1]-daty.wei2[22,1])

HP2$ASwarm<-c(daty.wei2[15,1]-daty.wei2[12,1],
             daty.wei2[15,1]-daty.wei2[13,1],
             daty.wei2[15,1]-daty.wei2[14,1],
             daty.wei2[15,1]-daty.wei2[15,1],
             daty.wei2[15,1]-daty.wei2[16,1],
             daty.wei2[15,1]-daty.wei2[17,1],
             daty.wei2[15,1]-daty.wei2[18,1],
             daty.wei2[15,1]-daty.wei2[19,1],
             daty.wei2[15,1]-daty.wei2[20,1],
             daty.wei2[15,1]-daty.wei2[21,1],
             daty.wei2[15,1]-daty.wei2[22,1])

HP2$SSwarm<-c(daty.wei2[16,1]-daty.wei2[12,1],
             daty.wei2[16,1]-daty.wei2[13,1],
             daty.wei2[16,1]-daty.wei2[14,1],
             daty.wei2[16,1]-daty.wei2[15,1],
             daty.wei2[16,1]-daty.wei2[16,1],
             daty.wei2[16,1]-daty.wei2[17,1],
             daty.wei2[16,1]-daty.wei2[18,1],
             daty.wei2[16,1]-daty.wei2[19,1],
             daty.wei2[16,1]-daty.wei2[20,1],
             daty.wei2[16,1]-daty.wei2[21,1],
             daty.wei2[16,1]-daty.wei2[22,1])

HP2$EDwarm<-c(daty.wei2[17,1]-daty.wei2[12,1],
             daty.wei2[17,1]-daty.wei2[13,1],
             daty.wei2[17,1]-daty.wei2[14,1],
             daty.wei2[17,1]-daty.wei2[15,1],
             daty.wei2[17,1]-daty.wei2[16,1],
             daty.wei2[17,1]-daty.wei2[17,1],
             daty.wei2[17,1]-daty.wei2[18,1],
             daty.wei2[17,1]-daty.wei2[19,1],
             daty.wei2[17,1]-daty.wei2[20,1],
             daty.wei2[17,1]-daty.wei2[21,1],
             daty.wei2[17,1]-daty.wei2[22,1])

HP2$AVwarm<-c(daty.wei2[18,1]-daty.wei2[12,1],
             daty.wei2[18,1]-daty.wei2[13,1],
             daty.wei2[18,1]-daty.wei2[14,1],
             daty.wei2[18,1]-daty.wei2[15,1],
             daty.wei2[18,1]-daty.wei2[16,1],
             daty.wei2[18,1]-daty.wei2[17,1],
             daty.wei2[18,1]-daty.wei2[18,1],
             daty.wei2[18,1]-daty.wei2[19,1],
             daty.wei2[18,1]-daty.wei2[20,1],
             daty.wei2[18,1]-daty.wei2[21,1],
             daty.wei2[18,1]-daty.wei2[22,1])


HP2$CCwarm<-c(daty.wei2[19,1]-daty.wei2[12,1],
             daty.wei2[19,1]-daty.wei2[13,1],
             daty.wei2[19,1]-daty.wei2[14,1],
             daty.wei2[19,1]-daty.wei2[15,1],
             daty.wei2[19,1]-daty.wei2[16,1],
             daty.wei2[19,1]-daty.wei2[17,1],
             daty.wei2[19,1]-daty.wei2[18,1],
             daty.wei2[19,1]-daty.wei2[19,1],
             daty.wei2[19,1]-daty.wei2[20,1],
             daty.wei2[19,1]-daty.wei2[21,1],
             daty.wei2[19,1]-daty.wei2[22,1])

HP2$CGwarm<-c(daty.wei2[20,1]-daty.wei2[12,1],
             daty.wei2[20,1]-daty.wei2[13,1],
             daty.wei2[20,1]-daty.wei2[14,1],
             daty.wei2[20,1]-daty.wei2[15,1],
             daty.wei2[20,1]-daty.wei2[16,1],
             daty.wei2[20,1]-daty.wei2[17,1],
             daty.wei2[20,1]-daty.wei2[18,1],
             daty.wei2[20,1]-daty.wei2[19,1],
             daty.wei2[20,1]-daty.wei2[20,1],
             daty.wei2[20,1]-daty.wei2[21,1],
             daty.wei2[20,1]-daty.wei2[22,1])

HP2$TDwarm<-c(daty.wei2[21,1]-daty.wei2[12,1],
             daty.wei2[21,1]-daty.wei2[13,1],
             daty.wei2[21,1]-daty.wei2[14,1],
             daty.wei2[21,1]-daty.wei2[15,1],
             daty.wei2[21,1]-daty.wei2[16,1],
             daty.wei2[21,1]-daty.wei2[17,1],
             daty.wei2[21,1]-daty.wei2[18,1],
             daty.wei2[21,1]-daty.wei2[19,1],
             daty.wei2[21,1]-daty.wei2[20,1],
             daty.wei2[21,1]-daty.wei2[21,1],
             daty.wei2[21,1]-daty.wei2[22,1])

HP2$PVwarm<-c(daty.wei2[22,1]-daty.wei2[12,1],
             daty.wei2[22,1]-daty.wei2[13,1],
             daty.wei2[22,1]-daty.wei2[14,1],
             daty.wei2[22,1]-daty.wei2[15,1],
             daty.wei2[22,1]-daty.wei2[16,1],
             daty.wei2[22,1]-daty.wei2[17,1],
             daty.wei2[22,1]-daty.wei2[18,1],
             daty.wei2[22,1]-daty.wei2[19,1],
             daty.wei2[22,1]-daty.wei2[20,1],
             daty.wei2[22,1]-daty.wei2[21,1],
             daty.wei2[22,1]-daty.wei2[22,1])




gooby<-data.frame(HP[1:1])
gooby<-cbind(gooby,round(HP2[2:12]-HP[2:12],digits=1))
colnames(gooby)<-c("Taxa",unique(HP$Taxa))
gooby[upper.tri(gooby)] <- ""

library(xtable)
xtable(gooby)
save.image("survmodel.Rda")
