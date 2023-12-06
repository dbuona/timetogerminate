###a script to amke my data survival formated
###need to add a way to elminate seeds that rotted.
rm(list=ls()) 
options(stringsAsFactors = FALSE)

library(survival)
library(ggplot2)
library(dplyr)
#library("ggfortify")


library(drc)
library(lme4)
library(brms)
library(ggstance)
library(tidybayes)
library(bayesplot)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
load("survmodel.Rda")
#save.image("survmodel.Rda")

d<-read.csv("..//survival_analysis/surival_dat_nointerval.csv")
d$DAY<-ifelse(d$DAY==0,0.00001,d$DAY)

unique(d$chill_time)
d$chillweeks<-d$chill_time/7 # make chilling weeks instead of days

d$force<-NA # make forcing numeric
d<- within(d, force[INC=="L"]<-0)
d<- within(d, force[INC=="H"]<-5)
###brms is backwords from other survival
d$censored<-ifelse(d$germinated==1,0,1)
d<- filter(d,!Taxa %in% c("Phlox cuspidata","Impatiens capensis","Carex grisea"))

d$censored<-ifelse(d$DAY==0.00001 & d$germinated==1,-1,d$censored)

###make histograms to help Lizzie and Megan
d.germ<-filter(d,germinated==1)
pdf("..//figures/histograms.pdf")
ggplot(d.germ,aes(DAY))+
  geom_density(aes(fill=as.factor(chillweeks),color=as.factor(chillweeks),alpha=0.3))+facet_grid(Taxa~as.factor(force),scales = "free")+
  theme_bw()
dev.off()
##try weibull

priorz.wei<-get_prior(DAY | cens(censored)~chillweeks+force+(chillweeks+force|Taxa),data=d,family= weibull())
fit.wei.all<- brm(DAY | cens(censored)~chillweeks+force+(chillweeks+force|Taxa), data=d, family =   weibull(),inits=0 ,prior=priorz.wei,iter=4000,warmup = 3000, control=list(adapt_delta=0.95),chains=4) 
fixef(fit.wei.all)

d.in<-dplyr::filter(d,Taxa %in% c("Hesperis matronalis","Cryptotaenia canadensis"))
colnames(d.in)
d.knb<-dplyr::select(d.in,-X,-INC,-COLD,-date,-warmT,-chill_time)
write.csv(d.knb,"germ_assays_knb.csv",row.names = FALSE)

priorz.wei.in<-get_prior(DAY | cens(censored)~chillweeks+force+Taxa+chillweeks:Taxa+force:Taxa,data=d.in,family= weibull())

fit.wei.in<-brm(DAY | cens(censored)~chillweeks+force+Taxa+chillweeks:Taxa+force:Taxa,data=d.in,family= weibull(),inits=0 ,prior=priorz.wei.in,iter=4000,warmup = 3000, chains=4) 



summary(fit.wei.all)
fixef(fit.wei.all,probs =c(0.05,0.25,.75,.95))
exp(fixef(fit.wei.all)[1])

exp(fixef(fit.wei.all)[1])-exp(fixef(fit.wei.all)[1]+fixef(fit.wei.all)[2]) #6.14 days decrease in t50 per week of chilling
exp(fixef(fit.wei.all)[1])-exp(fixef(fit.wei.all)[1]+fixef(fit.wei.all)[3]) #1.45 increase per degree C 
exp(fixef(fit.wei.all)[1]+fixef(fit.wei.all)[3])
yaya<-fit.wei.all%>%
  spread_draws(r_Taxa[Taxa,condition])
yaya<-filter(yaya,condition=="chillweeks")

yaya$species<-NA
yaya$species[which(yaya$Taxa=="Thalictrum.dioicum")]<-"Thalictrum dioicum"
yaya$species[which(yaya$Taxa=="Silene.vulgaris")]<-"Silene vulgaris"
yaya$species[which(yaya$Taxa=="Silene.stellata")]<-"Silene stellata"
yaya$species[which(yaya$Taxa=="Polygonum.virginiatum")]<-"Persicaria virginiana"
yaya$species[which(yaya$Taxa=="Eurbia.diviricata")]<-"Eurybia divaricata"
yaya$species[which(yaya$Taxa=="Oenethera.biennis")]<-"Oenethera biennis"
yaya$species[which(yaya$Taxa=="Hesperis.matronalis")]<-"Hesperis matronalis"
yaya$species[which(yaya$Taxa=="Cryptotaenia.canadensis")]<-"Cryptotaenia canadensis"
yaya$species[which(yaya$Taxa=="Carex.grayi")]<-"Carex grayi"
yaya$species[which(yaya$Taxa=="Asclepias.syriaca")]<-"Asclepias syriaca"
yaya$species[which(yaya$Taxa=="Anemone.virginana")]<-"Anemone virginiana"

yaya2<-fit.wei.all%>%
  spread_draws(b_chillweeks)
yaya$r_Taxa2<-yaya2$b_chillweeks+yaya$r_Taxa





yaya3<-fit.wei.all%>%
  spread_draws(r_Taxa[Taxa,condition])

yaya3<-filter(yaya3,condition=="force")

yaya3$species<-NA
yaya3$species[which(yaya3$Taxa=="Thalictrum.dioicum")]<-"Thalictrum dioicum"
yaya3$species[which(yaya3$Taxa=="Silene.vulgaris")]<-"Silene vulgaris"
yaya3$species[which(yaya3$Taxa=="Silene.stellata")]<-"Silene stellata"
yaya3$species[which(yaya3$Taxa=="Polygonum.virginiatum")]<-"Persicaria virginiana"
yaya3$species[which(yaya3$Taxa=="Eurbia.diviricata")]<-"Eurybia divaricata"
yaya3$species[which(yaya3$Taxa=="Oenethera.biennis")]<-"Oenethera biennis"
yaya3$species[which(yaya3$Taxa=="Hesperis.matronalis")]<-"Hesperis matronalis"
yaya3$species[which(yaya3$Taxa=="Cryptotaenia.canadensis")]<-"Cryptotaenia canadensis"
yaya3$species[which(yaya3$Taxa=="Carex.grayi")]<-"Carex grayi"
yaya3$species[which(yaya3$Taxa=="Asclepias.syriaca")]<-"Asclepias syriaca"
yaya3$species[which(yaya3$Taxa=="Anemone.virginana")]<-"Anemone virginiana"


yaya4<-fit.wei.all%>%
  spread_draws(b_force)
yaya3$r_Taxa2<-yaya4$b_force+yaya3$r_Taxa


yaya.n<-rbind(yaya,yaya3)
yaya.n$treatment<-ifelse(yaya.n$condition=="chillweeks","chilling","incubation")



mu1<-ggplot()+ stat_interval(data=yaya,aes(r_Taxa2,species),fill="skyblue1",.width = c(.5,.89,.975))+ggthemes::theme_few()+geom_vline(xintercept=0)+
  xlab("Estimated effect of chilling")+coord_cartesian(xlim=c(-.4,.25))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_text(face="italic"))+scale_y_discrete(limits=rev)
      

mu2<-ggplot()+ stat_interval(data=yaya3,aes(r_Taxa2,Taxa),fill="salmon",.width = c(.5,.89,.975))+ggthemes::theme_few()+geom_vline(xintercept=0)+
xlab("Estimated effect of incubation")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+scale_y_discrete(limits=rev)

jpeg("..//figures/mus_survival.jpeg",height=5,width=8, units="in",res=200)

ggplot()+ stat_interval(data=yaya.n,aes(r_Taxa2,reorder(species,-r_Taxa2)),.width = c(.5,.89,.975))+ggthemes::theme_few()+geom_vline(xintercept=0,color="firebrick4")+facet_wrap(~treatment)+
  xlab("estimated effect")+coord_cartesian(xlim=c(-.4,.4))+
  theme(axis.text.y=element_text(face="italic"),legend.position = "bottom")+ylab("")

#ggpubr::ggarrange(mu1,mu2,widths=c(.6,.4),common.legend = TRUE)#,labels = c("a)","b)"),hjust =c(-10, -1.5),vjust=-10)

dev.off()
brmsfamily("weibull")



####plot for all sps
pred.weeks.A<-c(6,12)
pred.force.A<-c(0,5)
rm(new.data.A)
new.data.A <- data.frame(Taxa=c(rep(c(unique(d$Taxa)),4)),
                        chillweeks = c(rep(pred.weeks.A,each=11)), force = c(rep(pred.force.A,each=22)))


pred<-fitted(fit.wei.all,newdata = new.data.A,probs = c(.055,.945))
new.data.A<-cbind(new.data.A,pred)



new.data.A$species<-NA
new.data.A$species[which(new.data.A$Taxa=="Thalictrum dioicum")]<-"Thalictrum dioicum"
new.data.A$species[which(new.data.A$Taxa=="Silene vulgaris")]<-"Silene vulgaris"
new.data.A$species[which(new.data.A$Taxa=="Silene stellata")]<-"Silene stellata"
new.data.A$species[which(new.data.A$Taxa=="Polygonum virginiatum")]<-"Persicaria virginiana"
new.data.A$species[which(new.data.A$Taxa=="Eurbia diviricata")]<-"Eurybia divaricata"
new.data.A$species[which(new.data.A$Taxa=="Oenethera biennis")]<-"Oenethera biennis"
new.data.A$species[which(new.data.A$Taxa=="Hesperis matronalis")]<-"Hesperis matronalis"
new.data.A$species[which(new.data.A$Taxa=="Cryptotaenia canadensis")]<-"Cryptotaenia canadensis"
new.data.A$species[which(new.data.A$Taxa=="Carex grayi")]<-"Carex grayi"
new.data.A$species[which(new.data.A$Taxa=="Asclepias syriaca")]<-"Asclepias syriaca"
new.data.A$species[which(new.data.A$Taxa=="Anemone virginana")]<-"Anemone virginiana"

smaps<-filter(new.data.A, species %in% c("Asclepias syriaca","Hesperis matronalis","Persicaria virginiana") )
smaps<-filter(smaps,force==0)
smaps$chilling<-ifelse(smaps$chillweeks==12, "12 weeks of chilling", "6 weeks of chilling")
smaps2<-filter(new.data.A, species %in% c("Persicaria virginiana","Hesperis matronalis") )
#smaps2<-filter(smaps2,force==0)

jpeg("..//figures/sps_case_examps.jpeg",height=5,width=8, units="in",res=200)
ggplot(smaps,aes(x=Estimate,y=0))+
  geom_point(aes(color=species),size=3)+geom_hline(yintercept=-.05)+
  geom_errorbarh(aes(xmin=`Q5.5`,xmax=`Q94.5`,color=species),height=0)+
  facet_wrap(~chilling,ncol=1)+scale_y_discrete()+
  scale_color_viridis_d(option="magma",begin = .1,end=.9)+ggthemes::theme_few()+xlab("predicted T50")+ylab("")+theme(legend.position = "bottom")+
  theme(legend.text = element_text(face = "italic"))
dev.off()

ggpubr::ggarrange(goob,exp,ncol=1)
ggplot(smaps2,aes(x=Estimate,y=0))+geom_point(aes(color=species),size=1,alpha=.8)+geom_errorbarh(aes(xmin=`Q5.5`,xmax=`Q94.5`,color=species),height=0)+facet_wrap(~chillweeks,ncol=1)+scale_y_discrete()


ggplot(new.data.A,aes(x = .epred,y=force))+stat_summary(aes(color=species,group=paste(species,force)))+facet_wrap(force~chillweeks)
 # geom_point(aes(y=.epred,group=paste(.draw,species,force)),size=.01)+
  stat_summary(fun=mean, geom="point", size = 1,aes(y=.epred,color=species,group=paste(Taxa,force)))+
  scale_color_viridis_d()+ggthemes::theme_few()+facet_wrap(~force)+
  labs(y="Model Estimated Days to 50% Germination",x="Weeks of cold stratification")+theme(legend.text = element_text(face = "italic"),legend.position = "top")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16))+
  scale_y_continuous(breaks = c(0,50,100,150,200,250,300))+theme(legend.title= element_blank())+ylim(0,50)
dev.off()



pred.weeks3<-c(0:16)
pred.force3<-c(0,5)
new.data3 <- data.frame(Taxa=c(rep(c(unique(d$Taxa)),each=17)),
                        chillweeks = c(rep(pred.weeks3,11)), force = c(rep(pred.force3,each=187)))


library(tidybayes)
new.data3<-new.data3 %>% add_epred_draws(fit.wei.all, ndraws = 100)


new.data3$species<-NA
new.data3$species[which(new.data3$Taxa=="Thalictrum dioicum")]<-"Thalictrum dioicum"
new.data3$species[which(new.data3$Taxa=="Silene vulgaris")]<-"Silene vulgaris"
new.data3$species[which(new.data3$Taxa=="Silene stellata")]<-"Silene stellata"
new.data3$species[which(new.data3$Taxa=="Polygonum virginiatum")]<-"Persicaria virginiana"
new.data3$species[which(new.data3$Taxa=="Eurbia diviricata")]<-"Eurybia divaricata"
new.data3$species[which(new.data3$Taxa=="Oenethera biennis")]<-"Oenethera biennis"
new.data3$species[which(new.data3$Taxa=="Hesperis matronalis")]<-"Hesperis matronalis"
new.data3$species[which(new.data3$Taxa=="Cryptotaenia canadensis")]<-"Cryptotaenia canadensis"
new.data3$species[which(new.data3$Taxa=="Carex grayi")]<-"Carex grayi"
new.data3$species[which(new.data3$Taxa=="Asclepias syriaca")]<-"Asclepias syriaca"
new.data3$species[which(new.data3$Taxa=="Anemone virginana")]<-"Anemone virginiana"
new.data3$incubation<-ifelse(new.data3$force==0,"low incubation","high incubation")

jpeg("..//figures/surv_prieff.jpeg",height=5,width=9, units="in",res=200)
ggplot(new.data3,aes(x = chillweeks, color = species))+
  geom_line(aes(y=.epred,group=paste(.draw,species,force)),alpha=.05)+
  stat_summary(fun=mean, geom="line", size = .75,aes(y=.epred,color=species,group=paste(species,force)))+
  scale_color_viridis_d(option = "turbo")+ggthemes::theme_few()+facet_wrap(~incubation)+
  labs(y="T50",x="Weeks of chilling")+theme(legend.text = element_text(face = "italic"),legend.position = "bottom")+
  scale_x_continuous(breaks = c(0,4,8,12,16))+coord_cartesian(ylim=c(0,365))+
  scale_y_continuous()+theme(legend.title= element_blank())
dev.off()

###########################
###### make plot 1############
pred.weeks<-c(5:16)
pred.force<-c(0)

new.data <- data.frame(Taxa=c("Cryptotaenia canadensis"),
  chillweeks = c(rep(pred.weeks,12)), force = c(rep(pred.force,each=12)))

pred.weeks2<-c(10:16)
pred.force2<-c(5)

new.data2 <- data.frame(Taxa=c("Cryptotaenia canadensis"),
                       chillweeks = c(rep(pred.weeks2,7)), force = c(rep(pred.force2,each=7)))



pred.weeks3<-c(0:16)
pred.force3<-c(0,5)
new.data3 <- data.frame(Taxa=c(rep(c(unique(yaya$Taxa)),each=17)),
                       chillweeks = c(rep(pred.weeks3,11)), force = c(rep(pred.force3,each=187)))


library(tidybayes)
new.data<-new.data %>% add_epred_draws(fit.wei.in, ndraws = 200)
new.data2<-new.data2 %>% add_epred_draws(fit.wei.in, ndraws = 200)
new.data4<-filter(new.data3,Taxa=="Hesperis matronalis")
new.data4<-new.data4 %>% add_epred_draws(fit.wei.in, ndraws = 200)

new.data$incubation<-"25/15"
new.data2$incubation<-"20/10"
new.data4$incubation<-ifelse(new.data4$force==0,"20/10","25/15")
new.data<-rbind(new.data,new.data2,new.data4)


new.data$Taxa2<-ifelse(new.data$Taxa=="Cryptotaenia canadensis","Cryptotaenia canadensis \n(native)","Hesperis matronalis \n(invasive)")

jpeg("..//figures/AFTsivansive.jpeg",height=8,width=8, units="in",res=200)
  
ggplot(new.data,aes(x = chillweeks, color = Taxa2))+
geom_line(aes(y=.epred,group=paste(.draw,Taxa2,incubation)),alpha=.05)+
  stat_summary(fun=mean, geom="line", size = .75,aes(y=.epred,color=Taxa2,group=paste(Taxa,incubation)))+
scale_color_viridis_d(begin=0,end=0.5)+ggthemes::theme_few()+
  labs(y="Model Estimated Days to 50% Germination",x="Weeks of cold stratification")+theme(legend.text = element_text(face = "italic"),legend.position = "top")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16))+
  scale_y_continuous(breaks = c(0,5,10,15,20,30,40,50))+theme(legend.title= element_blank())
dev.off()

  
  
  b<-new.data2 %>% add_epred_draws(fit.wei.in, ndraws = 200) %>%
  ggplot(aes(x = chillweeks, color = Taxa))+
  geom_line(aes(y=.epred,group=paste(.draw,Taxa)),alpha=.05)+
  stat_summary(fun=mean, geom="line", size = .75,aes(y=.epred,color=Taxa))+
  scale_color_viridis_d(begin=0,end=0.5)+ggthemes::theme_few()+
  labs(y="Model Estimated Days to 50% Germination",x="Weeks of cold stratification")+theme(legend.text = element_text(face = "italic"))+
  scale_x_continuous(breaks = c(10,12,14,16))+
  scale_y_continuous(breaks = c(5,10,15,20,30,40))

#scale_color_manual(values=c("#f0f921","#fdc527","#f89540","#440154","#e66c5c","#21918c","#cc4778","#aa2395","#7e03a8","#4c02a1", "#0d0887"))

new.data3$incubation<-ifelse(new.data3$force==0,"20/10","25/15")

c<-new.data3 %>% add_epred_draws(fit.wei.all, ndraws = 100) %>%
  ggplot(aes(x = chillweeks, color = Taxa))+
  scale_color_manual(values=c("#f0f921","#fdc527","#f89540","#440154","#e66c5c","#21918c","#cc4778","#aa2395","#7e03a8","#4c02a1", "#0d0887"))+
  geom_line(aes(y=.epred,group=paste(.draw,Taxa)),alpha=.05)+
  stat_summary(fun=mean, geom="line", size = .75,aes(y=.epred,color=Taxa))+
  ggthemes::theme_few()+facet_wrap(~incubation)+
  labs(y="",x="Weeks of cold stratification")+theme(legend.text = element_text(face = "italic"))

#fde725,#bddf26,#7ad151,#440154, #21918c,#2a788e,#355f8d,#414487,#482475,#2ab07f,#52c569

pred.weeks4<-c(6,12)
pred.force4<-c(0,5)
new.data4 <- data.frame(Taxa=c(rep(c(unique($Taxa)),each=2)),
                        chillweeks = c(rep(pred.weeks4,11)), force = c(rep(pred.force4,each=22)))

daty4<-fitted(fit.wei.all,probs =c(0.05,0.95),newdata=new.data4)### something is wrong with error
daty4<-cbind(daty4,new.data4)
daty4$incubation<-ifelse(daty4$force==0,"20/10","25/15")
daty4$stratification<-ifelse(daty4$chillweeks==6,"6 weeks","12 weeks")


daty4$invasive<-NA
daty4$invasive<-ifelse(new.data4$Taxa %in% c("Hesperis matronalis","Silene vulgaris"),"Invasive","Native")
daty4<-filter(daty4,!Taxa %in% c("Carex grayi","Thalictrum dioicum","Silene stellata"))
daty4$scenario<-NA

daty4$scenario[which(daty4$stratification=="12 weeks"& daty4$incubation=="20/10")]<-1
daty4$scenario[which(daty4$stratification=="12 weeks"& daty4$incubation=="25/15")]<-2
daty4$scenario[which(daty4$stratification=="6 weeks"& daty4$incubation=="25/15")]<-4
daty4$scenario[which(daty4$stratification=="6 weeks"& daty4$incubation=="20/10")]<-3

###clean species
daty4$Taxa[which(daty4$Taxa=="Eurbia diviricata")]<-"Eurybia divaricata"
daty4$Taxa[which(daty4$Taxa=="Polygonum virginiatum")]<-"Persicaria virginiana"
daty4$Taxa[which(daty4$Taxa=="Anemone virginana")]<-"Anemone virginiana"

ggplot(daty4,aes(Estimate,Taxa))+geom_point(aes(shape=invasive,group=chillweeks,color=scenario),size=3.5)+
scale_shape_manual(values = c(15,16))+facet_grid(incubation~stratification,scales = "free")+
geom_errorbarh(aes(xmin=Q5,xmax=Q95),height=0)+
  scale_color_viridis_c(option="turbo",direction = 1, begin=.2,end=.8)+
  ggthemes::theme_few(base_size = 11)+
  theme(axis.text.y = element_text(face = c("italic","italic","italic","italic","bold.italic","italic","italic","bold.italic"))
)+xlab("Day of Season")+ylab("")+xlim(0,20)+ guides(col = FALSE)
 
daty4a<-dplyr::filter(daty4,scenario %in% c(1,4))
daty4a$scenario2<-ifelse(daty4a$scenario==1,"historic","warming")
cc<-ggplot(daty4a,aes(Estimate,Taxa))+geom_point(aes(shape=invasive,group=chillweeks,color=scenario2),size=3,alpha=0.75)+
  scale_shape_manual(values = c(15,16))+
  geom_errorbarh(aes(xmin=Q5,xmax=Q95,color=scenario2),height=0)+
  scale_color_viridis_d(option="turbo",direction = 1, begin=.3,end=.95)+
  ggthemes::theme_few(base_size = 11)+
  theme(axis.text.y = element_text(face = c("bold.italic","italic","italic","bold.italic","italic","italic","italic","italic")))+
  xlab("Day of season for 50% germination")+ylab("")+guides(scale = "none")+geom_vline(xintercept=20, linetype="dashed")+theme(legend.title=element_blank(),legend.position = "top")+
  scale_y_discrete(limits=rev)

dev.off()


jpeg("..//figures/AFTall.jpeg",height=5,width=8, units="in",res=200)
c
dev.off()

jpeg("..//figures/commchange.jpeg",height=6,width=10, units="in",res=300)
cc
dev.off()

jpeg("..//figures/AFTsivansive.jpeg",height=8,width=8, units="in",res=200)
#ggpubr::ggarrange(a,b,common.legend = TRUE,legend="right",widths=c(.6,.3),labels = c("a)","b)"))
ggplot(new.data,aes(x = chillweeks, color = Taxa))+
  geom_line(aes(y=.epred,group=paste(.draw,Taxa,incubation)),alpha=.05)+
  stat_summary(fun=mean, geom="line", size = .75,aes(y=.epred,color=Taxa,group=paste(Taxa,incubation)))+
  scale_color_viridis_d(begin=0,end=0.5)+ggthemes::theme_few()+
  labs(y="Model Estimated Days to 50% Germination",x="Weeks of cold stratification")+theme(legend.text = element_text(face = "italic"),legend.position = "top")+
  scale_x_continuous(breaks = c(0,2,4,6,8,10,12,14,16))+
  scale_y_continuous(breaks = c(0,5,10,15,20,30,40,50))+theme(legend.title= element_blank())+scale_y_discrete()

dev.off()





stop()### below is scratch


daty.wei.all<-fitted(fit.wei.in,probs =c(0.055,0.25,.75,.945),newdata=new.data)### something is wrong with error
daty.wei<-cbind(daty.wei.all,new.data)
daty.wei$incubation<-ifelse(daty.wei$force==0,"20/10","25/15")



pd<-position_dodge(width=0.1)

library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 11
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
ggplot()+
  geom_line(data=daty.wei,aes(x=chillweeks,y=Estimate,color=Taxa),size=1)+

 #geom_line(data=daty.wei,aes(x=chillweeks,y=Q25,color=Taxa),size=.5,alpha=0.5)+
  #geom_line(data=daty.wei,aes(x=chillweeks,y=Q75,color=Taxa),size=.5,alpha=0.5)+
  geom_line(data=daty.wei,aes(x=chillweeks,y=Q5.5,color=Taxa),size=.5,alpha=0.5,linetype="dashed")+
  geom_line(data=daty.wei,aes(x=chillweeks,y=Q94.5,color=Taxa),size=.5,alpha=0.5,linetype="dashed")+##geom_point(aes(color=Taxa),position=pd)+
  #ggplot2::geom_errorbar(aes(ymax=Q75,ymin=Q25,color=Taxa),width=0,position=pd)+
  facet_wrap(~incubation,nrow = 2)+ggthemes::theme_few(base_size = 11)+scale_color_viridis_d(begin = 0,end=0.5)+ylim(0,40)+xlim(4,16)+
  labs(y="Model Estimated Days to 50% Germination",x="Weeks of cold stratification")+theme(legend.text = element_text(face = "italic"))
 
sp.use<-dplyr::filter(daty.wei,Taxa %in% c("Hesperis matronalis","Cryptotaenia canadensis"))

c<-ggplot()+
  geom_line(data=sp.use,aes(x=chillweeks,y=Estimate,color=Taxa),size=1)+
  
  geom_line(data=sp.use,aes(x=chillweeks,y=Q25,color=Taxa),size=.5,alpha=0.5)+
  geom_line(data=sp.use,aes(x=chillweeks,y=Q75,color=Taxa),size=.5,alpha=0.5)+
  geom_line(data=sp.use,aes(x=chillweeks,y=Q5.5,color=Taxa),size=.5,alpha=0.5,linetype="dashed")+
  geom_line(data=sp.use,aes(x=chillweeks,y=Q94.5,color=Taxa),size=.5,alpha=0.5,linetype="dashed")+##geom_point(aes(color=Taxa),position=pd)+
  #ggplot2::geom_errorbar(aes(ymax=Q75,ymin=Q25,color=Taxa),width=0,position=pd)+
  facet_wrap(~incubation,nrow = 2)+ggthemes::theme_base(base_size = 11)+scale_color_viridis_d()+ylim(0,60)+
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

jpeg("..//figures/AFTplots_spcomp.jpeg",height=8,width=10, units="in",res=300)
ggpubr::ggarrange(a,c,common.legend = TRUE,nrow=1,widths = c(.6,.45),legend="bottom",labels = c("a)","b)"))
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

