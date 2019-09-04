rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")


library(rstan)
library(brms)
#parameters alpha b1, b2, beta, sigma
default.parameter.values<-c(4.13,4.4,11.6,.5,0.1)
results<-array(NA)

model.fuction<-

Di_c<-c(5,4,3,2,1,0,10,8,6,4,2,0,20,16,12,8,4,0)
Di_h<-c(0,1,2,3,4,5,0,2,4,6,8,10,0,4,8,12,16,20)
Di_c<-rep(Di_c,2)
Di_h<-rep(Di_h,2)

plot(Di_c~Di_h)

Gt_h<-c(rep(2,18),rep(8,18))
Gt_c<-c(rep(10,18),rep(4,18))


#Rs<-1/(a_o+b1*ys+(b2((Wp/Ws)^B)*yp))

###this should be my eqution
#Rc<-1/A+Bdi_c*Di_c+(Bdi_h*((Gt_h/Gt_c)^Beta)*Di_h)

 ## parementers based on Table 5 in Connolly and Wayne 1996
Bdi_c<--4.4
Bdi_h<--11.6
A<-4.13
sigma<-1
#Beta<- 0.51 
hist(rnorm(1000,0,3))





repz<-1:5
df<-data.frame(C_density=numeric(),H_density=numeric(),Rc=numeric(),ID=numeric())
for(k in c(1:length(repz))){ 
  y <- c()
for (i in c(1:length(Di_c))){
  y<-c()
  for (j in c(1:length(Di_h))){
y<-A+Bdi_c*Di_c[i]+Bdi_h*Di_h[j]
dfhere <- data.frame(C_density=rep(Di_c[i],length(y)),H_density=rep(Di_h[i],length(y)),Rc=rnorm(length(y),y,sigma),ID=rep(repz[k],length(y)))
df <- rbind(df, dfhere)
}
}
}


data.list<-with(df,
                 list(Y=Rc,
                     Den_c=C_density,
                      Den_e=H_density,
                     #T_e=Gt_h,
                     #T_c=Gt_c,
                      N=nrow(df) # datalist
                 )
)

#brm(Rc~C_density+H_density,data=df)
#modq<-lm(Rc~(C_density+H_density),data=df)
summary(modq)

respmod1 = stan('responsemod_fake.stan', data = data.list, 
                          iter = 3000, warmup=2000, chain=4)

summary(respmod1)$summary[c("alpha","beta_c","beta_e","sigma"),]
