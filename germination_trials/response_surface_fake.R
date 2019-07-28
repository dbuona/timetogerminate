rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

if(length(grep("Lizzie", getwd())>0)) { 
  setwd("~/Documents/git/projects/misc/dan/timetogerminate/germination_trials") 
} else setwd("~/Documents/git/timetogerminate/germination_trials")


library(rstan)


D_c<-c(5,4,3,2,1,0,10,8,6,4,2,0,20,16,12,8,4,0,5,4,3,2,1,0,10,8,6,4,2,0,20,16,12,8,4,0) #densities I'd use
D_e<-c(0,1,2,3,4,5,0,2,4,6,8,10,0,4,8,12,16,20,0,1,2,3,4,5,0,2,4,6,8,10,0,4,8,12,16,20)#densities I'd use
T_c<-c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
T_e<-c(8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,8,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10)
#D_c<-c(10,20,18,26) ## desities from Connoley and Wayne
#D_e<-c(20,10,36,18) ## desities from Connoley and Wayne 
  
hist(rnorm(1000,8,5))

plot(D_c,D_e)
 ## parementers based on Table 5 in Connolly and Wayne 1996
b1<--11.6
b2<--4
a<-4.1
sigma<-.1
B<- 0.25 

repz<-1:100
df<-data.frame(C_density=numeric(),E_density=numeric(),T_c=numeric(),T_e=numeric(),RG_C=numeric(),ID=numeric())
for(k in c(1:length(repz))){ 
  y <- c()
for (i in c(1:length(D_c))){
#  y <- c()
#for (j in c(1:length(D_e))){ 
y<-1/(a+b1*D_c[i]+(b2*((T_c[i]/T_e[i])^B)*D_e[i]))
dfhere <- data.frame(C_density=rep(D_c[i],length(y)),E_density=rep(D_e[i],length(y)),T_c=rep(T_c[i],length(y)),T_e=rep(T_e[i],length(y)),RG_C=rnorm(length(y),y,sigma),ID=rep(repz[k],length(y)))
df <- rbind(df, dfhere)
}
}
#}

data.list<-with(df,
                 list(Y=RG_C,
                     Den_c=C_density,
                      Den_e=E_density,
                     T_e=T_e,
                     T_c=T_c,
                      N=nrow(df) # datalist
                 )
)


respmod1 = stan('responsemod_fake.stan', data = data.list, 
                          iter = 5000, warmup=4000, chain=4)

summary(respmod1)$summary[c("alpha","beta_c","beta_e","B","sigma"),]
