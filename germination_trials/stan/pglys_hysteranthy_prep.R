rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()
setwd("~/Documents/git/proterant/input")


#if you dont want to run the model: 
load("RData/zarchival/hystmodels.RData")

#1) need to get a VCV from the phylgeny

michVCV<-vcv(mich.tree.droughtprune,data=mich.data$pro2)


Lmat <- matrix(rep(1), nrow = nrow(michVCV), ncol = ncol(michVCV))
diag(Lmat) <- 0

datalist<- with(mich.data, 
                list(y=pro2,
                     pol=pol_cent,
                     flotime=flo_cent,
                     minP=precip_cent,
                     V=michVCV,
                     Lmat=Lmat,
                     N = nrow(mich.data)
                ))

z.funct.stan.phylo<- stan('..//stan/pgls_hyst2.stan', data = datalist,
             iter = 3000, warmup=2000) 
