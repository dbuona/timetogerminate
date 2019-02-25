## Started 15 February 2019 ##
## Working off fake_data_germination.R #

###fake data for germination mdodels
rm(list=ls())
options(stringsAsFactors = FALSE)
graphics.off()

# library(drc)
# library(tidyverse) # don't load this if you don't need it, overwrites a lot of useful stuff.
library(MCMCglmm) # rtnorm
library(shinystan)
library(drc)
library(rstan)
library(rstanarm)
library(bayesplot)
options(mc.cores = parallel::detectCores())

setwd("~/Desktop/meremoments/Dan")

##use lognoral to simulate the data
set.seed(613613)


###This make data folowing the log logistic function
germ<-function(t,d,b,t50){
y<- d/(1+((t/t50)^b))
return(data.frame(time=t, y=y))
  } 


## no chill, cool inc
# Not sure if needed, but an easier way to do letters...
A <- rep(1:10, each=9)
A <- as.character(A)

# reps
repz <- 10

##same parements in all functions
time <- seq(0,24,by=3)
d.low <- 0 
d.up <- 20

# code: parameter.chill(0-90).inc(0,1).

# sd values are pretty similar so just have one set
# If you think some parameters are very similar across treatments etc., ...
# then just ask the model to estimate one common parameter (this is a common choice for sigma)
# But (see below) I later noticed that these values are not used in the model ... 

d.sd <- 1
t50.sd <- 3
beta.sd <- 0.2


# An example of one way to think of working up the data, similar to how you were
    # I started this before realizing we may need to simplify down to build up
    # But I leave it here an example of a loop. 
    
# vectors of d.mu, t50, beta values
chilltreat <- c(0:8) # for ease, setting up as 8 levels of chilling ... this is the main driver of differences in this part of the code I think

# I need to better understand why we have multiple values for these, can it be modeled more simply?
# Do we expect the model to return these values? I think this is what the drc model may return, but
    # note that the crd model returns the CONTRASTS I believe, not the overall parameter estimates.
d.mus <- c(4, 5, 8, 9, 12, 13, 16, 17)
t50s <- c(20, 18, 15, 13, 11, 10, 8, 7)
betas <- seq(from=-6.5, to=-3, by=0.5)

# Checking
length(d.mus) 
length(t50s)
length(betas)

# I wasn't sure what all the apply commans were doing
    # I suggest we stick with loops so I can be more help
    # Here's an example one:
df <- data.frame(time=numeric(), y=numeric(), ID=numeric(), chilltreat=numeric())

for(i in c(1:length(d.mus))){ # i <- 1
    d.muhere <- d.mus[i]
    t50here <- t50s[i]
    betahere <- betas[i]
    chilltreathere <- chilltreat[i]
    dfhere <- data.frame(time=numeric(), y=numeric(), ID=numeric(), chilltreat=numeric())
    for(j in 1:repz){
        # not 100% sure that time works as you think here, would be good to double check.
        y <- rtnorm(1, d.muhere, d.sd, lower=d.low, upper=d.up)/
            (1+((time/rnorm(1, t50here, t50.sd))^rnorm(1, betahere, beta.sd)))
        dfrepadd <- data.frame(time=time, y=y, ID=rep(i, length(y)), chilltreat=chilltreathere)
        dfhere <- rbind(dfhere, dfrepadd)
    }
    df <- rbind(df, dfhere)
}


test.dat <- df
test.dat$y <- round(test.dat$y)
mod3<-drm(y~time,factor(chilltreat),fct=LL.3(),data=test.dat,type="continuous")
summary(mod3)
 # Need to plot the raw data, not just the drc output
plot(mod3,ylim=c(0,20),xlim=c(0,24),log="",pch=16,type="all")


    


# Easier example: Let's try to just model two levels of chilling
chilltreat <- c(0, 1) # We need to use this in the equation! Just the way you do in Stan
# Now, if we are starting simple, with just two levels of chill, we may need only one value of d, t50, beta ? Leaving it for now as I am not sure what this is doing...
d.mus <- c(4, 5, 8, 9, 12, 13, 16, 17)
t50s <- c(20, 18, 15, 13, 11, 10, 8, 7)
betas <- seq(from=-6.5, to=-3, by=0.5)
chill.d <- 0.2
chill.t50 <- 0.5
chill.beta <- 1
sigma_y <- 0.1


df <- data.frame(time=numeric(), y=numeric(), ID=numeric(), chilltreat=numeric())

for(i in c(1:length(chilltreat))){ # i <- 1
    chilltreathere <- chilltreat[i]
    for(j in c(1:length(d.mus))){ # j <- 1
    d.muhere <- d.mus[j]
    t50here <- t50s[j]
    betahere <- betas[j]
    dfhere <- data.frame(time=numeric(), y=numeric(), ID=numeric(), chilltreat=numeric())
    # This next equation needs to more closely match your stan code.
    # Here's the critcial bits of the stan model
    # y_hat[i] =(b_chill_d*chill[i]+a_d)/
    # (1+(((t[i])/(b_chill_t50*chill[i]+a_t50))^(b_chill_beta*chill[i]+a_beta)))
    # Also note: The only sigma you probably need in your Stan model is this one:
    # Y ~ normal(y_hat, sigma)
    # But you had three more below ... need to adjust .... here's my basic idea of how it should work
    for(k in 1:repz){
        y <- c()
        for (l in c(1:length(time))){ 
            yhere <- (chill.d*chilltreat[i]+d.muhere)/
                (1+((time[l]/(chill.t50*chilltreat[i]+t50here))^
                (chill.beta*chilltreat[i]*betahere)))
            y <- rbind(y, yhere)
        }
        # I add in sigma_y here...
        dfrepadd <- data.frame(time=time, y=rnorm(length(y), y, sigma_y), ID=rep(i, length(y)),
            chilltreat=chilltreathere)
        dfhere <- rbind(dfhere, dfrepadd)
    }
    df <- rbind(df, dfhere)
   }
}

test.dat <- df
test.dat$y <- round(test.dat$y)
mod3<-drm(y~time,factor(chilltreat),fct=LL.3(),data=test.dat,type="continuous")
summary(mod3)
plot(mod3,ylim=c(0,20),xlim=c(0,24),log="",pch=16)

plot(y~time, data=test.dat)


if(FALSE){
    data.list3 <- with(test.dat, 
                   list(Y=y, 
                        t = time,
                        chill=chilltreat,
                        N = nrow(test.dat)
                   )
)

germ.mod3 = stan('fakeseed_wchill_emw.stan', data = data.list3,
                 iter = 3000, warmup=2000) 
}
