
####multi species
library(drc)
data("chickweed")
data("germination")
germLL.2 <- drm(germinated ~ start + end, species:factor(temp), 
                data = germination[c(1:23, 25:61, 63:192), ], fct = LL.2(), type = "event")
plot(germLL.2, ylim=c(0, 1.5), legendPos=c(2.5,1.5))  # plotting the fitted curves and the data
summary(germLL.2)
b<-as.data.frame(coef(germLL.2))
b<-rownames_to_column(b,"garble")
b<-separate(b, garble, c("coef","species","treatment"))
b<-b[,c(2,3,1,4)]
colnames(b)<-c("species","treatment","coef","estimation")
b<-spread(b,coef,estimation)

a<-as.data.frame(ED(germLL.2,c(10,25,50,75)))
a<-rownames_to_column(a,"garble")
a<-separate(a, garble, c("extra","species", "treatment", "Tg"))
a$Estimate<-as.numeric(a$Estimate)
a$treatment<-as.numeric(a$treatment)
a$Tg<-as.character(a$Tg)

lmer(1/Estimate~treatment+(1|) data=a)

ggplot(a,aes(treatment,1/Estimate))+geom_jitter(aes(color=species,shape=Tg))+xlim(c(-40,40))+geom_smooth(method = "lm", se = FALSE,aes(color=species,linetype=Tg),fullrange=TRUE)+geom_hline(yintercept=0)+theme_classic()

###single
chickweed <- data.frame(start = c(0, chickweed0$time), end = c(chickweed0$time, Inf)) 
chickweed$count <- c(0, diff(chickweed0$count), 200 - tail(chickweed0$count, 1))
head(chickweed)  # showing top 6 lines of the dataset
tail(chickweed)  # showing bottom 6 lines

## Fitting the event-time model (by specifying the argument type explicitly)
chickweed.m1 <- drm(count~start+end, data = chickweed, fct = LL.3(), type = "event")
summary(chickweed.m1)  # showing a summmary of the model fit (including parameter estimates)

## Summary output with robust standard errors
## library(lmtest)
## library(sandwich)
## coeftest(chickweed.m1, vcov = sandwich)

## Calculating t10, t50, t90 for the distribution of viable seeds
ED(chickweed.m1, c(10, 50, 90))

### Do this in a loop maybe?