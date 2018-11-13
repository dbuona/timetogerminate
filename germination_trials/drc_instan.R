library(rstan)
library(standrc)
data(spinach)
head(spinach)

m <- standrm(SLOPE ~ DOSE, data=spinach, 
             fct=LL.5(fixed=c(NA, NA, NA, NA, 0)), 
             curveid=b + c + e ~ HERBICIDE, 
             random=c + e ~ CURVE, 
             iter=3000)
print(m)
ranef(m)
VarCorr(m)

ED(m, respLev=c(25, 50, 75))


###other drc models
#####this code will ultimately be moved to other scripts when modeling begins in earnest
unique(master.dat$Taxa)
goo<-filter(master.dat, Taxa %in% c("Oenethera biennis","Cryptotaenia canadensis", "Hesperis matronalis","Polygonum virginiatum", "Asclepias syriaca","Carex grayi","Silene stellata","Silene vulgaris","Eurbia diviricata","Thalictrum dioicum", "Anemone virginana"))

#Ft
d/(1+exp(b(logt-logt50)))





model <- brm(tru.daily~ start+COLD+INC +(1|Taxa/plate_num), data = goo, 
                 family =(),iter=2000) 
?brm()
### All the packages below kind of suck, we are going to need to code this shit in a nonlinear model, but first lets beild up aboe
library(lme4)
?glmer()

library(drc)

Hi<-filter(goo, INC=="H")
Lo<-filter(goo,INC=="L")
specieslist<-unique(Lo$Taxa)
listhere <- list()
TTG<-list()

for (sp in seq_along(specieslist)){
  dataonesp <- subset(goo, Taxa==specieslist[sp])
  model.1 <- drm(tru.daily ~ start + end,factor(COLD), 
                 data =dataonesp, fct = LL.3(), type = "event")
  listhere[[paste(sp, specieslist[sp])]] <- list(coef(model.1)) # adding species name and coefs for doy effect
  TTG[[paste(sp, specieslist[sp])]]<-list(ED(model.1,c(50))) 
}

View(drc::germination)

listhere
View(data.frame(TTG))



model.2 <- drm(tru.daily ~ start + end, data=goo, fct = LL.2(), type = "event")

summary(model.2)
gooby<-as.data.frame(ED(model.2,50))

summary(germLL.2)
?LL.2()
ED(germLL.2,0)

## Plotting the fitted curves and the data
## plot(germLL.3, ylim = c(0, 1.5), legendPos = c(2.5,1.5))






## first installing drc and drcData
devtools::install_github("DoseResponse/drcData")
#2devtools::install_github("DoseResponse/drc")
## then installing the development version of medrc
devtools::install_github("DoseResponse/medrc")
library(drcData)
library(medrc)
?broccoli
data(broccoli)

dater.no<-filter(chill.no,plate_num!=0)



?medrm()##3 The example doesnt even converge :()
bm <- medrm(tru.daily ~ start, curvid=b + c + e + d ~ COLD,
            random=b + c + e + d ~ 1|Taxa/plate_num,
            data=goob, fc=LL.5()) 

bm <- medrm(LeafLength ~ Day, curveid=b + c + e + d ~ Stress,
            random=c + e + d ~ 1|Genotype/ID,
            data=broccoli, fc=LL.5()) 

