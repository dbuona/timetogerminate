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
