master.dat<-filter(master.dat, Taxa!="Impatiens capensis")
master.dat<-filter(master.dat, Taxa!="Phlox cuspidata")
master.dat<-filter(master.dat, Taxa!="Carex grisea")
master.dat<-filter(master.dat, Taxa!="Oenethera biennis")
master.dat<-filter(master.dat,INC=="L")
class(master.dat$Taxa)
master.dat<-merge()

sp.list <- unique(master.dat$indentifyer)


listhere<-list()

for(sp in 1:length(sp.list)){
 
    
    subset.dat <- master.dat[master.dat$indentifyer == sp.list[sp],]
    mod <- drm(tru.daily~start+end, data = subset.dat, fct = LL.2(), type = "event")
    listhere[[paste(sp, specieslist[sp])]]<-c(ED(mod,50))
    
    }

listhere
