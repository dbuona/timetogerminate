rm(list=ls()) 
options(stringsAsFactors = FALSE)
graphics.off()

inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-cdr/632/8/009968b2ee17c7af3f039ea16758626a" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")
dt1 <-read.csv(infile1,header=F 
               ,skip=1
               ,sep="\t"  
               , col.names=c(
                 "Year",     
                 "Month",     
                 "Date",     
                 "Plot",     
                 "Species.richness",     
                 "Heat.treatment",     
                 "Depth",     
                 "meanTemp",     
                 "maximumTemp",     
                 "maximumTemp",     
                 "rangeTemp",     
                 "sdTemp",     
                 "meanHumidity",     
                 "maximumHumidity",     
                 "minimumHumidity",     
                 "rangeHumidity",     
                 "sdHumidity"    ), check.names=TRUE)

unlink(infile1)

# Fix any interval or ratio columns mistakenly read in as nominal and nominal columns read as numeric or dates read as strings

# attempting to convert dt1$Date dateTime string to R date structure (date or POSIXct)                                
tmpDateFormat<-"%m/%d/%Y"
tmp1Date<-as.Date(dt1$Date,format=tmpDateFormat)
# Keep the new dates only if they all converted correctly
if(length(tmp1Date) == length(tmp1Date[!is.na(tmp1Date)])){dt1$Date <- tmp1Date } else {print("Date conversion failed for dt1$Date. Please inspect the data and do the date conversion yourself.")}                                                                    
rm(tmpDateFormat,tmp1Date) 
if (class(dt1$Plot)!="factor") dt1$Plot<- as.factor(dt1$Plot)
if (class(dt1$Species.richness)=="factor") dt1$Species.richness <-as.numeric(levels(dt1$Species.richness))[as.integer(dt1$Species.richness) ]               
if (class(dt1$Species.richness)=="character") dt1$Species.richness <-as.numeric(dt1$Species.richness)
if (class(dt1$Heat.treatment)!="factor") dt1$Heat.treatment<- as.factor(dt1$Heat.treatment)
if (class(dt1$Depth)!="factor") dt1$Depth<- as.factor(dt1$Depth)
if (class(dt1$meanTemp)=="factor") dt1$meanTemp <-as.numeric(levels(dt1$meanTemp))[as.integer(dt1$meanTemp) ]               
if (class(dt1$meanTemp)=="character") dt1$meanTemp <-as.numeric(dt1$meanTemp)
if (class(dt1$maximumTemp)=="factor") dt1$maximumTemp <-as.numeric(levels(dt1$maximumTemp))[as.integer(dt1$maximumTemp) ]               
if (class(dt1$maximumTemp)=="character") dt1$maximumTemp <-as.numeric(dt1$maximumTemp)
if (class(dt1$maximumTemp)=="factor") dt1$maximumTemp <-as.numeric(levels(dt1$maximumTemp))[as.integer(dt1$maximumTemp) ]               
if (class(dt1$maximumTemp)=="character") dt1$maximumTemp <-as.numeric(dt1$maximumTemp)
if (class(dt1$rangeTemp)=="factor") dt1$rangeTemp <-as.numeric(levels(dt1$rangeTemp))[as.integer(dt1$rangeTemp) ]               
if (class(dt1$rangeTemp)=="character") dt1$rangeTemp <-as.numeric(dt1$rangeTemp)
if (class(dt1$sdTemp)=="factor") dt1$sdTemp <-as.numeric(levels(dt1$sdTemp))[as.integer(dt1$sdTemp) ]               
if (class(dt1$sdTemp)=="character") dt1$sdTemp <-as.numeric(dt1$sdTemp)
if (class(dt1$meanHumidity)=="factor") dt1$meanHumidity <-as.numeric(levels(dt1$meanHumidity))[as.integer(dt1$meanHumidity) ]               
if (class(dt1$meanHumidity)=="character") dt1$meanHumidity <-as.numeric(dt1$meanHumidity)
if (class(dt1$maximumHumidity)=="factor") dt1$maximumHumidity <-as.numeric(levels(dt1$maximumHumidity))[as.integer(dt1$maximumHumidity) ]               
if (class(dt1$maximumHumidity)=="character") dt1$maximumHumidity <-as.numeric(dt1$maximumHumidity)
if (class(dt1$minimumHumidity)=="factor") dt1$minimumHumidity <-as.numeric(levels(dt1$minimumHumidity))[as.integer(dt1$minimumHumidity) ]               
if (class(dt1$minimumHumidity)=="character") dt1$minimumHumidity <-as.numeric(dt1$minimumHumidity)
if (class(dt1$rangeHumidity)=="factor") dt1$rangeHumidity <-as.numeric(levels(dt1$rangeHumidity))[as.integer(dt1$rangeHumidity) ]               
if (class(dt1$rangeHumidity)=="character") dt1$rangeHumidity <-as.numeric(dt1$rangeHumidity)
if (class(dt1$sdHumidity)=="factor") dt1$sdHumidity <-as.numeric(levels(dt1$sdHumidity))[as.integer(dt1$sdHumidity) ]               
if (class(dt1$sdHumidity)=="character") dt1$sdHumidity <-as.numeric(dt1$sdHumidity)

# Convert Missing Values to NA for non-dates



# Here is the structure of the input data frame:
str(dt1)                            


library(dplyr)
unique(dt1$Depth)
dt1<-filter(dt1,Depth==10) ### measurements 10 M below soil


