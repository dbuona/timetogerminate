##This is the main cleaning file for germination work. Begin on november 12, 2018 by Dan B
###clean and combine data

rm(list=ls()) 
options(stringsAsFactors = FALSE)
setwd("~/Documents/git/timetogerminate/germination_trials/input")
library(tidyverse)
library(lubridate)
library("Hmisc")

chill.no<-read.csv("dat_CHILL_NO.CSV",header=TRUE)
colnames(chill.no)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","8/27/18","8/28/18","8/29/18",	"8/30/18",	"8/31/18",	"9/3/18",	"9/5/18",	"9/7/18",	"9/9/18",	"9/12/18",	"9/14/18",	"9/16/18",	"9/17/18","9/18/18","9/20/18","9/21/18","viable","chill_germ","mold","adj_total"))

chill.A<-read.csv("data_CHILL_A.csv",header=TRUE)
colnames(chill.A)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","9/10/18","9/12/18","9/14/18",	"9/16/18","9/17/18","9/18/18","9/20/18","9/21/18","9/23/18","9/26/18","9/28/18","9/30/18","10/3/18","10/5/18","viable","chill_germ","mold","adj_total"))

chill.B<-read.csv("dat_chill_B.csv",header=TRUE)
colnames(chill.B)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","9/24/18","9/26/18","9/28/18","9/30/18","10/3/18","10/5/18","10/7/18"  , "10/8/18" ,  "10/10/18","10/12/18" , "10/14/18" , "10/15/18", "10/17/18" , "10/19/18", "viable"     ,"chill_germ", "mold", "adj_total"))

chill.C<-read.csv("dat_chill_C.csv",header=TRUE)
colnames(chill.C)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/1/18","10/3/18","10/5/18","10/7/18"  , "10/8/18" ,  "10/10/18","10/12/18" , "10/14/18" , "10/15/18", "10/17/18" , "10/19/18","10/21/18",  "10/22/18" , "10/24/18",  "10/26/18","viable"     ,"chill_germ", "mold", "adj_total"))

chill.D<-read.csv("dat_chill_D.csv",header=TRUE)
colnames(chill.D)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num", "10/8/18" ,  "10/10/18","10/12/18" , "10/14/18" , "10/15/18", "10/17/18" , "10/19/18","10/21/18",  "10/22/18" , "10/24/18",  "10/26/18","10/28/18" ,"10/29/18",  "10/31/18" , "11/2/18","viable","chill_germ", "mold", "adj_total"))

chill.E<-read.csv("dat_chill_E.csv",header=TRUE)
colnames(chill.E)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/15/18", "10/17/18" , "10/19/18","10/21/18",  "10/22/18" , "10/24/18",  "10/26/18","10/28/18" ,"10/29/18",  "10/31/18" , "11/2/18","11/4/18",  "11/5/18",   "11/7/18" ,  "11/9/18","viable","chill_germ", "mold", "adj_total"))
chill.E$viable<-as.numeric(chill.E$viable)

chill.F<-read.csv("dat_chill_F.csv",header=TRUE)
colnames(chill.F)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num", "10/22/18" , "10/24/18",  "10/26/18","10/28/18" ,"10/29/18",  "10/31/18" , "11/2/18","11/4/18",  "11/5/18",   "11/7/18" ,  "11/9/18","11/11/18","11/12/18","11/14/18","11/16/18","viable","chill_germ", "mold", "adj_total"))
chill.F$COLD<-"f"

chill.G<-read.csv("dat_chill_G.csv",header=TRUE)
colnames(chill.G)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","10/29/18",  "10/31/18" , "11/2/18","11/4/18",  "11/5/18",   "11/7/18" ,  "11/9/18","11/11/18","11/12/18","11/14/18","11/16/18","11/18/18","11/19/18","11/21/18","11/23/18","viable","chill_germ", "mold", "adj_total"))

chill.H<-read.csv("dat_chill_H.csv",header=TRUE)
colnames(chill.H)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","11/12/18","11/14/18","11/16/18","11/18/18","11/19/18","11/21/18","11/23/18","11/26/18","11/28/18","11/30/18","12/2/18","12/3/18","12/5/18","12/7/18","viable","chill_germ", "mold", "adj_total"))

chill.i<-read.csv("dat_chill_i.csv",header=TRUE)
colnames(chill.i)<-(c("zero_day","tot_seed","Taxa","INC","COLD","plate_num","11/26/18","11/28/18","11/30/18","12/2/18","12/3/18","12/5/18","12/7/18","12/10/18","12/12/18","12/14/18","12/16/18","12/17/18","12/20/18","12/21/18","viable","chill_germ", "mold", "adj_total"))

##make all nas 0
na2zero<-function(x){x%>% replace(is.na(.), 0)}
chill.no<-na2zero(chill.no)
chill.A<-na2zero(chill.A)
chill.B<-na2zero(chill.B)
chill.C<-na2zero(chill.C)
chill.D<-na2zero(chill.D)
chill.E<-na2zero(chill.E)
chill.F<-na2zero(chill.F)
chill.G<-na2zero(chill.G)
chill.H<-na2zero(chill.H)
chill.i<-na2zero(chill.i)
##3adjust germination percentages

####1] total germinated/ totalviable+total germ
###not accounting for mold factor at the moment
chill.no$adj_total<-chill.no$"9/21/18"+chill.no$viable
chill.A$adj_total<-chill.A$"10/5/18"+chill.A$viable
chill.B$adj_total<-chill.B$"10/19/18"+chill.B$viable
chill.C$adj_total<-chill.C$"10/26/18"+chill.C$viable
chill.D$adj_total<-chill.D$"11/2/18"+chill.D$viable
chill.E$adj_total<-chill.E$"11/9/18"+chill.E$viable 
chill.F$adj_total<-chill.F$"11/16/18"+chill.F$viable
chill.G$adj_total<-chill.G$"11/23/18"+chill.G$viable
chill.H$adj_total<-chill.H$"12/7/18"+chill.H$viable
chill.i$adj_total<-chill.i$"12/21/18"+chill.i$viable

###This is a dummy varaible to be helpful when converting to drc format
chill.no$"1/1/18"<-chill.no$chill_germ
chill.A$"1/1/18"<-chill.A$chill_germ
chill.B$"1/1/18"<-chill.B$chill_germ
chill.C$"1/1/18"<-chill.C$chill_germ
chill.D$"1/1/18"<-chill.D$chill_germ
chill.E$"1/1/18"<-chill.E$chill_germ
chill.F$"1/1/18"<-chill.F$chill_germ
chill.G$"1/1/18"<-chill.G$chill_germ
chill.H$"1/1/18"<-chill.H$chill_germ
chill.i$"1/1/18"<-chill.i$chill_germ

colnames(chill.no)
chill.no <- chill.no[, c(1,2,3,4,5,6,27,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26)]
chill.A <- chill.A[, c(1,2,3,4,5,6,25,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)]
chill.B <- chill.B[, c(1,2,3,4,5,6,25,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)]
chill.C <- chill.C[, c(1,2,3,4,5,6,26,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
chill.D <- chill.D[, c(1,2,3,4,5,6,26,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
chill.E <- chill.E[, c(1,2,3,4,5,6,26,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
chill.F <- chill.F[, c(1,2,3,4,5,6,26,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
chill.G <- chill.G[, c(1,2,3,4,5,6,26,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25)]
chill.H <- chill.H[, c(1,2,3,4,5,6,25,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)]
chill.i <- chill.i[, c(1,2,3,4,5,6,25,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)]


##reformat data
chill.no<-gather(chill.no,"date","germ_num",7:23)
chill.no$germ_num<-as.numeric(chill.no$germ_num)
chill.A<-gather(chill.A,"date","germ_num",7:21)
chill.A$germ_num<-as.numeric(chill.A$germ_num)
chill.B<-gather(chill.B,"date","germ_num",7:21)
chill.B$germ_num<-as.numeric(chill.B$germ_num)
chill.C<-gather(chill.C,"date","germ_num",7:22)
chill.C$germ_num<-as.numeric(chill.C$germ_num)
chill.D<-gather(chill.D,"date","germ_num",7:22)
chill.D$germ_num<-as.numeric(chill.D$germ_num)
chill.E<-gather(chill.E,"date","germ_num",7:22)
chill.E$germ_num<-as.numeric(chill.E$germ_num)
chill.F<-gather(chill.F,"date","germ_num",7:22)
chill.F$germ_num<-as.numeric(chill.F$germ_num)
chill.G<-gather(chill.G,"date","germ_num",7:22)
chill.G$germ_num<-as.numeric(chill.G$germ_num)
chill.H<-gather(chill.H,"date","germ_num",7:21)
chill.H$germ_num<-as.numeric(chill.H$germ_num)
chill.i<-gather(chill.i,"date","germ_num",7:21)
chill.i$germ_num<-as.numeric(chill.i$germ_num)




####### calculate day of exp
chill.no$date<-as.Date(chill.no$date,format =  "%m/%d/%y")
chill.no$day<-yday(chill.no$date)
start<-yday("2018/08/27")
chill.no$DAY<-chill.no$day-start

chill.A$date<-as.Date(chill.A$date,format =  "%m/%d/%y")
chill.A$day<-yday(chill.A$date)
startA<-yday("2018/09/10")
chill.A$DAY<-chill.A$day-startA

chill.B$date<-as.Date(chill.B$date,format =  "%m/%d/%y")
chill.B$day<-yday(chill.B$date)
startB<-yday("2018/09/24")
chill.B$DAY<-chill.B$day-startB

chill.C$date<-as.Date(chill.C$date,format =  "%m/%d/%y")
chill.C$day<-yday(chill.C$date)
startC<-yday("2018/10/01")
chill.C$DAY<-chill.C$day-startC

chill.D$date<-as.Date(chill.D$date,format =  "%m/%d/%y")
chill.D$day<-yday(chill.D$date)
startD<-yday("2018/10/08")
chill.D$DAY<-chill.D$day-startD

chill.E$date<-as.Date(chill.E$date,format =  "%m/%d/%y")
chill.E$day<-yday(chill.E$date)
startE<-yday("2018/10/15")
chill.E$DAY<-chill.E$day-startE

chill.F$date<-as.Date(chill.F$date,format =  "%m/%d/%y")
chill.F$day<-yday(chill.F$date)
startF<-yday("2018/10/22")
chill.F$DAY<-chill.F$day-startF

chill.G$date<-as.Date(chill.G$date,format =  "%m/%d/%y")
chill.G$day<-yday(chill.G$date)
startG<-yday("2018/10/29")
chill.G$DAY<-chill.G$day-startG

chill.H$date<-as.Date(chill.H$date,format= "%m/%d/%y")
chill.H$day<-yday(chill.H$date)
startH<-yday("2018/11/12")
chill.H$DAY<-chill.H$day-startH

chill.i$date<-as.Date(chill.i$date,format= "%m/%d/%y")
chill.i$day<-yday(chill.i$date)
starti<-yday("2018/11/26")
chill.i$DAY<-chill.i$day-starti

cum_data_full<-rbind(chill.no,chill.A,chill.B,chill.C,chill.D,chill.E,chill.F,chill.G,chill.H,chill.i)

colnames(cum_data_full)
cum_data_full<-filter(cum_data_full,INC %in% c("H","L")) ##get rid of empty lines
cum_data_full<- within(cum_data_full, adj_total[Taxa=="Carex grayi"]<-20)
cum_data_full<- within(cum_data_full, adj_total[Taxa=="Impatiens capensis"]<-10)
cum_data_full<- within(cum_data_full, COLD[COLD==0]<-"O")


unique(cum_data_full$COLD)
#####
cum_data_A<-filter(cum_data_full,DAY>=0)
unique(cum_data_A$COLD)

cum_data_B<-filter(cum_data_full,DAY<0)
unique(cum_data_B$COLD)

cum_data_A$germ_num<-cum_data_A$germ_num+cum_data_A$chill_germ

cum_data_full<-bind_rows(cum_data_B,cum_data_A)
unique(cum_data_full$COLD)

cum_data_full$germ_perc<-cum_data_full$germ_num/cum_data_full$adj_total

cum_data_full$DAY<-ifelse(cum_data_full$DAY<0,-Inf,cum_data_full$DAY)

unique(cum_data_full$COLD)

write.csv(cum_data_full,"cumulative_data.csv")

chilli<-filter(cum_data_full,COLD=="i")
ggplot(chilli, aes(x = DAY, y = germ_num, shape=INC)) + stat_summary(alpha=0.7)+theme_bw()+geom_line(stat = "summary", fun.y = mean)+facet_wrap(~Taxa)
