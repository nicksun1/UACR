# Biomarker Regression Analysis
#
# Log transformed analysis of biomarker data from Study C comparing various doses of study drug against Metformin
# Outputs csv file including Treatment, Time Point, Number of Patients, Geometric Mean, SE for Geometric Mean, Mean of log(biomarker), SD of log(biomarker), and 95% confidence intervals
# Time points of interest include Baseline, Week 26 (midpoint), Week 52 (end of study), Change and Percent Change.
#
# Fit to model log(y) = log(y_b) + Treatment
# After model is finished transform back by LSM=exp(LSM) and SE=exp(LSM)*SE, CI for percent change given by [exp(L)-1, exp(U)-1]
##

library(dplyr)
library(zoo)
library(EnvStats)
library(reshape)
suppressMessages(library(coastr))


lab <- read.csv("labs.csv")

studyc <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studyc <- studyc %>% filter(LBTESTABR =="MAL/CR")
#bda <- gbda %>% filter(SAFFL=="Y")
trt_merge <- studyc %>% filter(VISID =="1" |VISID =="12"|VISID =="8"|VISID =="801")
#adsl_trt <- subjinfo %>% select(SUBJID ,TRT, DURDIABULNM)

#trt_merge <- merge(studyc, adsl_trt,all=TRUE)
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]
trt_merge <- trt_merge %>% filter (LBRUCD=="95")
trt_merge <- na.locf(trt_merge)
#trt_merge$LBBLVALTR[is.na(trt_merge$LBBLVALTR)] <- trt_merge$LBRN

#trt_merge <- trt_merge %>% filter(SAFFL=="Y")
trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change <- trt_merge$Change/trt_merge$LBBLVALTR
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)


detach("package:plyr", unload=TRUE)
visit_sum <- trt_merge %>% group_by(TRT, VISID) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))

asdf3 <- trt_merge%>% filter(VISID=="801")
asdf3$Change[asdf3$Change==0]<- NA
asdf3 <- asdf3[complete.cases(asdf3),]

#hold1 <- asdf2 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =log(geoMean(aval_unchanged)), geoSE= log(geoSD(aval_unchanged)))
hold2 <- asdf3 %>% group_by(TRT, VISID) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                       geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
hold3 <- asdf3 %>% group_by(TRT, VISID) %>% summarise( mean= mean(AVAL-BASE), n= n(), sd =(sd(AVAL-BASE)),
                                                       geoSE =(sd(AVAL-BASE)/sqrt(n())))
change_sum <- hold2
change_sum$VISID ="Change"


#something<- asdf3 %>% group_by(TRT, AVISIT) %>% summarise( mean= (exp(mean(log(aval_unchanged)-log(base_unchanged)))-1)*100, se =(exp(mean(log(aval_unchanged)-log(base_unchanged)))-1)*100 *(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())),n=n(),x=mean(Percent_change3))
#percent_change_sum <- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= (exp(mean(Percent_change3))-1)*100 , se = (exp(mean(Percent_change3))-1)*100*sd(Percent_change3),n=n() )
percent_change_sum2<- asdf3 %>% group_by(TRT, VISID) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                    geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )


#percent_change_sum3<- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= mean(Percent_change3), se = (sd(Percent_change3)/sqrt(n())),n=n() )

percent_change_sum2$VISID <- "Percent_Change"

check <- merge(visit_sum, change_sum, all=TRUE)
check <- merge(check, percent_change_sum2, all=TRUE)

diff <-check %>% filter(VISID=="Percent_Change")
diff$delta <- diff$geomean-diff$geomean[diff$TRT=="Insulin Glargine"]
diff$deltaSE <- sqrt((diff$geoSE[diff$TRT=="Insulin Glargine"])^2+(diff$geoSE)^2)
final_diff<-diff[1:2,c(1,2,8,9)]
final_diff$VISID <-"%Change vs Glargine"
final_diff$TRT <- as.character(final_diff$TRT)
final_diff$TRT[final_diff$TRT=="Drug 0.75"]<-"Drug_0.75"
final_diff$TRT[final_diff$TRT=="Drug 1.5"]<-"Drug_1.5"
names(final_diff)<-c("TRT","VISID","geomean","geoSE")
final_diff$Lower <- final_diff$geomean-qnorm(0.975)*final_diff$geoSE
final_diff$Upper <- final_diff$geomean+qnorm(0.975)*final_diff$geoSE
final_diff<-merge(final_diff, ci_merge,all=TRUE)
final_diff<-final_diff[order(final_diff$TRT),]
check <- merge(check, final_diff, all=TRUE)

#check <- check[c(1,4,2,3,5,8,6,7,9,12,10,11,13,14,16,17,15),]
check <- check[,c(1,2,3,5,4,6,7)]
check$VISID[check$VISID=="1"] <- "Baseline"
check$VISID[check$VISID=="801"] <- "Post-Baseline"
check$VISID[check$VISID=="8"] <- "Week 52"


#check[2,2] <-"Week 26"
#check[6,2] <-"Week 26"

names(check) <- c("Treatment","Time Point","N","Geometric Mean","SE for Geometric Mean","Mean of log(biomarker)","SD of log(biomarker)")


