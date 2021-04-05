# Biomarker Regression Analysis
#
# Log transformed analysis of biomarker data from Study G comparing various doses of study drug against competitor
# Outputs csv file including Treatment, Time Point, Number of Patients, Geometric Mean, SE for Geometric Mean, Mean of log(biomarker), SD of log(biomarker), and 95% confidence intervals
# Time points of interest include Baseline, Week 26 (midpoint), Week 52 (end of study), Change and Percent Change.
#
# Fit to model log(y) = log(y_b) + Treatment
# After model is finished transform back by LSM=exp(LSM) and SE=exp(LSM)*SE, CI for percent change given by [exp(L)-1, exp(U)-1]
##

library(dplyr)
library(EnvStats)
library(reshape)
suppressMessages(library(coastr))


lab <- read.csv("labs.csv")
adsl <- read.csv("adsl.csv")


studye <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
studye <- studye %>% filter(PARAMCD =="ALBCS49C")
studye <- studye %>% filter(SAFFL=="Y")
studye <- studye %>% filter(AVISITN =="10" |AVISITN =="100")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(studye, adsl_trt,all=TRUE)
trt_merge <- trt_merge %>% filter(SAFFL=="Y")
trt_merge$aval_unchanged <- trt_merge$AVAL
trt_merge$base_unchanged <- trt_merge$BASE
trt_merge$AVAL <- log(trt_merge$AVAL)
trt_merge$BASE <- log(trt_merge$BASE)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change <- trt_merge$Change/trt_merge$BASE
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)


detach("package:plyr", unload=TRUE)
visit_sum <- trt_merge %>% group_by(TRT01A, AVISIT) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))


hold_change <- trt_merge %>% filter(AVISITN=="10")
hold_change <- hold_change %>% select(USUBJID,Change, TRT01A)
hold_change <- hold_change[complete.cases(hold_change),]
hold_percentchange <- trt_merge %>% filter(AVISITN=="10")
hold_percentchange <-hold_percentchange %>% select(USUBJID,Percent_change,TRT01A)
hold_percentchange <- hold_percentchange[complete.cases(hold_percentchange),]
hold_percentchange$Percent_change <- hold_percentchange$Percent_change*100
hold_change_long <- melt(hold_change, id=c("USUBJID","TRT01A"))
oriscale <- trt_merge %>% filter(AVISITN=="10")
  oriscale <- oriscale %>% select(USUBJID,Change_oriscale, TRT01A)
  oriscale <- oriscale[complete.cases(oriscale),]
hold_change_long$aval_unchanged <- oriscale$Change_oriscale
hold_pchange_long <- melt(hold_percentchange, id=c("USUBJID","TRT01A"))
hold_pchange_long$aval_unchanged <- hold_percentchange$Percent_change
asdf <- merge(hold_pchange_long, hold_change_long, all=TRUE)
names(asdf) <- c("USUBJID","TRT01A","AVISIT", "AVAL","aval_unchanged")


asdf2 <- merge(asdf, trt_merge,all= TRUE)
asdf3 <- trt_merge%>% filter(AVISITN=="10")
asdf3$Change[asdf3$Change==0]<- NA
asdf3 <- asdf3[complete.cases(asdf3),]

hold2 <- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                           geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
change_sum <- hold2
change_sum$AVISIT ="Change"

something<- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= (exp(mean(log(aval_unchanged)-log(base_unchanged)))-1)*100, se =(exp(mean(log(aval_unchanged)-log(base_unchanged)))-1)*100 *(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())),n=n())
#percent_change_sum <- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= (exp(mean(Percent_change3))-1)*100 , se = (exp(mean(Percent_change3))-1)*100*sd(Percent_change3),n=n() )
percent_change_sum2<- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                        geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )

#percent_change_sum3<- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= mean(Percent_change3), se = (sd(Percent_change3)/sqrt(n())),n=n() )

percent_change_sum2$AVISIT <- "Percent_Change"

check <- merge(visit_sum, change_sum, all=TRUE)
check <- merge(check, percent_change_sum2, all=TRUE)


check <- check[c(1,4,2,3,5,8,6,7),]
check <- check[,c(1,2,3,5,4,6,7)]
check[2,2] <-"Week 26"
check[6,2] <-"Week 26"

names(check) <- c("Treatment","Time Point","N","Geometric Mean","SE for Geometric Mean","Mean of log(biomarker)","SD of log(biomarker)")


write.csv(check, "studye_biomarker.csv")










