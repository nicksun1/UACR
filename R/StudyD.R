# Biomarker Regression Analysis
#
# Log transformed analysis of biomarker data from Study D
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


## Load and clean data
lab <- read.csv("labs.csv")

studyd <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studyd <- studyd %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- studyd %>% filter( VISID =="13") ###################
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]
trt_merge <- trt_merge %>% filter (LBRUCD=="95")

trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change <- trt_merge$Change/trt_merge$LBBLVALTR
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)

### Geo mean summary
detach("package:plyr", unload=TRUE)
visit_sum_last <- trt_merge %>% group_by(TRT,VISID) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))
visit_sum_last$VISID <- "Week 52"
visit_sum_base <- trt_merge %>% group_by(TRT,VISID) %>% summarise( mean= mean(BASE),n=n(), sd =sd(BASE), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged))
visit_sum_base$VISID <- "Baseline"
visit_sum <- merge(visit_sum_last, visit_sum_base,all=TRUE)

#### Regression analysis and CI
trt_merge$TRT<-as.factor(trt_merge$TRT)
trt_merge$TRT = relevel(trt_merge$TRT, ref="Glargine")
lm_data <- trt_merge
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,13]),]
fit <- lm(Change~TRT, data=lm_data)
coefficients <-data.frame(summary(fit)$coefficients)
vsplacebo_glargine <- data.frame(coefficients[2:3,1:2])
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT","VISID")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT,VISID) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:3,1:2])
names(ci) <-c("Lower", "Upper","TRT","VISID")
ci_hold <- ci   %>% group_by(TRT,VISID) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:3,1:2])
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:3,1:2])
names(ci1) <-c("Lower", "Upper","TRT","VISID")

## Summary
ci_merge <- merge(trt_diff, ci1,all=TRUE)
trt_diff_merge <-merge(trt_diff,vsplacebo_glargine2,all=TRUE)
final_diff <- merge(ci_merge, trt_diff_merge,all=TRUE)
final_diff[c(1,3),3:6]<-final_diff[c(1,3),3:6]*100

asdf3 <- trt_merge%>% filter(VISID=="13")
asdf3$Change[asdf3$Change==0]<- NA
asdf3 <- asdf3[complete.cases(asdf3),]
hold2 <- asdf3 %>% group_by(TRT, VISID) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                       geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
change_sum <- hold2
change_sum$VISID ="Change"

## Summarize and add Percent Change and CI
percent_change_sum2<- asdf3 %>% group_by(TRT, VISID) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                    geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )
percent_change_sum2$VISID <- "Percent_Change"

check <- merge(visit_sum, change_sum, all=TRUE)
check <- merge(check, percent_change_sum2, all=TRUE)

diff <-check %>% filter(VISID=="Percent_Change")
diff$delta <- diff$geomean-diff$geomean[diff$TRT=="Glargine"]
diff$deltaSE <- sqrt((diff$geoSE[diff$TRT=="Glargine"])^2+(diff$geoSE)^2)
final_diff<-diff[1:2,c(1,2,8,9)]
final_diff$VISID <-"%Change vs Glargine"
final_diff$TRT <- as.character(final_diff$TRT)
names(final_diff)<-c("TRT","VISID","geomean","geoSE")
final_diff$Lower <- final_diff$geomean-qnorm(0.975)*final_diff$geoSE
final_diff$Upper <- final_diff$geomean+qnorm(0.975)*final_diff$geoSE
final_diff<-merge(final_diff, ci_merge,all=TRUE)
final_diff<-final_diff[order(final_diff$TRT),]
check <- merge(check, final_diff, all=TRUE)

## Reformat
check <- check[c(13,16,14,15,1,4,2,3,10,9,5,8,6,7,12,11),]
check <- check[,c(1,2,5,4,3,6,7,8,9)]

names(check) <- c("Treatment","Time Point","N","Geometric Mean","SE for Geometric Mean","Mean of log(biomarker)","SD of log(biomarker)","Lower","Upper")



