# Biomarker Regression Analysis
#
# Log transformed analysis of biomarker data from Study A comparing various doses of study drug against placebo
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

## Load and filter data
lab <- read.csv("labs.csv")

studya <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studya <- studya %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- studya %>% filter (LBRUCD=="95")
trt_merge <- trt_merge %>% filter(VISID =="1" |VISID =="10")
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]

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

## Run linear regression
trt_merge$TRT<-as.factor(trt_merge$TRT)
trt_merge$TRT = relevel(trt_merge$TRT, ref="Placebo")
lm_data <- trt_merge
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,13]),]
fit <- lm(Change~TRT, data=lm_data)

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:3,1:2])

ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:3,1:2])

ci_merge <- merge(trt_diff, ci1,all=TRUE)
asdf3 <- trt_merge%>% filter(VISID=="10")
asdf3$Change[asdf3$Change==0]<- NA
asdf3 <- asdf3[complete.cases(asdf3),]

# Summary
hold2 <- asdf3 %>% group_by(TRT, VISID) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                           geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
change_sum <- hold2
change_sum$VISID ="Change"

## Add and summarize percent change with CI
percent_change_sum2<- asdf3 %>% group_by(TRT, VISID) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                        geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )
percent_change_sum2$VISID <- "Percent_Change"

check <- merge(visit_sum, change_sum, all=TRUE)
check <- merge(check, percent_change_sum2, all=TRUE)

diff <-check %>% filter(VISID=="Percent_Change")
diff$delta <- diff$geomean-diff$geomean[diff$TRT=="Placebo"]
diff$deltaSE <- sqrt((diff$geoSE[diff$TRT=="Placebo"])^2+(diff$geoSE)^2)
final_diff<-diff[,c(1,2,8,9)]
names(final_diff)<-c("TRT","VISID","geomean","geoSE")
final_diff$Lower <- final_diff$geomean-qnorm(0.975)*final_diff$geoSE
final_diff$Upper <- final_diff$geomean+qnorm(0.975)*final_diff$geoSE
final_diff<-merge(final_diff, ci_merge,all=TRUE)
final_diff<-final_diff[order(final_diff$TRT),]

check <- merge(check, final_diff, all=TRUE)

## Reformat to specifications
check$VISID[check$VISID=="1"] <- "Baseline"
check$VISID[check$VISID=="10"] <- "Week 26"
check <- check[c(13,14,15,16,1,2,3,4,9,10,5,6,7,8,11,12),]
check <- check[,c(1,2,5,4,3,6,7,8,9)]

names(check) <- c("Treatment","Time Point","N","Geometric Mean","SE for Geometric Mean","Mean of log(biomaker)","SD of log(biomarker)","Lower","Upper")



