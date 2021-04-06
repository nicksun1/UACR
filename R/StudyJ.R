# Biomarker Regression Analysis
#
# Log transformed analysis of biomarker data from Study J
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
adsl <- read.csv("adsl.csv")


studyj <- dplyr::select(lab, USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL, VISIT,VISITNUM)
studyj <- studyj %>% filter(PARAMCD =="ALBCS49C")
studyj <- studyj %>% filter(SAFFL=="Y")

year2<- dplyr:: filter(studyj, grepl("12",AVISIT))
year2<- year2[!grepl("12",year2$AVISITN),]

studyj_new <- year2
studyj_new$VISIT <- "Visit 12 Months"
ttt <- studyj %>%filter(AVISIT=="BASELINE")
ttt$VISIT <-"Baseline"
studyj_new <- merge(studyj_new, ttt, all=TRUE)

adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(studyj_new, adsl_trt,all=TRUE)
trt_merge <- trt_merge[complete.cases(trt_merge),]
trt_merge <- trt_merge %>% filter(SAFFL=="Y")
trt_merge$aval_unchanged <- trt_merge$AVAL
trt_merge$base_unchanged <- trt_merge$BASE
trt_merge$AVAL <- log(trt_merge$AVAL)
trt_merge$BASE <- log(trt_merge$BASE)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change <- trt_merge$Change/trt_merge$BASE
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge<- trt_merge[-row(trt_merge)[trt_merge$aval_unchanged == 0],]
trt_merge<- trt_merge[-row(trt_merge)[trt_merge$base_unchanged == 0],]

detach("package:plyr", unload=TRUE)
visit_sum <- trt_merge %>% group_by(TRT01A, VISIT) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))


## Regression analysis
trt_merge$TRT01A<-as.factor(trt_merge$TRT01A)
trt_merge$TRT01A = relevel(trt_merge$TRT01A, ref="Placebo")
lm_data <- trt_merge
lm_data$Change[lm_data$Change==0]<- NA
lm_data$Change[lm_data$Change==Inf]<- NA
lm_data$Change[lm_data$Change==-Inf]<- NA
lm_data <- lm_data[complete.cases(lm_data[,13]),]

fit <- lm(Change~TRT01A, data=lm_data)
coefficients <-data.frame(summary(fit)$coefficients)
vsplacebo_glargine <- data.frame(coefficients[2,1:2])
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT01A","VISIT")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT01A,VISIT) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2,1:2])
names(ci) <-c("Lower", "Upper","TRT01A","VISIT")
ci_hold <- ci   %>% group_by(TRT01A,VISIT) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2,1:2])
names(trt_diff) <-c("geomean", "geoSE","TRT01A","VISIT")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2,1:2])
names(ci1) <-c("Lower", "Upper","TRT01A","VISIT")

ci_merge <- merge(trt_diff, ci1,all=TRUE)
names(ci_merge)<- c("TRT","VISID","geomean","geoSE","Lower","Upper")

## Summary
asdf3 <- trt_merge%>% filter(VISIT=="Visit 12 Months")
asdf3$Change[asdf3$Change==0]<- NA
asdf3 <- asdf3[complete.cases(asdf3),]
hold2 <- asdf3 %>% group_by(TRT01A, VISIT) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                           geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
change_sum <- hold2
change_sum$VISIT ="Change"
percent_change_sum2<- asdf3 %>% group_by(TRT01A, VISIT) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                        geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )

## Add Percent Change
percent_change_sum2$VISIT <- "Percent_Change"

check <- merge(visit_sum, change_sum, all=TRUE)
check <- merge(check, percent_change_sum2, all=TRUE)
names(check) <- c("TRT","VISID","n","geoSE","geomean","mean","sd")

diff <-check %>% filter(VISID=="Percent_Change")
diff$delta <- diff$geomean-diff$geomean[diff$TRT=="Placebo"]
diff$deltaSE <- sqrt((diff$geoSE[diff$TRT=="Placebo"])^2+(diff$geoSE)^2)
final_diff<-diff[,c(1,2,8,9)]
final_diff<-final_diff %>% filter(TRT!="Placebo")
final_diff$TRT <- as.character(final_diff$TRT)
names(final_diff)<-c("TRT","VISID","geomean","geoSE")
final_diff$Lower <- final_diff$geomean-qnorm(0.975)*final_diff$geoSE
final_diff$Upper <- final_diff$geomean+qnorm(0.975)*final_diff$geoSE
final_diff<-merge(final_diff, ci_merge,all=TRUE)
final_diff<-final_diff[order(final_diff$TRT),]
check <- merge(check, final_diff, all=TRUE)

check$VISID[check$VISID=="VISIT3"] <- "Month 0"
check$VISID[check$VISID=="FINAL VISIT"] <- "Month 72"
check <- check[c(7,10,8,9,2,6,5,3,1,4),]
check <- check[,c(1,2,5,4,3,6,7,8,9)]


names(check) <- c("Treatment","Time Point","N","Geometric Mean","SE for Geometric Mean","Mean of log(biomarker)","SD of log(biomarker)","Lower","Upper")


ten <- trt_merge %>% filter (base_unchanged<10)
thirty <- trt_merge %>% filter (base_unchanged>=10 & base_unchanged<30)
threehundred <- trt_merge %>% filter (base_unchanged>=30 & base_unchanged<300)
overhundred <- trt_merge %>% filter (base_unchanged>=300)

ten$Split <- "<10"
thirty$Split<- ">=10 and <30"
threehundred$Split<-  ">=30 and <300"
overhundred$Split<-  ">=300"


#### Summary Function
split_summary <-function( data, split){

  data$VISID<-"Visit 12 Months"
  visit_sum_last <- data %>% group_by(TRT01A,VISID,Split) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))
  visit_sum_base <- data %>% group_by(TRT01A,VISID,Split) %>% summarise( mean= mean(BASE),n=n(), sd =sd(BASE), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged))
  visit_sum_base$VISID <- "Baseline"
  visit_sum <- merge(visit_sum_last, visit_sum_base,all=TRUE)

  asdf3 <- data
  asdf3 <- asdf3[complete.cases(asdf3),]

  hold2 <- asdf3 %>% group_by(TRT01A, VISID,Split) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                               geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
  change_sum <- hold2
  change_sum$VISID ="Change"

  percent_change_sum2<- asdf3 %>% group_by(TRT01A, VISID, Split) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                             geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )
  percent_change_sum2$VISID <- "Percent_Change"

  trt_merge<-data
  trt_merge$TRT01A<-as.factor(trt_merge$TRT01A)
  trt_merge$TRT01A = relevel(trt_merge$TRT01A, ref="Placebo")
  lm_data <- trt_merge
  lm_data$Change[lm_data$Change==0]<- NA
  lm_data$Change[lm_data$Change==Inf]<- NA
  lm_data$Change[lm_data$Change==-Inf]<- NA
  lm_data <- lm_data[complete.cases(lm_data[,13]),]

  fit <- lm(Change~TRT01A, data=lm_data)
  coefficients <-data.frame(summary(fit)$coefficients)

  coeff_change <-data.frame(summary(fit)$coefficients)
  trt_diff <- data.frame(coeff_change[2,1:2])
  names(trt_diff) <-c("geomean", "geoSE","TRT01A","VISIT")
  ci1<-data.frame(confint(fit,level=0.95))
  ci1 <- data.frame(ci1[2,1:2])
  names(ci1) <-c("Lower", "Upper","TRT01A","VISIT")

  ci_merge <- merge(trt_diff, ci1,all=TRUE)
  names(ci_merge)<- c("TRT","VISID","geomean","geoSE","Lower","Upper")



  check <- merge(visit_sum, change_sum, all=TRUE)
  check <- merge(check, percent_change_sum2, all=TRUE)
  names(check) <- c("TRT","VISID","Split","n","geoSE","geomean","mean","sd")

  diff <-check %>% filter(VISID=="Percent_Change")
  diff$delta <- diff$geomean-diff$geomean[diff$TRT=="Placebo"]
  diff$deltaSE <- sqrt((diff$geoSE[diff$TRT=="Placebo"])^2+(diff$geoSE)^2)
  final_diff<-diff[,c(1,2,9,10)]
  final_diff<-final_diff %>% filter(TRT!="Placebo")
  final_diff$TRT <- as.character(final_diff$TRT)
  names(final_diff)<-c("TRT","VISID","geomean","geoSE")
  final_diff$Lower <- final_diff$geomean-qnorm(0.975)*final_diff$geoSE
  final_diff$Upper <- final_diff$geomean+qnorm(0.975)*final_diff$geoSE
  final_diff<-merge(final_diff, ci_merge,all=TRUE)
  final_diff<-final_diff[order(final_diff$TRT),]

  final_diff$Split <-"<10"
  if(split==30){ final_diff$Split <-">=10 and <30"}
  if(split==300){ final_diff$Split <-">=30 and <300"}
  if(split==3000){ final_diff$Split <-">=300"}
  names(final_diff)<-c("TRT","VISID","geomean","geoSE","Lower","Upper","Split")

  check <- merge(check, final_diff, all=TRUE)
  check <- check[,c(1,2,3,6,5,4,7,8,9,10)]
  check <- check[c(1,3,2,4,6,9,7,10,8,5),]
  asdf<<-check
}
split_summary(ten,0)
split_summary(thirty,30)
split_summary(threehundred,300)
split_summary(overhundred,3000)


asdf$TRT[asdf$TRT==i,1:2]





















