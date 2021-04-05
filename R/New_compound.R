## Biomarker Regression analysis
#
# Log transformed analysis of biomarker data for new compound for new compound comparing various doses of study drug against older study drug
# Outputs csv file including Treatment, Time Point, Number of Patients, Geometric Mean, SE for Geometric Mean, Mean of log(biomarker), SD of log(biomarker), and 95% confidence intervals
# Time points of interest include Baseline, Week 26 (midpoint), Week 52 (end of study), Change and Percent Change.
# Also outputs regression plots for a visual representation of regression analysis
#
# Fit to model log(y) = log(y_b) + Treatment
# After model is finished transform back by LSM=exp(LSM) and SE=exp(LSM)*SE, CI for percent change given by [exp(L)-1, exp(U)-1]
##



library(dplyr)
library(EnvStats)
library(reshape)
suppressMessages(library(coastr))
library(ggplot2)

lab <- read.csv("labs.csv")
adsl <- read.csv("adsl.csv")


new_study <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
new_study <- new_study %>% filter(PARAMCD =="ALBCS49C")
new_study <- new_study %>% filter(SAFFL=="Y")
new_study <- new_study %>% filter(AVISITN =="1001" |AVISITN =="1016")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(new_study, adsl_trt,all=TRUE)
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


detach("package:plyr", unload=TRUE)
visit_sum <- trt_merge %>% group_by(TRT01A, AVISIT) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))


trt_merge$TRT01A<-as.factor(trt_merge$TRT01A)
trt_merge$TRT01A = relevel(trt_merge$TRT01A, ref="Placebo")
lm_data <- trt_merge
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,13]),]
fit <- lm(Change~TRT01A, data=lm_data)
coefficients <-data.frame(summary(fit)$coefficients)
vsplacebo_glargine <- data.frame(coefficients[2:6,1:2])
vsplacebo_glargine$TRT01A <- "Drug 1.5mg"
vsplacebo_glargine[2,3] <-"NewDrug 10mg"
vsplacebo_glargine[3,3] <-"NewDrug 15mg"
vsplacebo_glargine[4,3] <-"NewDrug 1mg"
vsplacebo_glargine[5,3] <-"NewDrug 5mg"
vsplacebo_glargine$AVISIT <-"%Change vs Placebo"
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT01A","AVISIT")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT01A,AVISIT) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:6,1:2])
ci$TRT01A <-  "Drug 1.5mg"
ci[2,3] <-"NewDrug 10mg"
ci[3,3] <-"NewDrug 15mg"
ci[4,3] <-"NewDrug 1mg"
ci[5,3] <-"NewDrug 5mg"
ci$AVISIT <-"%Change vs Placebo"
names(ci) <-c("Lower", "Upper","TRT01A","AVISIT")
ci_hold <- ci   %>% group_by(TRT01A,AVISIT) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:6,1:2])
trt_diff$TRT01A <- "Drug 1.5mg"
trt_diff[2,3] <-"NewDrug 10mg"
trt_diff[3,3] <-"NewDrug 15mg"
trt_diff[4,3] <-"NewDrug 1mg"
trt_diff[5,3] <-"NewDrug 5mg"
trt_diff$AVISIT <-"Change in log(biomarker) vs Placebo"
names(trt_diff) <-c("geomean", "geoSE","TRT01A","AVISIT")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:6,1:2])
ci1$TRT01A <- "Drug 1.5mg"
ci1[2,3] <-"NewDrug 10mg"
ci1[3,3] <-"NewDrug 15mg"
ci1[4,3] <-"NewDrug 1mg"
ci1[5,3] <-"NewDrug 5mg"
ci1$AVISIT <-"Change in log(biomarker) vs Placebo"
names(ci1) <-c("Lower", "Upper","TRT01A","AVISIT")

ci_merge <- merge(trt_diff, ci1,all=TRUE)
names(ci_merge)<- c("TRT","VISID","geomean","geoSE","Lower","Upper")

#trt_diff_merge <-merge(trt_diff,vsplacebo_glargine2,all=TRUE)
#final_diff <- merge(ci_merge, trt_diff_merge,all=TRUE)
#final_diff[c(1,3,5,7,9),3:6]<-final_diff[c(1,3,5,7,9),3:6]*100

#final_diff<- final_diff[,c(1,2,5,6,3,4)]


asdf3 <- trt_merge%>% filter(AVISITN=="1016")
asdf3$Change[asdf3$Change==0]<- NA
asdf3 <- asdf3[complete.cases(asdf3),]

hold2 <- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                           geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
change_sum <- hold2
change_sum$AVISIT ="Change"
percent_change_sum2<- asdf3 %>% group_by(TRT01A, AVISIT) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                        geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )

percent_change_sum2$AVISIT <- "Percent_Change"

check <- merge(visit_sum, change_sum, all=TRUE)
check <- merge(check, percent_change_sum2, all=TRUE)
names(check) <- c("TRT","VISID","n","geoSE","geomean","mean","sd")

diff <-check %>% filter(VISID=="Percent_Change")
diff$delta <- diff$geomean-diff$geomean[diff$TRT=="Placebo"]
diff$deltaSE <- sqrt((diff$geoSE[diff$TRT=="Placebo"])^2+(diff$geoSE)^2)
final_diff<-diff[,c(1,2,8,9)]
final_diff<-final_diff %>% filter(TRT!="Placebo")
final_diff$VISID <-"%Change vs Placebo"
final_diff$TRT <- as.character(final_diff$TRT)
#final_diff$TRT[final_diff$TRT=="Drug 0.75"]<-"Drug_0.75"
#final_diff$TRT[final_diff$TRT=="Drug 1.5"]<-"Drug_1.5"
names(final_diff)<-c("TRT","VISID","geomean","geoSE")
final_diff$Lower <- final_diff$geomean-qnorm(0.975)*final_diff$geoSE
final_diff$Upper <- final_diff$geomean+qnorm(0.975)*final_diff$geoSE
final_diff<-merge(final_diff, ci_merge,all=TRUE)
final_diff<-final_diff[order(final_diff$TRT),]

check <- merge(check, final_diff, all=TRUE)


check <- check[c(33,34,31,32,5,6,2,4,3,1,23,24,20,22,21,19,29,20,26,28,27,25,11,12,8,10,9,7,17,18,14,16,15,13),]
check <- check[,c(1,2,5,4,3,6,7,8,9)]
#check[2,2] <-"Week 26"
#check[6,2] <-"Week 26"

names(check) <- c("Treatment","Time Point","N","Geometric Mean","SE for Geometric Mean","Mean of log(biomarker)","SD of log(biomarker)","Lower","Upper")


write.csv(check, "new_study_biomarker.csv")


## Regression analysis

#fit_placebo <- lm(Change~BASE + TRT, data=all)
#regressiondata_placebo <- all %>% filter(TRT=="Placebo")

pdf("new_study_biomarker_regression_plots.pdf")

plot <-ggplot(trt_merge, aes(BASE, Change,shape=TRT01A, colour=TRT01A, fill=TRT01A)) +geom_smooth(method="lm", se=FALSE) +
      geom_vline(xintercept=log(30))+ geom_vline(xintercept=log(10) )
print(plot)
dev.off()









