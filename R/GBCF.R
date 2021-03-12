# UACR (Urine Albumin/Creatinine (g/kg))
#
# Log transformed analysis of UACR data from study GBCF comparing various doses of study drug against placebo.
# Outputs csv file including Treatment, Time Point, Number of Patients, Geometric Mean, SE for Geometric Mean, Mean of log(UACR), SD of log(UACR), and 95% confidence intervals
# Time points of interest include Baseline, Week 26 (midpoint), Week 52 (end of study), Change and Percent Change.
#
# Fit to model log(y) = log(y_b) + Treatment
# After model is finished transform back by LSM=exp(LSM) and SE=exp(LSM)*SE, CI for percent change given by [exp(L)-1, exp(U)-1]
#
# *Library coastr is an internal package with the sole purpose of increasing the efficiency and speed of pulling data from the internal server
#


library(dplyr)
library(zoo)
library(EnvStats)
library(reshape)
suppressMessages(library(coastr))


lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbcf/final/data/analysis/shared",data_file="labs.sas7bdat")
#subjinfo <- import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbcf/final/data/analysis/shared",data_file="subjinfo.sas7bdat")

gbcf <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
gbcf <- gbcf %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- gbcf %>% filter (LBRUCD=="95")
trt_merge <- trt_merge  %>% filter( VISID =="9") ###################
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]

trt_merge$TRT[trt_merge$TRT=="Placebo/Sitagliptin"] <-"Placebo"
trt_merge$TRT[trt_merge$TRT=="LY 1.5mg"] <-"Dula 1.5"
trt_merge$TRT[trt_merge$TRT=="LY 0.75mg"] <-"Dula 0.75"
trt_merge <- trt_merge %>% filter(TRT=="Placebo" | TRT=="Dula 1.5" | TRT=="Dula 0.75"|TRT=="Sitagliptin")

trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change <- trt_merge$Change/trt_merge$LBBLVALTR
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[complete.cases(trt_merge[,14]),]
trt_merge <- trt_merge %>% filter(TRT!="Sitagliptin")

detach("package:plyr", unload=TRUE)
visit_sum_last <- trt_merge %>% group_by(TRT,VISID) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))
visit_sum_last$VISID <- "Week 26"
visit_sum_base <- trt_merge %>% group_by(TRT,VISID) %>% summarise( mean= mean(BASE),n=n(), sd =sd(BASE), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged))
visit_sum_base$VISID <- "Baseline"
visit_sum <- merge(visit_sum_last, visit_sum_base,all=TRUE)

trt_merge$TRT<-as.factor(trt_merge$TRT)
trt_merge$TRT = relevel(trt_merge$TRT, ref="Placebo")
lm_data <- trt_merge
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,13]),]
fit <- lm(Change~TRT, data=lm_data)

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:3,1:2])
trt_diff$TRT <- "Dula_0.75"
trt_diff[2,3] <-"Dula_1.5"
trt_diff$VISID <-"Change in log(UACR) vs Placebo"
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:3,1:2])
ci1$TRT <- "Dula_0.75"
ci1[2,3] <-"Dula_1.5"
ci1$VISID <-"Change in log(UACR) vs Placebo"
names(ci1) <-c("Lower", "Upper","TRT","VISID")

ci_merge <- merge(trt_diff, ci1,all=TRUE)
temp3 <- trt_merge%>% filter(VISID=="9")
temp3$Change[temp3$Change==0]<- NA
temp3 <- temp3[complete.cases(temp3),]

hold2 <- temp3 %>% group_by(TRT, VISID) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                       geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
change_sum <- hold2
change_sum$VISID ="Change"

percent_change_sum2<- temp3 %>% group_by(TRT, VISID) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                    geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )
percent_change_sum2$VISID <- "Percent_Change"


check <- merge(visit_sum, change_sum, all=TRUE)
check <- merge(check, percent_change_sum2, all=TRUE)

diff <-check %>% filter(VISID=="Percent_Change")
diff$delta <- diff$geomean-diff$geomean[diff$TRT=="Placebo"]
diff$deltaSE <- sqrt((diff$geoSE[diff$TRT=="Placebo"])^2+(diff$geoSE)^2)
final_diff<-diff[1:2,c(1,2,8,9)]
final_diff$VISID <-"%Change vs Placebo"
final_diff$TRT <- as.character(final_diff$TRT)
final_diff$TRT[final_diff$TRT=="Dula 0.75"]<-"Dula_0.75"
final_diff$TRT[final_diff$TRT=="Dula 1.5"]<-"Dula_1.5"
names(final_diff)<-c("TRT","VISID","geomean","geoSE")
final_diff$Lower <- final_diff$geomean-qnorm(0.975)*final_diff$geoSE
final_diff$Upper <- final_diff$geomean+qnorm(0.975)*final_diff$geoSE
final_diff<-merge(final_diff, ci_merge,all=TRUE)
final_diff<-final_diff[order(final_diff$TRT),]

check <- merge(check, final_diff, all=TRUE)

check <- check[c(13,16,14,15,4,2,3,10,9,5,8,6,7,12,11),]
check <- check[,c(1,2,5,4,3,6,7,8,9)]

names(check) <- c("Treatment","Time Point","N","Geometric Mean","SE for Geometric Mean","Mean of log(UACR)","SD of log(UACR)")
write.csv(check, "GBCF_UACR.csv")


