# Biomarker Regression Analysis
#
# Log transformed analysis of combined biomarker data from studies B, D and X comparing various combined doses of study drug against an insulin glargine comparator
# Outputs csv file including Treatment,Number of Patients,  Mean of log(biomarker), SE of log(biomarker), and 95% confidence intervals
# Results of interest include %Change vs Glargine and Change in log(biomarker) vs Glargine for separate doses of study drug.
#
# Fit to model Change in log(biomarker) = Treatment
#
#



library(dplyr)
library(EnvStats)
library(reshape)
suppressMessages(library(coastr))

## studyd

lab <- read.csv("labsd.csv")

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
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
studyd_comb <- trt_merge
studyd_comb$SUBJID <- paste("studyd",studyd_comb$SUBJID, sep="_")


## studyb

lab <- read.csv("labsb.csv")

studyb <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studyb <- studyb %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- studyb %>% filter( VISID =="16") ###################
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]
trt_merge <- trt_merge %>% filter (LBRUCD=="95")
#trt_merge$LBBLVALTR[is.na(trt_merge$LBBLVALTR)] <- trt_merge$LBRN

trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
studyb_comb <- trt_merge
studyb_comb$SUBJID <- paste("studyb",studyb_comb$SUBJID, sep="_")


## studyx

lab <- read.csv("labsx.csv")
adsl <- read.csv("adslx.csv")

studyx <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
studyx <- studyx %>% filter(PARAMCD =="ALBCS49C")
studyx <- studyx %>% filter(SAFFL=="Y")
studyx <- studyx %>% filter(AVISITN =="0" |AVISITN =="25")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(studyx, adsl_trt,all=TRUE)
trt_merge <- trt_merge[complete.cases(trt_merge),]
trt_merge <- trt_merge %>% filter(SAFFL=="Y")
trt_merge$aval_unchanged <- trt_merge$AVAL
trt_merge$base_unchanged <- trt_merge$BASE
trt_merge$AVAL <- log(trt_merge$AVAL)
trt_merge$BASE <- log(trt_merge$BASE)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)

trt_merge <- trt_merge[,c(-2,-5,-10,-8)]
trt_merge <- trt_merge[,c(1,2,6,3,7,8,4,5,9,10,11)]
names(trt_merge) <- c(colnames(studyd_comb))
studyx_comb <- trt_merge

combined <- merge(studyd_comb, studyb_comb, all=TRUE)
combined <- merge(combined, studyx_comb, all=TRUE)
combined$TRT[combined$TRT=="Drug 1.5"] <-"Drug_1.5"
combined$TRT[combined$TRT=="Drug 0.75"] <-"Drug_0.75"
combined$TRT[combined$TRT=="Insulin Glargine"] <-"Glargine"
combined$Threshold <- ifelse(combined$base_unchanged>30,  "Above 30 Baseline", "Below 30 Baseline")

combined_above <- combined %>% filter (base_unchanged>30)
## Above 30

combined_above$TRT<-as.factor(combined_above$TRT)
combined_above$TRT = relevel(combined_above$TRT, ref="Glargine")
lm_data <- combined_above
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,9]),]
fit <- lm(Change~TRT, data=lm_data)
coefficients <-data.frame(summary(fit)$coefficients)

vsplacebo_glargine <- data.frame(coefficients[2:3,1:2])
vsplacebo_glargine$TRT <- "Drug_0.75"
vsplacebo_glargine[2,3] <-"Drug_1.5"
vsplacebo_glargine$VISID <-"%Change vs Glargine"
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT","VISID")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT,VISID) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
#vsplacebo_glargine2$TRT <- "Drug_0.75"
#vsplacebo_glargine2[2,3] <-"Drug_1.5"
#vsplacebo_glargine2$VISID <-"%Change vs Glargine"
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:3,1:2])
ci$TRT <- "Drug_0.75"
ci[2,3] <-"Drug_1.5"
ci$VISID <-"%Change vs Glargine"
names(ci) <-c("Lower", "Upper","TRT","VISID")
ci_hold <- ci   %>% group_by(TRT,VISID) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)
#ci_hold$TRT <- "Drug_0.75"
#ci_hold[2,3] <-"Drug_1.5"
#ci_hold$VISID <-"%Change vs Glargine"

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:3,1:2])
trt_diff$TRT <- "Drug_0.75"
trt_diff[2,3] <-"Drug_1.5"
trt_diff$VISID <-"Change in log(biomarker) vs Glargine"
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:3,1:2])
ci1$TRT <- "Drug_0.75"
ci1[2,3] <-"Drug_1.5"
ci1$VISID <-"Change in log(biomarker) vs Glargine"
names(ci1) <-c("Lower", "Upper","TRT","VISID")

ci_merge <- merge(ci_hold, ci1,all=TRUE)
trt_diff_merge <-merge(trt_diff,vsplacebo_glargine2,all=TRUE)
final_diff <- merge(ci_merge, trt_diff_merge,all=TRUE)
final_diff[c(1,3),3:6]<-final_diff[c(1,3),3:6]*100

final_diff$Threshold <- "Above 30 Baseline"
final_diff<- final_diff[,c(7,1,2,5,6,3,4)]
names(final_diff) <- c("Threshold","Treatment", "Comparison", "Mean", "SE","Lower","Upper")

#write.csv(final_diff, "Glargine-comparator_combined.csv")


## Below 30 Baseline
combined_below <- combined %>% filter (base_unchanged<30)
combined_below$TRT<-as.factor(combined_below$TRT)
combined_below$TRT = relevel(combined_below$TRT, ref="Glargine")
lm_data <- combined_below
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,9]),]
fit <- lm(Change~TRT, data=lm_data)
coefficients <-data.frame(summary(fit)$coefficients)

vsplacebo_glargine <- data.frame(coefficients[2:3,1:2])
vsplacebo_glargine$TRT <- "Drug_0.75"
vsplacebo_glargine[2,3] <-"Drug_1.5"
vsplacebo_glargine$VISID <-"%Change vs Glargine"
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT","VISID")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT,VISID) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
#vsplacebo_glargine2$TRT <- "Drug_0.75"
#vsplacebo_glargine2[2,3] <-"Drug_1.5"
#vsplacebo_glargine2$VISID <-"%Change vs Glargine"
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:3,1:2])
ci$TRT <- "Drug_0.75"
ci[2,3] <-"Drug_1.5"
ci$VISID <-"%Change vs Glargine"
names(ci) <-c("Lower", "Upper","TRT","VISID")
ci_hold <- ci   %>% group_by(TRT,VISID) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)
#ci_hold$TRT <- "Drug_0.75"
#ci_hold[2,3] <-"Drug_1.5"
#ci_hold$VISID <-"%Change vs Glargine"

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:3,1:2])
trt_diff$TRT <- "Drug_0.75"
trt_diff[2,3] <-"Drug_1.5"
trt_diff$VISID <-"Change in log(biomarker) vs Glargine"
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:3,1:2])
ci1$TRT <- "Drug_0.75"
ci1[2,3] <-"Drug_1.5"
ci1$VISID <-"Change in log(biomarker) vs Glargine"
names(ci1) <-c("Lower", "Upper","TRT","VISID")

ci_merge <- merge(ci_hold, ci1,all=TRUE)
trt_diff_merge <-merge(trt_diff,vsplacebo_glargine2,all=TRUE)
final_diff_below <- merge(ci_merge, trt_diff_merge,all=TRUE)
final_diff_below[c(1,3),3:6]<-final_diff_below[c(1,3),3:6]*100

final_diff_below$Threshold <- "Below 30 Baseline"
final_diff_below<- final_diff_below[,c(7,1,2,5,6,3,4)]
names(final_diff_below) <- c("Threshold","Treatment", "Comparison", "Mean", "SE","Lower","Upper")

## Combine
output <- merge(final_diff,final_diff_below,all=TRUE)

write.csv(output, "Glargine_compare_threshold30.csv")
















