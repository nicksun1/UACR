# Biomarker Regression Analysis
#
# Log transformed analysis of combined biomarker data from studies A, CF, I, G comparing various combined doses of study drug (and competitor drugs) against placebo
# Outputs csv file including Treatment, Number of Patients, Mean of log(biomarker), SE of log(biomarker), and 95% confidence intervals
# Results of interest include %Change vs Placebo and Change in log(biomarker) vs Placebo for separate doses of study drug as well as competitor drugs.
#
# Fit to model Change in log(biomarker) = Treatment
##


library(dplyr)
library(zoo)
library(EnvStats)
library(reshape)
suppressMessages(library(coastr))


## studycf

lab <- read.csv("labscf.csv")
studycf <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studycf <- studycf %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- studycf %>% filter (LBRUCD=="95")
trt_merge <- trt_merge  %>% filter( VISID =="9") ###################
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]

trt_merge$TRT[trt_merge$TRT=="Placebo/Sitagliptin"] <-"Placebo"
trt_merge$TRT[trt_merge$TRT=="LY 1.5mg"] <-"Drug_1.5"
trt_merge$TRT[trt_merge$TRT=="LY 0.75mg"] <-"Drug_0.75"
trt_merge <- trt_merge %>% filter(TRT=="Placebo" | TRT=="Drug_1.5" | TRT=="Drug_0.75"|TRT=="Sitagliptin")

trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[complete.cases(trt_merge[,14]),]
trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
studycf_comb <- trt_merge
studycf_comb$SUBJID <- paste("studycf",studycf_comb$SUBJID, sep="_")


## studya

lab <- read.csv("labsa.csv")
studya <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studya <- studya %>% filter(LBTESTABR =="MAL/CR")
#bda <- studya %>% filter(SAFFL=="Y")
trt_merge <- studya %>% filter (LBRUCD=="95")
trt_merge <- trt_merge %>% filter(VISID =="1" |VISID =="10")
#adsl_trt <- subjinfo %>% select(USUBJID ,TRT, TRTSORT)

#trt_merge <- merge(studya, adsl_trt,all=TRUE)
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]

trt_merge$TRT[trt_merge$TRT=="Placebo/LY2189265 1.5 mg"] <-"Placebo"
trt_merge$TRT[trt_merge$TRT=="Placebo/LY2189265 0.75 mg"] <-"Placebo"
trt_merge$TRT[trt_merge$TRT=="LY2189265 0.75 mg"] <-"Drug_0.75"
trt_merge$TRT[trt_merge$TRT=="LY2189265 1.5 mg"] <-"Drug_1.5"
trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
studya_comb <- trt_merge
studya_comb$SUBJID <- paste("studya",studya_comb$SUBJID, sep="_")


## studyg
lab <- read.csv("labsg.csv")
adsl <- read.csv("adslg.csv")

studyg <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
studyg <- studyg %>% filter(PARAMCD =="ALBCS49C")
studyg <- studyg %>% filter(SAFFL=="Y")
studyg <- studyg %>% filter(AVISITN =="9" |AVISITN =="100")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(studyg, adsl_trt,all=TRUE)
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
names(trt_merge) <- c(colnames(studycf_comb))
#names(trt_merge)<- c("SUBJID", "VISID", "TRT", TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studyg_comb <- trt_merge


## studyi


lab <- read.csv("labsi.csv")
adsl <- read.csv("adsli.csv")
studyi <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
studyi <- studyi %>% filter(PARAMCD =="ALBCS49C")
studyi <- studyi %>% filter(SAFFL=="Y")
studyi <- studyi %>% filter(AVISITN =="100" |AVISITN =="204")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(studyi, adsl_trt,all=TRUE)
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
names(trt_merge) <- c(colnames(studycf_comb))
studyi_comb <- trt_merge

combined <- merge(studycf_comb, studya_comb, all=TRUE)
combined <- merge(combined, studyi_comb, all=TRUE)
combined <- merge(combined, studyg_comb, all=TRUE)
combined$TRT[combined$TRT=="Drug 1.5"] <-"Drug_1.5"
combined$TRT[combined$TRT=="Drug 0.75"] <-"Drug_0.75"
combined$Threshold <- ifelse(combined$base_unchanged>30,  "Above 30 Baseline", "Below 30 Baseline")

combined_above <- combined %>% filter (base_unchanged>30)
## Above 30

combined_above$TRT<-as.factor(combined_above$TRT)
combined_above$TRT = relevel(combined_above$TRT, ref="Placebo")
lm_data <- combined_above
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,9]),]
fit <- lm(Change~TRT, data=lm_data)
coefficients <-data.frame(summary(fit)$coefficients)
vsplacebo_glargine <- data.frame(coefficients[2:5,1:2])
vsplacebo_glargine$TRT <- "Drug_0.75"
vsplacebo_glargine[2,3] <-"Drug_1.5"
vsplacebo_glargine[3,3] <-"Exenatide"
vsplacebo_glargine[4,3] <-"Sitagliptin"
vsplacebo_glargine$VISID <-"%Change vs Placebo"
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT","VISID")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT,VISID) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
#vsplacebo_glargine2$TRT <- "Drug_0.75"
#vsplacebo_glargine2[2,3] <-"Drug_1.5"
#vsplacebo_glargine2[3,3] <-"Exenatide"
#vsplacebo_glargine2[4,3] <-"Sitagliptin"
#vsplacebo_glargine2$VISID <-"%Change vs Placebo"
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:5,1:2])
ci$TRT <- "Drug_0.75"
ci[2,3] <-"Drug_1.5"
ci[3,3] <-"Exenatide"
ci[4,3] <-"Sitagliptin"
ci$VISID <-"%Change vs Placebo"
names(ci) <-c("Lower", "Upper","TRT","VISID")
ci_hold <- ci   %>% group_by(TRT,VISID) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)
#ci_hold$TRT <- "Drug_0.75"
#ci_hold[2,3] <-"Drug_1.5"
#ci_hold[3,3] <-"Exenatide"
#ci_hold[4,3] <-"Sitagliptin"
#ci_hold$VISID <-"%Change vs Placebo"

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:5,1:2])
trt_diff$TRT <- "Drug_0.75"
trt_diff[2,3] <-"Drug_1.5"
trt_diff[3,3] <-"Exenatide"
trt_diff[4,3] <-"Sitagliptin"
trt_diff$VISID <-"Change in log(biomarker) vs Placebo"
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:5,1:2])
ci1$TRT <- "Drug_0.75"
ci1[2,3] <-"Drug_1.5"
ci1[3,3] <-"Exenatide"
ci1[4,3] <-"Sitagliptin"
ci1$VISID <-"Change in log(biomarker) vs Placebo"
names(ci1) <-c("Lower", "Upper","TRT","VISID")

ci_merge <- merge(ci_hold, ci1,all=TRUE)
trt_diff_merge <-merge(trt_diff,vsplacebo_glargine2,all=TRUE)
final_diff <- merge(ci_merge, trt_diff_merge,all=TRUE)
final_diff[c(1,3,5,7),3:6]<-final_diff[c(1,3,5,7),3:6]*100

final_diff$Threshold <- "Above 30 Baseline"
final_diff<- final_diff[,c(7,1,2,5,6,3,4)]
names(final_diff) <- c("Threshold","Treatment", "Comparison", "Mean", "SE","Lower","Upper")

#write.csv(final_diff, "Placebo-comparator_combined.csv")


## Below 30 Baseline
combined_below <- combined %>% filter (base_unchanged<30)

combined_below$TRT<-as.factor(combined_below$TRT)
combined_below$TRT = relevel(combined_below$TRT, ref="Placebo")
lm_data <- combined_below
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,9]),]
fit <- lm(Change~TRT, data=lm_data)
coefficients <-data.frame(summary(fit)$coefficients)
vsplacebo_glargine <- data.frame(coefficients[2:5,1:2])
vsplacebo_glargine$TRT <- "Drug_0.75"
vsplacebo_glargine[2,3] <-"Drug_1.5"
vsplacebo_glargine[3,3] <-"Exenatide"
vsplacebo_glargine[4,3] <-"Sitagliptin"
vsplacebo_glargine$VISID <-"%Change vs Placebo"
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT","VISID")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT,VISID) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
#vsplacebo_glargine2$TRT <- "Drug_0.75"
#vsplacebo_glargine2[2,3] <-"Drug_1.5"
#vsplacebo_glargine2[3,3] <-"Exenatide"
#vsplacebo_glargine2[4,3] <-"Sitagliptin"
#vsplacebo_glargine2$VISID <-"%Change vs Placebo"
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:5,1:2])
ci$TRT <- "Drug_0.75"
ci[2,3] <-"Drug_1.5"
ci[3,3] <-"Exenatide"
ci[4,3] <-"Sitagliptin"
ci$VISID <-"%Change vs Placebo"
names(ci) <-c("Lower", "Upper","TRT","VISID")
ci_hold <- ci   %>% group_by(TRT,VISID) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)
#ci_hold$TRT <- "Drug_0.75"
#ci_hold[2,3] <-"Drug_1.5"
#ci_hold[3,3] <-"Exenatide"
#ci_hold[4,3] <-"Sitagliptin"
#ci_hold$VISID <-"%Change vs Placebo"

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:5,1:2])
trt_diff$TRT <- "Drug_0.75"
trt_diff[2,3] <-"Drug_1.5"
trt_diff[3,3] <-"Exenatide"
trt_diff[4,3] <-"Sitagliptin"
trt_diff$VISID <-"Change in log(biomarker) vs Placebo"
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:5,1:2])
ci1$TRT <- "Drug_0.75"
ci1[2,3] <-"Drug_1.5"
ci1[3,3] <-"Exenatide"
ci1[4,3] <-"Sitagliptin"
ci1$VISID <-"Change in log(biomarker) vs Placebo"
names(ci1) <-c("Lower", "Upper","TRT","VISID")

ci_merge <- merge(ci_hold, ci1,all=TRUE)
trt_diff_merge <-merge(trt_diff,vsplacebo_glargine2,all=TRUE)
final_diff_below <- merge(ci_merge, trt_diff_merge,all=TRUE)
final_diff_below[c(1,3,5,7),3:6]<-final_diff_below[c(1,3,5,7),3:6]*100

final_diff_below$Threshold <- "Below 30 Baseline"
final_diff_below<- final_diff_below[,c(7,1,2,5,6,3,4)]
names(final_diff_below) <- c("Threshold","Treatment", "Comparison", "Mean", "SE","Lower","Upper")


## Combine
output_placebo <- merge(final_diff,final_diff_below,all=TRUE)

write.csv(output_placebo, "Placebo_compare_threshold30.csv")




