# UACR (Urine Albumin/Creatinine (g/kg))
#
# Log transformed analysis of combined UACR data from studies GBDA, GBCF, GBDI, GBDG comparing various combined doses of study drug (and competitor drugs) against placebo
# Outputs csv file including Treatment, Number of Patients, Mean of log(UACR), SE of log(UACR), and 95% confidence intervals
# Results of interest include %Change vs Placebo and Change in log(UACR) vs Placebo for separate doses of study drug as well as competitor drugs.
#
# Fit to model Change in log(UACR) = Treatment
#
# *Library coastr is an internal package with the sole purpose of increasing the efficiency and speed of pulling data from the internal server
#


library(dplyr)
library(zoo)
library(EnvStats)
library(reshape)
suppressMessages(library(coastr))


## GBCF

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbcf/final/data/analysis/shared",data_file="labs.sas7bdat")
gbcf <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
gbcf <- gbcf %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- gbcf %>% filter (LBRUCD=="95")
trt_merge <- trt_merge  %>% filter( VISID =="9") ###################
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]

trt_merge$TRT[trt_merge$TRT=="Placebo/Sitagliptin"] <-"Placebo"
trt_merge$TRT[trt_merge$TRT=="LY 1.5mg"] <-"Dula_1.5"
trt_merge$TRT[trt_merge$TRT=="LY 0.75mg"] <-"Dula_0.75"
trt_merge <- trt_merge %>% filter(TRT=="Placebo" | TRT=="Dula_1.5" | TRT=="Dula_0.75"|TRT=="Sitagliptin")

trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[complete.cases(trt_merge[,14]),]
trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
gbcf_comb <- trt_merge
gbcf_comb$SUBJID <- paste("GBCF",gbcf_comb$SUBJID, sep="_")


## GBDA

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbda/final/data/analysis/shared",data_file="labs.sas7bdat")
gbda <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
gbda <- gbda %>% filter(LBTESTABR =="MAL/CR")
#bda <- gbda %>% filter(SAFFL=="Y")
trt_merge <- gbda %>% filter (LBRUCD=="95")
trt_merge <- trt_merge %>% filter(VISID =="1" |VISID =="10")
#adsl_trt <- subjinfo %>% select(USUBJID ,TRT, TRTSORT)

#trt_merge <- merge(gbda, adsl_trt,all=TRUE)
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]

trt_merge$TRT[trt_merge$TRT=="Placebo/LY2189265 1.5 mg"] <-"Placebo"
trt_merge$TRT[trt_merge$TRT=="Placebo/LY2189265 0.75 mg"] <-"Placebo"
trt_merge$TRT[trt_merge$TRT=="LY2189265 0.75 mg"] <-"Dula_0.75"
trt_merge$TRT[trt_merge$TRT=="LY2189265 1.5 mg"] <-"Dula_1.5"
trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
gbda_comb <- trt_merge
gbda_comb$SUBJID <- paste("GBDA",gbda_comb$SUBJID, sep="_")


## GBDG

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdg/final/data/analysis/shared/adam",data_file="adlbcn.sas7bdat")
adsl <- import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdg/final/data/analysis/shared/adam",data_file="adsl.sas7bdat")

gbdg <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
gbdg <- gbdg %>% filter(PARAMCD =="ALBCS49C")
gbdg <- gbdg %>% filter(SAFFL=="Y")
gbdg <- gbdg %>% filter(AVISITN =="9" |AVISITN =="100")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(gbdg, adsl_trt,all=TRUE)
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
names(trt_merge) <- c(colnames(gbcf_comb))
#names(trt_merge)<- c("SUBJID", "VISID", "TRT", TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
gbdg_comb <- trt_merge


## GBDI


lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdi/final/data/analysis/shared/adam",data_file="adlbcn.sas7bdat")
adsl <- import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdi/final/data/analysis/shared/adam",data_file="adsl.sas7bdat")

gbdi <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
gbdi <- gbdi %>% filter(PARAMCD =="ALBCS49C")
gbdi <- gbdi %>% filter(SAFFL=="Y")
gbdi <- gbdi %>% filter(AVISITN =="100" |AVISITN =="204")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(gbdi, adsl_trt,all=TRUE)
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
names(trt_merge) <- c(colnames(gbcf_comb))
gbdi_comb <- trt_merge

combined <- merge(gbcf_comb, gbda_comb, all=TRUE)
combined <- merge(combined, gbdi_comb, all=TRUE)
combined <- merge(combined, gbdg_comb, all=TRUE)
combined$TRT[combined$TRT=="Dula 1.5"] <-"Dula_1.5"
combined$TRT[combined$TRT=="Dula 0.75"] <-"Dula_0.75"
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
vsplacebo_glargine$TRT <- "Dula_0.75"
vsplacebo_glargine[2,3] <-"Dula_1.5"
vsplacebo_glargine[3,3] <-"Exenatide"
vsplacebo_glargine[4,3] <-"Sitagliptin"
vsplacebo_glargine$VISID <-"%Change vs Placebo"
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT","VISID")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT,VISID) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
#vsplacebo_glargine2$TRT <- "Dula_0.75"
#vsplacebo_glargine2[2,3] <-"Dula_1.5"
#vsplacebo_glargine2[3,3] <-"Exenatide"
#vsplacebo_glargine2[4,3] <-"Sitagliptin"
#vsplacebo_glargine2$VISID <-"%Change vs Placebo"
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:5,1:2])
ci$TRT <- "Dula_0.75"
ci[2,3] <-"Dula_1.5"
ci[3,3] <-"Exenatide"
ci[4,3] <-"Sitagliptin"
ci$VISID <-"%Change vs Placebo"
names(ci) <-c("Lower", "Upper","TRT","VISID")
ci_hold <- ci   %>% group_by(TRT,VISID) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)
#ci_hold$TRT <- "Dula_0.75"
#ci_hold[2,3] <-"Dula_1.5"
#ci_hold[3,3] <-"Exenatide"
#ci_hold[4,3] <-"Sitagliptin"
#ci_hold$VISID <-"%Change vs Placebo"

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:5,1:2])
trt_diff$TRT <- "Dula_0.75"
trt_diff[2,3] <-"Dula_1.5"
trt_diff[3,3] <-"Exenatide"
trt_diff[4,3] <-"Sitagliptin"
trt_diff$VISID <-"Change in log(UACR) vs Placebo"
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:5,1:2])
ci1$TRT <- "Dula_0.75"
ci1[2,3] <-"Dula_1.5"
ci1[3,3] <-"Exenatide"
ci1[4,3] <-"Sitagliptin"
ci1$VISID <-"Change in log(UACR) vs Placebo"
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
vsplacebo_glargine$TRT <- "Dula_0.75"
vsplacebo_glargine[2,3] <-"Dula_1.5"
vsplacebo_glargine[3,3] <-"Exenatide"
vsplacebo_glargine[4,3] <-"Sitagliptin"
vsplacebo_glargine$VISID <-"%Change vs Placebo"
names(vsplacebo_glargine) <-c("AVAL", "SE","TRT","VISID")
vsplacebo_glargine2 <- vsplacebo_glargine %>% group_by(TRT,VISID) %>% summarise( geomean =exp(AVAL)-1, geoSE= exp(AVAL)*SE)
#vsplacebo_glargine2$TRT <- "Dula_0.75"
#vsplacebo_glargine2[2,3] <-"Dula_1.5"
#vsplacebo_glargine2[3,3] <-"Exenatide"
#vsplacebo_glargine2[4,3] <-"Sitagliptin"
#vsplacebo_glargine2$VISID <-"%Change vs Placebo"
ci<-data.frame(confint(fit,level=0.95))
ci <- data.frame(ci[2:5,1:2])
ci$TRT <- "Dula_0.75"
ci[2,3] <-"Dula_1.5"
ci[3,3] <-"Exenatide"
ci[4,3] <-"Sitagliptin"
ci$VISID <-"%Change vs Placebo"
names(ci) <-c("Lower", "Upper","TRT","VISID")
ci_hold <- ci   %>% group_by(TRT,VISID) %>% summarise( Lower =exp(Lower)-1, Upper=exp(Upper)-1)
#ci_hold$TRT <- "Dula_0.75"
#ci_hold[2,3] <-"Dula_1.5"
#ci_hold[3,3] <-"Exenatide"
#ci_hold[4,3] <-"Sitagliptin"
#ci_hold$VISID <-"%Change vs Placebo"

coeff_change <-data.frame(summary(fit)$coefficients)
trt_diff <- data.frame(coeff_change[2:5,1:2])
trt_diff$TRT <- "Dula_0.75"
trt_diff[2,3] <-"Dula_1.5"
trt_diff[3,3] <-"Exenatide"
trt_diff[4,3] <-"Sitagliptin"
trt_diff$VISID <-"Change in log(UACR) vs Placebo"
names(trt_diff) <-c("geomean", "geoSE","TRT","VISID")
ci1<-data.frame(confint(fit,level=0.95))
ci1 <- data.frame(ci1[2:5,1:2])
ci1$TRT <- "Dula_0.75"
ci1[2,3] <-"Dula_1.5"
ci1[3,3] <-"Exenatide"
ci1[4,3] <-"Sitagliptin"
ci1$VISID <-"Change in log(UACR) vs Placebo"
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




