# UACR (Urine Albumin/Creatinine (g/kg))
#
# Log transformed analysis of combined UACR data from studies GBDA, GBCF, GBDI, GBDG comparing various combined doses of study drug (and competitor drugs) against placebo
# Separates combined data in 4 different summary sets based on baseline UACR value.
# Outputs csv file that includes all data split by thresholds where results of interest include %Change vs Placebo and Change in log(UACR) vs Placebo
# Also outputs separate analysis and summary files for each individual threshold in format similar to individual study analysis output (GBCF).
#
# Fit to model Change in log(UACR) = Treatment and Fit to model log(y) = log(y_b) + Treatment
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
trt_merge <- trt_merge %>% filter(VISID =="10")
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
gbdg <- gbdg %>% filter(AVISITN =="9")
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
gbdi <- gbdi %>% filter(AVISITN =="204")
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
combined <- combined %>% filter(TRT =="Dula_1.5" | TRT =="Dula_0.75" | TRT =="Placebo")
#combined$Threshold <- ifelse(combined$base_unchanged>30,  "Above 30 Baseline", "Below 30 Baseline")

ten <- combined %>% filter (base_unchanged<10)
thirty <- combined %>% filter (base_unchanged>=10 & base_unchanged<30)
threehundred <- combined %>% filter (base_unchanged>=30 & base_unchanged<300)
overhundred <- combined %>% filter (base_unchanged>=300)

ten$Split <- "<10"
thirty$Split<- ">=10 and <30"
threehundred$Split<-  ">=30 and <300"
overhundred$Split<-  ">=300"
num_patients <- merge(ten,thirty,all=TRUE)
num_patients <- merge(num_patients,threehundred,all=TRUE)
num_patients <- merge(num_patients,overhundred,all=TRUE)
num <- num_patients %>% group_by(TRT,Split) %>% summarise( n=n() )

#combined_above <- ten
#combined_above <- thirty
#combined_above <- threehundred
combined_above <- overhundred


hold_percentchange <- combined_above %>% filter(VISID=="9" | VISID=="10" | VISID=="204")
hold_percentchange$Change[hold_percentchange$Change==0]<- NA
hold_percentchange <- hold_percentchange[complete.cases(hold_percentchange),]
percent_change_sum<- hold_percentchange %>% group_by(TRT) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                         geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )

percent_change_sum$delta <- percent_change_sum$geomean-percent_change_sum$geomean[percent_change_sum$TRT=="Placebo"]
percent_change_sum$deltaSE <- sqrt((percent_change_sum$geoSE[percent_change_sum$TRT=="Placebo"])^2+(percent_change_sum$geoSE)^2)
percent_change_sum$VISID <- "%Change vs Placebo"
percent_change_sum<-percent_change_sum[1:2,c(1,7,5,6)]
names(percent_change_sum)<-c("TRT","VISID","geomean","geoSE")
percent_change_sum$Lower <- percent_change_sum$geomean-qnorm(0.975)*percent_change_sum$geoSE
percent_change_sum$Upper <- percent_change_sum$geomean+qnorm(0.975)*percent_change_sum$geoSE

combined_above$TRT<-as.factor(combined_above$TRT)
combined_above$TRT = relevel(combined_above$TRT, ref="Placebo")
lm_data <- combined_above
lm_data$Change[lm_data$Change==0]<- NA
lm_data <- lm_data[complete.cases(lm_data[,9]),]
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

trt_diff_merge <-merge(trt_diff,ci1,all=TRUE)
final_diff <- merge(percent_change_sum, trt_diff_merge,all=TRUE)
final_diff$Threshold <- "<10"
final_diff<- final_diff[,c(7,1,2,3,4,5,6)]
names(final_diff) <- c("Threshold","Treatment", "Comparison", "Mean", "SE","Lower","Upper")

#final_diff_10 <- final_diff
#final_diff_30 <- final_diff
#final_diff_300 <- final_diff
final_diff_max <- final_diff

final_diff_30$Threshold <- ">=10 and <30"
final_diff_300$Threshold <- ">=30 and <300"
final_diff_max$Threshold <- ">=300"

final_merge<-merge(final_diff_10,final_diff_30,all=TRUE)
final_merge<-merge(final_merge,final_diff_300,all=TRUE)
final_merge<-merge(final_merge,final_diff_max,all=TRUE)

write.csv(final_merge, "Placebo_comparison_4_thresholds.csv")



####
split_summary <-function( data, split){

  data$VISID<-"Last Visit"
  visit_sum_last <- data %>% group_by(TRT,VISID,Split) %>% summarise( mean= mean(AVAL),n=n(), sd =sd(AVAL), geomean =geoMean(aval_unchanged), geoSE= geoSD(aval_unchanged))
  visit_sum_base <- data %>% group_by(TRT,VISID,Split) %>% summarise( mean= mean(BASE),n=n(), sd =sd(BASE), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged))
  visit_sum_base$VISID <- "Baseline"
  visit_sum <- merge(visit_sum_last, visit_sum_base,all=TRUE)

  temp3 <- data
  temp3$Change[temp3$Change==0]<- NA
  temp3 <- temp3[complete.cases(temp3),]

  hold2 <- temp3 %>% group_by(TRT, VISID,Split) %>% summarise( mean= mean(log(aval_unchanged)-log(base_unchanged)), n= n(), sd =(sd(log(aval_unchanged)-log(base_unchanged))),
                                                               geoSE =(sd(log(aval_unchanged)-log(base_unchanged))/sqrt(n())))
  change_sum <- hold2
  change_sum$VISID ="Change"

  percent_change_sum2<- temp3 %>% group_by(TRT, VISID, Split) %>% summarise( geomean= (exp(mean(Percent_change3))-1)*100 ,
                                                                             geoSE =( exp(mean(Percent_change3))*100*(sd(Percent_change3)/sqrt(n()))),n=n() )
  percent_change_sum2$VISID <- "Percent_Change"

  final_diff<-final_diff_10
  if(split==30){final_diff<-final_diff_30}
  if(split==300){final_diff<-final_diff_300}
  if(split==3000){final_diff<-final_diff_max}
  names(final_diff)<-c("Split","TRT","VISID","geomean","geoSE","Lower","Upper")

  check <- merge(visit_sum, change_sum, all=TRUE)
  check <- merge(check, percent_change_sum2, all=TRUE)
  check <- merge(check, final_diff, all=TRUE)
  check <- check[,c(1,2,3,6,5,4,7,8,9,10)]
  check <- check[c(13,15,14,16,2,5,3,6,4,1,8,11,9,12,10,7),]
  output<<-check
}
split_summary(ten,0)
write.csv(output, "Placebo_combined_<10_summary.csv")
split_summary(thirty,30)
write.csv(output, "Placebo_combined_<30_>10_summary.csv")
split_summary(threehundred,300)
write.csv(output, "Placebo_combined_<300_>30_summary.csv")
split_summary(overhundred,3000)
write.csv(output, "Placebo_combined_>300_summary.csv")


