# UACR (Urine Albumin/Creatinine (g/kg))
#
# Log transformed analysis of combined UACR data from all studies for baseline UACR values
# Summarizes combined UACR results from all studies where baseline UACR <30 into 3 outputs (placebo <30, glargine <30 and combined <30)
# Separates combined data in 4 different summary sets based on baseline UACR value (<10, 10-30, 30-300, >300).
# Resulting csv files include Treatment, Time Point, Number of Patients, Geometric Mean, SE for Geometric Mean, Mean of log(UACR), SD of log(UACR), Max and Min
#
# Additionally outputs regression plots in pdf file for placebo <30 and glargine <30 based on LSM.
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
gbcf_comb$Study <- "GBCF"


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
gbda_comb$Study <- "GBDA"


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
gbdg_comb$Study <- "GBDG"


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
gbdi_comb$Study <- "GBDI"


## GBDD

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdd/final/data/analysis/shared",data_file="labs.sas7bdat")

gbdd <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
gbdd <- gbdd %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- gbdd %>% filter( VISID =="13") ###################
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
gbdd_comb <- trt_merge
gbdd_comb$SUBJID <- paste("GBDD",gbdd_comb$SUBJID, sep="_")
gbdd_comb$Study <- "GBDD"


## GBDB

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdb/final/data/analysis/shared",data_file="labs.sas7bdat")
#subjinfo <- import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdx/final/data/analysis/shared/adam",data_file="adsl.sas7bdat")


gbdb <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
gbdb <- gbdb %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- gbdb %>% filter( VISID =="16") ###################
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
gbdb_comb <- trt_merge
gbdb_comb$SUBJID <- paste("GBDB",gbdb_comb$SUBJID, sep="_")
gbdb_comb$Study <- "GBDB"


## GBDX

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdx/final/data/analysis/shared/adam",data_file="adlbcn.sas7bdat")
adsl <- import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdx/final/data/analysis/shared/adam",data_file="adsl.sas7bdat")

gbdx <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
gbdx <- gbdx %>% filter(PARAMCD =="ALBCS49C")
gbdx <- gbdx %>% filter(SAFFL=="Y")
gbdx <- gbdx %>% filter(AVISITN =="0" |AVISITN =="25")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(gbdx, adsl_trt,all=TRUE)
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
names(trt_merge) <- c(colnames(gbdd_comb))
gbdx_comb <- trt_merge
gbdx_comb$Study <- "GBDX"


## GBDC

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbdc/final/data/analysis/shared",data_file="labs.sas7bdat")
gbdc <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
gbdc <- gbdc %>% filter(LBTESTABR =="MAL/CR")
#bda <- gbda %>% filter(SAFFL=="Y")
trt_merge <- gbdc %>% filter(VISID =="1" |VISID =="12"|VISID =="8"|VISID =="801")
#adsl_trt <- subjinfo %>% select(SUBJID ,TRT, DURDIABULNM)

#trt_merge <- merge(gbdc, adsl_trt,all=TRUE)
trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]
trt_merge <- trt_merge %>% filter (LBRUCD=="95")
trt_merge <- na.locf(trt_merge)
#trt_merge$LBBLVALTR[is.na(trt_merge$LBBLVALTR)] <- trt_merge$LBRN

#trt_merge <- trt_merge %>% filter(SAFFL=="Y")
trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)

trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
gbdc_comb <- trt_merge
gbdc_comb$SUBJID <- paste("GBDC",gbdc_comb$SUBJID, sep="_")
gbdc_comb$Study <- "GBDC"


## GBDE

lab <-import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbde/final/data/analysis/shared/adam",data_file="adlbcn.sas7bdat")
adsl <- import_cluwe_data(source_path="/lillyce/prd/ly2189265/h9x_mc_gbde/final/data/analysis/shared/adam",data_file="adsl.sas7bdat")

gbde <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
gbde <- gbde %>% filter(PARAMCD =="ALBCS49C")
gbde <- gbde %>% filter(SAFFL=="Y")
gbde <- gbde %>% filter(AVISITN =="10" |AVISITN =="100")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(gbde, adsl_trt,all=TRUE)
trt_merge <- trt_merge %>% filter(SAFFL=="Y")
trt_merge$aval_unchanged <- trt_merge$AVAL
trt_merge$base_unchanged <- trt_merge$BASE
trt_merge$AVAL <- log(trt_merge$AVAL)
trt_merge$BASE <- log(trt_merge$BASE)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change <- trt_merge$Change/trt_merge$BASE
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)

trt_merge <- trt_merge[,c(-2,-5,-10,-8)]
trt_merge <- trt_merge[,c(1,2,6,3,7,8,4,5,9,10,11)]
names(trt_merge) <- c(colnames(gbdd_comb))
gbde_comb <- trt_merge
gbde_comb$Study <- "GBDE"

gbda_comb <- gbda_comb[complete.cases(gbda_comb),]
gbde_comb <- gbde_comb[complete.cases(gbde_comb),]
names(gbda_comb) <- c(colnames(gbdd_comb))
names(gbdb_comb) <- c(colnames(gbdd_comb))
names(gbdg_comb) <- c(colnames(gbdd_comb))
names(gbdi_comb) <- c(colnames(gbdd_comb))
names(gbdx_comb) <- c(colnames(gbdd_comb))
names(gbdc_comb) <- c(colnames(gbdd_comb))
names(gbde_comb) <- c(colnames(gbdd_comb))

all <- merge(gbcf_comb, gbda_comb, all=TRUE)
all <- merge(all, gbdb_comb, all=TRUE)
all <- merge(all, gbdd_comb, all=TRUE)
all <- merge(all, gbdg_comb, all=TRUE)
all <- merge(all, gbdi_comb, all=TRUE)
all <- merge(all, gbdx_comb, all=TRUE)
all <- merge(all, gbdc_comb, all=TRUE)
all <- merge(all, gbde_comb, all=TRUE)

all <- all%>% mutate(TRT=replace(TRT,TRT=="Dula 1.5","Dula_1.5"))
all <- all%>% mutate(TRT=replace(TRT,TRT=="Dula 0.75","Dula_0.75"))
all <- all%>% mutate(TRT=replace(TRT,TRT=="Insulin Glargine","Glargine"))

all <- all%>% filter(TRT=="Dula_1.5" |TRT=="Dula_0.75"|TRT=="Placebo"|TRT=="Glargine")
placebo_studies <- all %>% filter(Study=="GBCF"|Study=="GBDA"|Study=="GBDG"|Study=="GBDI")
glargine_studies <- all %>% filter(Study=="GBDX"|Study=="GBDD"|Study=="GBDB")

all_above30 <- all %>% filter(base_unchanged>=30)
all_below30 <- all %>% filter(base_unchanged<30)




## Summary

visit_sum_base <- all_above30 %>% group_by(TRT) %>% summarise( n=n(),mean= mean(base_unchanged),sd =sd(base_unchanged), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged),
                                                         median=median(base_unchanged),max = max(base_unchanged), min = min(base_unchanged))
visit_sum_base <-cbind(Threshold=">=30 for Baseline", visit_sum_base)
names(visit_sum_base) <- c("Threshold","Treatment","N","Mean","Standard Deviation","Geometric Mean","SE for Geometric Mean","Median","Maximum","Minimum")
## Below 30 baseline
visit_sum_base_below <- all_below30 %>% group_by(TRT) %>% summarise( n=n(),mean= mean(base_unchanged),sd =sd(base_unchanged), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged),
                                                               median=median(base_unchanged),max = max(base_unchanged), min = min(base_unchanged))
visit_sum_base_below <-cbind(Threshold="<30 for Baseline", visit_sum_base_below)
names(visit_sum_base_below) <- c("Threshold","Treatment","N","Mean","Standard Deviation","Geometric Mean","SE for Geometric Mean","Median","Maximum","Minimum")
threshold <-merge(visit_sum_base, visit_sum_base_below,all=TRUE)

write.csv(threshold, "Baseline Summary All_30.csv")

## Placebo Studies
all_placebo_above30 <- placebo_studies %>% filter(base_unchanged>=30)
all_placebo_below30 <- placebo_studies %>% filter(base_unchanged<30)

placebo_summary_above <- all_placebo_above30 %>% group_by(TRT) %>% summarise( n=n(),mean= mean(base_unchanged),sd =sd(base_unchanged), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged),
                                                               median=median(base_unchanged),max = max(base_unchanged), min = min(base_unchanged))
placebo_summary_above <-cbind(Threshold=">=30 for Baseline", placebo_summary_above)
names(placebo_summary_above) <- c("Threshold","Treatment","N","Mean","Standard Deviation","Geometric Mean","SE for Geometric Mean","Median","Maximum","Minimum")
## Below 30 baseline
placebo_summary_below <- all_placebo_below30 %>% group_by(TRT) %>% summarise( n=n(),mean= mean(base_unchanged),sd =sd(base_unchanged), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged),
                                                                     median=median(base_unchanged),max = max(base_unchanged), min = min(base_unchanged))
placebo_summary_below <-cbind(Threshold="<30 for Baseline", placebo_summary_below)
names(placebo_summary_below) <- c("Threshold","Treatment","N","Mean","Standard Deviation","Geometric Mean","SE for Geometric Mean","Median","Maximum","Minimum")
threshold_placebo <-merge(placebo_summary_above, placebo_summary_below,all=TRUE)

write.csv(threshold_placebo, "Baseline Summary Placebo_30.csv")



## Glargine Studies
all_glargine_above30 <- glargine_studies %>% filter(base_unchanged>=30)
all_glargine_below30 <- glargine_studies %>% filter(base_unchanged<30)

glargine_summary_above <- all_glargine_above30 %>% group_by(TRT) %>% summarise( n=n(),mean= mean(base_unchanged),sd =sd(base_unchanged), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged),
                                                                              median=median(base_unchanged),max = max(base_unchanged), min = min(base_unchanged))
glargine_summary_above <-cbind(Threshold=">=30 for Baseline", glargine_summary_above)
names(glargine_summary_above) <- c("Threshold","Treatment","N","Mean","Standard Deviation","Geometric Mean","SE for Geometric Mean","Median","Maximum","Minimum")
## Below 30 baseline
glargine_summary_below <- all_glargine_below30 %>% group_by(TRT) %>% summarise( n=n(),mean= mean(base_unchanged),sd =sd(base_unchanged), geomean =geoMean(base_unchanged), geoSE= geoSD(base_unchanged),
                                                                              median=median(base_unchanged),max = max(base_unchanged), min = min(base_unchanged))
glargine_summary_below <-cbind(Threshold="<30 for Baseline", glargine_summary_below)
names(glargine_summary_below) <- c("Threshold","Treatment","N","Mean","Standard Deviation","Geometric Mean","SE for Geometric Mean","Median","Maximum","Minimum")
threshold_glargine <-merge(glargine_summary_above, glargine_summary_below,all=TRUE)

write.csv(threshold_glargine, "Baseline Summary Glargine_30.csv")


## Regression analysis

#fit_placebo <- lm(Change~BASE + TRT, data=all)
#regressiondata_placebo <- all %>% filter(TRT=="Placebo")

pdf("UACR_regression_plots.pdf")

placebo <-ggplot(placebo_studies, aes(BASE, Change,shape=TRT, colour=TRT, fill=TRT)) +geom_smooth(method="lm", se=FALSE) + geom_vline(xintercept=log(30))
glargine<-ggplot(glargine_studies, aes(BASE, Change,shape=TRT, colour=TRT, fill=TRT)) +geom_smooth(method="lm", se=FALSE) + geom_vline(xintercept=log(30))
print(placebo)
print(glargine)
dev.off()












