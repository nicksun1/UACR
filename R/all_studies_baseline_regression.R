# Biomarker Regression Analysis
#
# Log transformed analysis of combined biomarker data from all studies for baseline biomarker values
# Summarizes combined biomarker results from all studies where baseline biomarker <30 into 3 outputs (placebo <30, glargine <30 and combined <30)
# Separates combined data in 4 different summary sets based on baseline biomarker value (<10, 10-30, 30-300, >300).
# Resulting csv files include Treatment, Time Point, Number of Patients, Geometric Mean, SE for Geometric Mean, Mean of log(biomarker), SD of log(biomarker), Max and Min
#
# Additionally outputs regression plots in pdf file for placebo <30 and glargine <30 based on LSM.
#
# Fit to model Change in log(biomarker) = Treatment and Fit to model log(y) = log(y_b) + Treatment
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
studycf_comb$Study <- "studycf"


## studya

lab <- read.csv("labsa.csv")
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
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)
trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
studya_comb <- trt_merge
studya_comb$SUBJID <- paste("studya",studya_comb$SUBJID, sep="_")
studya_comb$Study <- "studya"


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
studyg_comb <- trt_merge
studyg_comb$Study <- "studyg"


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
studyi_comb$Study <- "studyi"


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
studyd_comb$Study <- "studyd"


## studyb
lab <- read.csv("labsb.csv")

studyb <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studyb <- studyb %>% filter(LBTESTABR =="MAL/CR")
trt_merge <- studyb %>% filter( VISID =="16") ###################
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
studyb_comb <- trt_merge
studyb_comb$SUBJID <- paste("studyb",studyb_comb$SUBJID, sep="_")
studyb_comb$Study <- "studyb"


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
studyx_comb$Study <- "studyx"


## studyc

lab <- read.csv("labsc.csv")
studyc <- lab %>% select(SUBJID, VISID, TRT, TRTSORT, LBTESTABR, LBTEST, LBRN, LBBLVALTR,LBRUCD)
studyc <- studyc %>% filter(LBTESTABR =="MAL/CR")
#bda <- studya %>% filter(SAFFL=="Y")
trt_merge <- studyc %>% filter(VISID =="1" |VISID =="12"|VISID =="8"|VISID =="801")

trt_merge <- trt_merge[complete.cases(trt_merge[,4]),]
trt_merge <- trt_merge %>% filter (LBRUCD=="95")
trt_merge <- na.locf(trt_merge)

trt_merge$aval_unchanged <- trt_merge$LBRN
trt_merge$base_unchanged <- trt_merge$LBBLVALTR
trt_merge$AVAL <- log(trt_merge$LBRN)
trt_merge$BASE <- log(trt_merge$LBBLVALTR)
trt_merge$Change <- trt_merge$AVAL - trt_merge$BASE
trt_merge$Change_oriscale <- trt_merge$aval_unchanged-trt_merge$base_unchanged
trt_merge$Percent_change3 <- log(trt_merge$aval_unchanged/trt_merge$base_unchanged)

trt_merge <- trt_merge[,c(-4,-5,-7,-8,-9)]
studyc_comb <- trt_merge
studyc_comb$SUBJID <- paste("studyc",studyc_comb$SUBJID, sep="_")
studyc_comb$Study <- "studyc"


## studye
lab <- read.csv("labse.csv")
adsl <- read.csv("adsle.csv")


studye <- lab %>% select(USUBJID, AVISIT, AVISITN, PARAM, PARAMCD, AVAL, BASE, SAFFL)
studye <- studye %>% filter(PARAMCD =="ALBCS49C")
studye <- studye %>% filter(SAFFL=="Y")
studye <- studye %>% filter(AVISITN =="10" |AVISITN =="100")
adsl_trt <- adsl %>% select(USUBJID ,TRT01A, TRT01AN)

trt_merge <- merge(studye, adsl_trt,all=TRUE)
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
names(trt_merge) <- c(colnames(studyd_comb))
studye_comb <- trt_merge
studye_comb$Study <- "studye"


###### Combine studies and Filter by biomarker value

studya_comb <- studya_comb[complete.cases(studya_comb),]
studye_comb <- studye_comb[complete.cases(studye_comb),]
names(studya_comb) <- c(colnames(studyd_comb))
names(studyb_comb) <- c(colnames(studyd_comb))
names(studyg_comb) <- c(colnames(studyd_comb))
names(studyi_comb) <- c(colnames(studyd_comb))
names(studyx_comb) <- c(colnames(studyd_comb))
names(studyc_comb) <- c(colnames(studyd_comb))
names(studye_comb) <- c(colnames(studyd_comb))

all <- merge(studycf_comb, studya_comb, all=TRUE)
all <- merge(all, studyb_comb, all=TRUE)
all <- merge(all, studyd_comb, all=TRUE)
all <- merge(all, studyg_comb, all=TRUE)
all <- merge(all, studyi_comb, all=TRUE)
all <- merge(all, studyx_comb, all=TRUE)
all <- merge(all, studyc_comb, all=TRUE)
all <- merge(all, studye_comb, all=TRUE)

all <- all%>% mutate(TRT=replace(TRT,TRT=="Drug 1.5","Drug_1.5"))
all <- all%>% mutate(TRT=replace(TRT,TRT=="Drug 0.75","Drug_0.75"))
all <- all%>% mutate(TRT=replace(TRT,TRT=="Insulin Glargine","Glargine"))

all <- all%>% filter(TRT=="Drug_1.5" |TRT=="Drug_0.75"|TRT=="Placebo"|TRT=="Glargine")
placebo_studies <- all %>% filter(Study=="studycf"|Study=="studya"|Study=="studyg"|Study=="studyi")
glargine_studies <- all %>% filter(Study=="studyx"|Study=="studyd"|Study=="studyb")

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


pdf("biomarker_regression_plots.pdf")

placebo <-ggplot(placebo_studies, aes(BASE, Change,shape=TRT, colour=TRT, fill=TRT)) +geom_smooth(method="lm", se=FALSE) + geom_vline(xintercept=log(30))
glargine<-ggplot(glargine_studies, aes(BASE, Change,shape=TRT, colour=TRT, fill=TRT)) +geom_smooth(method="lm", se=FALSE) + geom_vline(xintercept=log(30))
print(placebo)
print(glargine)
dev.off()












