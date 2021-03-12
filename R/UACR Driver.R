# UACR (Urine Albumin/Creatinine (g/kg))
#
# Summary and Analysis of UACR including percent change from baseline via log-transformed analysis MMRM from baseline through 52 weeks
# Least square mean percent change and standard error between study drug and placebo
#
# Fit to model log(y) = log(y_b) + trt
# After model is finished transform back by LSM=exp(LSM) and SE=exp(LSM)*SE
#
# Also find geometric mean and SE
# (Find n, lsm, geometric mean, se for baseline, post-baseline (week 52), change from baseline and percent change from baseline)
#
# Placebo controlled (GBCF, CBDG, GBDA, GBDI), Insulin Glargine (GBDB, GBDD, GBDX)
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

setwd('H:/R/UACR/R')
UACR_analysis <- function(study_name)
{
  if (study_name == "GBCF"){
    source("check.R")
  }
  else if (study_name == "GBDA"){
    source("GBDA.R")
  }
  else if (study_name == "GBDB"){
    source("GBDB.R")
  }
  else if (study_name == "GBDC"){
    source("GBDC.R")
  }
  else if (study_name == "GBDD"){
    source("GBDD.R")
  }

  else if (study_name == "GBDG"){
    source("GBDG.R")
  }
  else if (study_name == "GBDI"){
    source("GBDI.R")
  }
  else if (study_name == "GBDJ"){
    source("GBDJ.R")
  }
  else if (study_name == "GBDX"){
    source("GBDX.R")
  }
  else if (study_name == "Glargine_combined"){
    source("Glargine-combined_UACR.R")
  }
  else if (study_name == "Placebo_combined"){
    source("Placebo-combined_UACR.R")
  }
  else if (study_name == "Glargine_split"){
    source("UACR_4split_Glargine_combined.R")
  }
  else if (study_name == "Placebo_split"){
    source("UACR_4split_Placebo_combined.R")
  }
  else if (study_name == "Baseline_analysis"){
    source("UACR_all_studies_baseline.R")
  }
  else if (study_name == "TZP_GPGB"){
    source("TZP_GPGB.R")
  }
  else{
    print("No studies by that name, valid outputs include GBCF, GBDA, GBDB, GBDC, GBDD, GBDE, GBDG, GBDI, GBDJ, GBDX,
              Glargine_combined, Placebo_combined, Glargine_split, Placebo_split, Baseline_analysis, TZP_GPGB")
  }
}





