# Biomarker Regression Analysis
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


UACR_analysis <- function(study_name)
{
 if (study_name == "StudyA"){
    source("StudyA.R")
  }
  else if (study_name == "StudyB"){
    source("StudyB.R")
  }
  else if (study_name == "StudyC"){
    source("StudyC.R")
  }
  else if (study_name == "StudyD"){
    source("StudyD.R")
  }
  else if (study_name == "StudyE"){
    source("StudyE.R")
  }
  else if (study_name == "StudyCF"){
    source("StudyCF.R")
  }
  else if (study_name == "StudyG"){
    source("StudyG.R")
  }
  else if (study_name == "StudyI"){
    source("StudyI.R")
  }
  else if (study_name == "StudyJ"){
    source("StudyJ.R")
  }
  else if (study_name == "StudyX"){
    source("StudyX.R")
  }
  else if (study_name == "Glargine_combined"){
    source("Glargine-combined.R")
  }
  else if (study_name == "Placebo_combined"){
    source("Placebo-combined.R")
  }
  else if (study_name == "Glargine_split"){
    source("4split_Glargine_combined.R")
  }
  else if (study_name == "Placebo_split"){
    source("4split_Placebo_combined.R")
  }
  else if (study_name == "Baseline_analysis"){
    source("all_studies_baseline_regression.R")
  }
  else if (study_name == "New_compound"){
    source("New_compound.R")
  }
  else{
    print("No studies by that name, valid outputs include StudyA, StudyB, StudyC, StudyD, StudyE, StudyCF, StudyG, StudyI, StudyJ, StudyX,
              Glargine_combined, Placebo_combined, Glargine_split, Placebo_split, Baseline_analysis, New_compound")
  }
}





