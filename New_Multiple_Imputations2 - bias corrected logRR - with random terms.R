#################################################################
# title: 'Missing SD and SS in meta-analysis: Run the Simulations 
# by: Stephan Kambach
# Date: 21.12.2018
#################################################################

# libraries and data
rm(list=ls())
gc()

library(mice)
library(ggplot2)
library(metafor)
library(ROCR)
library(data.table)
library(randomForest)
library(cowplot)
library(extrafont)
library(mi)
library(Amelia)
library(missForest)
library(Hmisc)
library(betareg)

setwd("C:/Users/kambach/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("C:/Users/Agando/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("D:/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("C:/Users/localadmin/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("C:/Users/sk85xupa/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
source("additional functions2 - bias corrected logRR.R")

# 0. notes ----------------------------------------------------------------

# 3 effect sizes:
#     - Hedge's d
#     - log RR
#     - Fisher's z

# 4 dataset sizes:
#     - 20
#     - 50
#     - 100
#     - 1000

# dataset properties:

# 2 weighting schemes:
#     - based on 1/effect size variance
#     - based on sample sizes

# missingness patterns:
#     - missing completely at random (MCAR)
#     - missing with 1/effect size (MAR, conditional on other variable)
#     - missing with effect size variance (MNAR, depending on the target variable)

# correlation/dataset patterns: 
#     - not correlations
#     - SD ~ 1 / mean and SS ~ mean

# imputation methods:
#         - complete - complete case analysis
#         - unweighted - unweighted analysis
#         - mean_value - mean - mean value imputation
#         - median_value - median imputaion - median value imputation
#         - random_sample - sample - random sample imputation
#         - bayes_predict - mice::norm:predict - Predicted vbalues from linear regression
#         - pmm - mice::pmm - predictive mean matching
#         - cart - mice::cart - classification and regression trees
#         - random_forest - mice::rf - random forest
#         - bayes_pmm - mi::mi - bayes framework (pmm), chained equation bootstrap algorithm to model missing variables based on observed variable (bayesglm)
#         - bootstrap_EMB - Amelia-package: bootstrap EMB algorithm - needs multivariate normal distribution (for SDs?) and data missing at random
#         - missForest - missForest - random forest
#         - Hmisc - Hmisc - uses bootstrapping, then additive regression, then pmm, assumes linearity in the imputed variables


# 1.a - function to create logRR/SMDH datasets ----------------------------

create.logRR.SMDH.dataset = function(SD.correlation, SS.correlation, cov.correlation, dataset.size, SD.max, SS.max){
  # arguments:
  #   - SD.correlation can be one of the following:
  #       "no correlation","SD ~ mean", "SD ~ 1/mean"
  #   - SS.correlation can be one of the following:
  #       "no correlation","SS_treat = SS_contr", "SS ~ mean", "SS ~ 1/mean", "SS ~ 1/SD"
  #   - cov.correlation can be one of the following:
  #       "no correlation", "cov ~ SD", "cov ~ SS"
  #   - dataset.size can be one of the following:
  #       20, 50, 500
  #   - SD.max can be one of the following:
  #       - 1, 5, 10
  #   - SS.max can be one of the following:
  #       - 10, 50, 100

  dataset.temp = data.frame("study_Id" = factor(1:dataset.size),
                            "contr_mean" = rnorm(n = dataset.size, mean = 1, sd = 0.1))
  
  dataset.temp$treat_mean = dataset.temp$contr_mean * seq(from = 0.5, to  = 2.5, length.out = dataset.size) + rnorm(dataset.size, mean = 0 , sd = 0.01)
  
  # add SD values
  if(SD.correlation == "no correlation"){
    dataset.temp$contr_SD = sample(seq(from = 0.1, to  = SD.max, length.out = dataset.size), size = dataset.size, replace = F)
    dataset.temp$treat_SD = sample(seq(from = 0.1, to  = SD.max, length.out = dataset.size), size = dataset.size, replace = F)}
  
  if(SD.correlation == "SD ~ mean"){
    dataset.temp$contr_SD = seq(from = 0.1, to  = SD.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SD.max / 100)
    dataset.temp$treat_SD = seq(from = 0.1, to  = SD.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SD.max / 100)}
  
  if(SD.correlation == "SD ~ 1/mean"){
    dataset.temp$contr_SD = rev(seq(from = 0.1, to  = SD.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SD.max / 100))
    dataset.temp$treat_SD = rev(seq(from = 0.1, to  = SD.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SD.max / 100))}
  
  # add SS values
  if(SS.correlation == "no correlation"){
    dataset.temp$contr_SS = round(sample(seq(from = 3, to  = SS.max, length.out = dataset.size), size = dataset.size, replace = F) + rnorm(dataset.size, mean = 0, sd = SS.max / 100))
    dataset.temp$treat_SS = round(sample(seq(from = 3, to  = SS.max, length.out = dataset.size), size = dataset.size, replace = F) + rnorm(dataset.size, mean = 0, sd = SS.max / 100))}
  
  if(SS.correlation == "SS_treat = SS_contr"){
    dataset.temp$contr_SS = round(sample(seq(from = 3, to  = SS.max, length.out = dataset.size), size = dataset.size, replace = F) + rnorm(dataset.size, mean = 0, sd = SS.max / 100))
    dataset.temp$treat_SS = round(dataset.temp$contr_SS + rnorm(dataset.size, mean = 0, sd = SS.max / 100))}
  
  if(SS.correlation == "SS ~ mean"){
    dataset.temp$contr_SS = round(seq(from = 3, to  = SS.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SS.max / 100))
    dataset.temp$treat_SS = round(seq(from = 3, to  = SS.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SS.max / 100))}
  
  if(SS.correlation == "SS ~ 1/mean"){
    dataset.temp$contr_SS = round(rev(seq(from = 3, to  = SS.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SS.max / 100)))
    dataset.temp$treat_SS = round(rev(seq(from = 3, to  = SS.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SS.max / 100)))}
  
  if(SS.correlation == "SS ~ 1/SD"){
    contr_n.order = match(seq(1:dataset.size), order(dataset.temp$contr_SD))
    treat_n.order = match(seq(1:dataset.size), order(dataset.temp$treat_SD))
    dataset.temp$contr_SS = rev(round(seq(from = 3, to  = SS.max, length.out = dataset.size))[contr_n.order]  + round(rnorm(dataset.size, mean = 0, sd = SS.max / 100)))
    dataset.temp$treat_SS = rev(round(seq(from = 3, to  = SS.max, length.out = dataset.size))[treat_n.order]  + round(rnorm(dataset.size, mean = 0, sd = SS.max / 100)))}
  
  # SD should be between 0.01 and 10
  dataset.temp$contr_SD[which(dataset.temp$contr_SD < 0.01)] = 0.01
  dataset.temp$treat_SD[which(dataset.temp$treat_SD < 0.01)] = 0.01
  dataset.temp$contr_SD[which(dataset.temp$contr_SD > 10)] = 10
  dataset.temp$treat_SD[which(dataset.temp$treat_SD > 10)] = 10
  
  # SS should be between 3 and 100
  dataset.temp$contr_SS[which(dataset.temp$contr_SS < 3)] = 3
  dataset.temp$treat_SS[which(dataset.temp$treat_SS < 3)] = 3
  dataset.temp$contr_SS[which(dataset.temp$contr_SS > 100)] = 100
  dataset.temp$treat_SS[which(dataset.temp$treat_SS > 100)] = 100
  
  
  # add cov values
  if(cov.correlation == "no correlation"){
    dataset.temp$cov =  rnorm(dataset.size, mean = 0, sd = 10)}
  
  if(cov.correlation == "cov ~ SD"){
    cov.order = match(seq(1:dataset.size), order(dataset.temp$contr_SD))
    dataset.temp$cov = seq(from = 1, to  = 100, length.out = dataset.size)[cov.order]  + rnorm(dataset.size, mean = 0, sd = 1)}

  if(cov.correlation == "cov ~ SS"){
    cov.order = match(seq(1:dataset.size), order(dataset.temp$contr_SS))
    dataset.temp$cov = seq(from = 1, to  = 100, length.out = dataset.size)[cov.order]  + rnorm(dataset.size, mean = 0, sd = 1)}
  
  return(dataset.temp)
}
# 1b. function to create Fisher's z datasets ------------------------------

create.ZCOR.dataset = function(dataset.size, SS.max, SS.correlation, cov.correlation){
  # arguments:
  #   - SS.correlation can be one of the following:
  #       "no correlation","SS ~ cor_coef", "SS ~ 1/cor_coef"
  #   - cov.correlation can be one of the following:
  #       "no correlation", "cov ~ SS"
  #   - dataset.size can be one of the following:
  #       20, 50, 100, 500
  #   - SS.max can be one of the following: (SS.min = 5)
  #       - 10, 50, 100
  
  dataset.temp = data.frame("study_Id" = factor(1:dataset.size),
                            "cor_coef" = seq(from = 0.01, to  = 0.99, length.out = dataset.size))
  
  # add SS values
  if(SS.correlation == "no correlation"){
    dataset.temp$SS = round(sample(seq(from = 5, to  = SS.max, length.out = dataset.size), size = dataset.size, replace = F))}
  
  if(SS.correlation == "SS ~ cor_coef"){
    dataset.temp$SS = round(seq(from = 5, to  = SS.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SS.max / 100))}
  
  if(SS.correlation == "SS ~ 1/cor_coef"){
    dataset.temp$SS = round(rev(seq(from = 5, to  = SS.max, length.out = dataset.size) + rnorm(dataset.size, mean = 0, sd = SS.max / 100)))}
  
  # SS should be between 5 and 100
  dataset.temp$SS[dataset.temp$SS < 5] = 5
  dataset.temp$SS[dataset.temp$SS > 100] = 100
  
  # add cov values
  if(cov.correlation == "no correlation"){
    dataset.temp$cov =  rnorm(dataset.size, mean = 0, sd = 10)}
  
  if(cov.correlation == "cov ~ SS"){
    cov.order = match(seq(1:dataset.size), order(dataset.temp$SS))
    dataset.temp$cov = seq(from = 1, to  = 100, length.out = dataset.size)[cov.order]  + rnorm(dataset.size, mean = 0, sd = 1)}
  
  return(dataset.temp)
}
# 1c. test the functions to create the datasets ---------------------------

test = create.logRR.SMDH.dataset(dataset.size = 1000, 
                                 SD.correlation = "SD ~ 1/mean", 
                                 SS.correlation = "SS ~ 1/mean", 
                                 SD.max = 10, 
                                 SS.max = 100, 
                                 cov.correlation = "cov ~ SD")  
plot(treat_SD ~ treat_mean, data = test) # works

test = create.ZCOR.dataset(SS.correlation = "SS ~ 1/cor_coef", 
                           cov.correlation = "cov ~ SS", 
                           dataset.size = 100, 
                           SS.max = 100)  
plot(SS ~ cor_coef, data = test) # works
# 2. create all the datasets ----------------------------------------------
all.logRR.SMDH.datasets = list()
for(SD.correlation in c("no correlation", "SD ~ 1/mean")){
  for(SS.correlation in c("SS_treat = SS_contr", "SS ~ mean", "SS ~ 1/SD")){
    for(cov.correlation in c("no correlation", "cov ~ SD", "cov ~ SS")){
      for(dataset.size in c("50", "100","500")){
        for(SD.max in c("1", "10")){
          for(SS.max in c("10", "100")){
            all.logRR.SMDH.datasets[[SD.correlation]][[SS.correlation]][[cov.correlation]][[dataset.size]][[SD.max]][[SS.max]] =
              create.logRR.SMDH.dataset(SD.correlation = SD.correlation, 
                                        SS.correlation = SS.correlation,
                                        cov.correlation = cov.correlation,
                                        dataset.size = as.numeric(dataset.size),
                                        SD.max = as.numeric(SD.max),
                                        SS.max = as.numeric(SS.max))
          }
        }
      }
    }
  }
}

all.ZCOR.datasets = list()
for(SS.correlation in c("no correlation","SS ~ cor_coef", "SS ~ 1/cor_coef")){
  for(cov.correlation in c("no correlation", "cov ~ SS")){
    for(dataset.size in c("50","100", "500")){
      for(SS.max in c("10", "100")){
        all.ZCOR.datasets[[SS.correlation]][[cov.correlation]][[dataset.size]][[SS.max]] =
          create.ZCOR.dataset(SS.correlation = SS.correlation,
                              cov.correlation = cov.correlation,
                              dataset.size = as.numeric(dataset.size),
                              SS.max = as.numeric(SS.max))
      }
    }
  }
}
# 3a. save the full data sets ----------------------------------------------

saveRDS(all.logRR.SMDH.datasets, "intermediate results2 - bias corrected logRR/all_logRR_SMDH_datasets.rds")
saveRDS(all.ZCOR.datasets, "intermediate results2 - bias corrected logRR/all_ZCOR_datasets.rds")
# 3b. load the full data sets ----------------------------------------------

all.logRR.SMDH.datasets = readRDS("intermediate results2 - bias corrected logRR/all_logRR_SMDH_datasets.rds")
all.ZCOR.datasets = readRDS("intermediate results2 - bias corrected logRR/all_ZCOR_datasets.rds")

# 4. create Figure 2 ( = Fig. 3 in the manuscript) - MCAR -----------------------------------------

# Figure 2: Effect of MCAR with a sample size of 100 and no further correlation
logRR.SMDH.data.fig2 = all.logRR.SMDH.datasets[["no correlation"]][["SS_treat = SS_contr"]][["no correlation"]][["100"]][["10"]][["100"]]
ZCOR.data.fig2 = all.ZCOR.datasets[["no correlation"]][["no correlation"]][["100"]][["100"]]

# delete study_ID because some of the algorithm do not work otherwise
logRR.SMDH.data.fig2$study_Id = NULL
ZCOR.data.fig2$study_Id = NULL
# 4a. - full and unweighted analysis + delete and impute SDs in the logRR and SMDH dataset ---------------

# create dataframe to store results
fig2.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

# grand mean from full datasets (SMDH and log RR)
fig2.results.temp = imp.full.dataset.logRR.SMDH(logRR.SMDH.data.fig2, what.is.missing = "nothing")
fig2.results = rbind(fig2.results, fig2.results.temp)
  
# grand mean from full datasets (ZCOR)
fig2.results.temp = imp.full.dataset.ZCOR(ZCOR.data.fig2, what.is.missing = "nothing")
fig2.results = rbind(fig2.results, fig2.results.temp)

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig2

  # delete random SDs
  rows.to.del = sample(1:100, size = 100 * del.rate, replace = F)
  data.del.logRR.SMDH$contr_SD[rows.to.del] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del] = NA

  # 1. complete case analysis
  fig2.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 2. unweighted analysis
  fig2.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # x. sample-size-weighted analysis
  fig2.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 3. mean value imputation
  fig2.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 4. median value imputation
  fig2.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 5. random sample imputation
  fig2.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 6. linear prediction
  fig2.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 7. predictive mean matching
  fig2.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 8. classification and regression trees
  fig2.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 9. random forest
  fig2.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 10. bayes pmm - mi
  fig2.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})  
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig2.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)}
  # 12. random forest (missForest)
  fig2.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig2.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)}, 
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # give status back
  print(paste("Deletion rate computed:", del.rate))

  # save and create new temporary table
  write.table(fig2.results, paste(c("intermediate results2 - bias corrected logRR/fig2_data_part1_",del.rate*100,".txt"), collapse = ""))
  fig2.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
}

# 4b. - delete and impute SDs and SS for the logRR and SMDH dataset --------
# create dataframe to store results
fig2.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig2
  
  # delete random SSs
  rows.to.del.SD = sample(1:100, size = 100 * del.rate, replace = F)
  rows.to.del.SS = sample(1:100, size = 100 * del.rate, replace = F)
  
  data.del.logRR.SMDH$contr_SD[rows.to.del.SD] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del.SD] = NA
  data.del.logRR.SMDH$contr_SS[rows.to.del.SS] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del.SS] = NA
  
  # 1. complete case analysis
  fig2.results.temp = tryCatch({
    imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 2. unweighted analysis
  fig2.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # x. sample-size-weighted analysis
  fig2.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 3. mean value imputation
  fig2.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 4. median value imputation
  fig2.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 5. random sample imputation
  fig2.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 6. linear prediction
  fig2.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 7. predictive mean matching
  fig2.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 8. classification and regression trees
  fig2.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 9. random forest
  fig2.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 10. bayes pmm - mi
  fig2.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig2.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)}
  # 12. random forest (missForest)
  fig2.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig2.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))

  # save and create new temporary table
  write.table(fig2.results, paste(c("intermediate results2 - bias corrected logRR/fig2_data_part2_",del.rate*100,".txt"), collapse = ""))
  fig2.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
}  

# 4c. - delete and impute SS (for the logRR, SMDH and ZCOR - datasets) --------
# create dataframe to store results
fig2.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig2
  data.del.ZCOR = ZCOR.data.fig2
  
  # delete random rows with SS
  rows.to.del = sample(1:100, size = 100 * del.rate, replace = F)
  data.del.logRR.SMDH$contr_SS[rows.to.del] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del] = NA
  data.del.ZCOR$SS[rows.to.del] = NA
  
  # 1. complete case analysis
  fig2.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = imp.complete.case.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 2. unweighted analysis
  fig2.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = imp.unweighted.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # x. sample-size-weighted analysis
  fig2.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 3. mean value imputation
  fig2.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = imp.mean.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 4. median value imputation
  fig2.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = imp.median.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 5. random sample imputation
  fig2.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    imp.random.sample.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 6. linear prediction
  fig2.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    imp.linear.predict.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 7. predictive mean matching
  fig2.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    imp.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 8. classification and regression trees
  fig2.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    imp.cart.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 9. random forest
  fig2.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    imp.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 10. bayes pmm - mi
  fig2.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    imp.mi.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig2.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    boot.ebm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)}
  # 12. random forest (missForest)
  fig2.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    miss.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig2.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  fig2.results.temp = tryCatch({
    boot.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig2.results = rbind(fig2.results, fig2.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig2.results, paste(c("intermediate results2 - bias corrected logRR/fig2_data_part3_",del.rate*100,".txt"), collapse = ""))
  fig2.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
}

# 5. create Figure 3 ( = Fig. 4 in the manuscript) - MAR ------------------------------------------

# Figure 3: Effect of MAR with a sample size of 100 and no further correlation
logRR.SMDH.data.fig3 = all.logRR.SMDH.datasets[["no correlation"]][["SS_treat = SS_contr"]][["no correlation"]][["100"]][["10"]][["100"]]
ZCOR.data.fig3 = all.ZCOR.datasets[["no correlation"]][["no correlation"]][["100"]][["100"]]

# delete study_ID because some of the algorithm do not work otherwise
logRR.SMDH.data.fig3$study_Id = NULL
ZCOR.data.fig3$study_Id = NULL
# 5a. - full and unweighted analysis + delete and impute SDs in the logRR and SMDH dataset ---------------

# create dataframe to store results
fig3.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

# grand mean from full datasets (SMDH and log RR)
fig3.results.temp = imp.full.dataset.logRR.SMDH(logRR.SMDH.data.fig3, what.is.missing = "nothing")
fig3.results = rbind(fig3.results, fig3.results.temp)

# grand mean from full datasets (ZCOR)
fig3.results.temp = imp.full.dataset.ZCOR(ZCOR.data.fig3, what.is.missing = "nothing")
fig3.results = rbind(fig3.results, fig3.results.temp)

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig3
  
  # deletion probability of SDs corrlates with 1/mean
  rows.to.del = sample(1:100, size = 100 * del.rate, prob = 100:1, replace = F)
  data.del.logRR.SMDH$contr_SD[rows.to.del] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del] = NA
  
  # 1. complete case analysis
  fig3.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 2. unweighted analysis
  fig3.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # x. sample-size-weighted analysis
  fig3.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 3. mean value imputation
  fig3.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 4. median value imputation
  fig3.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 5. random sample imputation
  fig3.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 6. linear prediction
  fig3.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 7. predictive mean matching
  fig3.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 8. classification and regression trees
  fig3.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 9. random forest
  fig3.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 10. bayes pmm - mi
  fig3.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})  
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig3.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)}
  # 12. random forest (missForest)
  fig3.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig3.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)}, 
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig3.results, paste(c("intermediate results2 - bias corrected logRR/fig3_data_part1_",del.rate*100,".txt"), collapse = ""))
  fig3.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
}

# 5b. - delete and impute SDs and SS for the logRR and SMDH datase --------
# create dataframe to store results
fig3.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig3
  
  # deletion probability of SDs and SSs corrlates with 1/mean^2
  rows.to.del.SD = sample(1:100, size = 100 * del.rate, prob = 100:1, replace = F)
  rows.to.del.SS = sample(1:100, size = 100 * del.rate, prob = 100:1, replace = F)
  
  data.del.logRR.SMDH$contr_SD[rows.to.del.SD] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del.SD] = NA
  data.del.logRR.SMDH$contr_SS[rows.to.del.SS] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del.SS] = NA
  
  # 1. complete case analysis
  fig3.results.temp = tryCatch({
    imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  
  # 2. unweighted analysis
  fig3.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # x. sample-size-weighted analysis
  fig3.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 3. mean value imputation
  fig3.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 4. median value imputation
  fig3.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 5. random sample imputation
  fig3.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 6. linear prediction
  fig3.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 7. predictive mean matching
  fig3.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 8. classification and regression trees
  fig3.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 9. random forest
  fig3.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 10. bayes pmm - mi
  fig3.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig3.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)}
  # 12. random forest (missForest)
  fig3.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig3.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig3.results, paste(c("intermediate results2 - bias corrected logRR/fig3_data_part2_",del.rate*100,".txt"), collapse = ""))
  fig3.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
} 

# 5c. - delete and impute SS (for the logRR, SMDH and ZCOR - datasets) --------
# create dataframe to store results
fig3.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig3
  data.del.ZCOR = ZCOR.data.fig3
  
  # deletion probability of SS corrlates with 1/mean^2
  rows.to.del = sample(1:100, size = 100 * del.rate, prob = 100:1, replace = F)
  data.del.logRR.SMDH$contr_SS[rows.to.del] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del] = NA
  data.del.ZCOR$SS[rows.to.del] = NA
  
  # 1. complete case analysis
  fig3.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = imp.complete.case.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 2. unweighted analysis
  fig3.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = imp.unweighted.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # x. sample-size-weighted analysis
  fig3.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 3. mean value imputation
  fig3.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = imp.mean.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 4. median value imputation
  fig3.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = imp.median.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 5. random sample imputation
  fig3.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    imp.random.sample.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 6. linear prediction
  fig3.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    imp.linear.predict.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 7. predictive mean matching
  fig3.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    imp.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 8. classification and regression trees
  fig3.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    imp.cart.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 9. classification and regression trees
  fig3.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    imp.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 10. bayes pmm - mi
  fig3.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    imp.mi.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){ fig3.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    boot.ebm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)}
  # 12. random forest (Forestmiss)
  fig3.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    miss.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig3.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  fig3.results.temp = tryCatch({
    boot.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig3.results = rbind(fig3.results, fig3.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig3.results, paste(c("intermediate results2 - bias corrected logRR/fig3_data_part3_",del.rate*100,".txt"), collapse = ""))
  fig3.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
} 

# 6. Figure 4 ( = Fig. 5 in the manuscript) - MNAR ------------------------------------------

# Figure 4: Effect of MNAR with a sample size of 100 and no further correlation
logRR.SMDH.data.fig4 = all.logRR.SMDH.datasets[["no correlation"]][["SS_treat = SS_contr"]][["no correlation"]][["100"]][["10"]][["100"]]
ZCOR.data.fig4 = all.ZCOR.datasets[["no correlation"]][["no correlation"]][["100"]][["100"]]

# delete study_ID because some of the algorithm do not work otherwise
logRR.SMDH.data.fig4$study_Id = NULL
ZCOR.data.fig4$study_Id = NULL
# 6a. - SS/SD ~ full and unweighted analysis + delete and impute SDs in the logRR and SMDH dataset ---------------

# create dataframe to store results
fig4.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

# grand mean from full datasets (SMDH and log RR)
fig4.results.temp = imp.full.dataset.logRR.SMDH(logRR.SMDH.data.fig4, what.is.missing = "nothing")
fig4.results = rbind(fig4.results, fig4.results.temp)

# grand mean from full datasets (ZCOR)
fig4.results.temp = imp.full.dataset.ZCOR(ZCOR.data.fig4, what.is.missing = "nothing")
fig4.results = rbind(fig4.results, fig4.results.temp)

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig4
  # deletion probability of SDs corrlates with SD
  rows.to.del.SD.logRR.SMDH = sample(1:100, size = 100 * del.rate, 
                                     prob = rank(data.del.logRR.SMDH$contr_SD +  data.del.logRR.SMDH$treat_SD), replace = F)
  data.del.logRR.SMDH$contr_SD[rows.to.del.SD.logRR.SMDH] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del.SD.logRR.SMDH] = NA
  
  # 1. complete case analysis
  fig4.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 2. unweighted analysis
  fig4.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # x. sample-size-weighted analysis
  fig4.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 3. mean value imputation
  fig4.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 4. median value imputation
  fig4.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 5. random sample imputation
  fig4.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 6. linear prediction
  fig4.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 7. predictive mean matching
  fig4.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 8. classification and regression trees
  fig4.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 9. random forest
  fig4.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 10. bayes pmm - mi
  fig4.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})  
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){ fig4.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)}
  # 12. random forest (missForest)
  fig4.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig4.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)}, 
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig4.results, paste(c("intermediate results2 - bias corrected logRR/fig4_data_part1_",del.rate*100,".txt"), collapse = ""))
  fig4.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
}

# 6b. - delete and impute SDs and SS for the logRR and SMDH datase --------
# create dataframe to store results
fig4.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig4
  
  # deletion probability of SDs and SSs corrlates with 1/mean^2
  rows.to.del.SD.logRR.SMDH = sample(1:100, size = 100 * del.rate, 
                                     prob = rank(data.del.logRR.SMDH$treat_SD + data.del.logRR.SMDH$contr_SD), replace = F)
  rows.to.del.SS.logRR.SMDH = sample(1:100, size = 100 * del.rate, 
                                     prob = rank(-(data.del.logRR.SMDH$treat_SS + data.del.logRR.SMDH$contr_SS)), replace = F)
  
  data.del.logRR.SMDH$contr_SD[rows.to.del.SD.logRR.SMDH] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del.SD.logRR.SMDH] = NA
  data.del.logRR.SMDH$contr_SS[rows.to.del.SS.logRR.SMDH] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del.SS.logRR.SMDH] = NA
  
  # 1. complete case analysis
  fig4.results.temp = tryCatch({
    imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 2. unweighted analysis
  fig4.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # x. sample-size-weighted analysis
  fig4.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 3. mean value imputation
  fig4.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 4. median value imputation
  fig4.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 5. random sample imputation
  fig4.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 6. linear prediction
  fig4.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 7. predictive mean matching
  fig4.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 8. classification and regression trees
  fig4.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 9. random forest
  fig4.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 10. bayes pmm - mi
  fig4.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig4.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)}
  # 12. random forest (missForest)
  fig4.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig4.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig4.results, paste(c("intermediate results2 - bias corrected logRR/fig4_data_part2_",del.rate*100,".txt"), collapse = ""))
  fig4.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
} 

# 6c. - delete and impute SS (for the logRR, SMDH and ZCOR - datasets) --------
# create dataframe to store results
fig4.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig4
  data.del.ZCOR = ZCOR.data.fig4
  
  # deletion probability of SS corrlates with 1/mean^2
  rows.to.del.SS.logRR.SMDH = sample(1:100, size = 100 * del.rate, prob = rank(-(data.del.logRR.SMDH$treat_SS + data.del.logRR.SMDH$contr_SS)), replace = F)
  rows.to.del.SS.ZCOR = sample(1:100, size = 100 * del.rate, prob = rank(-(data.del.ZCOR$SS)), replace = F)
  
  data.del.logRR.SMDH$contr_SS[rows.to.del.SS.logRR.SMDH] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del.SS.logRR.SMDH] = NA
  data.del.ZCOR$SS[rows.to.del.SS.ZCOR] = NA
  
  # 1. complete case analysis
  fig4.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = imp.complete.case.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 2. unweighted analysis
  fig4.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = imp.unweighted.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # x. sample-size-weighted analysis
  fig4.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 3. mean value imputation
  fig4.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = imp.mean.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 4. median value imputation
  fig4.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = imp.median.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 5. random sample imputation
  fig4.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    imp.random.sample.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 6. linear prediction
  fig4.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    imp.linear.predict.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 7. predictive mean matching
  fig4.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    imp.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 8. classification and regression trees
  fig4.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    imp.cart.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 9. classification and regression trees
  fig4.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    imp.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 10. bayes pmm - mi
  fig4.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    imp.mi.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig4.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    boot.ebm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)}
  # 12. bootstrap EMB
  fig4.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    miss.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig4.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  fig4.results.temp = tryCatch({
    boot.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig4.results = rbind(fig4.results, fig4.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig4.results, paste(c("intermediate results2 - bias corrected logRR/fig4_data_part3_",del.rate*100,".txt"), collapse = ""))
  fig4.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
} 

# 7. create Figure 4 ( = Fig. 6 in the manuscript) - MCAR - SD ~ 1/mean, & SS ~ mean-----------------------------------------

# Figure 5: Effect of SD ~ mean & SS ~ 1/mean, with MCAR, sample size of 100 and no further correlation
logRR.SMDH.data.fig5 = all.logRR.SMDH.datasets[["SD ~ 1/mean"]][["SS ~ mean"]][["no correlation"]][["100"]][["10"]][["100"]]
ZCOR.data.fig5 = all.ZCOR.datasets[["SS ~ cor_coef"]][["no correlation"]][["100"]][["100"]]

# delete study_ID because some of the algorithm do not work otherwise
logRR.SMDH.data.fig5$study_Id = NULL
ZCOR.data.fig5$study_Id = NULL
# 7a. - full and unweighted analysis + delete and impute SDs in the logRR and SMDH dataset ---------------

# create dataframe to store results
fig5.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

# grand mean from full datasets (SMDH and log RR)
fig5.results.temp = imp.full.dataset.logRR.SMDH(logRR.SMDH.data.fig5, what.is.missing = "nothing")
fig5.results = rbind(fig5.results, fig5.results.temp)

# grand mean from full datasets (ZCOR)
fig5.results.temp = imp.full.dataset.ZCOR(ZCOR.data.fig5, what.is.missing = "nothing")
fig5.results = rbind(fig5.results, fig5.results.temp)

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig5
  # deletion probability of SDs corrlates with 1/mean
  rows.to.del = sample(1:100, size = 100 * del.rate, replace = F)
  data.del.logRR.SMDH$contr_SD[rows.to.del] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del] = NA
  
  # 1. complete case analysis
  fig5.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 2. unweighted analysis
  fig5.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # x. sample-size-weighted analysis
  fig5.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 3. mean value imputation
  fig5.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 4. median value imputation
  fig5.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 5. random sample imputation
  fig5.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 6. linear prediction
  fig5.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 7. predictive mean matching
  fig5.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 8. classification and regression trees
  fig5.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 9. random forest
  fig5.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 10. bayes pmm - mi
  fig5.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})  
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig5.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)}
  # 12. random forest (missForest)
  fig5.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig5.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)}, 
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)

  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig5.results, paste(c("intermediate results2 - bias corrected logRR/fig5_data_part1_",del.rate*100,".txt"), collapse = ""))
  fig5.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
}

# 7b. - delete and impute SDs and SS for the logRR and SMDH datase --------
# create dataframe to store results
fig5.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig5
  
  # deletion probability of SDs and SSs corrlates with 1/mean^2
  rows.to.del.SD = sample(1:100, size = 100 * del.rate, replace = F)
  rows.to.del.SS = sample(1:100, size = 100 * del.rate, replace = F)
  
  data.del.logRR.SMDH$contr_SD[rows.to.del.SD] = NA
  data.del.logRR.SMDH$treat_SD[rows.to.del.SD] = NA
  data.del.logRR.SMDH$contr_SS[rows.to.del.SS] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del.SS] = NA
  
  # 1. complete case analysis
  fig5.results.temp = tryCatch({
    imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 2. unweighted analysis
  fig5.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # x. sample-size-weighted analysis
  fig5.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 3. mean value imputation
  fig5.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 4. median value imputation
  fig5.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 5. random sample imputation
  fig5.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 6. linear prediction
  fig5.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 7. predictive mean matching
  fig5.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 8. classification and regression trees
  fig5.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 9. random forest
  fig5.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 10. bayes pmm - mi
  fig5.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig5.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)}
  # 12. random forest (missForest)
  fig5.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig5.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD & SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig5.results, paste(c("intermediate results2 - bias corrected logRR/fig5_data_part2_",del.rate*100,".txt"), collapse = ""))
  fig5.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
} 

# 7c. - delete and impute SS (for the logRR, SMDH and ZCOR - datasets) --------
# create dataframe to store results
fig5.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                           "imputation_method" = NA, "missing_percent" = NA,
                           "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]

for(del.rate in seq(from = 0.1, to = 0.9, by = 0.05)){
  
  # create temporary datasets
  data.del.logRR.SMDH = logRR.SMDH.data.fig5
  data.del.ZCOR = ZCOR.data.fig5
  
  # deletion probability of SS corrlates with 1/mean^2
  rows.to.del = sample(1:100, size = 100 * del.rate, replace = F)
  data.del.logRR.SMDH$contr_SS[rows.to.del] = NA
  data.del.logRR.SMDH$treat_SS[rows.to.del] = NA
  data.del.ZCOR$SS[rows.to.del] = NA
  
  # 1. complete case analysis
  fig5.results.temp = imp.complete.case.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = imp.complete.case.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 2. unweighted analysis
  fig5.results.temp = imp.unweighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = imp.unweighted.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # x. sample-size-weighted analysis
  fig5.results.temp = imp.sample.size.weighted.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SD", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 3. mean value imputation
  fig5.results.temp = imp.mean.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = imp.mean.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 4. median value imputation
  fig5.results.temp = imp.median.value.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = imp.median.value.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 5. random sample imputation
  fig5.results.temp = tryCatch({
    imp.random.sample.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    imp.random.sample.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 6. linear prediction
  fig5.results.temp = tryCatch({
    imp.linear.predict.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    imp.linear.predict.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 7. predictive mean matching
  fig5.results.temp = tryCatch({
    imp.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    imp.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 8. classification and regression trees
  fig5.results.temp = tryCatch({
    imp.cart.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    imp.cart.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 9. classification and regression trees
  fig5.results.temp = tryCatch({
    imp.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    imp.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 10. bayes pmm - mi
  fig5.results.temp = tryCatch({
    imp.mi.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    imp.mi.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 11. bootstrap EMB
  if(del.rate < 0.65){fig5.results.temp = tryCatch({
    boot.ebm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    boot.ebm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)}
  # 12. bootstrap EMB
  fig5.results.temp = tryCatch({
    miss.rf.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    miss.rf.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  # 13. bootstrap pmm - Hmisc
  fig5.results.temp = tryCatch({
    boot.pmm.logRR.SMDH(data.del.logRR.SMDH, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  fig5.results.temp = tryCatch({
    boot.pmm.ZCOR(data.del.ZCOR, what.is.missing = "SS", del.rate = del.rate)},
    error = function(e){create.empty.result.data.frame()})
  fig5.results = rbind(fig5.results, fig5.results.temp)
  
  # give status back
  print(paste("Deletion rate computed:", del.rate))
  
  # save and create new temporary table
  write.table(fig5.results, paste(c("intermediate results2 - bias corrected logRR/fig5_data_part3_",del.rate*100,".txt"), collapse = ""))
  fig5.results =  data.frame("effect_size" = NA, "weighting_by" = NA, "what_is_missing" = NA,
                             "imputation_method" = NA, "missing_percent" = NA,
                             "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, "tau" = NA, "i2" = NA, "h2" = NA)[0,]
} 
