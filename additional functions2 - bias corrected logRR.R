
# y - calculate bias-corrected logRR --------------------------------------
calculate.logRR.bias.corrected = function(dat.temp){
  yi = dat.temp$yi + ( 0.5 * 
                         (((dat.temp$treat_SD^2)/(dat.temp$treat_SS * (dat.temp$treat_mean^2))) -
                            ((dat.temp$contr_SD^2)/(dat.temp$contr_SS * (dat.temp$contr_mean^2)))))
  vi = dat.temp$vi + ( 0.5 * 
                         (((dat.temp$treat_SD^4)/((dat.temp$treat_SS^2) * (dat.temp$treat_mean^4))) +
                            ((dat.temp$contr_SD^4)/((dat.temp$contr_SS^2) * (dat.temp$contr_mean^4)))))
  
  dat.temp$yi = yi
  dat.temp$vi = vi
  
  return(dat.temp)

}


# z. create empty data frame when algorithms fail -------------------------

create.empty.result.data.frame = function(){
  results.df = data.frame("effect_size" = NA, 
                          "weighting_by" = NA, 
                          "what_is_missing" = NA,
                          "imputation_method" = NA, 
                          "missing_percent" = NA,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[0,]
  return(results.df)
}

# 0a. analysis of full dataset logRR/SMDH ---------------------------------------------
imp.full.dataset.logRR.SMDH = function(data.full, what.is.missing){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR", "SMDH"), 
                          "weighting_by" = c("vi", "vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = "full dataset analysis", 
                          "missing_percent" = rep(0,2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # calculate effect sizes # now the bias corrected logRR
  data.full.temp.log.RR = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                            m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                            data = data.full, measure="ROM"))
  data.full.temp.log.RR = calculate.logRR.bias.corrected(data.full.temp.log.RR)
  
  data.full.temp.SMDH = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                          m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                          data = data.full, measure="SMDH"))
  
  # run rma
  rma.imp.logRR = rma(yi = data.full.temp.log.RR$yi, vi = data.full.temp.log.RR$vi, measure = "ROM",
                      control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
  rma.imp.SMDH = rma(yi = data.full.temp.SMDH$yi, vi = data.full.temp.SMDH$vi, measure = "SMDH",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = c(rma.imp.logRR$beta, rma.imp.SMDH$beta)
  results.df$ci_lb = c(rma.imp.logRR$ci.lb, rma.imp.SMDH$ci.lb)
  results.df$ci_ub = c(rma.imp.logRR$ci.ub, rma.imp.SMDH$ci.ub)
  results.df$tau = c(rma.imp.logRR$tau2, rma.imp.SMDH$tau2)
  results.df$i2 = c(rma.imp.logRR$I2, rma.imp.SMDH$I2)
  results.df$h2 = c(rma.imp.logRR$H2, rma.imp.SMDH$H2)
  
  # return results
  return(results.df)
  
}

# 0b. analysis of full dataset ZCOR---------------------------------------------
imp.full.dataset.ZCOR = function(data.full, what.is.missing){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = what.is.missing,
                          "imputation_method" = "full dataset analysis", 
                          "missing_percent" = 0,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # calculate effect sizes
  data.full.temp.log.RR = data.frame(escalc(ri = cor_coef, ni = SS, data = data.full, measure = "ZCOR"))
  
  # run rma
  rma.imp.ZCOR = rma(yi = data.full.temp.log.RR$yi, vi = data.full.temp.log.RR$vi, measure = "ZCOR",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = rma.imp.ZCOR$beta
  results.df$ci_lb = rma.imp.ZCOR$ci.lb
  results.df$ci_ub = rma.imp.ZCOR$ci.ub
  results.df$tau = rma.imp.ZCOR$tau2
  results.df$i2 = rma.imp.ZCOR$I2
  results.df$h2 = rma.imp.ZCOR$H2
  
  # return results
  return(results.df)
}

# 1a. complete case analysis logRR/SMDH -----------------------------------------------

imp.complete.case.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR", "SMDH"), 
                          "weighting_by" = c("vi", "vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("complete case analysis",2), 
                          "missing_percent" = rep(del.rate,2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # calculate effect sizes
  data.del.temp.log.RR = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                           m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                           data = data.del, measure="ROM"))
  data.del.temp.log.RR = calculate.logRR.bias.corrected(data.del.temp.log.RR)
  data.del.temp.SMDH = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                         m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                         data = data.del, measure="SMDH"))
  
  # run rma
  rma.imp.logRR = rma(yi = data.del.temp.log.RR$yi, vi = data.del.temp.log.RR$vi, measure = "ROM",
                      control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
  rma.imp.SMDH = rma(yi = data.del.temp.SMDH$yi, vi = data.del.temp.SMDH$vi, measure = "SMDH",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = c(rma.imp.logRR$beta, rma.imp.SMDH$beta)
  results.df$ci_lb = c(rma.imp.logRR$ci.lb, rma.imp.SMDH$ci.lb)
  results.df$ci_ub = c(rma.imp.logRR$ci.ub, rma.imp.SMDH$ci.ub)
  results.df$tau = c(rma.imp.logRR$tau2, rma.imp.SMDH$tau2)
  results.df$i2 = c(rma.imp.logRR$I2, rma.imp.SMDH$I2)
  results.df$h2 = c(rma.imp.logRR$H2, rma.imp.SMDH$H2)
  
  # return results
  return(results.df)
}


# 1b. complete case analysis ZCOR -----------------------------------------------

imp.complete.case.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "complete case analysis", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # calculate effect sizes
  data.del.temp.ZCOR = data.frame(escalc(ri = cor_coef, ni = SS,
                                         data = data.del, measure="ZCOR"))
  
  # run rma
  rma.imp.ZCOR = rma(yi = data.del.temp.ZCOR$yi, vi = data.del.temp.ZCOR$vi, measure = "ROM",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = rma.imp.ZCOR$beta
  results.df$ci_lb = rma.imp.ZCOR$ci.lb
  results.df$ci_ub = rma.imp.ZCOR$ci.ub
  results.df$tau = rma.imp.ZCOR$tau2
  results.df$i2 = rma.imp.ZCOR$I2
  results.df$h2 = rma.imp.ZCOR$H2
  
  # return results
  return(results.df)
}

# 2a. unweighted analysis logRR/SMDH --------------------------------------------------

imp.unweighted.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR", "SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = "unweighted analysis", 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # calculate effect sizes <- Hedge's d cannot be calculated without SDs
  data.del.temp.log.RR = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                           m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                           data = data.del, measure="ROM"))
  #data.del.temp.log.RR = calculate.logRR.bias.corrected(data.del.temp.log.RR)
  
  
  # run rma
  rma.imp.logRR = rma(yi = data.del.temp.log.RR$yi, vi = rep(1, nrow(data.del.temp.log.RR)), measure = "ROM",
                      control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
  
  # extract output
  results.df$grand_mean[1] = rma.imp.logRR$beta
  results.df$ci_lb[1] = rma.imp.logRR$ci.lb
  results.df$ci_ub[1] = rma.imp.logRR$ci.ub
  results.df$tau[1] = rma.imp.logRR$tau2
  results.df$i2[1] = rma.imp.logRR$I2
  results.df$h2[1] = rma.imp.logRR$H2
  
  # return results
  return(results.df)
}
# 2b. unweighted analysis ZCOR --------------------------------------------------

imp.unweighted.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "unweighted analysis", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # calculate effect sizes
  data.del.temp.ZCOR = data.frame(escalc(ri = cor_coef, ni = SS,
                                         data = data.del, measure="ZCOR"))
  
  # run rma
  rma.imp.ZCOR = rma(yi = data.del.temp.ZCOR$yi, vi = rep(1, nrow(data.del.temp.ZCOR)), measure = "ZCOR",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = rma.imp.ZCOR$beta
  results.df$ci_lb = rma.imp.ZCOR$ci.lb
  results.df$ci_ub = rma.imp.ZCOR$ci.ub
  results.df$tau = rma.imp.ZCOR$tau2
  results.df$i2 = rma.imp.ZCOR$I2
  results.df$h2 = rma.imp.ZCOR$H2
  
  # return results
  return(results.df)
}


# xa. sample-size-weighted analysis logRR/SMDH --------------------------------------------------

imp.sample.size.weighted.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("SS","SS"), 
                          "what_is_missing" = rep(what.is.missing,2),
                          "imputation_method" = rep("sample-size-weighted analysis",2), 
                          "missing_percent" = rep(del.rate,2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # calculate effect sizes <- Hedge's d cannot be calculated without SDs
  data.del.temp.log.RR = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                           m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                           data = data.del, measure="ROM"))
  data.del.temp.log.RR = calculate.logRR.bias.corrected(data.del.temp.log.RR)
  
  # add sample-size-weights
  data.del.temp.log.RR$SS_var = as.vector(1 / (data.del.temp.log.RR$treat_SS * data.del.temp.log.RR$contr_SS) / 
                                            (data.del.temp.log.RR$treat_SS + data.del.temp.log.RR$contr_SS))
  
  # run rma
  rma.imp.logRR = rma(yi = data.del.temp.log.RR$yi, vi = data.del.temp.log.RR$SS_var, measure = "ROM",
                      control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
  
  # extract output
  results.df$grand_mean[1] = rma.imp.logRR$beta
  results.df$ci_lb[1] = rma.imp.logRR$ci.lb
  results.df$ci_ub[1] = rma.imp.logRR$ci.ub
  results.df$tau[1] = rma.imp.logRR$tau2
  results.df$i2[1] = rma.imp.logRR$I2
  results.df$h2[1] = rma.imp.logRR$H2
  
  return(results.df)
}

# 3a. mean value imputation logRR/SMDH------------------------------------------------

imp.mean.value.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("mean value imputation",2), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # add mean values
  data.del$contr_SD[is.na(data.del$contr_SD)] = mean(data.del$contr_SD, na.rm = T)
  data.del$treat_SD[is.na(data.del$treat_SD)] = mean(data.del$treat_SD, na.rm = T)
  data.del$contr_SS[is.na(data.del$contr_SS)] = mean(data.del$contr_SS, na.rm = T)
  data.del$treat_SS[is.na(data.del$treat_SS)] = mean(data.del$treat_SS, na.rm = T)
  
  # calculate effect sizes
  data.del.temp.log.RR = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                           m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                           data = data.del, measure="ROM"))
  data.del.temp.log.RR = calculate.logRR.bias.corrected(data.del.temp.log.RR)
  
  data.del.temp.SMDH = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                         m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                         data = data.del, measure="SMDH"))
  
  # run rma
  rma.imp.logRR = rma(yi = data.del.temp.log.RR$yi, vi = data.del.temp.log.RR$vi, measure = "ROM",
                      control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
  rma.imp.SMDH = rma(yi = data.del.temp.SMDH$yi, vi = data.del.temp.SMDH$vi, measure = "SMDH",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = c(rma.imp.logRR$beta, rma.imp.SMDH$beta)
  results.df$ci_lb = c(rma.imp.logRR$ci.lb, rma.imp.SMDH$ci.lb)
  results.df$ci_ub = c(rma.imp.logRR$ci.ub, rma.imp.SMDH$ci.ub)
  results.df$tau = c(rma.imp.logRR$tau2, rma.imp.SMDH$tau2)
  results.df$i2 = c(rma.imp.logRR$I2, rma.imp.SMDH$I2)
  results.df$h2 = c(rma.imp.logRR$H2, rma.imp.SMDH$H2)
  
  # return results
  return(results.df)
}

# 3b. mean value imputation ZCOR ------------------------------------------------

imp.mean.value.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "mean value imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # add mean values
  data.del$SS[is.na(data.del$SS)] = mean(data.del$SS, na.rm = T)
  
  # calculate effect sizes
  data.del.temp.ZCOR = data.frame(escalc(ri = cor_coef, ni = SS, data = data.del, measure="ZCOR"))
  
  # run rma
  rma.imp.ZCOR = rma(yi = data.del.temp.ZCOR$yi, vi = data.del.temp.ZCOR$vi, measure = "ZCOR",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = rma.imp.ZCOR$beta
  results.df$ci_lb = rma.imp.ZCOR$ci.lb
  results.df$ci_ub = rma.imp.ZCOR$ci.ub
  results.df$tau = rma.imp.ZCOR$tau2
  results.df$i2 = rma.imp.ZCOR$I2
  results.df$h2 = rma.imp.ZCOR$H2
  
  # return results
  return(results.df)
}

# 4a. median value imputation logRR/SMDH----------------------------------------------

imp.median.value.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("median value imputation",2), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # add mean values
  data.del$contr_SD[is.na(data.del$contr_SD)] = median(data.del$contr_SD, na.rm = T)
  data.del$treat_SD[is.na(data.del$treat_SD)] = median(data.del$treat_SD, na.rm = T)
  data.del$contr_SS[is.na(data.del$contr_SS)] = median(data.del$contr_SS, na.rm = T)
  data.del$treat_SS[is.na(data.del$treat_SS)] = median(data.del$treat_SS, na.rm = T)
  
  # calculate effect sizes
  data.del.temp.log.RR = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                           m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                           data = data.del, measure="ROM"))
  data.del.temp.log.RR = calculate.logRR.bias.corrected(data.del.temp.log.RR)
  
  data.del.temp.SMDH = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                         m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                         data = data.del, measure="SMDH"))
  
  # run rma
  rma.imp.logRR = rma(yi = data.del.temp.log.RR$yi, vi = data.del.temp.log.RR$vi, measure = "ROM",
                      control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
  rma.imp.SMDH = rma(yi = data.del.temp.SMDH$yi, vi = data.del.temp.SMDH$vi, measure = "SMDH",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = c(rma.imp.logRR$beta, rma.imp.SMDH$beta)
  results.df$ci_lb = c(rma.imp.logRR$ci.lb, rma.imp.SMDH$ci.lb)
  results.df$ci_ub = c(rma.imp.logRR$ci.ub, rma.imp.SMDH$ci.ub)
  results.df$tau = c(rma.imp.logRR$tau2, rma.imp.SMDH$tau2)
  results.df$i2 = c(rma.imp.logRR$I2, rma.imp.SMDH$I2)
  results.df$h2 = c(rma.imp.logRR$H2, rma.imp.SMDH$H2)
  
  # return results
  return(results.df)
}

# 4b. mean value imputation ZCOR------------------------------------------------

imp.median.value.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "median value imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # add median values
  data.del$SS[is.na(data.del$SS)] = median(data.del$SS, na.rm = T)
  
  # calculate effect sizes
  data.del.temp.ZCOR = data.frame(escalc(ri = cor_coef, ni = SS, data = data.del, measure="ZCOR"))
  
  # run rma
  rma.imp.ZCOR = rma(yi = data.del.temp.ZCOR$yi, vi = data.del.temp.ZCOR$vi, measure = "ZCOR",
                     control=list(maxiter = 600, stepadj = 0.2))
  
  # extract output
  results.df$grand_mean = rma.imp.ZCOR$beta
  results.df$ci_lb = rma.imp.ZCOR$ci.lb
  results.df$ci_ub = rma.imp.ZCOR$ci.ub
  results.df$tau = rma.imp.ZCOR$tau2
  results.df$i2 = rma.imp.ZCOR$I2
  results.df$h2 = rma.imp.ZCOR$H2
  
  # return results
  return(results.df)
}

# 5a. random sample imputation logRR/SMDH---------------------------------------------

imp.random.sample.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("mice::random sample imputation",2), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # impute missing values, post-process values to be positive
  imputed.list = list()
  
  # impute missing values, in order sd_control, sd_treat, n_control, n_treat
  columns_with_missing_data = colnames(data.del)[colSums(is.na(data.del)) > 0]
  columns_with_missing_data = columns_with_missing_data[order(match(columns_with_missing_data,
                                                                    c("contr_SD", "treat_SD", "contr_SS", "treat_SS")))]
  data.temp = data.del[,-which(names(data.del) %in% columns_with_missing_data)]
  
  for(i in 1:100){
    
    dat.imp = data.temp
    
    for(imp_step in columns_with_missing_data){
      dat.imp[,imp_step] = data.del[,which(names(data.del) == imp_step)]
      dat.imp = mice::complete(mice::mice(dat.imp, method = "sample", m = 1, printFlag=F))}
    
    imputed.list[[i]] = dat.imp
    }

  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the tau, I� and H쾤alues
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 5b. random sample imputation ZCOR---------------------------------------------

imp.random.sample.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "mice::random sample imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # impute missing values, post-process values to be positive
  data.imp.ZCOR = list()
  
  imp = mice(data.del, method = "sample", m = 100, printFlag=F)
  
  for(i in 1:100){
    data.imp.ZCOR[[i]] = mice::complete(imp, i)}
  
  # calculate effect sizes
  for(i in 1:100){
    data.imp.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = data.imp.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = data.imp.ZCOR[[i]]$yi, vi = data.imp.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the tau, I� and H쾤alues
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}

# 6a. linear prediction logRR/SMDH----------------------------------------------------

imp.linear.predict.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"),
                          "weighting_by" = c("vi","vi"),
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("mice:: linear prediction imputation", 2),
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # imputing with predictor matrix and visitSequence (contr_SD, treat_SD, contr_SS, treat_SS)
  columns_with_missing_data = colnames(data.del)[colSums(is.na(data.del)) > 0]
  columns_with_missing_data = columns_with_missing_data[order(match(columns_with_missing_data,
                                                                    c("contr_SD", "treat_SD", "contr_SS", "treat_SS")))]
  
  predictorMatrix = matrix(1, nrow= ncol(data.del), ncol = ncol(data.del))
  predictorMatrix[which(names(data.del) %in% columns_with_missing_data),
                  which(names(data.del) %in% columns_with_missing_data)] = 0
  predictorMatrix[c(4,5,6), c(3,4,5)] = 1
  
  visitSequence = match(columns_with_missing_data, names(data.del))
  
  # initial imputations
  imp.ini = mice::mice(data = data.del, method = "norm.predict", m = 100, printFlag = F,
                    predictorMatrix = predictorMatrix, 
                    visitSequence = visitSequence)
  
  # restrict imputed values to only positive
  post = imp.ini$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  
  # restricted imputations
  imp = mice::mice(data = data.del, method = "norm.predict", m = 100, printFlag = F,
                   predictorMatrix = predictorMatrix,
                   visitSequence = visitSequence,
                   post = post)
  
  # fill imputation list
  imputed.list = list()
  for(i in 1:100){
    imputed.list[[i]] = mice::complete(imp, i)}
  
    
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the tau, I� and H쾤alues
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 6b. linear prediction ZCOR ----------------------------------------------------
imp.linear.predict.ZCOR = function(data.del, what.is.missing, del.rate){
  
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "mice:: linear prediction imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # impute missing values, post-process values to be positive
  data.imp.ZCOR = list()
  
  ini.imp = mice(data.del, method = "norm.predict", m = 100, printFlag=F)
  
  # restrict imputed values to only positive
  post = ini.imp$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  
  imp = mice(data.del, method = "norm.predict", post = post, m = 100, printFlag=F)
  
  for(i in 1:100){
    data.imp.ZCOR[[i]] = mice::complete(imp, i)}
  
  # calculate effect sizes
  for(i in 1:100){
    data.imp.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = data.imp.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = data.imp.ZCOR[[i]]$yi, vi = data.imp.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the tau, I� and H쾤alues
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}

# 7a. predictive mean matching logRR/SMDH---------------------------------------------

imp.pmm.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"),
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("mice::pmm imputation",2), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # imputing with predictor matrix and visitSequence (contr_SD, treat_SD, contr_SS, treat_SS)
  columns_with_missing_data = colnames(data.del)[colSums(is.na(data.del)) > 0]
  columns_with_missing_data = columns_with_missing_data[order(match(columns_with_missing_data,
                                                                    c("contr_SD", "treat_SD", "contr_SS", "treat_SS")))]
  
  predictorMatrix = matrix(1, nrow= ncol(data.del), ncol = ncol(data.del))
  predictorMatrix[which(names(data.del) %in% columns_with_missing_data),
                  which(names(data.del) %in% columns_with_missing_data)] = 0
  predictorMatrix[c(4,5,6), c(3,4,5)] = 1
  
  visitSequence = match(columns_with_missing_data, names(data.del))
  
  # initial imputations
  imp.ini = mice::mice(data = data.del, method = "pmm", m = 100, printFlag = F,
                       predictorMatrix = predictorMatrix, 
                       visitSequence = visitSequence)
  
  # restrict imputed values to only positive
  post = imp.ini$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  
  # restricted imputations
  imp = mice::mice(data = data.del, method = "pmm", m = 100, printFlag = F,
                   predictorMatrix = predictorMatrix,
                   visitSequence = visitSequence,
                   post = post)
  
  # fill imputation list
  imputed.list = list()
  for(i in 1:100){
    imputed.list[[i]] = mice::complete(imp, i)}
  
  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the tau, I� and H쾤alues
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 7b. predictive mean matching ZCOR----------------------------------------------------
imp.pmm.ZCOR = function(data.del, what.is.missing, del.rate){
  
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "mice::pmm imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # impute missing values, post-process values to be positive
  data.imp.ZCOR = list()
  
  ini.imp = mice(data.del, method = "pmm", m = 100, printFlag=F)
  
  # restrict imputed values to only positive
  post = ini.imp$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  
  imp = mice(data.del, method = "pmm", post = post, m = 100, printFlag=F)
  
  for(i in 1:100){
    data.imp.ZCOR[[i]] = mice::complete(imp, i)}
  
  # calculate effect sizes
  for(i in 1:100){
    data.imp.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = data.imp.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = data.imp.ZCOR[[i]]$yi, vi = data.imp.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the tau, I� and H쾤alues
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}

# 8a. classification and regression trees logRR/SMDH----------------------------------

imp.cart.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"),
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("mice::cart imputation",2), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # imputing with predictor matrix and visitSequence (contr_SD, treat_SD, contr_SS, treat_SS)
  columns_with_missing_data = colnames(data.del)[colSums(is.na(data.del)) > 0]
  columns_with_missing_data = columns_with_missing_data[order(match(columns_with_missing_data,
                                                                    c("contr_SD", "treat_SD", "contr_SS", "treat_SS")))]
  
  predictorMatrix = matrix(1, nrow= ncol(data.del), ncol = ncol(data.del))
  predictorMatrix[which(names(data.del) %in% columns_with_missing_data),
                  which(names(data.del) %in% columns_with_missing_data)] = 0
  predictorMatrix[c(4,5,6), c(3,4,5)] = 1
  
  visitSequence = match(columns_with_missing_data, names(data.del))
  
  # initial imputations
  imp.ini = mice::mice(data = data.del, method = "cart", m = 100, printFlag = F,
                       predictorMatrix = predictorMatrix, 
                       visitSequence = visitSequence)
  
  # restrict imputed values to only positive
  post = imp.ini$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  
  # restricted imputations
  imp = mice::mice(data = data.del, method = "cart", m = 100, printFlag = F,
                   predictorMatrix = predictorMatrix,
                   visitSequence = visitSequence,
                   post = post)
  
  # fill imputation list
  imputed.list = list()
  for(i in 1:100){
    imputed.list[[i]] = mice::complete(imp, i)}
  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the tau, I� and H쾤alues
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 8b. classification and regression trees ZCOR----------------------------------

imp.cart.ZCOR = function(data.del, what.is.missing, del.rate){
  
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "mice::cart imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # impute missing values, post-process values to be positive
  data.imp.ZCOR = list()
  
  ini.imp = mice(data.del, method = "cart", m = 100, printFlag=F)
  
  # restrict imputed values to only positive
  post = ini.imp$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  
  imp = mice(data.del, method = "cart", post = post, m = 100, printFlag=F)
  
  for(i in 1:100){
    data.imp.ZCOR[[i]] = mice::complete(imp, i)}
  
  # calculate effect sizes
  for(i in 1:100){
    data.imp.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = data.imp.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = data.imp.ZCOR[[i]]$yi, vi = data.imp.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the tau, I� and H쾤alues
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}


# 9a. random forest -------------------------------------------------------

imp.rf.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("mice::rf imputation",2), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # imputing with predictor matrix and visitSequence (contr_SD, treat_SD, contr_SS, treat_SS)
  columns_with_missing_data = colnames(data.del)[colSums(is.na(data.del)) > 0]
  columns_with_missing_data = columns_with_missing_data[order(match(columns_with_missing_data,
                                                                    c("contr_SD", "treat_SD", "contr_SS", "treat_SS")))]
  
  predictorMatrix = matrix(1, nrow= ncol(data.del), ncol = ncol(data.del))
  predictorMatrix[which(names(data.del) %in% columns_with_missing_data),
                  which(names(data.del) %in% columns_with_missing_data)] = 0
  predictorMatrix[c(4,5,6), c(3,4,5)] = 1
  
  visitSequence = match(columns_with_missing_data, names(data.del))
  
  # initial imputations
  imp.ini = mice::mice(data = data.del, method = "rf", m = 100, printFlag = F,
                       predictorMatrix = predictorMatrix, 
                       visitSequence = visitSequence)
  
  # restrict imputed values to only positive
  post = imp.ini$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
  
  # restricted imputations
  imp = mice::mice(data = data.del, method = "rf", m = 100, printFlag = F,
                   predictorMatrix = predictorMatrix,
                   visitSequence = visitSequence,
                   post = post)
  
  # fill imputation list
  imputed.list = list()
  for(i in 1:100){
    imputed.list[[i]] = mice::complete(imp, i)}
  
  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the tau, I� and H쾤alues
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 9b. random forest -------------------------------------------------------

imp.rf.ZCOR = function(data.del, what.is.missing, del.rate){
  
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "mice::rf imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # impute missing values, post-process values to be positive
  data.imp.ZCOR = list()
  
  ini.imp = mice(data.del, method = "rf", m = 100, printFlag=F)
  
  # restrict imputed values to only positive
  post = ini.imp$post
  post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
  post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(5, Inf))"
  
  imp = mice(data.del, method = "rf", post = post, m = 100, printFlag=F)
  
  for(i in 1:100){
    data.imp.ZCOR[[i]] = mice::complete(imp, i)}
  
  # calculate effect sizes
  for(i in 1:100){
    data.imp.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = data.imp.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = data.imp.ZCOR[[i]]$yi, vi = data.imp.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the tau, I� and H쾤alues
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}

# 10a. bayes pmm - mi ------------------------------------------------------

imp.mi.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("mi::bayes pmm imputation", 2),
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  
  # create a mi::missing data frame
  data.del.temp = missing_data.frame(data.del, favor_positive = T)  
  
  # define the SDs and SSs should be positive
  #data.del.temp = change(data.del.temp, y = "contr_SD", what = "type", to = "positive-continuous")
  #data.del.temp = change(data.del.temp, y = "treat_SD", what = "type", to = "positive-continuous")
  #data.del.temp = change(data.del.temp, y = "contr_SS", what = "type", to = "positive-continuous")
  #data.del.temp = change(data.del.temp, y = "treat_SS", what = "type", to = "positive-continuous")
  
  # impute missing values
  imp = mi(data.del.temp, n.chains = 100)
  
  # collect imputed datasets
  imputed.list = mi::complete(imp, m = 100)
  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 10b. bayes pmm - mi ------------------------------------------------------

imp.mi.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "mi::bayes pmm imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # create a mi::missing data frame
  data.del.temp = missing_data.frame(data.del)
  
  # define the SDs and SSs should be positive
  data.del.temp = change(data.del.temp, y = "SS", what = "type", to = "positive-continuous")
  
  # impute missing values
  imp = mi(data.del.temp, n.chains = 100)
  
  # collect imputed datasets
  imputed.list.ZCOR = mi::complete(imp, m = 100)
  
  for(i in 1:100){
    # calculate effect sizes
    imputed.list.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = imputed.list.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = imputed.list.ZCOR[[i]]$yi, vi = imputed.list.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control = list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}


# 11a. bootstrap EBM -------------------------------------------------------

boot.ebm.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("Amelia::boot ebm imputation",2), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # determine if contr_SS and treat_SS are perfectly correlated
  cor.SS = ifelse(cor(data.del$contr_SS, data.del$treat_SS, use= "pairwise.complete.obs") == 1, "yes", "no")
  
  # create list to store imputation results
  imputed.list = list()
  
  # impute columns in order contr_SD, treat_SD, contr_SS, treat_SS
  missing_data_cols = colnames(data.del)[colSums(is.na(data.del)) > 0]
  
  # run 100 imputations
  for(i in 1:100){
    
    # create temporary data subset, exclude treat_SS if this is perfectly correlated with contr_SS
    data.temp = data.del[, -(which(names(data.del) %in% missing_data_cols))]
    if("treat_SS" %in% names(data.temp) & cor.SS == "yes"){data.temp$treat_SS = NULL}
    
    for(to.imput in missing_data_cols){
      
      data.temp[,to.imput] = data.del[,to.imput]
      
      if(to.imput %in% c("contr_SD", "treat_SD")){
        bounds = matrix(c(which(names(data.temp) == to.imput), 0.01, 10), 
                        nrow = 1, ncol = 3)
        }else{
        bounds = matrix(c(which(names(data.temp) == to.imput), 3, 100), 
                        nrow = 1, ncol = 3)}
      
      data.temp =  amelia(data.temp, m = 1,  bounds = bounds, p2s = 0)$imputations[[1]]
      }
    
    # add treat_SS if this was deleted and not imputed
    if(!("treat_SS" %in% names(data.temp))){
      data.temp$treat_SS = data.del$treat_SS}
    
    imputed.list[[i]] = data.temp
  }
  
  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}


# 11b. bootstrap EBM -------------------------------------------------------

boot.ebm.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "Amelia::boot ebm imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  
  # define the SDs and SSs should be positive
  bounds = matrix(c(which(names(data.del) == "SS"), 4, 100), 
                  nrow = 1, ncol = 3)
  
  # impute missing values
  imp =  amelia(data.del, m = 100,  bounds = bounds, p2s = 0)
  
  # collect imputed datasets
  imputed.list.ZCOR = list()
  for(i in 1:100){
    imputed.list.ZCOR[[i]] = imp$imputations[[i]]}
  
  for(i in 1:100){
    # calculate effect sizes
    imputed.list.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = imputed.list.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = imputed.list.ZCOR[[i]]$yi, vi = imputed.list.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}

# 12a. random forest (missForest) ------------------------------------------

miss.rf.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("missForest::rf imputation", 2),
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # impute and collect imputed datasets
  imputed.list = list()
  for(i in 1:100){
    imputed.list[[i]] = missForest(data.del)$ximp}
  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 12b. random forest (missForest) ------------------------------------------

miss.rf.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "missForest::rf imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # impute and collect imputed datasets
  imputed.list.ZCOR = list()
  for(i in 1:100){
    imputed.list.ZCOR[[i]] = missForest(data.del)$ximp}
  
  for(i in 1:100){
    # calculate effect sizes
    imputed.list.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = imputed.list.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = imputed.list.ZCOR[[i]]$yi, vi = imputed.list.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}

# 13a. bootstrap pmm - Hmisc -------------------------------------------------------

boot.pmm.logRR.SMDH = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("logRR","SMDH"), 
                          "weighting_by" = c("vi","vi"), 
                          "what_is_missing" = rep(what.is.missing, 2),
                          "imputation_method" = rep("Hmisc::boot pmm imputation"), 
                          "missing_percent" = rep(del.rate, 2),
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1:2,]
  
  # impute and datasets
  imp = aregImpute(~ contr_mean + treat_mean + contr_SD + treat_SD + contr_SS + treat_SS + cov, 
                   data = data.del, n.impute = 100, pr = F )
  
  # store imputed datasets
  imputed.list = list()
  for(i in 1:100){
    imputed.list[[i]] = as.data.frame(impute.transcan(imp, imputation = i, data = data.del, list.out=TRUE, pr=FALSE, check=FALSE))}
  
  # calculate effect sizes
  data.del.temp.log.RR = imputed.list
  data.del.temp.SMDH = imputed.list
  
  for(i in 1:100){
    # calculate effect sizes
    data.del.temp.log.RR[[i]] =  data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                   m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                   data = data.del.temp.log.RR[[i]], measure="ROM"))
    data.del.temp.log.RR[[i]] = calculate.logRR.bias.corrected(data.del.temp.log.RR[[i]])
    
    data.del.temp.SMDH[[i]] = data.frame(escalc(m2i = contr_mean, sd2i = contr_SD, n2i = contr_SS,
                                                m1i = treat_mean, sd1i = treat_SD, n1i = treat_SS,
                                                data = data.del.temp.SMDH[[i]], measure="SMDH"))
  }
  
  # create lists to store the rma models
  rma.imp.logRR = list()
  rma.imp.SMDH = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.logRR[[i]] = rma(yi = data.del.temp.log.RR[[i]]$yi, vi = data.del.temp.log.RR[[i]]$vi, 
                             measure = "ROM", control=list(maxiter = 600, stepadj = 0.2, tol = 1e-10))
    rma.imp.SMDH[[i]] = rma(yi = data.del.temp.SMDH[[i]]$yi, vi = data.del.temp.SMDH[[i]]$vi, 
                            measure = "SMDH", control=list(maxiter = 600, stepadj = 0.2))
  }
  
  # pool models
  rma.imp.logRR.pooled = mice::pool(as.mira(rma.imp.logRR))
  rma.imp.SMDH.pooled = mice::pool(as.mira(rma.imp.SMDH))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.logRR = c(); tau.imp.SMDH = c()
  i2.imp.logRR = c(); i2.imp.SMDH = c()
  h2.imp.logRR = c(); h2.imp.SMDH = c()
  
  for(i in 1:100){
    tau.imp.logRR = c(tau.imp.logRR, rma.imp.logRR[[i]]$tau2)
    tau.imp.SMDH = c(tau.imp.SMDH, rma.imp.SMDH[[i]]$tau2)
    
    i2.imp.logRR = c(i2.imp.logRR, rma.imp.logRR[[i]]$I2)
    i2.imp.SMDH = c(i2.imp.SMDH, rma.imp.SMDH[[i]]$I2)
    
    h2.imp.logRR = c(h2.imp.logRR, rma.imp.logRR[[i]]$H2)
    h2.imp.SMDH = c(h2.imp.SMDH, rma.imp.SMDH[[i]]$H2)
  }
  
  # extract output
  results.df$grand_mean = c(summary(rma.imp.logRR.pooled)[,"est"], summary(rma.imp.SMDH.pooled)[,"est"])
  results.df$ci_lb = c(summary(rma.imp.logRR.pooled)[,"lo 95"], summary(rma.imp.SMDH.pooled)[,"lo 95"])
  results.df$ci_ub = c(summary(rma.imp.logRR.pooled)[,"hi 95"], summary(rma.imp.SMDH.pooled)[,"hi 95"])
  results.df$tau = c(mean(tau.imp.logRR), mean(tau.imp.SMDH))
  results.df$i2 = c(mean(i2.imp.logRR), mean(i2.imp.SMDH))
  results.df$h2 = c(mean(h2.imp.logRR), mean(h2.imp.SMDH))
  
  # return results
  return(results.df)
}

# 13b. bootstrap pmm - Hmisc -------------------------------------------------------

boot.pmm.ZCOR = function(data.del, what.is.missing, del.rate){
  
  # create dataframe to store results
  results.df = data.frame("effect_size" = c("ZCOR"), 
                          "weighting_by" = c("SS"), 
                          "what_is_missing" = "SS",
                          "imputation_method" = "Hmisc::boot pmm imputation", 
                          "missing_percent" = del.rate,
                          "grand_mean" = NA, "ci_lb" = NA, "ci_ub" = NA, 
                          "tau" = NA, "i2" = NA, "h2" = NA)[1,]
  
  # impute and datasets
  imp = aregImpute(~ cor_coef + SS + cov, 
                   data = data.del, n.impute = 100)
  
  # store imputed datasets
  imputed.list.ZCOR = list()
  for(i in 1:100){
    imputed.list.ZCOR[[i]] = as.data.frame(impute.transcan(imp, imputation = i, data = data.del, list.out=TRUE, pr=FALSE, check=FALSE))}
  
  # calculate effect sizes
  for(i in 1:100){
    imputed.list.ZCOR[[i]] =  data.frame(escalc(ri = cor_coef, ni = SS, data = imputed.list.ZCOR[[i]], measure="ZCOR"))}
  
  # create lists to store the rma models
  rma.imp.ZCOR = list()
  
  # run rmas
  for(i in 1:100){
    rma.imp.ZCOR[[i]] = rma(yi = imputed.list.ZCOR[[i]]$yi, vi = imputed.list.ZCOR[[i]]$vi, 
                            measure = "ZCOR", control=list(maxiter = 600, stepadj = 0.2))}
  
  # pool models
  rma.imp.ZCOR.pooled = mice::pool(as.mira(rma.imp.ZCOR))
  
  # get the mean tau, I� and H쾤alues from 100 imputations
  tau.imp.ZCOR = c()
  i2.imp.ZCOR = c()
  h2.imp.ZCOR = c()
  
  for(i in 1:100){
    tau.imp.ZCOR = c(tau.imp.ZCOR, rma.imp.ZCOR[[i]]$tau2)
    i2.imp.ZCOR = c(i2.imp.ZCOR, rma.imp.ZCOR[[i]]$I2)
    h2.imp.ZCOR = c(h2.imp.ZCOR, rma.imp.ZCOR[[i]]$H2)}
  
  # extract output
  results.df$grand_mean = summary(rma.imp.ZCOR.pooled)[,"est"]
  results.df$ci_lb = summary(rma.imp.ZCOR.pooled)[,"lo 95"]
  results.df$ci_ub = summary(rma.imp.ZCOR.pooled)[,"hi 95"]
  results.df$tau = mean(tau.imp.ZCOR)
  results.df$i2 = mean(i2.imp.ZCOR)
  results.df$h2 = mean(h2.imp.ZCOR)
  
  # return results
  return(results.df)
}


###########################
# resterampe


# former functions used for imputations - now replaced with shorter visitSequence - version
# # impute missing values, in case that variables are correlated do the imputation step by step
# imp.ini = tryCatch({
#   mice::mice(data = data.del, method = "norm.predict", m = 100, printFlag=F)},
#   error = function(e){
#     predictorMatrix = matrix(1, nrow= ncol(dat.del), ncol = ncol(dat.del))
#     missing_columns = c(which(colSums(is.na(data.del)) > 0))
#     predictorMatrix[missing_columns, missing_columns] = 0
#     imp.ini = mice::mice(data = data.del, method = "norm.predict", m = 100, printFlag=F, predictorMatrix = predictorMatrix)
#     return(imp.ini)})
# 
# # restrict imputed values to only positive
# post = imp.ini$post
# post["contr_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
# post["treat_SD"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(0.1, 10))"
# post["contr_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
# post["treat_SS"] = "imp[[j]][,i] <- squeeze(imp[[j]][,i],c(3, Inf))"
# 
# # now the imputations with restrictions
# imp = tryCatch({
#   mice::mice(data = data.del, method = "norm.predict", m = 100, printFlag=F, post = post)},
#   error = function(e){
#     predictorMatrix = matrix(1, nrow= ncol(dat.del), ncol = ncol(dat.del))
#     missing_columns = c(which(colSums(is.na(data.del)) > 0))
#     predictorMatrix[missing_columns, missing_columns] = 0
#     imp = mice::mice(data = data.del, method = "norm.predict", m = 100, post = post, printFlag=F, predictorMatrix = predictorMatrix)
#     return(imp)})
# 
# 
# for(i in 1:100){
#   imputed.list[[i]] = mice::complete(imp, i)}
# 
# # check if imputation was sucessfull
# if(length(which(is.na(imputed.list[[i]]$contr_SD))) > 0 | length(which(is.na(imputed.list[[i]]$contr_SS))) > 0){
#   data.temp = data.del
#   columns_with_missing_data = names(imp$nmis)[which(imp$nmis >0)]
#   data.temp = data.temp[,-which(names(data.temp) %in% columns_with_missing_data)]
#   
#   for(i in 1:100){
#     
#     dat.imp = data.temp
#     
#     for(imp_step in columns_with_missing_data){
#       dat.imp[,imp_step] = data.del[,which(names(data.del) == imp_step)]
#       dat.imp = mice::complete(mice::mice(dat.imp, method = "norm.predict",post = post, m = 1, printFlag=F))}
#     
#     imputed.list[[i]] = dat.imp
#   }
# }





# ebm
# # define the SDs and SSs should be positive
# bounds = matrix(c(which(names(data.del) == "contr_SD"), which(names(data.del) == "treat_SD"),
#                   which(names(data.del) == "contr_SS"), which(names(data.del) == "treat_SS"),
#                   0.01, 0.01, 0.01, 0.01, 10, 10, 10, 10), 
#                 nrow = 4, ncol = 3)
# 
# 
# # impute missing values
# imp =  amelia(data.del, m = 100,  bounds = bounds, p2s = 0)
# 
# # collect imputed datasets
# imputed.list = list()
# for(i in 1:100){
#   imputed.list[[i]] = imp$imputations[[i]]}
