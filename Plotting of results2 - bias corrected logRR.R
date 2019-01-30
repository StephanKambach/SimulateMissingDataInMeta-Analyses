##########################################################
# title: 'Missing SD and SS in meta-analysis: Plot results
# by: Stephan Kambach
# Date: 04.2018
##########################################################

# libraries and data
rm(list=ls())
gc()

library(ggplot2)
library(data.table)
library(cowplot)
library(extrafont)
library(XLConnect)  
library(scales)
font_import(pattern ="[T/t]imes")
y
loadfonts(device="win")
fonts()

setwd("C:/Users/kambach/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("C:/Users/Agando/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("D:/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("C:/Users/localadmin/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")
setwd("C:/Users/sk85xupa/Dropbox/current_tasks/SESYNC_multiple_imputation/new analyses after MEE revision/Simulations")


# Figure 1 - literature review --------------------------------------------
wk = loadWorkbook("data\\Review - coding table.xlsx") 
dat.review = data.table(readWorksheet(wk, sheet=1))
dat.review = data.table(read.table("data\\Review - coding table.txt", header = T, sep = "\t", stringsAsFactors = F, dec= ","))
str(dat.review)
dat.review


# 1a. Number of studies with missing data per effect size
dat.fig1.a = dat.review

table(dat.fig1.a$TYPE_EFFECT_SIZE) # check that all cases of multiple effect sizes are dealt with in the below section

# reformat in case studies looked at multiple effect sizes
for(i in 1:nrow(dat.fig1.a)){
  if(dat.fig1.a$TYPE_EFFECT_SIZE[i] == "mean difference and correlation coefficient"){
    dat.fig1.a = rbind(dat.fig1.a,dat.fig1.a[i,])
    dat.fig1.a$TYPE_EFFECT_SIZE[i] = "standardized mean difference"
    dat.fig1.a$TYPE_EFFECT_SIZE[nrow(dat.fig1.a)] = "correlation coef"}
  
  if(dat.fig1.a$TYPE_EFFECT_SIZE[i] == "response ratio and mean difference"){
    dat.fig1.a = rbind(dat.fig1.a,dat.fig1.a[i,])
    dat.fig1.a$TYPE_EFFECT_SIZE[i] = "response ratio"
    dat.fig1.a$TYPE_EFFECT_SIZE[nrow(dat.fig1.a)] = "standardized mean difference"}
  
  if(dat.fig1.a$TYPE_EFFECT_SIZE[i] == "response ratio and correlation coef"){
    dat.fig1.a = rbind(dat.fig1.a,dat.fig1.a[i,])
    dat.fig1.a$TYPE_EFFECT_SIZE[i] = "response ratio"
    dat.fig1.a$TYPE_EFFECT_SIZE[nrow(dat.fig1.a)] = "correlation coef"}
  
}

table(dat.fig1.a$TYPE_EFFECT_SIZE) # should only have the three basic effect sizes

# replace maybes with yes
table(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED)
dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED[which(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED == "maybe - only complete cases included")] = "yes"
dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED[which(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED == "yes - both")] = "yes"
dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED[which(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED == "yes - sample sizes")] = "yes"
dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED[which(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED == "yes - variances")] = "yes"
table(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED)

# replace NA with "not clear"
dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED[which(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED == "NA")] = "unclear"
dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED[which(is.na(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED))] = "unclear"
dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED = factor(dat.fig1.a$MISSING_VARIANCES_ENCOUNTERED,
                                                  levels = c("no", "unclear", "yes"))

# replace effect size names
dat.fig1.a$TYPE_EFFECT_SIZE[which(dat.fig1.a$TYPE_EFFECT_SIZE == "response ratio")] = "response\nratio"
dat.fig1.a$TYPE_EFFECT_SIZE[which(dat.fig1.a$TYPE_EFFECT_SIZE == "standardized mean difference")] = "standardized\nmean difference"
dat.fig1.a$TYPE_EFFECT_SIZE[which(dat.fig1.a$TYPE_EFFECT_SIZE == "correlation coef")] = "correlation\ncoefficient"

# check all conversions
levels(factor(dat.fig1.a$TYPE_EFFECT_SIZE))
levels(factor(dat.fig1.a$TYPE_EFFECT_SIZE))


svg(filename = "results/fig1a.svg", height = 4, width = 6)
ggplot(data = dat.fig1.a) +
  geom_bar(aes(x = TYPE_EFFECT_SIZE,fill = MISSING_VARIANCES_ENCOUNTERED), color = "black") +
  scale_fill_manual(name = "", values = c("#6FB76F","white", "#fd8d3c")) +
  ggtitle("Did ecological meta-analyses encounter \nmissing SDs or sample sizes?") +
  xlab("Effect size") + ylab("Number of studies") +
  theme(legend.position = "left", text = element_text(family = "Times New Roman")) + 
  ylim(c(0, max(table(dat.fig1.a$TYPE_EFFECT_SIZE))))
graphics.off()

# 1b. treatment of missing information
dat.fig1.b = dat.fig1.a
dat.fig1.b = dat.fig1.b[MISSING_VARIANCES_ENCOUNTERED %in% "yes"]

# renaming
table(dat.fig1.b$TREATMENT_OF_MISSING_VALUES)
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES %in% c("bayesian imputation","multiple imputation"))] = "Bayes or multiple imputation"
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES  %in% c("maybe - only complete cases included","complete case analysis","complete case and sensitivity analysis" ))] = "Complete case analysis"
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES  %in% c("imputation + sensitivity analysis", "imputation"))] = "Imputation"
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES  %in% c("approximated variances from sample sizes"))] = "Variances approximated from sample sizes"
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES  %in% c("complete case + unweighted analysis"))] = "Complete case + unweighted analysis"
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES  %in% c("different effect size used"))] = "Different effect size applied"
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES  %in% c("different weighting scheme"))] = "Different weighting scheme"
dat.fig1.b$TREATMENT_OF_MISSING_VALUES[which(dat.fig1.b$TREATMENT_OF_MISSING_VALUES  %in% c("unweighted analysis","unweighted + sensitivity analysis"))] = "Unweighted analysis"

table(dat.fig1.b$TREATMENT_OF_MISSING_VALUES)

dat.fig1.b$TREATMENT_OF_MISSING_VALUES = factor(dat.fig1.b$TREATMENT_OF_MISSING_VALUES,
                                                levels = c("Bayes or multiple imputation",
                                                           "Imputation",
                                                           "Variances approximated from sample sizes",
                                                           "Different weighting scheme",
                                                           "Different effect size applied",
                                                           "Unweighted analysis",
                                                           "Complete case + unweighted analysis",
                                                           "Complete case analysis"))

svg(filename = "results/fig1b.svg", height = 4, width = 8)
ggplot(data = dat.fig1.b) +
  geom_bar(aes(x = TYPE_EFFECT_SIZE,fill = TREATMENT_OF_MISSING_VALUES), color = "black") +
  scale_fill_manual(name = "", values = c('#B9C8DF','#7291BF','#F7FFD7','#FBDD8C','#CCA664','#FDA500','#FF0000','#CF1A1A')) +
  ggtitle("How did ecological meta-analyses treat\nmissing SD or sample size information?") +
  xlab("Effect size") + ylab("Number of studies") +
  theme(legend.position = "right", text = element_text(family = "Times New Roman")) + 
  ylim(c(0, max(table(dat.fig1.a$TYPE_EFFECT_SIZE))))
graphics.off()

# 1c. timeline
dat.fig1.c = dat.review
dat.fig1.c = dat.fig1.c[,.(prop.complete.cases.per.year = length(which(TREATMENT_OF_MISSING_VALUES %in% c("complete case analysis","complete case + unweighted analysis", "maybe - only complete cases included"))) ,
                           prop.imputation.per.year = length(which(TREATMENT_OF_MISSING_VALUES %in% c("imputation","imputation + sensitivity analysis", "multiple imputation", "bayesian imputation"))),
                           nr.studies = length(ID)),
                        by = (YEAR)]

svg(filename = "results/fig1c.svg", height = 3, width = 11.5)
ggplot(data = dat.fig1.c) +
  geom_line(aes(x = YEAR, y = prop.complete.cases.per.year), color = "#a63603", size= 1) +
  geom_line(aes(x = YEAR, y = prop.imputation.per.year), color = "#feedde", size= 1) +
  geom_line(aes(x = YEAR, y = nr.studies), color = "grey", size= 1) +
  
  geom_point(aes(x = YEAR, y = prop.complete.cases.per.year), color = "black", size= 7) +
  geom_point(aes(x = YEAR, y = prop.imputation.per.year), color = "black", size= 7) +
  geom_point(aes(x = YEAR, y = nr.studies), color = "black", size= 7) +
  
  geom_point(aes(x = YEAR, y = prop.complete.cases.per.year), color = "#a63603", size= 5) +
  geom_point(aes(x = YEAR, y = prop.imputation.per.year), color = "#feedde", size= 5) +
  geom_point(aes(x = YEAR, y = nr.studies), color = "grey", size= 5) +
  
  scale_x_continuous(breaks = c(2000:2018))+ 
  xlab("Year") + ylab("Number of studies") +
  theme(legend.position = "right", text = element_text(family = "Times New Roman"),
        axis.text.x=element_text(angle = -90, vjust = 0.5)) 
graphics.off()


#--------------------
# Figure 2 - read data, MCAR,. no correlations ---------------------------------------
getwd()
list.files("intermediate results2 - bias corrected logRR")
pos = grep('fig2', list.files("intermediate results2 - bias corrected logRR"))
dat.fig2 = read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[1]], sep = ""))
for(i in 2:length(pos)){
  dat.fig2 = rbind(dat.fig2,
                   read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[i]], sep = "")))
}
dat.fig2 = data.table(dat.fig2)

# Figure 2 - arrange the correct estimates from the full model -----------------------

# grand mean
dat.fig2$grand_mean_full_data = NA

dat.fig2$grand_mean_full_data[which(dat.fig2$effect_size == "logRR")] = 
  dat.fig2$grand_mean[which(dat.fig2$effect_size == "logRR" & 
                              dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$grand_mean_full_data[which(dat.fig2$effect_size == "SMDH")] = 
  dat.fig2$grand_mean[which(dat.fig2$effect_size == "SMDH" & 
                              dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$grand_mean_full_data[which(dat.fig2$effect_size == "ZCOR")] = 
  dat.fig2$grand_mean[which(dat.fig2$effect_size == "ZCOR" & 
                              dat.fig2$imputation_method == "full dataset analysis")]

# ci_lb
dat.fig2$ci_lb_full_data = NA

dat.fig2$ci_lb_full_data[which(dat.fig2$effect_size == "logRR")] = 
  dat.fig2$ci_lb[which(dat.fig2$effect_size == "logRR" & 
                         dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$ci_lb_full_data[which(dat.fig2$effect_size == "SMDH")] = 
  dat.fig2$ci_lb[which(dat.fig2$effect_size == "SMDH" & 
                         dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$ci_lb_full_data[which(dat.fig2$effect_size == "ZCOR")] = 
  dat.fig2$ci_lb[which(dat.fig2$effect_size == "ZCOR" & 
                         dat.fig2$imputation_method == "full dataset analysis")]

# ci_ub
dat.fig2$ci_ub_full_data = NA

dat.fig2$ci_ub_full_data[which(dat.fig2$effect_size == "logRR")] = 
  dat.fig2$ci_ub[which(dat.fig2$effect_size == "logRR" & 
                         dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$ci_ub_full_data[which(dat.fig2$effect_size == "SMDH")] = 
  dat.fig2$ci_ub[which(dat.fig2$effect_size == "SMDH" & 
                         dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$ci_ub_full_data[which(dat.fig2$effect_size == "ZCOR")] = 
  dat.fig2$ci_ub[which(dat.fig2$effect_size == "ZCOR" & 
                         dat.fig2$imputation_method == "full dataset analysis")]

# tau
dat.fig2$tau_full_data = NA

dat.fig2$tau_full_data[which(dat.fig2$effect_size == "logRR")] = 
  dat.fig2$tau[which(dat.fig2$effect_size == "logRR" & 
                       dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$tau_full_data[which(dat.fig2$effect_size == "SMDH")] = 
  dat.fig2$tau[which(dat.fig2$effect_size == "SMDH" & 
                       dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$tau_full_data[which(dat.fig2$effect_size == "ZCOR")] = 
  dat.fig2$tau[which(dat.fig2$effect_size == "ZCOR" & 
                       dat.fig2$imputation_method == "full dataset analysis")]

# i2
dat.fig2$i2_full_data = NA

dat.fig2$i2_full_data[which(dat.fig2$effect_size == "logRR")] = 
  dat.fig2$i2[which(dat.fig2$effect_size == "logRR" & 
                      dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$i2_full_data[which(dat.fig2$effect_size == "SMDH")] = 
  dat.fig2$i2[which(dat.fig2$effect_size == "SMDH" & 
                      dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$i2_full_data[which(dat.fig2$effect_size == "ZCOR")] = 
  dat.fig2$i2[which(dat.fig2$effect_size == "ZCOR" & 
                      dat.fig2$imputation_method == "full dataset analysis")]

# h2
dat.fig2$h2_full_data = NA

dat.fig2$h2_full_data[which(dat.fig2$effect_size == "logRR")] = 
  dat.fig2$h2[which(dat.fig2$effect_size == "logRR" & 
                      dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$h2_full_data[which(dat.fig2$effect_size == "SMDH")] = 
  dat.fig2$h2[which(dat.fig2$effect_size == "SMDH" & 
                      dat.fig2$imputation_method == "full dataset analysis")]
dat.fig2$h2_full_data[which(dat.fig2$effect_size == "ZCOR")] = 
  dat.fig2$h2[which(dat.fig2$effect_size == "ZCOR" & 
                      dat.fig2$imputation_method == "full dataset analysis")]

# Figure 2 - further data preparation ------------------------------------------------

# remove the full data as an impuation methods
dat.fig2 = dat.fig2[imputation_method != "full dataset analysis"]

# re-name imputation methods
dat.fig2$imputation_method = as.character(dat.fig2$imputation_method)
dat.fig2$imputation_method[dat.fig2$imputation_method == "complete case analysis"] = "Complete-case\nanalysis"
dat.fig2$imputation_method[dat.fig2$imputation_method == "unweighted analysis"] = "Unweighted\nanalysis"
dat.fig2$imputation_method[dat.fig2$imputation_method == "sample-size-weighted analysis"] = "Sample-size-weighted\nanalysis"
dat.fig2$imputation_method[dat.fig2$imputation_method == "mean value imputation"] = "Mean value\nimputation"
dat.fig2$imputation_method[dat.fig2$imputation_method == "median value imputation"] = "Median value\nimputation"
dat.fig2$imputation_method[dat.fig2$imputation_method == "mice::random sample imputation"] = "mice: random sample\nimputation"
dat.fig2$imputation_method[dat.fig2$imputation_method == "mice:: linear prediction imputation"] = "mice: prediction from\n linear regression"
dat.fig2$imputation_method[dat.fig2$imputation_method == "mice::pmm imputation"] = "mice: predictive\nmean matching"
dat.fig2$imputation_method[dat.fig2$imputation_method == "mice::cart imputation"] = "mice: classification \nand regression trees"
dat.fig2$imputation_method[dat.fig2$imputation_method == "mice::rf imputation"] = "mice: random forest"
dat.fig2$imputation_method[dat.fig2$imputation_method == "mi::bayes pmm imputation"] = "mi: Bayes predictive\n mean matching"
dat.fig2$imputation_method[dat.fig2$imputation_method == "Amelia::boot ebm imputation"] = "Amelia: bootstrap \nexpectation maximization"
dat.fig2$imputation_method[dat.fig2$imputation_method == "missForest::rf imputation"] = "missForest: non-parametric\n random forest"
dat.fig2$imputation_method[dat.fig2$imputation_method == "Hmisc::boot pmm imputation"] = "Hmisc: additive regression and\nbootstrap predictive mean matching"

# re-order the factor levels of imputation_method for plotting
unique(dat.fig2$imputation_method)

dat.fig2$imputation_method = factor(dat.fig2$imputation_method, 
                                    levels = unique(dat.fig2$imputation_method))

levels(dat.fig2$imputation_method)

# add rows with sample-size-weighted analyses in order get 14 plots for every effect size
dat.fig2 = rbind(dat.fig2,data.frame("effect_size" = c("logRR", "logRR", "SMDH", "SMDH", "SMDH", "ZCOR"),
                                     "what_is_missing" = c("SD & SS", "SS", "SD", "SD & SS", "SS", "SS"),
                                     "imputation_method" = "Sample-size-weighted\nanalysis"), fill = T)

# Figure 2 - plotting -----------------------------------------------------

# I need to do individual plots for every effect size
create.plot.for.fig2 = function(data.temp){
  
  dist.mean.and.ci = data.temp$grand_mean_full_data[1] - data.temp$ci_lb_full_data[1]
  plot2 = ggplot(data = data.temp) +
    
    geom_ribbon(aes(x = missing_percent, ymin = ci_lb, ymax = ci_ub, fill = imputation_method, color = imputation_method), alpha = 0.35) +
    geom_line(aes(x = missing_percent, y = grand_mean, color = imputation_method), size = 1) + 
    
    #geom_segment(aes(x = missing_percent, y = ci_lb, color = imputation_method ), size = 1, alpha = 0.35)+
    #geom_line(aes(x = missing_percent, y = ci_ub, color = imputation_method ), size = 1, alpha = 0.35)+
    
    geom_line(aes(x = missing_percent, y = grand_mean_full_data), color = "black") + 
    geom_line(aes(x = missing_percent, y = ci_lb_full_data), color = "black", linetype = "twodash") + 
    geom_line(aes(x = missing_percent, y = ci_ub_full_data), color = "black", linetype = "twodash") +
    
    facet_grid(imputation_method ~ .) +
    scale_x_reverse(breaks = c(0.1, 0.9)) +
    coord_flip(ylim = c(data.temp$grand_mean_full_data[1] - (3*dist.mean.and.ci),
                        data.temp$grand_mean_full_data[1] + (3*dist.mean.and.ci)),
               xlim = c(0.1, 0.9)) +
    
    xlab("") + ylab("") + ggtitle("") +
    
    theme(strip.text.y = element_text(angle = 0), 
          legend.position = "none",
          text = element_text(family = "Times New Roman"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  return(plot2)
}

# do the actual plotting and combine the figures later by hand

svg(filename = "results2 - bias corrected logRR/Figure2a.svg", height = 16, width = 5)
create.plot.for.fig2(data.temp = dat.fig2[effect_size == "logRR" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure2b.svg", height = 16, width = 5)
create.plot.for.fig2(data.temp = dat.fig2[effect_size == "logRR" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure2c.svg", height = 16, width = 5)
create.plot.for.fig2(data.temp = dat.fig2[effect_size == "logRR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure2d.svg", height = 16, width = 5)
create.plot.for.fig2(data.temp = dat.fig2[effect_size == "SMDH" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure2e.svg", height = 16, width = 5)
create.plot.for.fig2(data.temp = dat.fig2[effect_size == "SMDH" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure2f.svg", height = 16, width = 5)
create.plot.for.fig2(data.temp = dat.fig2[effect_size == "SMDH" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure2g.svg", height = 16, width = 5)
create.plot.for.fig2(data.temp = dat.fig2[effect_size == "ZCOR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()


#--------------------
# Figure 3 - read data, MCAR,. no correlations ---------------------------------------
getwd()
list.files("intermediate results2 - bias corrected logRR")
pos = grep('fig3', list.files("intermediate results2 - bias corrected logRR"))
dat.fig3 = read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[1]], sep = ""))
for(i in 2:length(pos)){
  dat.fig3 = rbind(dat.fig3,
                   read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[i]], sep = "")))
}
dat.fig3 = data.table(dat.fig3)

# Figure 3- arrange the correct estimates from the full model -----------------------

# grand mean
dat.fig3$grand_mean_full_data = NA

dat.fig3$grand_mean_full_data[which(dat.fig3$effect_size == "logRR")] = 
  dat.fig3$grand_mean[which(dat.fig3$effect_size == "logRR" & 
                              dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$grand_mean_full_data[which(dat.fig3$effect_size == "SMDH")] = 
  dat.fig3$grand_mean[which(dat.fig3$effect_size == "SMDH" & 
                              dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$grand_mean_full_data[which(dat.fig3$effect_size == "ZCOR")] = 
  dat.fig3$grand_mean[which(dat.fig3$effect_size == "ZCOR" & 
                              dat.fig3$imputation_method == "full dataset analysis")]

# ci_lb
dat.fig3$ci_lb_full_data = NA

dat.fig3$ci_lb_full_data[which(dat.fig3$effect_size == "logRR")] = 
  dat.fig3$ci_lb[which(dat.fig3$effect_size == "logRR" & 
                         dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$ci_lb_full_data[which(dat.fig3$effect_size == "SMDH")] = 
  dat.fig3$ci_lb[which(dat.fig3$effect_size == "SMDH" & 
                         dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$ci_lb_full_data[which(dat.fig3$effect_size == "ZCOR")] = 
  dat.fig3$ci_lb[which(dat.fig3$effect_size == "ZCOR" & 
                         dat.fig3$imputation_method == "full dataset analysis")]

# ci_ub
dat.fig3$ci_ub_full_data = NA

dat.fig3$ci_ub_full_data[which(dat.fig3$effect_size == "logRR")] = 
  dat.fig3$ci_ub[which(dat.fig3$effect_size == "logRR" & 
                         dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$ci_ub_full_data[which(dat.fig3$effect_size == "SMDH")] = 
  dat.fig3$ci_ub[which(dat.fig3$effect_size == "SMDH" & 
                         dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$ci_ub_full_data[which(dat.fig3$effect_size == "ZCOR")] = 
  dat.fig3$ci_ub[which(dat.fig3$effect_size == "ZCOR" & 
                         dat.fig3$imputation_method == "full dataset analysis")]

# tau
dat.fig3$tau_full_data = NA

dat.fig3$tau_full_data[which(dat.fig3$effect_size == "logRR")] = 
  dat.fig3$tau[which(dat.fig3$effect_size == "logRR" & 
                       dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$tau_full_data[which(dat.fig3$effect_size == "SMDH")] = 
  dat.fig3$tau[which(dat.fig3$effect_size == "SMDH" & 
                       dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$tau_full_data[which(dat.fig3$effect_size == "ZCOR")] = 
  dat.fig3$tau[which(dat.fig3$effect_size == "ZCOR" & 
                       dat.fig3$imputation_method == "full dataset analysis")]

# i2
dat.fig3$i2_full_data = NA

dat.fig3$i2_full_data[which(dat.fig3$effect_size == "logRR")] = 
  dat.fig3$i2[which(dat.fig3$effect_size == "logRR" & 
                      dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$i2_full_data[which(dat.fig3$effect_size == "SMDH")] = 
  dat.fig3$i2[which(dat.fig3$effect_size == "SMDH" & 
                      dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$i2_full_data[which(dat.fig3$effect_size == "ZCOR")] = 
  dat.fig3$i2[which(dat.fig3$effect_size == "ZCOR" & 
                      dat.fig3$imputation_method == "full dataset analysis")]

# h2
dat.fig3$h2_full_data = NA

dat.fig3$h2_full_data[which(dat.fig3$effect_size == "logRR")] = 
  dat.fig3$h2[which(dat.fig3$effect_size == "logRR" & 
                      dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$h2_full_data[which(dat.fig3$effect_size == "SMDH")] = 
  dat.fig3$h2[which(dat.fig3$effect_size == "SMDH" & 
                      dat.fig3$imputation_method == "full dataset analysis")]
dat.fig3$h2_full_data[which(dat.fig3$effect_size == "ZCOR")] = 
  dat.fig3$h2[which(dat.fig3$effect_size == "ZCOR" & 
                      dat.fig3$imputation_method == "full dataset analysis")]

# Figure 3 - further data preparation ------------------------------------------------

# remove the full data as an impuation methods
dat.fig3 = dat.fig3[imputation_method != "full dataset analysis"]

# re-name imputation methods
dat.fig3$imputation_method = as.character(dat.fig3$imputation_method)
dat.fig3$imputation_method[dat.fig3$imputation_method == "complete case analysis"] = "Complete-case\nanalysis"
dat.fig3$imputation_method[dat.fig3$imputation_method == "unweighted analysis"] = "Unweighted\nanalysis"
dat.fig3$imputation_method[dat.fig3$imputation_method == "sample-size-weighted analysis"] = "Sample-size-weighted\nanalysis"
dat.fig3$imputation_method[dat.fig3$imputation_method == "mean value imputation"] = "Mean value\nimputation"
dat.fig3$imputation_method[dat.fig3$imputation_method == "median value imputation"] = "Median value\nimputation"
dat.fig3$imputation_method[dat.fig3$imputation_method == "mice::random sample imputation"] = "mice: random sample\nimputation"
dat.fig3$imputation_method[dat.fig3$imputation_method == "mice:: linear prediction imputation"] = "mice: prediction from\n linear regression"
dat.fig3$imputation_method[dat.fig3$imputation_method == "mice::pmm imputation"] = "mice: predictive\nmean matching"
dat.fig3$imputation_method[dat.fig3$imputation_method == "mice::cart imputation"] = "mice: classification \nand regression trees"
dat.fig3$imputation_method[dat.fig3$imputation_method == "mice::rf imputation"] = "mice: random forest"
dat.fig3$imputation_method[dat.fig3$imputation_method == "mi::bayes pmm imputation"] = "mi: Bayes predictive\n mean matching"
dat.fig3$imputation_method[dat.fig3$imputation_method == "Amelia::boot ebm imputation"] = "Amelia: bootstrap \nexpectation maximization"
dat.fig3$imputation_method[dat.fig3$imputation_method == "missForest::rf imputation"] = "missForest: non-parametric\n random forest"
dat.fig3$imputation_method[dat.fig3$imputation_method == "Hmisc::boot pmm imputation"] = "Hmisc: additive regression and\nbootstrap predictive mean matching"

# re-order the factor levels of imputation_method for plotting
unique(dat.fig3$imputation_method)

dat.fig3$imputation_method = factor(dat.fig3$imputation_method, 
                                    levels = unique(dat.fig3$imputation_method))

levels(dat.fig3$imputation_method)

# add rows with sample-size-weighted analyses in order get 14 plots for every effect size
dat.fig3 = rbind(dat.fig3,data.frame("effect_size" = c("logRR", "logRR", "SMDH", "SMDH", "SMDH", "ZCOR"),
                                     "what_is_missing" = c("SD & SS", "SS", "SD", "SD & SS", "SS", "SS"),
                                     "imputation_method" = "Sample-size-weighted\nanalysis"), fill = T)

# Figure 3 - plotting -----------------------------------------------------

# I need to do individual plots for every effect size
create.plot.for.fig3 = function(data.temp){
  
  dist.mean.and.ci = data.temp$grand_mean_full_data[1] - data.temp$ci_lb_full_data[1]
  plot2 = ggplot(data = data.temp) +
    
    geom_ribbon(aes(x = missing_percent, ymin = ci_lb, ymax = ci_ub, fill = imputation_method, color = imputation_method), alpha = 0.35) +
    geom_line(aes(x = missing_percent, y = grand_mean, color = imputation_method), size = 1) + 
    
    #geom_segment(aes(x = missing_percent, y = ci_lb, color = imputation_method ), size = 1, alpha = 0.35)+
    #geom_line(aes(x = missing_percent, y = ci_ub, color = imputation_method ), size = 1, alpha = 0.35)+
    
    geom_line(aes(x = missing_percent, y = grand_mean_full_data), color = "black") + 
    geom_line(aes(x = missing_percent, y = ci_lb_full_data), color = "black", linetype = "twodash") + 
    geom_line(aes(x = missing_percent, y = ci_ub_full_data), color = "black", linetype = "twodash") +
    
    facet_grid(imputation_method ~ .) +
    scale_x_reverse(breaks = c(0.1, 0.9)) +
    coord_flip(ylim = c(data.temp$grand_mean_full_data[1] - (3*dist.mean.and.ci),
                        data.temp$grand_mean_full_data[1] + (3*dist.mean.and.ci)),
               xlim = c(0.1, 0.9)) +
    
    xlab("") + ylab("") + ggtitle("") +
    
    theme(strip.text.y = element_text(angle = 0), 
          legend.position = "none",
          text = element_text(family = "Times New Roman"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  return(plot2)
}

# do the actual plotting and combine the figures later by hand

svg(filename = "results2 - bias corrected logRR/Figure3a.svg", height = 16, width = 5)
create.plot.for.fig3(data.temp = dat.fig3[effect_size == "logRR" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure3b.svg", height = 16, width = 5)
create.plot.for.fig3(data.temp = dat.fig3[effect_size == "logRR" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure3c.svg", height = 16, width = 5)
create.plot.for.fig3(data.temp = dat.fig3[effect_size == "logRR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure3d.svg", height = 16, width = 5)
create.plot.for.fig3(data.temp = dat.fig3[effect_size == "SMDH" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure3e.svg", height = 16, width = 5)
create.plot.for.fig3(data.temp = dat.fig3[effect_size == "SMDH" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure3f.svg", height = 16, width = 5)
create.plot.for.fig3(data.temp = dat.fig3[effect_size == "SMDH" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure3g.svg", height = 16, width = 5)
create.plot.for.fig3(data.temp = dat.fig3[effect_size == "ZCOR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()


#--------------------
# Figure 4 - read data, MCAR,. no correlations ---------------------------------------
getwd()
list.files("intermediate results2 - bias corrected logRR")
pos = grep('fig4', list.files("intermediate results2 - bias corrected logRR"))
dat.fig4 = read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[1]], sep = ""))
for(i in 2:length(pos)){
  dat.fig4 = rbind(dat.fig4,
                   read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[i]], sep = "")))
}
dat.fig4 = data.table(dat.fig4)

# Figure 4- arrange the correct estimates from the full model -----------------------

# grand mean
dat.fig4$grand_mean_full_data = NA

dat.fig4$grand_mean_full_data[which(dat.fig4$effect_size == "logRR")] = 
  dat.fig4$grand_mean[which(dat.fig4$effect_size == "logRR" & 
                              dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$grand_mean_full_data[which(dat.fig4$effect_size == "SMDH")] = 
  dat.fig4$grand_mean[which(dat.fig4$effect_size == "SMDH" & 
                              dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$grand_mean_full_data[which(dat.fig4$effect_size == "ZCOR")] = 
  dat.fig4$grand_mean[which(dat.fig4$effect_size == "ZCOR" & 
                              dat.fig4$imputation_method == "full dataset analysis")]

# ci_lb
dat.fig4$ci_lb_full_data = NA

dat.fig4$ci_lb_full_data[which(dat.fig4$effect_size == "logRR")] = 
  dat.fig4$ci_lb[which(dat.fig4$effect_size == "logRR" & 
                         dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$ci_lb_full_data[which(dat.fig4$effect_size == "SMDH")] = 
  dat.fig4$ci_lb[which(dat.fig4$effect_size == "SMDH" & 
                         dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$ci_lb_full_data[which(dat.fig4$effect_size == "ZCOR")] = 
  dat.fig4$ci_lb[which(dat.fig4$effect_size == "ZCOR" & 
                         dat.fig4$imputation_method == "full dataset analysis")]

# ci_ub
dat.fig4$ci_ub_full_data = NA

dat.fig4$ci_ub_full_data[which(dat.fig4$effect_size == "logRR")] = 
  dat.fig4$ci_ub[which(dat.fig4$effect_size == "logRR" & 
                         dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$ci_ub_full_data[which(dat.fig4$effect_size == "SMDH")] = 
  dat.fig4$ci_ub[which(dat.fig4$effect_size == "SMDH" & 
                         dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$ci_ub_full_data[which(dat.fig4$effect_size == "ZCOR")] = 
  dat.fig4$ci_ub[which(dat.fig4$effect_size == "ZCOR" & 
                         dat.fig4$imputation_method == "full dataset analysis")]

# tau
dat.fig4$tau_full_data = NA

dat.fig4$tau_full_data[which(dat.fig4$effect_size == "logRR")] = 
  dat.fig4$tau[which(dat.fig4$effect_size == "logRR" & 
                       dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$tau_full_data[which(dat.fig4$effect_size == "SMDH")] = 
  dat.fig4$tau[which(dat.fig4$effect_size == "SMDH" & 
                       dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$tau_full_data[which(dat.fig4$effect_size == "ZCOR")] = 
  dat.fig4$tau[which(dat.fig4$effect_size == "ZCOR" & 
                       dat.fig4$imputation_method == "full dataset analysis")]

# i2
dat.fig4$i2_full_data = NA

dat.fig4$i2_full_data[which(dat.fig4$effect_size == "logRR")] = 
  dat.fig4$i2[which(dat.fig4$effect_size == "logRR" & 
                      dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$i2_full_data[which(dat.fig4$effect_size == "SMDH")] = 
  dat.fig4$i2[which(dat.fig4$effect_size == "SMDH" & 
                      dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$i2_full_data[which(dat.fig4$effect_size == "ZCOR")] = 
  dat.fig4$i2[which(dat.fig4$effect_size == "ZCOR" & 
                      dat.fig4$imputation_method == "full dataset analysis")]

# h2
dat.fig4$h2_full_data = NA

dat.fig4$h2_full_data[which(dat.fig4$effect_size == "logRR")] = 
  dat.fig4$h2[which(dat.fig4$effect_size == "logRR" & 
                      dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$h2_full_data[which(dat.fig4$effect_size == "SMDH")] = 
  dat.fig4$h2[which(dat.fig4$effect_size == "SMDH" & 
                      dat.fig4$imputation_method == "full dataset analysis")]
dat.fig4$h2_full_data[which(dat.fig4$effect_size == "ZCOR")] = 
  dat.fig4$h2[which(dat.fig4$effect_size == "ZCOR" & 
                      dat.fig4$imputation_method == "full dataset analysis")]

# Figure 4 - further data preparation ------------------------------------------------

# remove the full data as an impuation methods
dat.fig4 = dat.fig4[imputation_method != "full dataset analysis"]

# re-name imputation methods
dat.fig4$imputation_method = as.character(dat.fig4$imputation_method)
dat.fig4$imputation_method[dat.fig4$imputation_method == "complete case analysis"] = "Complete-case\nanalysis"
dat.fig4$imputation_method[dat.fig4$imputation_method == "unweighted analysis"] = "Unweighted\nanalysis"
dat.fig4$imputation_method[dat.fig4$imputation_method == "sample-size-weighted analysis"] = "Sample-size-weighted\nanalysis"
dat.fig4$imputation_method[dat.fig4$imputation_method == "mean value imputation"] = "Mean value\nimputation"
dat.fig4$imputation_method[dat.fig4$imputation_method == "median value imputation"] = "Median value\nimputation"
dat.fig4$imputation_method[dat.fig4$imputation_method == "mice::random sample imputation"] = "mice: random sample\nimputation"
dat.fig4$imputation_method[dat.fig4$imputation_method == "mice:: linear prediction imputation"] = "mice: prediction from\n linear regression"
dat.fig4$imputation_method[dat.fig4$imputation_method == "mice::pmm imputation"] = "mice: predictive\nmean matching"
dat.fig4$imputation_method[dat.fig4$imputation_method == "mice::cart imputation"] = "mice: classification \nand regression trees"
dat.fig4$imputation_method[dat.fig4$imputation_method == "mice::rf imputation"] = "mice: random forest"
dat.fig4$imputation_method[dat.fig4$imputation_method == "mi::bayes pmm imputation"] = "mi: Bayes predictive\n mean matching"
dat.fig4$imputation_method[dat.fig4$imputation_method == "Amelia::boot ebm imputation"] = "Amelia: bootstrap \nexpectation maximization"
dat.fig4$imputation_method[dat.fig4$imputation_method == "missForest::rf imputation"] = "missForest: non-parametric\n random forest"
dat.fig4$imputation_method[dat.fig4$imputation_method == "Hmisc::boot pmm imputation"] = "Hmisc: additive regression and\nbootstrap predictive mean matching"

# re-order the factor levels of imputation_method for plotting
unique(dat.fig4$imputation_method)

dat.fig4$imputation_method = factor(dat.fig4$imputation_method, 
                                    levels = unique(dat.fig4$imputation_method))

levels(dat.fig4$imputation_method)

# add rows with sample-size-weighted analyses in order get 14 plots for every effect size
dat.fig4 = rbind(dat.fig4,data.frame("effect_size" = c("logRR", "logRR", "SMDH", "SMDH", "SMDH", "ZCOR"),
                                     "what_is_missing" = c("SD & SS", "SS", "SD", "SD & SS", "SS", "SS"),
                                     "imputation_method" = "Sample-size-weighted\nanalysis"), fill = T)

# Figure 4 - plotting -----------------------------------------------------

# I need to do individual plots for every effect size
create.plot.for.fig4 = function(data.temp){
  
  dist.mean.and.ci = data.temp$grand_mean_full_data[1] - data.temp$ci_lb_full_data[1]
  plot2 = ggplot(data = data.temp) +
    
    geom_ribbon(aes(x = missing_percent, ymin = ci_lb, ymax = ci_ub, fill = imputation_method, color = imputation_method), alpha = 0.35) +
    geom_line(aes(x = missing_percent, y = grand_mean, color = imputation_method), size = 1) + 
    
    #geom_segment(aes(x = missing_percent, y = ci_lb, color = imputation_method ), size = 1, alpha = 0.35)+
    #geom_line(aes(x = missing_percent, y = ci_ub, color = imputation_method ), size = 1, alpha = 0.35)+
    
    geom_line(aes(x = missing_percent, y = grand_mean_full_data), color = "black") + 
    geom_line(aes(x = missing_percent, y = ci_lb_full_data), color = "black", linetype = "twodash") + 
    geom_line(aes(x = missing_percent, y = ci_ub_full_data), color = "black", linetype = "twodash") +
    
    facet_grid(imputation_method ~ .) +
    scale_x_reverse(breaks = c(0.1, 0.9)) +
    coord_flip(ylim = c(data.temp$grand_mean_full_data[1] - (3*dist.mean.and.ci),
                        data.temp$grand_mean_full_data[1] + (3*dist.mean.and.ci)),
               xlim = c(0.1, 0.9)) +
    
    xlab("") + ylab("") + ggtitle("") +
    
    theme(strip.text.y = element_text(angle = 0), 
          legend.position = "none",
          text = element_text(family = "Times New Roman"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  return(plot2)
}

# do the actual plotting and combine the figures later by hand

svg(filename = "results2 - bias corrected logRR/Figure4a.svg", height = 16, width = 5)
create.plot.for.fig4(data.temp = dat.fig4[effect_size == "logRR" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure4b.svg", height = 16, width = 5)
create.plot.for.fig4(data.temp = dat.fig4[effect_size == "logRR" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure4c.svg", height = 16, width = 5)
create.plot.for.fig4(data.temp = dat.fig4[effect_size == "logRR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure4d.svg", height = 16, width = 5)
create.plot.for.fig4(data.temp = dat.fig4[effect_size == "SMDH" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure4e.svg", height = 16, width = 5)
create.plot.for.fig4(data.temp = dat.fig4[effect_size == "SMDH" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure4f.svg", height = 16, width = 5)
create.plot.for.fig4(data.temp = dat.fig4[effect_size == "SMDH" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure4g.svg", height = 16, width = 5)
create.plot.for.fig4(data.temp = dat.fig4[effect_size == "ZCOR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()



#--------------------
# Figure 5 - read data, MCAR,. no correlations ---------------------------------------
getwd()
list.files("intermediate results2 - bias corrected logRR")
pos = grep('fig5', list.files("intermediate results2 - bias corrected logRR"))
dat.fig5 = read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[1]], sep = ""))
for(i in 2:length(pos)){
  dat.fig5 = rbind(dat.fig5,
                   read.table(paste("intermediate results2 - bias corrected logRR/", list.files("intermediate results2 - bias corrected logRR")[pos[i]], sep = "")))
}
dat.fig5 = data.table(dat.fig5)

# Figure 5 - arrange the correct estimates from the full model -----------------------

# grand mean
dat.fig5$grand_mean_full_data = NA

dat.fig5$grand_mean_full_data[which(dat.fig5$effect_size == "logRR")] = 
  dat.fig5$grand_mean[which(dat.fig5$effect_size == "logRR" & 
                              dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$grand_mean_full_data[which(dat.fig5$effect_size == "SMDH")] = 
  dat.fig5$grand_mean[which(dat.fig5$effect_size == "SMDH" & 
                              dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$grand_mean_full_data[which(dat.fig5$effect_size == "ZCOR")] = 
  dat.fig5$grand_mean[which(dat.fig5$effect_size == "ZCOR" & 
                              dat.fig5$imputation_method == "full dataset analysis")]

# ci_lb
dat.fig5$ci_lb_full_data = NA

dat.fig5$ci_lb_full_data[which(dat.fig5$effect_size == "logRR")] = 
  dat.fig5$ci_lb[which(dat.fig5$effect_size == "logRR" & 
                         dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$ci_lb_full_data[which(dat.fig5$effect_size == "SMDH")] = 
  dat.fig5$ci_lb[which(dat.fig5$effect_size == "SMDH" & 
                         dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$ci_lb_full_data[which(dat.fig5$effect_size == "ZCOR")] = 
  dat.fig5$ci_lb[which(dat.fig5$effect_size == "ZCOR" & 
                         dat.fig5$imputation_method == "full dataset analysis")]

# ci_ub
dat.fig5$ci_ub_full_data = NA

dat.fig5$ci_ub_full_data[which(dat.fig5$effect_size == "logRR")] = 
  dat.fig5$ci_ub[which(dat.fig5$effect_size == "logRR" & 
                         dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$ci_ub_full_data[which(dat.fig5$effect_size == "SMDH")] = 
  dat.fig5$ci_ub[which(dat.fig5$effect_size == "SMDH" & 
                         dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$ci_ub_full_data[which(dat.fig5$effect_size == "ZCOR")] = 
  dat.fig5$ci_ub[which(dat.fig5$effect_size == "ZCOR" & 
                         dat.fig5$imputation_method == "full dataset analysis")]

# tau
dat.fig5$tau_full_data = NA

dat.fig5$tau_full_data[which(dat.fig5$effect_size == "logRR")] = 
  dat.fig5$tau[which(dat.fig5$effect_size == "logRR" & 
                       dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$tau_full_data[which(dat.fig5$effect_size == "SMDH")] = 
  dat.fig5$tau[which(dat.fig5$effect_size == "SMDH" & 
                       dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$tau_full_data[which(dat.fig5$effect_size == "ZCOR")] = 
  dat.fig5$tau[which(dat.fig5$effect_size == "ZCOR" & 
                       dat.fig5$imputation_method == "full dataset analysis")]

# i2
dat.fig5$i2_full_data = NA

dat.fig5$i2_full_data[which(dat.fig5$effect_size == "logRR")] = 
  dat.fig5$i2[which(dat.fig5$effect_size == "logRR" & 
                      dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$i2_full_data[which(dat.fig5$effect_size == "SMDH")] = 
  dat.fig5$i2[which(dat.fig5$effect_size == "SMDH" & 
                      dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$i2_full_data[which(dat.fig5$effect_size == "ZCOR")] = 
  dat.fig5$i2[which(dat.fig5$effect_size == "ZCOR" & 
                      dat.fig5$imputation_method == "full dataset analysis")]

# h2
dat.fig5$h2_full_data = NA

dat.fig5$h2_full_data[which(dat.fig5$effect_size == "logRR")] = 
  dat.fig5$h2[which(dat.fig5$effect_size == "logRR" & 
                      dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$h2_full_data[which(dat.fig5$effect_size == "SMDH")] = 
  dat.fig5$h2[which(dat.fig5$effect_size == "SMDH" & 
                      dat.fig5$imputation_method == "full dataset analysis")]
dat.fig5$h2_full_data[which(dat.fig5$effect_size == "ZCOR")] = 
  dat.fig5$h2[which(dat.fig5$effect_size == "ZCOR" & 
                      dat.fig5$imputation_method == "full dataset analysis")]

# Figure 5 - further data preparation ------------------------------------------------

# remove the full data as an impuation methods
dat.fig5 = dat.fig5[imputation_method != "full dataset analysis"]

# re-name imputation methods
dat.fig5$imputation_method = as.character(dat.fig5$imputation_method)
dat.fig5$imputation_method[dat.fig5$imputation_method == "complete case analysis"] = "Complete-case\nanalysis"
dat.fig5$imputation_method[dat.fig5$imputation_method == "unweighted analysis"] = "Unweighted\nanalysis"
dat.fig5$imputation_method[dat.fig5$imputation_method == "sample-size-weighted analysis"] = "Sample-size-weighted\nanalysis"
dat.fig5$imputation_method[dat.fig5$imputation_method == "mean value imputation"] = "Mean value\nimputation"
dat.fig5$imputation_method[dat.fig5$imputation_method == "median value imputation"] = "Median value\nimputation"
dat.fig5$imputation_method[dat.fig5$imputation_method == "mice::random sample imputation"] = "mice: random sample\nimputation"
dat.fig5$imputation_method[dat.fig5$imputation_method == "mice:: linear prediction imputation"] = "mice: prediction from\n linear regression"
dat.fig5$imputation_method[dat.fig5$imputation_method == "mice::pmm imputation"] = "mice: predictive\nmean matching"
dat.fig5$imputation_method[dat.fig5$imputation_method == "mice::cart imputation"] = "mice: classification \nand regression trees"
dat.fig5$imputation_method[dat.fig5$imputation_method == "mice::rf imputation"] = "mice: random forest"
dat.fig5$imputation_method[dat.fig5$imputation_method == "mi::bayes pmm imputation"] = "mi: Bayes predictive\n mean matching"
dat.fig5$imputation_method[dat.fig5$imputation_method == "Amelia::boot ebm imputation"] = "Amelia: bootstrap \nexpectation maximization"
dat.fig5$imputation_method[dat.fig5$imputation_method == "missForest::rf imputation"] = "missForest: non-parametric\n random forest"
dat.fig5$imputation_method[dat.fig5$imputation_method == "Hmisc::boot pmm imputation"] = "Hmisc: additive regression and\nbootstrap predictive mean matching"

# re-order the factor levels of imputation_method for plotting
unique(dat.fig5$imputation_method)

dat.fig5$imputation_method = factor(dat.fig5$imputation_method, 
                                    levels = unique(dat.fig5$imputation_method))

levels(dat.fig5$imputation_method)

# add rows with sample-size-weighted analyses in order get 14 plots for every effect size
dat.fig5 = rbind(dat.fig5,data.frame("effect_size" = c("logRR", "logRR", "SMDH", "SMDH", "SMDH", "ZCOR"),
                                     "what_is_missing" = c("SD & SS", "SS", "SD", "SD & SS", "SS", "SS"),
                                     "imputation_method" = "Sample-size-weighted\nanalysis"), fill = T)

# Figure 5 - plotting -----------------------------------------------------

# I need to do individual plots for every effect size
create.plot.for.fig5 = function(data.temp){
  
  dist.mean.and.ci = data.temp$grand_mean_full_data[1] - data.temp$ci_lb_full_data[1]
  plot2 = ggplot(data = data.temp) +
    
    geom_ribbon(aes(x = missing_percent, ymin = ci_lb, ymax = ci_ub, fill = imputation_method, color = imputation_method), alpha = 0.35) +
    geom_line(aes(x = missing_percent, y = grand_mean, color = imputation_method), size = 1) + 
    
    #geom_segment(aes(x = missing_percent, y = ci_lb, color = imputation_method ), size = 1, alpha = 0.35)+
    #geom_line(aes(x = missing_percent, y = ci_ub, color = imputation_method ), size = 1, alpha = 0.35)+
    
    geom_line(aes(x = missing_percent, y = grand_mean_full_data), color = "black") + 
    geom_line(aes(x = missing_percent, y = ci_lb_full_data), color = "black", linetype = "twodash") + 
    geom_line(aes(x = missing_percent, y = ci_ub_full_data), color = "black", linetype = "twodash") +
    
    facet_grid(imputation_method ~ .) +
    scale_x_reverse(breaks = c(0.1, 0.9)) +
    coord_flip(ylim = c(data.temp$grand_mean_full_data[1] - (3*dist.mean.and.ci),
                        data.temp$grand_mean_full_data[1] + (3*dist.mean.and.ci)),
               xlim = c(0.1, 0.9)) +
    
    xlab("") + ylab("") + ggtitle("") +
    
    theme(strip.text.y = element_text(angle = 0), 
          legend.position = "none",
          text = element_text(family = "Times New Roman"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  return(plot2)
}

# do the actual plotting and combine the figures later by hand

svg(filename = "results2 - bias corrected logRR/Figure5a.svg", height = 16, width = 5)
create.plot.for.fig5(data.temp = dat.fig5[effect_size == "logRR" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure5b.svg", height = 16, width = 5)
create.plot.for.fig5(data.temp = dat.fig5[effect_size == "logRR" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure5c.svg", height = 16, width = 5)
create.plot.for.fig5(data.temp = dat.fig5[effect_size == "logRR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure5d.svg", height = 16, width = 5)
create.plot.for.fig5(data.temp = dat.fig5[effect_size == "SMDH" & what_is_missing == "SD"]) + 
  ggtitle("logRR, weighted by vi, missing is SD")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure5e.svg", height = 16, width = 5)
create.plot.for.fig5(data.temp = dat.fig5[effect_size == "SMDH" & what_is_missing == "SD & SS"]) +
  ggtitle("logRR, weighted by vi, missing is SD & SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure5f.svg", height = 16, width = 5)
create.plot.for.fig5(data.temp = dat.fig5[effect_size == "SMDH" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()

svg(filename = "results2 - bias corrected logRR/Figure5g.svg", height = 16, width = 5)
create.plot.for.fig5(data.temp = dat.fig5[effect_size == "ZCOR" & what_is_missing == "SS"]) +
  ggtitle("logRR, weighted by SS, missing is SS")
graphics.off()


