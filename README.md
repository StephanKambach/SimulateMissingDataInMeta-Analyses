# SimulateMissingDataInMeta-Analyses

Ecological meta-analyses often encounter incompletely reported studies that have missing variance measures (e.g. standard deviations, SDs) or sample sizes (SSs) that are necessary for incorporating the study in a weighted meta-analysis.
Here, we evaluated 14 different options to treat missing information in simulated data sets of known structures and between 10% and 90% of missing standard deviations and/or sample sizes.

Script 1 contains the functions used to treat / impute missing SDs and/or SSs.
Script 2 contains functions to 
  - create meta-analysis data sets (group comparisons and correlation analyses) with varying size, range of SD and SS values and correlation structures
  - create between 10 and 90% of missing SDs and/or SSs either MCAR, MAR or MNAR
  - apply the 14 approaches from script 1 to treat/impute those missing values
 Script 3 contains functions to plot the grand mean effect size and approximated 95% confidence intervals obtained from the treatment of missing values versus the estimates obtained from complete data sets. 
