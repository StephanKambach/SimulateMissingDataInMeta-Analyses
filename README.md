# SimulateMissingDataInMeta-Analyses

Ecological meta-analyses often encounter incompletely reported studies that have missing variance measures (e.g. standard deviations, SDs) or sample sizes (SSs) that are necessary for incorporating the study in a weighted meta-analysis.
Here, we evaluated 14 different options to treat missing information in simulated data sets of known structures and between 10% and 90% of missing standard deviations and/or sample sizes.

Script 1 (Mainscript) 
 - creates meta-analysis data sets (group comparisons and correlation analyses) with varying size, range of SD and SS values and correlation structures
  - deletes between 10 and 90% of SDs and/or SSs completely at random (MCAR), at random (MAR) or not at random (MNAR)
  - applies 14 approaches to treat/impute those deleted values
  - saves the obtained grand mean estimates in a folder that must be named "intermediate output"

Script 2(additional functions2 - bias corrected logRR) contains the functions used in script 1
  
Script 3 (Plotting of results2 - bias corrected logRR) creates figures 
- for the literature review
- for comparing the grand mean and confidence intervals obtained from the treatment of missing data with the grand mean and CIs of the complete data sets
