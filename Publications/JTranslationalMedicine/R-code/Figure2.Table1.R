# This script saves the data for each comparison in the 'comp' variable
source("./setComparisons.R")
library(RSvgDevice)

################################################################################
#                  Figure 2B and Table 1                                      #
################################################################################

# Pretreatment LMR values versus Post-combination treatment LMR values


# Variables needed for the analysis
time = comp$s1$unpaired.time
patient = comp$s1$unpaired.patient
lmr.values = comp$s1$lmr.values
ave.dlmrs = comp$s1$ave.dlmrs

# paired t-test
foo = function(x) ttest.paired(x,time,patient)
pvals = apply(lmr.values,1,foo)
topmiRNAs = pvals < 0.01 & abs(ave.dlmrs) > 0.5 # cutoffs
topTable(pvals,ave.dlmrs,topmiRNAs)

devSVG(file='Figure2B.svg', height=4, width=6, onefile=TRUE)
volcanoPlot(pvals,ave.dlmrs,topmiRNAs,
            xlim=c(-1.4,1.4),ylim=c(0,4.2))
opendots = c("hsa-miR-10b","hsa-miR-140-3p","hsa-miR-4328")
points(ave.dlmrs[opendots], -log10(pvals[opendots]), col="blue" )
dev.off()


################################################################################
#                               Figure 2A                                      #
################################################################################

# Variables needed for the analysis
time = comp$s3$unpaired.time
patient = comp$s3$unpaired.patient
lmr.values = comp$s3$lmr.values
ave.dlmrs = comp$s3$ave.dlmrs

# paired t-test
foo = function(x) ttest.paired(x,time,patient)
pvals = apply(lmr.values,1,foo)

# miRNAs with p<0.01
topmiRNAs = pvals < 0.01 & abs(ave.dlmrs) > 0 # cutoffs
topTable(pvals,ave.dlmrs,topmiRNAs)

# miRNAs with effect size > 0.5
topmiRNAs = pvals < 2 & abs(ave.dlmrs) > 0.5 # cutoffs
topTable(pvals,ave.dlmrs,topmiRNAs)

#miRNAs with p<0.01 and effect size>0.5
topmiRNAs = pvals < 0.01 & abs(ave.dlmrs) > 0.5 # cutoffs
topTable(pvals,ave.dlmrs,topmiRNAs)

devSVG(file='Figure2A.svg', height=4, width=6, onefile=TRUE)
volcanoPlot(pvals,ave.dlmrs,topmiRNAs,
            xlim=c(-1.4,1.4),ylim=c(0,4.2))
dev.off()






