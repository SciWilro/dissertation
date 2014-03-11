# This script saves the data for each comparison in the 'comp' variable
source("./setComparisons.R")

# Type of clustering
distance.function = function(x) as.dist(1 - cor(t(x)))
cluster.function = function(x) hclust( x, method="ward" )

################################################################################
#                                Figure 4A                                     #
################################################################################

############## Pre-treatment, PD versus SD/PR ####################
# Variables needed for the analysis
classes = comp[[9]]$classes
lmr.values = comp[[9]]$lmr.values
ave.dlmrs = comp[[9]]$effect.size

# unpaired t-test
foo = function(x) ttest(x,classes,c("sd","pr"),"pd")
pvals = apply(lmr.values,1,foo)
topmiRNAs = pvals<0.01 & abs(ave.dlmrs)>0.5 #cutoffs

# clustering analysis
lmr.values.clus = lmr.values[,classes %in% c("pd","sd","pr") ]
heatmap.2( lmr.values.clus[topmiRNAs,], col=bluered,trace="none",
           distfun=distance.function,hclustfun=cluster.function,
           mar=c(5,15), breaks=100, symbreaks=T,density.info="none")

################################################################################
#                                Figure 4B                                     #
################################################################################

############### dLMR from com-pre, PD versus SD/PR ####################
# Variables needed for the analysis
dlmr.values = comp[[7]]$dlmr.values
classes = comp[[7]]$classes
ddLMR = comp[[7]]$effect.size # effect size

# unpaired t-test on dLMR values
foo = function(x) ttest(x,classes,"sd","pd")
pvals = apply(dlmr.values,1,foo)
topmiRNAs = pvals<0.04 & abs(ddLMR)>0.5

# clustering analysis
dlmr.values.clus = dlmr.values[,classes %in% c("pd","sd","pr") ]
heatmap.2( dlmr.values.clus[topmiRNAs,], col=bluered,trace="none",
           distfun=distance.function,hclustfun=cluster.function,
           mar=c(5,15), breaks=100, symbreaks=T, density.info="none")


################################################################################
#                                Figure 4C                                     #
################################################################################

############### Pre-treatment, braf mutant versus braf wt ######################
# Variables needed for the analysis
classes = comp[[14]]$classes
lmr.values = comp[[14]]$lmr.values
ave.dlmrs = comp[[14]]$effect.size

# unpaired t-test
foo = function(x) ttest(x,classes,"wt","m")
pvals = apply(lmr.values,1,foo)
topmiRNAs = pvals<0.02 & abs(ave.dlmrs)>0.5

# clustering analysis
lmr.values.clus = lmr.values[,classes %in% c("m","wt") ]
heatmap.2( lmr.values.clus[topmiRNAs,], col=bluered,trace="none",
           distfun=distance.function,hclustfun=cluster.function,
           mar=c(5,15), breaks=100, symbreaks=T, density.info="none")


