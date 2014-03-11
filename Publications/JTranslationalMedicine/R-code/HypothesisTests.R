# This script saves the data for each comparison in the 'comp' variable
source("./setComparisons.R")

################################################################################
#           T-tests and Significant Analysis of Microarrays                    #
################################################################################

# Save the top miRNAs from each hypothesis test to be used in the cluster
# analyses
topHits = list()

########################     S1 and S2    ###################################
# Paired, Pretreatment LMR values and Post-combination treatment LMR values

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
volcanoPlot(pvals,ave.dlmrs,topmiRNAs)
topHits$s1 = names(which(topmiRNAs))
topHits$s2 = names(which(topmiRNAs))

###########################    S3 and S4    ####################################
# Paired, Pretreatment LMR values and Post-temsirolimus treatment LMR values

# Variables needed for the analysis
time = comp$s3$unpaired.time
patient = comp$s3$unpaired.patient
lmr.values = comp$s3$lmr.values
ave.dlmrs = comp$s3$ave.dlmrs

# paired t-test
foo = function(x) ttest.paired(x,time,patient)
pvals = apply(lmr.values,1,foo)
topmiRNAs = pvals < 0.05 & abs(ave.dlmrs) > 0.1 # cutoffs
topTable(pvals,ave.dlmrs,topmiRNAs)
volcanoPlot(pvals,ave.dlmrs,topmiRNAs)
topHits$s3 = names(which(topmiRNAs))
topHits$s4 = names(which(topmiRNAs))


#########################          S5          #################################
# Unpaired, Pretreatment LMR values and Post-combination treatment LMR values

lmr.values = comp$s5$lmr.values
classes = comp$s5$classes
effect.size = comp$s5$effect.size

foo = function(x) ttest(x,classes,"com","pre")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.1 & abs(effect.size)>0.5 #cutoffs
topTable(pvals,effect.size,topmiRNAs)
volcanoPlot(pvals,effect.size,topmiRNAs)
topHits$s5 = names(which(topmiRNAs))

#########################          S6          #################################
# Unpaired, Pretreatment LMR values and Post-temsirolimus treatment LMR values

lmr.values = comp$s6$lmr.values
classes = comp$s6$classes
effect.size = comp$s6$effect.size

foo = function(x) ttest(x,classes,"sin","pre")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.1  #cutoffs
topTable(pvals,effect.size,topmiRNAs)
volcanoPlot(pvals,effect.size,topmiRNAs)
topHits$s6 = names(which(topmiRNAs))



###############################    s7    #######################################
# dLMR combination-pretreatment. SD/PR versus PD

dlmr.values = comp$s7$dlmr.values
classes = comp$s7$classes
effect = comp$s7$effect.size

foo = function(x) ttest(x,classes,c("sd","pr"),"pd")
pvals = apply(dlmr.values,1,foo)

topmiRNAs = pvals<0.04 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s7 = names(which(topmiRNAs))

##############################    s8      ######################################

dlmr.values = comp$s8$dlmr.values
classes = comp$s8$classes
effect = comp$s8$effect.size

foo = function(x) ttest(x,classes,c("sd","pr"),"pd")
pvals = apply(dlmr.values,1,foo)

topmiRNAs = pvals<0.05 & abs(effect)>0.2 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s8 = names(which(topmiRNAs))

##########################       s9       ######################################

lmr.values = comp$s9$lmr.values
classes = comp$s9$classes
effect = comp$s9$effect.size

foo = function(x) ttest(x,classes,c("sd","pr"),"pd")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.01 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s9 = names(which(topmiRNAs))

##########################       s10       ######################################

lmr.values = comp$s10$lmr.values
classes = comp$s10$classes
effect = comp$s10$effect.size

foo = function(x) ttest(x,classes,c("sd","pr"),"pd")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.03 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s10 = names(which(topmiRNAs))

##########################       s11       ######################################

lmr.values = comp$s11$lmr.values
classes = comp$s11$classes
effect = comp$s11$effect.size

foo = function(x) ttest(x,classes,c("sd","pr"),"pd")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.03 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s11 = names(which(topmiRNAs))

###############################    S12    #######################################
# dLMR combination-pretreatment. SD/PR versus PD

dlmr.values = comp$s12$dlmr.values
classes = comp$s12$classes
effect = comp$s12$effect.size

foo = function(x) ttest(x,classes,"wt","m")
pvals = apply(dlmr.values,1,foo)

topmiRNAs = pvals<0.05 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s12 = names(which(topmiRNAs))

###############################    S13    #######################################
# dLMR combination-pretreatment. SD/PR versus PD

dlmr.values = comp$s13$dlmr.values
classes = comp$s13$classes
effect = comp$s13$effect.size

foo = function(x) ttest(x,classes,"wt","m")
pvals = apply(dlmr.values,1,foo)

topmiRNAs = pvals<0.05 & abs(effect)>0.2 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s13 = names(which(topmiRNAs))

##########################       S14       ######################################

lmr.values = comp$s14$lmr.values
classes = comp$s14$classes
effect = comp$s14$effect.size

foo = function(x) ttest(x,classes,"wt","m")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.02 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s14 = names(which(topmiRNAs))

##########################       S15       ######################################

lmr.values = comp$s15$lmr.values
classes = comp$s15$classes
effect = comp$s15$effect.size

foo = function(x) ttest(x,classes,"wt","m")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.05 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s15 = names(which(topmiRNAs))

##########################       S16       ######################################

lmr.values = comp$s16$lmr.values
classes = comp$s16$classes
effect = comp$s16$effect.size

foo = function(x) ttest(x,classes,"wt","m")
pvals = apply(lmr.values,1,foo)

topmiRNAs = pvals<0.05 & abs(effect)>0.5 #cutoffs
topTable(pvals,effect,topmiRNAs)
volcanoPlot(pvals,effect,topmiRNAs)
topHits$s16 = names(which(topmiRNAs))







