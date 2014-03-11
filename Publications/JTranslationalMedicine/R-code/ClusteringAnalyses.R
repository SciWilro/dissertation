# Clustering analyses for the five comparisons in HypthesisTests.SAM.R and
# HypothesisTests.t-tests.R

# Parse and format the "full" miniML file from NCBI GEO
source("./HypothesisTests.R")

################################################################################
#                         Clustering Analyses                                  #
################################################################################

########################        S1       ###################################

# Variables needed for the analysis
dlmr.values = comp$s1$dlmr.values
paired.classes = comp$s1$paired.classes
topmiRNAs = topHits$s1

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( paired.classes )
levels(dlmr.colors.clus) =  c("gray","pink","lightgreen")

heatmap.3(dlmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3(dlmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )

########################         S2       ####################################

# Variables needed for the analysis
lmr.values = comp$s2$lmr.values
unpaired.classes = comp$s2$unpaired.classes
topmiRNAs = topHits$s2

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( unpaired.classes )
levels(dlmr.colors.clus) =  c("gray","pink","lightgreen")

heatmap.3( lmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )

########################         S3       ####################################

# Variables needed for the analysis
dlmr.values = comp$s3$dlmr.values
paired.classes = comp$s3$paired.classes
topmiRNAs = topHits$s3

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( paired.classes )
levels(dlmr.colors.clus) =  c("gray","pink","lightgreen","lightgreen")

heatmap.3(dlmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3(dlmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )


########################         S4       ####################################

# Variables needed for the analysis
lmr.values = comp$s4$lmr.values
unpaired.classes = comp$s4$unpaired.classes
topmiRNAs = topHits$s4

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( unpaired.classes )
levels(dlmr.colors.clus) =  c("gray","pink","lightgreen","lightgreen")

heatmap.3( lmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )


#########################          S5          #################################
# Unaired, Pretreatment LMR values and Post-combination treatment LMR values

lmr.values = comp$s5$lmr.values
topmiRNAs = topHits$s5
classes = comp$s5$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("pink","lightgreen")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )


#########################          S6          #################################
# Unaired, Pretreatment LMR values and Post-combination treatment LMR values

lmr.values = comp$s6$lmr.values
topmiRNAs = topHits$s6
classes = comp$s6$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("pink","lightgreen")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )



#######################       S7        #######################################

dlmr.values = comp$s7$dlmr.values
topmiRNAs = topHits$s7
classes = comp$s7$classes

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( classes )
levels(dlmr.colors.clus) =  c("pink","lightgreen")

heatmap.3( dlmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3( dlmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )


#######################       S8        #######################################

dlmr.values = comp$s8$dlmr.values
topmiRNAs = topHits$s8
classes = comp$s8$classes

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( classes )
levels(dlmr.colors.clus) =  c("pink","lightgreen","lightgreen")

heatmap.3( dlmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3( dlmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )


###########################     S9       #######################################

lmr.values = comp$s9$lmr.values
topmiRNAs = topHits$s9
classes = comp$s9$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("pink","lightgreen","lightgreen")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )


###########################     S10       #######################################

lmr.values = comp$s10$lmr.values
topmiRNAs = topHits$s10
classes = comp$s10$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("pink","lightgreen","lightgreen")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )


###########################     S11       #######################################

lmr.values = comp$s11$lmr.values
topmiRNAs = topHits$s11
classes = comp$s11$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("pink","lightgreen","lightgreen")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )


#######################       S12        #######################################

dlmr.values = comp$s12$dlmr.values
topmiRNAs = topHits$s12
classes = comp$s12$classes

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( classes )
levels(dlmr.colors.clus) =  c("orange","lightblue")

heatmap.3( dlmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3( dlmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )


#######################       S13        #######################################

dlmr.values = comp$s13$dlmr.values
topmiRNAs = topHits$s13
classes = comp$s13$classes

# Setup color bars for cluster analysis
dlmr.colors.clus = factor( classes )
levels(dlmr.colors.clus) =  c("orange","lightblue")

heatmap.3( dlmr.values, ColSideColors=dlmr.colors.clus)
heatmap.3( dlmr.values[topmiRNAs,], ColSideColors=dlmr.colors.clus )


###########################     S14       #######################################

lmr.values = comp$s14$lmr.values
topmiRNAs = topHits$s14
classes = comp$s14$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("orange","lightblue")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )


###########################     S15       #######################################

lmr.values = comp$s15$lmr.values
topmiRNAs = topHits$s15
classes = comp$s15$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("orange","lightblue")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )


###########################     S16       #######################################

lmr.values = comp$s16$lmr.values
topmiRNAs = topHits$s16
classes = comp$s16$classes

# Setup color bars for cluster analysis
lmr.colors.clus = factor( classes )
levels(lmr.colors.clus) =  c("orange","lightblue")

heatmap.3( lmr.values, ColSideColors=lmr.colors.clus)
heatmap.3( lmr.values[topmiRNAs,], ColSideColors=lmr.colors.clus )



