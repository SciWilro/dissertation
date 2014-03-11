# Parse and format the "full" miniML file from NCBI GEO
source("./helperFunctions.R")
source("./loadData.R")
eset = loadData(miniML.file="./GSE37131.xml")

################################################################################
#                                Figure 1                                      #
################################################################################

# Automatically download the miRNA annotations for this platform from NCBI
id.to.miRNA = getGEO.ids("GPL11434")
id.to.miRNA = id.to.miRNA[ featureNames(eset) ]

# remove any row with a missing value
lmr.values = exprs(eset)[ apply(!is.na(exprs(eset)),1,all), ]

# select the 50 most variable miRNAs
top50 = names( sort( apply(lmr.values,1,sd), decreasing=TRUE ) )[1:50]
lmr.values = lmr.values[top50,]
colnames(lmr.values) = 1:ncol(lmr.values) # must be done for set.dend.order function

# A cluster analysis
distance.function = function(x) as.dist(1 - cor(t(x)))
cluster.function = function(x) hclust( x, method="ward" )
dend = as.dendrogram( cluster.function( distance.function(t(lmr.values)) ) )

# choose a color to display for each patient
pat = as.numeric(pData(eset)$patient)
names(pat) = colnames(lmr.values)
cc = timPalette(n=length(unique(pat)))

# Reorder the columns and draw the final heatmap
final.order = c(10,33,11,12,24,25,26,36,27,28,7,8,9,14,13,34,15,17,35,16,
                3,32,4,21,22,1,5,6,20,19,18,29,30,31,2,23)
dend2 = set.dend.order(dend,final.order)
colnames(lmr.values) = paste( pData(eset)$patient, pData(eset)$treatment, sep="." )
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
heatmap.2(lmr.values,Colv=dend2,col=bluered,trace="none",
          ColSideColors=cc[pat],mar=c(5,15),breaks=100,symbreaks=T,density.info="none")


