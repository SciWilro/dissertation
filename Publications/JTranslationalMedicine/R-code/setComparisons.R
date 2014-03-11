# Set up all the comparisons so that these results can be displayed in several
# other separate scripts

# Parse and format the "full" miniML file from NCBI GEO
source("./helperFunctions.R")
source("./loadData.R")
eset = loadData(miniML.file="./GSE37131.xml")
#eset = removeDuplicates(eset)
#exprs(eset) = omitNA(eset)

# Automatically download the miRNA annotations for this platform from NCBI
id.to.miRNA = getGEO.ids("GPL11434")
id.to.miRNA = id.to.miRNA[ featureNames(eset) ]

# Make a list to store the information for each comparison
comp = list()

###########################    S1 and S2    ####################################
# Paired, Pretreatment LMR values and Post-combination treatment LMR values

treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"com"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
  pData(eset)$treatment %in% c("pre","com")
eset.1 = eset[ , combo.samples ]
eset.1 = removeDuplicates(eset.1)
lmr.values = omitNA(eset.1)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.1)$patient,pData(eset.1)$treatment)

# dLMR
time = pData(eset.1)$treatment
patient = factor(pData(eset.1)$patient)
dlmr.values = dLMR.paired(lmr.values,time,patient)
ave.dlmrs = rowMeans(dlmr.values)

# sd/pr and pd classes
paired.classes = pData(eset.1[ , pData(eset.1)$treatment=="pre" ])$response
unpaired.classes = pData(eset.1)$response

# patients for paired and unpaired data sets
paired.patient = pData(eset.1[ , pData(eset.1)$treatment=="pre" ])$patient
paired.patient = factor(paired.patient,levels=unique(paired.patient))


comp$s1 = list( lmr.values = lmr.values, 
                  dlmr.values=dlmr.values,
                  eset=eset.1, 
                  unpaired.time=time, 
                  unpaired.patient=patient,
                  unpaired.classes=pData(eset.1)$response,
                  paired.classes=paired.classes, 
                  paired.patient=paired.patient,
                  ave.dlmrs=ave.dlmrs)
comp$s2 = comp$s1

###########################    S3 and S4    ####################################
# Paired, Pretreatment LMR values and Post-temsirolimus treatment LMR values

treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"sin"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
  pData(eset)$treatment %in% c("pre","sin")
eset.3 = eset[ , combo.samples ]
eset.3 = removeDuplicates(eset.3)
lmr.values = omitNA(eset.3)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.3)$patient,pData(eset.3)$treatment)

# dLMR
time = pData(eset.3)$treatment
patient = factor(pData(eset.3)$patient)
dlmr.values = dLMR.paired(lmr.values,time,patient)
ave.dlmrs = rowMeans(dlmr.values)

# sd/pr and pd classes
paired.classes = pData(eset.3[ , pData(eset.1)$treatment=="pre" ])$response
unpaired.classes = pData(eset.1)$response

# patients for paired and unpaired data sets
paired.patient = pData(eset.3[ , pData(eset.1)$treatment=="pre" ])$patient
paired.patient = factor(paired.patient,levels=unique(paired.patient))


comp$s3 = list( lmr.values = lmr.values, 
                dlmr.values=dlmr.values,
                eset=eset.3, 
                unpaired.time=time, 
                unpaired.patient=patient,
                unpaired.classes=pData(eset.3)$response,
                paired.classes=paired.classes, 
                paired.patient=paired.patient,
                ave.dlmrs=ave.dlmrs)
comp$s4 = comp$s3


##############################     S5       ##################################
# Unaired, Pretreatment LMR values and Post-combination treatment LMR values

treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"com"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
  pData(eset)$treatment %in% c("pre","com")
eset.5 = eset[ , combo.samples ]
eset.5 = averageDuplicates(eset.5)
lmr.values = omitNA.knn(eset.5,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.5)$patient,pData(eset.5)$treatment)

# dLMR
classes = pData(eset.5)$treatment
patient = factor(pData(eset.5)$patient)

effect.size = rowMeans(lmr.values[,classes == "com"]) -
              rowMeans(lmr.values[,classes == "pre"])

comp$s5 = list( lmr.values = lmr.values, 
                eset=eset.5, 
                classes=classes, 
                patient=patient,
                response=pData(eset.1)$response,
                effect.size = effect.size
                )

##############################     S6       ##################################
# Unpaired, Pretreatment LMR values and Post-temsirolimus treatment LMR values

treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"sin"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
                pData(eset)$treatment %in% c("pre","sin")
eset.6 = eset[ , combo.samples ]
eset.6 = averageDuplicates(eset.6)
lmr.values = omitNA.knn(eset.6,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.6)$patient,pData(eset.6)$treatment)

# dLMR
classes = pData(eset.6)$treatment
patient = factor(pData(eset.6)$patient)

effect.size = rowMeans(lmr.values[,classes == "sin"]) -
  rowMeans(lmr.values[,classes == "pre"])

comp$s6 = list( lmr.values = lmr.values, 
                eset=eset.6, 
                classes=classes, 
                patient=patient,
                response=pData(eset.6)$response,
                effect.size = effect.size
                )


###############################    S7    #######################################
# dLMR combination-pretreatment. SD/PR versus PD


treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"com"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
  pData(eset)$treatment %in% c("pre","com")
eset.7 = eset[ , combo.samples ]
eset.7 = eset.7[ ,pData(eset.7)$response != "ne" ]

eset.7 = averageDuplicates(eset.7)
lmr.values = omitNA.knn(eset.7,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.7)$patient,pData(eset.7)$treatment)

time = pData(eset.7)$treatment
patient = factor(pData(eset.7)$patient)
dlmr.values = dLMR.paired(lmr.values,time,patient)

classes = pData( eset.7[,pData(eset.7)$treatment=="pre"] )$response

effect.size = rowMeans(dlmr.values[,classes == "pd"]) -
              rowMeans(dlmr.values[,classes %in% c("pr","sd")])

comp$s7 = list( eset = eset.7,
                dlmr.values = dlmr.values,
                classes = classes,
                effect.size = effect.size )

###############################    S8    #######################################
# dLMR temsirolimus-pretreatment. SD/PR versus PD

treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"sin"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
  pData(eset)$treatment %in% c("pre","sin")
eset.8 = eset[ , combo.samples ]
eset.8 = eset.8[ ,pData(eset.8)$response != "ne" ]

eset.8 = averageDuplicates(eset.8)
lmr.values = omitNA.knn(eset.8,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.8)$patient,pData(eset.8)$treatment)

time = pData(eset.8)$treatment
patient = factor(pData(eset.8)$patient)
dlmr.values = dLMR.paired(lmr.values,time,patient)

classes = pData( eset.8[,pData(eset.8)$treatment=="pre"] )$response

effect.size = rowMeans(dlmr.values[,classes == "pd"]) -
  rowMeans(dlmr.values[,classes %in% c("pr","sd")])

comp$s8 = list( eset = eset.8,
                dlmr.values = dlmr.values,
                classes = classes,
                effect.size = effect.size )


###############################    S9    #######################################
# Preatreatment LMR values. SD/PR versus PD

eset.9 = eset[ , pData(eset)$treatment == "pre" &
                 pData(eset)$response != "ne"]
eset.9 = averageDuplicates(eset.9)
lmr.values = omitNA.knn(eset.9,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = pData(eset.9)$patient

# for an unpaired t-test
classes = pData(eset.9)$response

# The effect size between the two groups
effect.size = rowMeans(lmr.values[,classes=="pd"]) - 
             rowMeans(lmr.values[,classes %in% c("sd","pr")])

comp$s9 = list( eset = eset.9,
                classes = classes,
                lmr.values = lmr.values,
                effect.size = effect.size )


###############################    S10    #######################################
# Post-temsirolimus LMR values. SD/PR versus PD

eset.10 = eset[ , pData(eset)$treatment == "sin" &
  pData(eset)$response != "ne"]
eset.10 = averageDuplicates(eset.10)
lmr.values = omitNA.knn(eset.10,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = pData(eset.10)$patient

# for an unpaired t-test
classes = pData(eset.10)$response

# The effect size between the two groups
effect.size = rowMeans(lmr.values[,classes=="pd"]) - 
             rowMeans(lmr.values[,classes %in% c("sd","pr")])

comp$s10 = list( eset = eset.10,
                classes = classes,
                lmr.values = lmr.values,
                effect.size = effect.size )

###############################    S11    #######################################
# Post-combination LMR values. SD/PR versus PD

eset.11 = eset[ , pData(eset)$treatment == "com" &
                 pData(eset)$response != "ne"]
eset.11 = averageDuplicates(eset.11)
lmr.values = omitNA.knn(eset.11,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = pData(eset.11)$patient

# for an unpaired t-test
classes = pData(eset.11)$response

# The effect size between the two groups
effect.size = rowMeans(lmr.values[,classes=="pd"]) - 
              rowMeans(lmr.values[,classes %in% c("sd","pr")])

comp$s11 = list( eset = eset.11,
                classes = classes,
                lmr.values = lmr.values,
                effect.size = effect.size )


###############################    S12    ######################################
# dLMR combination-pretreatment. BRAF m versus BRAF wt

treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"com"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
  pData(eset)$treatment %in% c("pre","com")
eset.12 = eset[ , combo.samples ]

eset.12 = averageDuplicates(eset.12)
lmr.values = omitNA.knn(eset.12,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.12)$patient,pData(eset.12)$treatment)

time = pData(eset.12)$treatment
patient = factor(pData(eset.12)$patient)
dlmr.values = dLMR.paired(lmr.values,time,patient)

classes = pData( eset.12[,pData(eset.12)$treatment=="pre"] )$braf

effect.size = rowMeans(dlmr.values[,classes == "m"]) -
  rowMeans(dlmr.values[,classes=="wt"])

comp$s12 = list( eset = eset.12,
                dlmr.values = dlmr.values,
                classes = classes,
                effect.size = effect.size )


###############################    S13    ######################################
# dLMR temsirolimus-pretreatment. BRAF m versus BRAF wt

treatments = table(pData(eset)[,c("patient","treatment")])
combo.patients = names( which( treatments[,"pre"] & treatments[,"sin"] ) )
combo.samples = pData(eset)$patient %in% combo.patients &
  pData(eset)$treatment %in% c("pre","sin")
eset.13 = eset[ , combo.samples ]

eset.13 = averageDuplicates(eset.13)
lmr.values = omitNA.knn(eset.13,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = paste(pData(eset.13)$patient,pData(eset.13)$treatment)

time = pData(eset.13)$treatment
patient = factor(pData(eset.13)$patient)
dlmr.values = dLMR.paired(lmr.values,time,patient)

classes = pData( eset.13[,pData(eset.13)$treatment=="pre"] )$braf

effect.size = rowMeans(dlmr.values[,classes == "m"]) -
  rowMeans(dlmr.values[,classes=="wt"])

comp$s13 = list( eset = eset.13,
                 dlmr.values = dlmr.values,
                 classes = classes,
                 effect.size = effect.size )


###############################    S14    #######################################
# Preatreatment LMR values. BRAF m versus BRAF wt

eset.14 = eset[ , pData(eset)$treatment == "pre"]
eset.14 = averageDuplicates(eset.14)
lmr.values = omitNA.knn(eset.14,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = pData(eset.14)$patient

# for an unpaired t-test
classes = pData(eset.14)$braf

# The effect size between the two groups
effect.size = rowMeans(lmr.values[,classes=="m"]) - 
              rowMeans(lmr.values[,classes=="wt"])

comp$s14 = list( eset = eset.14,
                classes = classes,
                lmr.values = lmr.values,
                effect.size = effect.size )

###############################    S15    #######################################
# Preatreatment LMR values. BRAF m versus BRAF wt

eset.15 = eset[ , pData(eset)$treatment == "sin"]
eset.15 = averageDuplicates(eset.15)
lmr.values = omitNA.knn(eset.15,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = pData(eset.15)$patient

# for an unpaired t-test
classes = pData(eset.15)$braf

# The effect size between the two groups
effect.size = rowMeans(lmr.values[,classes=="m"]) - 
              rowMeans(lmr.values[,classes=="wt"])

comp$s15 = list( eset = eset.15,
                 classes = classes,
                 lmr.values = lmr.values,
                 effect.size = effect.size )

###############################    S16    #######################################
# Preatreatment LMR values. BRAF m versus BRAF wt

eset.16 = eset[ , pData(eset)$treatment == "com"]
eset.16 = averageDuplicates(eset.16)
lmr.values = omitNA.knn(eset.16,1/2)
rownames(lmr.values) = id.to.miRNA[rownames(lmr.values)]
colnames(lmr.values) = pData(eset.16)$patient

# for an unpaired t-test
classes = pData(eset.16)$braf

# The effect size between the two groups
effect.size = rowMeans(lmr.values[,classes=="m"]) - 
              rowMeans(lmr.values[,classes=="wt"])

comp$s16 = list( eset = eset.16,
                 classes = classes,
                 lmr.values = lmr.values,
                 effect.size = effect.size )









