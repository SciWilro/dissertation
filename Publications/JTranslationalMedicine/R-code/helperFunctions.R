require(gplots)
require(fBasics)
require(impute)

# A wrapper function to change the defaults for heatmap.2
heatmap.3 = function(x, ColSideColors=NULL) {
  
  distance.function = function(x) as.dist(1 - cor(t(x)))
  cluster.function = function(x) hclust( x, method="ward" )
  
  if(is.null(ColSideColors)) {
    heatmap.2(x, trace="none", col=bluered, breaks=100,
              hclustfun=cluster.function, distfun=distance.function,
              density.info="none",mar=c(5,15),symbreaks=TRUE )
  } else {
    heatmap.2(x, trace="none", col=bluered, breaks=100,
              hclustfun=cluster.function, distfun=distance.function,
              density.info="none", 
              ColSideColors=as.character(ColSideColors), mar=c(5,15),
              symbreaks=TRUE)
  }
}

# Input:  p-values, effect sizes, and a vector of logicals indicating the 
# "significant" miRNAs. Prints the significant miRNAs
topTable = function( pvals, effectsizes, sigs ) {
  table1 = data.frame(p.val=pvals[sigs],effect=effectsizes[sigs])
  table1 = table1[order(as.numeric(table1$p.val),decreasing=FALSE),]
  return(table1)
}

# Make a volcano plot (effect size versus log10(p-values) ). The second argument
# is a vector of logicals indicating the "significant" miRNAs
volcanoPlot = function( pvals, effectsizes, sigs, xlim=NULL, ylim=NULL ) {
  
  if( is.null(ylim) | is.null(xlim) ) {
    plot(effectsizes[!sigs],
         -log10(pvals[!sigs]),
         pch=16,cex=0.5,
         xlim=range(effectsizes)*1.05,
         ylim=c(0,max(-log10(pvals))*1.05))
    points(effectsizes[sigs],
           -log10(pvals[sigs]),
           pch=16,cex=1,col="red")
  } else {
    plot(effectsizes[!sigs],
         -log10(pvals[!sigs]),
         pch=16,cex=0.5,xlim=xlim, ylim=ylim)
    points(effectsizes[sigs],
           -log10(pvals[sigs]),
           pch=16,cex=1,col="red")
  }
}

# Given an hclust object, makes a dendrogram with a specified leaf order
# The order must be provided using the hclust.object$labels
# All the values for hclust.object$labels must be unique
set.dend.order.labels = function( hclust.object, new.order ) {
  
  dend = as.dendrogram(hclust.object)
  
  leaf.names = hclust.object$labels
  names(leaf.names) = 1:ncol(data)
  
  sort.vector = names(leaf.names)
  names(sort.vector) = leaf.names
  
  new.order2 = as.character(sort.vector[as.character(new.order)])
  s.reorder = sort(as.numeric(new.order2),index.return=T)$ix
  hclust.object$labels = sort.vector
  
  dend2 = reorder(dend,s.reorder,mean)
  return(dend2)
}

# Given an allowable leaf order of the leaves of a 2D dendrogram, 
# output a dendrogram with that order of leaves. If the order cannot
# possible be drawn in 2D, then a dendrogram that's "fairly" close is returned
# the leafs must be integers from 1 to the number of leaves. Hence, the 
# 'leaf.order' given must correspond to these integers. Honestly, a very messy 
# function
set.dend.order = function( dend, leaf.order ) {
  names(leaf.order) = 1:length(labels(dend))
  reorder.vector = as.numeric(names(leaf.order))
  names(reorder.vector) = leaf.order
  s.reorder = reorder.vector[ as.character(1:length(labels(dend))) ]
  dend.out = reorder(x=dend,wts=s.reorder,agglo.FUN=mean)
  return(dend.out)
}

# calculate the dLMR for a paired data set
dLMR.paired = function( values, time, patient ) {
  time = as.character(time)
  patient = as.character(patient)
  tb = data.frame(time,patient)
  tb1 = tb[ tb$time == unique(time)[1], ]
  tb2 = tb[ tb$time == unique(time)[2], ]
  tb1$order = 1:nrow(tb1)
  tb1 = tb1[ order(tb1$patient), ]
  tb2 = tb2[ order(tb2$patient), ]
  
  group1.values = values[,as.numeric(rownames(tb1))]
  group2.values = values[,as.numeric(rownames(tb2))]

  dlmrs = group2.values - group1.values
  colnames(dlmrs) = tb1$patient
  dlmrs = dlmrs[,unique(patient)]
  return(dlmrs)
}

# A wrapper function for the t.test function customized for different inputs
# Written specifically for a paired t test
ttest.paired = function( values, time, patient ) {
  tb = data.frame(time,patient,values)
  tb1 = tb[ tb$time == unique(time)[1], ]
  tb2 = tb[ tb$time == unique(time)[2], ]
  tb1 = tb1[ order(tb1$patient), ]
  tb2 = tb2[ order(tb2$patient), ]
  group1.values = tb1$values
  group2.values = tb2$values
  p = t.test(group1.values,group2.values,paired=TRUE,var.equal=FALSE)$p.value
  return(p)
}

# A wrapper function for the t.test function customized for different inputs
# this is only intended for unpaired test, thought it could be used for paired tests
ttest = function(values, classes, group1.classes, group2.classes) {
  group1.values = values[ classes %in% group1.classes ] 
  group2.values = values[ classes %in% group2.classes ]
  p = t.test(group1.values,group2.values)$p.value
  return(p)
}

# onit any rows with a missing value
omitNA = function( eset ) {
  good.rows = rowSums(is.na(exprs(eset))) < 1
  exprs(eset) = exprs(eset)[good.rows,]
  return(exprs(eset))
}


# omit rows missing over specified percentage of values and then do knn
# imputation for the other rows
omitNA.knn = function( eset,fraction,k=10 ) {
  good.rows = rowSums(is.na(exprs(eset))) < ncol(eset)*fraction
  exprs(eset) = exprs(eset)[good.rows,]
  
  values = impute.knn(exprs(eset),k)$data
  return(values)
}

# take the average of any replicates, thus reducing the number of samples
# Duplicates determined according to 'treatment' and 'patient' slots in pData
averageDuplicates = function( eset ) {
  
  dups = which(duplicated( pData(eset)[,c("patient","treatment")] ))
  samples.to.remove = numeric()
  for( i in dups ) {
    samp.type = pData(eset)[i,c("patient","treatment")]
    samples = which( pData(eset)$patient == samp.type$patient & 
      pData(eset)$treatment == samp.type$treatment )
    values = exprs(eset)[,samples]
    foo = function(x) {
      if( sum(is.na(x)) == 1 )
        return(x[!is.na(x)])
      else if( all(is.na(x)) )
        return(NA)
      else if( !any(is.na(x)) )
        return(mean(x))
    }
    exprs(eset)[,samples[1]] = apply(values,1,foo)
    samples.to.remove = c(samples.to.remove,samples[-1])
  }
  eset = eset[, setdiff(1:ncol(eset),samples.to.remove)]
  return(eset)
  
}

# remove any replicates
# Duplicates determined according to 'treatment' and 'patient' slots in pData
removeDuplicates = function( eset ) {
  
  dups = which(duplicated( pData(eset)[,c("patient","treatment")] ))
  non.dups = setdiff(1:ncol(eset),dups)
  return( eset[,non.dups])
  
}


