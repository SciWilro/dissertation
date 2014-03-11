require(XML)
require(GEOquery)

################################################################################
#                  Downloading and formatting the data                         #
################################################################################

# the function to convert miniML file to and ExpressionSet object
loadData = function( miniML.file ) {
  tree = xmlTreeParse(miniML.file)
  eset = miniML.to.eSet(tree)
  
  # rename classes as single treatment, combination treatment, or pre-treatment
  levels(pData(eset)$treatment) = c("sin","com","pre")
  
  # temporary fix because one of the microarrays is mislabeled
  pData(eset)["GSM911902","patient"] = 7
  
  # Mark the patients response to treatment: stable disease (sd), 
  # progressive disease(pd), partial response (pr), no data (nd)
  responses = c("sd","sd","sd","sd","pd","pr","pd","pd","sd","ne","sd","sd")
  names(responses) = c(1,2,3,4,5,6,7,8,9,10,11,12)
  pData(eset)$response = responses[ as.character(pData(eset)$patient) ]
  
  # Mark the patients with BRAF mutation
  braf = c("wt","m","m","m","m","wt","wt","m","wt","wt","wt","wt")
  names(braf) = c(1,2,3,4,5,6,7,8,9,10,11,12)
  pData(eset)$braf = braf[ as.character(pData(eset)$patient) ]
  
  # remove columns not needed for this analysis
  pData(eset)[["days after treatment"]] = NULL
  pData(eset)[["tumor"]] = NULL
  pData(eset)[["tissue"]] = NULL
  pData(eset)
  
  # Remove rows where all of the samples have missing values
  exprs(eset) = exprs(eset)[ apply( !is.na(exprs(eset)), 1, any ), ]
  return(eset)
}

################################################################################
#                          Helper Functions                                    #
################################################################################

# Take a MiniML file and return a list where each element contains 
# (1) a data matrix and (2) a vector of characteristics
miniML.to.List = function( xmlTree ) {
  
  root = xmlRoot(xmlTree)
  samples = root[ sapply(xmlChildren(root),xmlName) == "Sample" ]
  
  sample.Rlist = list()
  for( i in seq_along(samples) ) {
    sample = samples[[i]]
    sample.kids = xmlChildren(sample)
    
    # get the characteristics of each sample
    channel = sample.kids$Channel
    characteristics = channel[ sapply(xmlChildren(channel),xmlName) == "Characteristics" ]
    iid = xmlAttrs(sample,"iid")
    sample.Rlist[[iid]] = list()
    tags = sapply(characteristics, function(x) xmlAttrs(x,"tag") )
    values = sapply(characteristics, xmlValue )
    names(values) = tags
    sample.Rlist[[iid]]$characteristics = values
    
    # get the data of each sample
    data = sample[["Data-Table"]]
    columns = data[ sapply(xmlChildren(data),xmlName) == "Column" ]
    col.names = sapply(columns,function(x) xmlValue(xmlChildren(x)$Name) )
    data.string = xmlValue( data[["Internal-Data"]] )
    data.rows = strsplit(data.string,"\n")
    data.vector = unlist( sapply(data.rows,function(x) strsplit(x,"\t")) )
    data = t( matrix(data.vector,nrow=length(columns)) )
    colnames(data) = col.names
    sample.Rlist[[iid]]$data = data
  }
  return(sample.Rlist)  
}

# Take a MiniML file and return an eSet object
miniML.to.eSet = function( xmlTree ) {
  
  root = xmlRoot(xmlTree)
  samples = root[ sapply(xmlChildren(root),xmlName) == "Sample" ]
  xmlList = miniML.to.List(xmlTree)
  
  characteristics = as.data.frame(t(sapply(xmlList,function(x) x$characteristics)))
  
  # Remove the ID_REF column and make it the rownames of the matrix
  foo = function(x) {
    key.col = "ID_REF"
    row.names = x[,key.col]
    col.names = colnames(x)[colnames(x) != key.col]
    x = matrix( x[,col.names], ncol=length(col.names) )
    rownames(x) = row.names
    colnames(x) = col.names
    return(x)
  }
  dataList = lapply(xmlList,function(x) x$data)
  named.dataList = lapply( dataList, foo )
  probeNames = rownames(named.dataList[[1]])
  dataMatrix = sapply(named.dataList,function(x) x[,1] )
  dataMatrix[dataMatrix=="null"] = NA
  data = matrix(as.numeric(dataMatrix),ncol=ncol(dataMatrix),
                dimnames=list(probeNames,colnames(dataMatrix)))
  
  out = new("ExpressionSet",exprs=data)
  pData(out) = characteristics
  return(out)
}

# Download ID to miRNA mappings directly from GEO
getGEO.ids = function( GPL.ID ) {
  
  platform = getGEO(GPL.ID)@dataTable@table
  platform.ids = platform[,1]
  
  # If there are multiple names (i.e. hsa/mmu for human/mouse) just take the human
  platform.names = as.character(platform[,2])
  
  foo = function(x) {
    ids = strsplit(as.character(x),"/")[[1]] 
    if( length(ids) > 1 && length(grep("hsa",ids))>0 ) {
      id = ids[grep("hsa",ids)]
      id = paste(id,collapse="/")
      return(id)
    }
    return(as.character(x))
  }
  platform.names = unlist(sapply(platform.names,foo))
  names(platform.names) = platform.ids
  return(platform.names)
}




