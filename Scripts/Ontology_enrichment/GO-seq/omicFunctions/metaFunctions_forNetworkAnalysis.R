#######
#######
makeBinaryTable <- function(listDFs) {
  cat ("-- makeBinaryTable function called \n")
  # First make an empty table
  uniqRows <- Reduce(union,lapply(listDFs,rownames))
  zeroTable <- as.data.frame(matrix(0,
                                    nrow = length(uniqRows),
                                    length(names(listDFs))),
                             row.names =uniqRows)
  colnames(zeroTable) = names(listDFs)
  # Then Fill it
  for (each in names(listDFs)){
    cat ("Filling binary table:", each,"\n")
    zeroTable[rownames(listDFs[[each]]),each] <- 1
  }
  return(zeroTable)
  cat ("-- binary table done \n")
  cat ("\n")
}
## binary table function end.


######################################
######################################
#### Finding number of clusters
# Requires a table with rows as genes and columns of values (binary or)

#Use PCA for clustering:
#http://ranger.uta.edu/~chqding/papers/KmeansPCA1.pdf

idealK <- function(mydata,maxFind){
  data <- mydata
  ##
  require(cluster)
  ## Convert data
  distance = dist(mydata,method = "euclidean")
  pca <- princomp(distance, cor=T) # principal components analysis using correlation matrix
  pc.comp <- pca$scores
  pc.comp1 <- -1*pc.comp[,1] # principal component 1 scores (negated for convenience)
  pc.comp2 <- -1*pc.comp[,2] # principal component 2 scores (negated for convenience)
  xy = cbind(pc.comp1, pc.comp2)
  
  # Check if K is ok
  maxPoints <- table(!duplicated(rowSums(xy)))[2]
  if (!is.na(maxPoints) & maxFind > maxPoints){ #!is.na() in case there are no duplicates
    cat ("K choice is greater than data points available, choosing another K:",maxPoints, "\n")
    maxFind <- maxPoints
  } else {
    cat ("K choice is",maxFind, "\n")
  }
  ##
  wss <- sapply(1:maxFind, 
                function(k){kmeans(data, k)$tot.withinss})
  par(mfrow=c(2,1))
  plot(1:maxFind, wss,
       type="b", pch = 19, frame = F, 
       xlim = c(0,maxFind+1),
       
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares")
  plot.new()
}




exploreK <- function(mydata,maxFind) {  
  
  require(cluster)
  ## Convert data
  distance = dist(mydata,method = "euclidean")
  pca <- princomp(distance, cor=T) # principal components analysis using correlation matrix
  pc.comp <- pca$scores
  pc.comp1 <- -1*pc.comp[,1] # principal component 1 scores (negated for convenience)
  pc.comp2 <- -1*pc.comp[,2] # principal component 2 scores (negated for convenience)
  xy = cbind(pc.comp1, pc.comp2)
  
  # Check if K is ok
  maxPoints <- table(!duplicated(rowSums(xy)))[2]
  if (!is.na(maxPoints) & maxFind > maxPoints){ #!is.na() in case there are no duplicates
    cat ("K choice is greater than data points available, choosing another K:",maxPoints, "\n")
    maxFind <- maxPoints
  } else {
    cat ("K choice is",maxFind, "\n")
  }
  # --
  par(mfrow=c(3,4))
  #n = ncol(mydata)
  n = maxFind
  # --
  for (i in 1:n){
    ChooseK <- i
    fit = kmeans(xy, ChooseK, iter.max=10000, nstart=10)
    clusplot(xy, fit$cluster, color=TRUE, shade=TRUE,  lines=0, plotchar=F, 
             main = paste(ChooseK,'Clusters',sep=" "))
    barplot(table(fit$cluster))
  }
  par(mfrow=c(1,1))
}


### kFinding end//


obtainClusters <- function (mydata,decidedK){
  require(cluster)
  ## Convert data
  distance = dist(mydata,method = "euclidean")
  pca <- princomp(distance, cor=T) # principal components analysis using correlation matrix
  pc.comp <- pca$scores
  pc.comp1 <- -1*pc.comp[,1] # principal component 1 scores (negated for convenience)
  pc.comp2 <- -1*pc.comp[,2] # principal component 2 scores (negated for convenience)
  xy = cbind(pc.comp1, pc.comp2)
  # --
  par(mfrow=c(1,1))
  decidedfit = kmeans(xy, decidedK, iter.max=100, nstart=1000)
  clusplot(xy, decidedfit$cluster, color=TRUE, shade=TRUE,  lines=0, plotchar=F,
           main = paste(decidedK,'Clusters',sep=" "))
  clusters <- lapply(split(decidedfit$cluster,decidedfit$cluster),names)
  return(clusters)
}





######################################
######################################
## GO term enrichment using goseq
findGOTerms <- function(GeneLists,Universe,GFFfile,GOFile) {
  cat ("-- findGOTerms function called \n")
  
  #
  require(goseq)
  library(rtracklayer)
  library(Rsamtools)
  require(GenomicRanges)
  ## Read and reduce GFF file to get length of sequences
  GFF <- import.gff(GFFfile,version="3",feature.type="gene")
  grl <- reduce(split(GFF, mcols(GFF)$Name))
  reducedGTF <- unlist(grl, use.names=T)
  mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
  reducedGTF
  AllLengths <- width(reducedGTF)
  names(AllLengths) <- mcols(reducedGTF)$Name
  # --
  head(reducedGTF)
  head(AllLengths)
  # --
  
  ## Read file with GO terms
  go <- read.delim(GOFile,header = F)
  
  colnames(go) <- c(
    "locusName", "TAIRAccession", "AGI", "GOrelationship", "GOterm",
    "GOID", "TAIRKeywordID", "Aspect", "GOslimTerm", "EvidenceCode",
    "EvidenceDesc", "EvidenceWith", "Reference", "Annotator", "Date"
  )
  
  ## -- Select only AGI and GO id
  go.goseq <- go[,c("AGI", "GOID")]
  head(go.goseq)
  
  # -- Loop on each cluster
  GOList <- list()
  for (each in names(GeneLists)){
    print (each) 
    
    ##
    
    de.genes <- GeneLists[[each]]
    length(Universe)
    length(de.genes)
    
    # Create a binary vector of genes present/absent in the set vs the universe.
    gene.vector=as.integer(Universe%in%de.genes)
    names(gene.vector)=Universe
    head(gene.vector)
    table(gene.vector)
    
    # Gene lengths for the list
    GeneLengths <- AllLengths[names(gene.vector)]
    cbind(GeneLengths,gene.vector)
    ##
    
    all.genes <- rownames(gene.vector)
    pwf=nullp(gene.vector, "tair10", id=all.genes, bias.data=GeneLengths)
    head(pwf)
    
    # 
    go.all <- goseq(pwf, "tair10", gene2cat=go.goseq)
    go.all <- go.all[!is.na(go.all$term),] #remove NAs
    
    head(go.all)
    dim(go.all)
    
    table(go.all$over_represented_pvalue < 0.05)
    
    
    
    
    #View(go.sign)
    GOList[[each]] <- go.all
  }
  return(GOList)
}


#############################################
### Create a larger table with all columns 

condenseListTables <- function(listDFs) {
  cat ("-- make condenseListTables function called \n")
  # First make an empty table
  uniqRows <- Reduce(union,lapply(listDFs,rownames))
  uniqCols <- unlist(sapply(listDFs, colnames))
  zeroTable <- as.data.frame(matrix(NA,
                                    nrow = length(uniqRows),
                                    length(uniqCols)),
                             row.names =uniqRows)
  colnames(zeroTable) = uniqCols
  # Then Fill it
  for (each in names(listDFs)){
    cat ("Filling binary table:", each,"\n")
    zeroTable[rownames(listDFs[[each]]),
              colnames(listDFs[[each]])] <- listDFs[[each]]
  }
  #zeroTable <- apply(zeroTable,c(1,2),as.numeric)
  #zeroTable[is.na(zeroTable)] <- 0
  return(zeroTable)
  cat ("-- binary table done \n")
  cat ("\n")
}


#############################################
###### Draw heatmap function
# Some palettes
hmCols <- colorRampPalette(c("darkgoldenrod1","white", "deepskyblue3"))#-colorRampPalette(c("white","royalblue"))
hmColsTwo <- colorRampPalette(c("white", "steelblue"))#-colorRampPalette(c("white","royalblue"))
hm.OGWCYP <- colorRampPalette(
  c("gold", "orange",
    "white",
    "cyan","royalblue"),
  space="Lab")
#
hm.WCRGO <- colorRampPalette(
  c("white","cyan","royalblue","gold", "orange"),
  space="Lab")
## Heatmap function
hmFunction <- function(hmMatrix,clrPalette,relSize,hmRows,hmName) {
  
  # load
  require(gplots)
  # 
  hmMatrix <- as.matrix(hmMatrix)
  mylmat = rbind(c(4,3,0),
                 c(2,1,0),
                 c(0,5,0)) # creates 3x3 table with location of heatmap elements defined
  mylwid = c(1,4,0.1)
  mylhei = c(0.15,4,0.4)
  par(#mar=c(4,2,4,2),
    oma=c(2,4,1,0),
    cex.main=1)
  heatmap.2(hmMatrix, dendrogram ='row',
            Colv=F, col=clrPalette, 
            key=T, keysize=0.9, symkey=T, density.info='none',
            trace='none', 
            colsep=rep(1:ncol(hmMatrix)),
            
            #rowsep=rep(1:nrow(hmData)),
            sepcolor='white', sepwidth=c(0.025),
            scale="none",
            cexRow=relSize,#srtRow = 45,
            cexCol=2,srtCol=45,
            labRow = hmRows,
            #
            margins = c(2,8),
            #
            #hclustfun=function(c){hclust(c, method='mcquitty')},
            lmat=mylmat, 
            lhei = mylhei, 
            lwid = mylwid,
            main=hmName)
}


################
######## Write genes to clusters

condenseGeneList_toMatrix <- function(listOfClusters){
  ## Takes in a list of clusters, each element has genes, converts it into an uneven matrix, each column is a cluster, filled with the genes corresponding to that cluster. 
  #Clusters don't have the same number of elements (hence, uneven)
  
  # To avoid having errors, output might be an empty matrix of 1 row. 
  totElements <- max(sapply(listOfClusters, length))
  if(totElements==0){totElements=1}
  
  totClusters <- length(listOfClusters)
  matrixClusters <- matrix("",nrow = totElements,ncol = totClusters,dimnames = list(seq(1,totElements),names(listOfClusters)))
  
  for (each in names(listOfClusters)){
    geneList <- listOfClusters[[each]]
    if(length(geneList) == 0) {geneList=""}
    matrixClusters[1:length(geneList),each] <- geneList
  }
  return(matrixClusters)
}

######## Custom palettes
colorPalette <- c(
  "#74909f",
  "#8baea8",
  "#a79b7d",
  "#b48080",
  "#7d7f96",
  "#6b6d5a",
  "#6a836f",
  "#83b7db",
  "#6897bb",
  "#f37054",
  "#0000cd",
  "#ccccff",
  "#30d5c8",
  "#ec9a7e",
  "#7bceea",
  "#fff498",
  "#981010",
  "#5e7143",
  "#06357a",
  "#ec891d",
  "#ffc125",
  "#005555")

#testPalette <- sample(colorPalette,6,replace = F)
#pie(rep(1,length(testPalette)),col = testPalette,)


######## Custom palettes
colorPalette <- c(
  "#74909f",
  "#8baea8",
  "#a79b7d",
  "#b48080",
  "#7d7f96",
  "#6b6d5a",
  "#6a836f",
  "#83b7db",
  "#6897bb",
  "#f37054",
  "#0000cd",
  "#ccccff",
  "#30d5c8",
  "#ec9a7e",
  "#7bceea",
  "#fff498",
  "#981010",
  "#5e7143",
  "#06357a",
  "#ec891d",
  "#ffc125",
  "#005555")

#testPalette <- sample(colorPalette,6,replace = F)
#pie(rep(1,length(testPalette)),col = testPalette,)
