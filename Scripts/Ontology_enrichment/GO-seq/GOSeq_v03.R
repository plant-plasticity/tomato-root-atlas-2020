
### SET YOUR WD HERE
setwd("~/mywd/")
## INSTALL THESE PACKAGES IF NOT ALREADY INSTALLED

## https://pathwaycommons.github.io/guide/primers/statistics/fishers_exact_test/

## Load libraries and Functions
{
  ########
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(c("GenomicRanges"))
#  biocLite(c("rtracklayer"))
#  biocLite(c("goseq"))
  library(goseq)
  library(rtracklayer)
  library(GenomicRanges)
  library(Rsamtools)
  library(purrrlyr) #install.packages("purrrlyr")
}

source("GOseq/omicFunctions/metaFunctions_forNetworkAnalysis.R")
source("GOseq/omicFunctions/generic_functions.R")
source("GOseq/omicFunctions/GOseq_functions.R")

differentGenomicBuilds <- T
######


  ##########################################################################################
  ## UPLOAD ENRICHED GENE LISTS (TABLE S1, WORKSHEET "Cell type lists")
  {  
  GeneList <- read.table("Cell type lists.csv",sep=",",header = T,as.is = T)
  GenesPerCategory <- apply(GeneList,2,function(x){length(x[x!=""])})
  names(GenesPerCategory) <- gsub(".FALSE","",names(GenesPerCategory)) #modify the column names
  }
  #########################
  
barplot(GenesPerCategory,ylim = c(0,max(GenesPerCategory)*1.25),las=2)


######## Create your universe: i.e. all genes that are expressed.
allgenes <- read.csv("GOseq/meta/20180920_ATLAS_ITAG3.2_Kallisto_quantile_norm_median_tpm.csv", sep=",",header = T, row.names = 1)
if (differentGenomicBuilds) row.names(allgenes) <- removeDotGene(row.names(allgenes))

myuniverse <- allgenes[which(rowSums(allgenes)>0),] # consider only genes that the sum of row is > 0. 

## All genes
 assayed.genes <- rownames(myuniverse) # create a list of genes out of the table of expressed genes
  ###

 ######### Extract TRANSCRIPT LENGTH (i.e. effective length) from Kallisto counts
 tmp_table <- read.csv("GOseq/meta/KallistoCounts_G01abundance.tsv", header = T, sep = "\t") #row.names = 1, 
 
 efflength <- tmp_table$eff_length
 names(efflength) <- tmp_table$target_id 
 if (differentGenomicBuilds) names(efflength) <- removeDotGene(names(efflength)) #transcript length based on Kallisto counts
 
######### Prepare list of GO terms
####################################

{ # If reading the 3.2 version
  go <- read.delim("GOseq/meta/ITAG3.2_protein_go_Split_jrm.tsv",header = F,stringsAsFactors = F)
  colnames(go) <- c("ITAG",
                    "GOID")
  
  go.goseq <- go[,c("ITAG", "GOID")]
  head(go.goseq)
}

######### GO terms for each cluster
####################################################################################################

### This script assumes you've created a matrix of genes of interest (GeneList) using the condenseGeneList_toMatrix() function.
# This function essentially creates a table with each column contains genes that were significant or of interest for a particular experiment
#(ie, a contrast or the result of clustering)
# The function takes a list of genes and creates an uneven (not symmetrical) matrix. 

GOList <- list()
maxN <- ncol(GeneList)
for( each in seq(1,maxN)){
  print (each) 
  #tmp <- GeneList[GenePatterns[,each]!="",each]
  
  #assayed.genes <- rownames(tmp)
  ### Since the GeneList matrix is uneven, the function fills in blanks with empty quotes "". 
  de.genes <- GeneList[GeneList[,each]!="",each] 

  # goseq needs a binary vector with 0 corresponding to genes not DE/interest and 1 being a DE gene (or gene of interest)
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  table(gene.vector)
  
  ## Get the gene lengths only of the genes assayed 
  # (the AllLengths contains information for ALL the genes in the genome but not all of them are *normally* used)
  GeneLengths <- efflength[names(gene.vector)] #transcript length based on Kallisto counts
  
  ##
  #cbind(GeneLengths,gene.vector)
  
  all.genes <- names(gene.vector)
  
  ########## GOseq functions
  ##################################################
  ## nullp creates a weight probability of a gene based DE based on it's length
  pwf=nullp(gene.vector, bias.data=GeneLengths)
  head(pwf)
  
  if (differentGenomicBuilds) rownames(pwf) <- removeDotGene( rownames(pwf) )
  if (differentGenomicBuilds) go.goseq$ITAG <- removeDotGene ( go.goseq$ITAG)  
  
  
  go.all <- goseq(pwf, "ITAG3.1", gene2cat=go.goseq)
  ################################################################################
  
  ### Filter significant categories
  #table(go.all$over_represented_pvalue < 0.05)
  #go.sign <- go.all[go.all$over_represented_pvalue < 0.05,]
  
  # Save to a list
  #View(go.sign)
  GOList[[as.character(each)]] <- go.all
}


######### Filter by Significant Categories
####################################################################################################
sapply(GOList, dim)
signGOList <- lapply(GOList, function(x){  x[x$over_represented_pvalue < 0.05,] })
sapply(signGOList, nrow)
names(signGOList) <- c("EN","EP","EXO","gCOR","iCOR","MiCO","MZ","PH","QC","V","XY")

# EN   EP  EXO gCOR iCOR  MCO   MZ   PH   QC    V   XY
# 1  2  3  4  5  6  7  8  9 10 11 
# 34 61 55 58 43 21 29 58 48 31 34

##### Get DE genes on each category
#########################
library(tidyverse)
getGenesPerCategory <- lapply(seq_along(signGOList),function(DP){
  cat("Genes per category in pattern", DP,"\n")
  
  ## Get DE genes
  de.genes <- GeneList[GeneList[,DP]!="",DP]
  
  ### Get categories
  categories <- signGOList[[DP]][,"category"]
  
  ## From the full GO seq list, get the genes that overlap with the DE
  DEgenesInCategory <- go.goseq %>%
    filter(.$GOID %in% categories) %>% group_by(GOID) %>% filter(ITAG %in% removeDotGene(de.genes)) %>% nest()
  ### Convert the table with significant GO categories to a tibble and join the list-columns table
  test <- as.tibble(signGOList[[DP]])
  test <- test %>% 
    left_join(DEgenesInCategory,by=c("category"="GOID")) %>% mutate("Set"=DP) #%>% rename("DEGenesinCategory"="data") 
  return(test)
})
names(getGenesPerCategory) <- names(signGOList)
names(getGenesPerCategory) 

### Make a single table
#####################
SlsignificantGOtibble <- bind_rows(getGenesPerCategory) %>% 
  mutate(uniqueID = row_number())

nrow(SlsignificantGOtibble) == sum(sapply(getGenesPerCategory, nrow)) # Should be true


## Fix terms with no name associated
SlsignificantGOtibble <- SlsignificantGOtibble %>%  
  mutate("term"=ifelse( is.na(term),category,term  ) )

## Save
save(SlsignificantGOtibble,file = "FILENAME.RData")
SlysigGOterm <- subset(SlsignificantGOtibble, select = c("over_represented_pvalue", "numDEInCat", "numInCat", "category", "term", "Set", "ontology")) #,"data"
write.csv(SlysigGOterm, "FILENAME.csv")

################################

######### Write genes from GO categories to files
####################################################################################################
##### Create directories to save outputs
{
  resultsPath <- "GenesPerCategory/"
  itagPath <- paste0(resultsPath,"ITAG_byCategory/")
  #orthoPath <- paste0(resultsPath,"AthOrthologues_byCategory/")
  dir.create(itagPath,recursive = T) 
  #dir.create(orthoPath,recursive = T)
}

##############################
unite_(SlsignificantGOtibble,"tmpName",c("term","ontology","Set"),sep="_",remove=F) %>%  ## Create new names for the
  mutate("itagFile"=paste0(itagPath,"/",gsub("[[:space:]]|/|>","-",tmpName),".txt")) %>%
  by_row(~write.table(.$DE_ITAGsinCats, file = .$itagFile,sep = "\t",row.names = F,quote = F)) %>%
  mutate("orthoFile"=paste0(orthoPath,"/",gsub("[[:space:]]|/|>","-",tmpName),".txt")) %>%
  by_row(~write.table(.$Sly2AGI_orthologues, file = .$orthoFile,sep = "\t",row.names = F,quote = F))

################################################################################################################################


########### Subset the p values of overrepresented categories and split by ontology
############################################
{    
  ## Get only two columns
  go2 <- lapply(signGOList, function(go){
    #go <- signGOList[[1]]
    go[is.na(go$term),"term"] <- go[is.na(go$term),"category"]
    rownames(go) <- go$term
    go <- go[,c("ontology","over_represented_pvalue")]
  })
  
  ## Split by category
  splitted <- lapply(seq_along(go2), function(x){
    tmpName <- names(go2)[[x]]
    go3 <- go2[[x]]
    go3
    if (nrow(go3) != 0){
      go3 <- split(go3,go3$ontology)
      names(go3) <- paste0(tmpName,".",names(go3))
      return(go3)
    }
  })
  
  splitted <- unlist(splitted,recursive = F)
  names(splitted)
  
  BioProcess <- splitted[grep("BP",names(splitted))]
  MolFunction <- splitted[grep("MF",names(splitted))]
  CellComponent <- splitted[grep("CC",names(splitted))]
  
  ##
  BioProcess <- lapply(seq_along(BioProcess), function(x){
    DF <- BioProcess[[x]][,2,drop=F]
    colnames(DF) <- gsub("\\...","",names(BioProcess)[[x]])
    return(DF)
  })
  names(BioProcess) <- lapply(BioProcess,colnames)
  ##
  MolFunction <- lapply(seq_along(MolFunction), function(x){
    DF <- MolFunction[[x]][,2,drop=F]
    colnames(DF) <- gsub("\\...","",names(MolFunction)[[x]])
    return(DF)
  })
  names(MolFunction) <- lapply(MolFunction,colnames)
  ##
  CellComponent <- lapply(seq_along(CellComponent), function(x){
    DF <- CellComponent[[x]][,2,drop=F]
    colnames(DF) <- gsub("\\...","",names(CellComponent)[[x]])
    return(DF)
  })
  names(CellComponent) <- lapply(CellComponent,colnames)
  
  ##### Convert to tables
  ###############
  BioProcessDF <- condenseListTables(BioProcess)
  MolFunctionDF <- condenseListTables(MolFunction)
  CellComponentDF <- condenseListTables(CellComponent)
}
dir.create("Results",showWarnings = T)
write.csv(BioProcessDF,file = "FILENAME/GOenrich_BioProcess.csv",quote = F,row.names = T)
write.csv(MolFunctionDF,file = "FILENAME/GOenrich_MolFunction.csv",quote = F,row.names = T)
write.csv(CellComponentDF,file = "FILENAME/GOenrich_CellComponent.csv",quote = F,row.names = T)

