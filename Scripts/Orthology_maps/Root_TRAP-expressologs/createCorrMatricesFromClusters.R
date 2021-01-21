# Expressolog Analysis, step 1

# This script was adapted from https://github.com/VinLau/expressolog-guide 

# It uses expression values of homologous cell/tissue types and creates all possible between species correlation metrices 
# for each gene family (i.e., cluster), based on Phytozome v12 gene family file that was updated with ITAG3.2 orphan genes


#NB### SET YOUR WD HERE
setwd("~/mywd/")
#NB### 

#NB### INSTALL THESE PACKAGES IF NOT ALREADY INSTALLED
library("dplyr")
library("tidyverse")
library(data.table)
#NB###
#REMOVE TOMATO GENOME BUILD VERSION

removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}

#NB### NAME YOUR OUTPUT FILE NAME HERE (WARNING, THIS SCRIPT WILL DELETE A FILE NAMED THIS IN YOUR DIR)
flname <- "Sly-Ath2-expressologCorMatrices_phyto12ITAG3.2.csv"
#flname <- "Sly-Osa-expressologCorMatrices_phyto12ITAG3.2.csv"
#flname <- "Ath-Osa-expressologCorMatrices_phyto12ITAG3.2.csv"

#NB### 
if (length(commandArgs(trailingOnly = TRUE)) >= 1){ flname <- commandArgs(trailingOnly = TRUE)[1][1] }

if (file.exists(flname)){ file.remove(flname) }

#TOMATO
tomatoCSVDataTPM <- read.csv(file="Data/20180920_Sl_ITAG3.2_Kallisto_quantile_norm_median_tpm_CompCT.csv",
                             header=TRUE, sep=",") #35768

#NB### THESE LINES ARE FOR CLEANING TPM DATA, COMMENT THESE OUT IF YOU DONT NEED TO FILTER
zeroClndTomCSV <- tomatoCSVDataTPM[rowSums(tomatoCSVDataTPM[,-1]) > 0,] #23055
fullyFltredTomCSV <- filter_all(zeroClndTomCSV, any_vars(. > 2)) #retains any rows that has any value >2, 19034 row
#NB### 

#NB### THIS LINE REMOVES ROWNAMES, PARTICUARLY THE FIRST CELL
tomatoGenesExpr <- fullyFltredTomCSV %>% remove_rownames() %>% column_to_rownames(var="SlGeneID")
row.names(tomatoGenesExpr) <- removeDotGene(row.names(tomatoGenesExpr))


#NB### THIS LINE THEN REPLACES THE OLDER 'TISSUE' COLUMN NAMES WITH A UNIVERSAL ONE SO THAT IT IS EASIER FOR MERGING IN THE FUTURE
colnames(tomatoGenesExpr) <- c("COR", "EN", "EP", "MCO", "MZ", "PH", "V", "QC", "x35S")
#NB### 

#NB### SAME LOGIC ABOVE APPLIED TO OTHER SPECIES' EXPR DATA

#RICE
riceCSVDataTPM <- read.csv(file="Data/20200313_Os_Atlas_IRGSP-1.0_TRANSCRIPT_Kallisto_qnorm_MeanMedian_tpm_genes_wQC.csv",
                             header=TRUE, sep=",") #TRANSCRIPT 37849 (CDS 35667 genes)
zeroClndriceCSV <- riceCSVDataTPM[rowSums(riceCSVDataTPM[,-1]) > 0,] #TRANSCRIPT 33818 (CDS 30998 genes)
fullyFltredriceCSV <- filter_all(zeroClndriceCSV, any_vars(. > 2)) #TRANSCRIPT 23368 (CDS 23527 genes)

# #Convert RAP to MSU's LOC_Os IDs using the conversion table in https://rapdb.dna.affrc.go.jp/download/irgsp1.html
# Osconvert <- read.csv("Data/RAP-MSU_2019-08-29.txt", header = F, sep = "\t",stringsAsFactors = F)
# RAP_fullyFltredriceCSV<-merge(fullyFltredriceCSV, Osconvert, by.x = "X", by.y = "V1")
# RAP_fullyFltredriceCSV$MUSuniq <- removeDotGene(RAP_fullyFltredriceCSV$V2)
# SaveRAP_fullyFltredrice <- RAP_fullyFltredriceCSV[,-8]
# write.csv(SaveRAP_fullyFltredrice, "20200313_Os_RAP-MUS_filtered_meanMedian_qnorm_tpm_genes_TRANSCRIPT.csv", quote = F, row.names = F)
# 
# #MUS = none (3236rows) & duplicates (774 genes) were removed (19359 genes left)
# RAPdrop_FltredriceCSV <- RAP_fullyFltredriceCSV[!RAP_fullyFltredriceCSV$MUSuniq=="None",c(2:7,9)]#20132 (3236 genes don't have a MUS ID)
# #
# uniqRAPdrop_FltredriceCSV <-RAPdrop_FltredriceCSV[which(!duplicated(RAPdrop_FltredriceCSV$MUSuniq)),] #19920 (335 duplicated genes were excluded)

uniqRAPdrop_FltredriceCSV <- read.csv("Data/20200313_Os_RAP-MUS_filterCurated_meanMedian_qnorm_tpm_genes_TRANSCRIPT.csv", header = T)

riceGenesExpr <- uniqRAPdrop_FltredriceCSV[,-1] %>% remove_rownames() %>% column_to_rownames(var="MUSuniq")
colnames(riceGenesExpr) <- c("MCO","EN","V","MZ","QC","x35S")

#ARABIDOPSIS
athCSVData <- read.csv(file="Data/At_Mustroph_RMA_mean_unlogged_COMP_CT.csv", header=TRUE, sep=",") #17418 genes
athGenesExpr <- athCSVData %>% remove_rownames() %>% column_to_rownames(var="AtGeneID")
colnames(athGenesExpr) <- c("x35S", "V", "EN", "MZ", "EP", "PH", "MCO", "COR")
### 

### IF YOU NEED TO DELETE A COLUMN (BECAUSE AN EQUILVALENT TISSUE DOES NOT EXIST IN THE OTHER EXPR DATA)
#Comparing tomato & Arabidopsis 
tomatoGenesExpr$QC = NULL

#Comparing tomato & rice 
tomatoGenesExpr[,c(1,3,6)] = NULL

#Comparing Arabidopsis & rice 
 riceGenesExpr[,5] = NULL
 athGenesExpr[,c(5,6,8)] = NULL

#phyto12virid with ITAG3.2 orphans- only Sl, At and Os orthologs
phyto <- read.table("Data/phyto12update_add3.2_Evalue0.01_SUB_Sl-At-Os.txt", header = T, stringsAsFactors = F) 
#stringsAsFactors = F: importent to keep genes as character and not factor (for sapply)


#NB### NAME YOUR CELL NAMES, WHAT YOU HAD NOTED ABOVE
#FOR TOMATO AND ARABIDOPSIS
cellNames <- c("COR", "EN", "MZ", "x35S", "V", "EP", "PH", "MCO")

#FOR TOMATO AND RICE
cellNames <- c("EN", "MZ", "x35S", "V", "MCO", "QC")

#FOR ARABIDOPSIS AND RICE
cellNames <- c("EN", "MZ", "x35S", "V", "MCO")

#Assign genes to clusters
# uniqu_cluster <- unique(phyto$ClusterID) #create a character of unique clusters 54,709
# uniqu_cluster <-uniqu_cluster[-140] #54,708 (one of the clusters is NA!)
# 
# cluster_df <- data.frame(uniqu_cluster=uniqu_cluster) # create a dataframe with the unique cluster names
# cluster_df[,"Genes"] <- NA #add a column for the genes

############################################################
#Match each uniqu_cluster with phyto$ClusterID and append the corresponding genes to column "Genes"
#commented out because can take forever to run...
# cluster_df$Genes <- sapply(cluster_df$uniqu_cluster,function(x){
#   phyto[grep(x,phyto$ClusterID),"Genes"]})
############################################################
#saveRDS(cluster_df,"20190925_Sl-At-Os_byClusters_phytoV12_ITAG3.2orphan_Eval01.RDS")


# Desired matrix for correlation matrix output
#       GeneA GeneA` GeneA```
# 35S  100  432    322
# EN   25    231    213
# I.e. the cell types are collated together to be the observations whilst Genes are variables

### If using phytozome gene families data, use this block
#Loop through every gene
#############################################
# tmp <- cluster_df[1:6,,drop=F]
# View(tmp$Genes[[4]])

cluster_df <- readRDS("Data/20190925_Sl-At-Os_byClusters_phytoV12_ITAG3.2orphan_Eval01.RDS")

by(cluster_df, 1:nrow(cluster_df), function(row) { 
  #row <- tmp$Genes[[3]]
  genes <- unlist(row)
  
  # For each pairwise comparison delete one species: LOC_Os/ Solyc/ AT
  #NB### Below line deletes Rice containing genes in the Phytozome cluster, add more of these if you need to delete more species
  genes <- genes[lapply(genes, function(x) length(grep("LOC_Os", x, value=FALSE))) == 0]
  
  #NB### 
  
  # Create a empty matrix for correlation calculation
  geneCellMtx <- matrix(nrow=8,ncol=0) #8:Sl&At, 6:Sl&Os, 5:At&Os
  rownames(geneCellMtx) <- cellNames
  
  # loop through every gene in a cluster
  for(gene in genes){
    tempDF <- matrix()
    #NB###  Make sure to change the regex below to the species you want:  if (grepl("Solyc"/"LOC_Os"/"AT"
    if (grepl("Solyc", gene)) { 
      gene <- removeDotGene(gene) #Use this line for tomato!
      tempDF <- data.frame(t(tomatoGenesExpr[gene,])) #tomatoGenesExpr/athGenesExpr
    }
    else if (grepl("AT", gene)){
      tempDF <- data.frame(t(athGenesExpr[gene,])) #riceGenesExpr/athGenesExpr
    }
    #NB### 
    
    if (anyNA(tempDF)) {next} # Check if we even have expr values for the gene
    
    # Check if column of expression values has already been added (e.g. Ath expr for single gene for multiple transcripts)
    if (!colnames(tempDF)[1] %in% colnames(geneCellMtx)) { 
      geneCellMtx <- transform(merge(geneCellMtx, tempDF, by="row.names", all = TRUE ), row.names=Row.names, Row.names=NULL)
    }
    
    # print(tempDF)
  }
  # print(geneCellMtx)
  # print(cor(geneCellMtx))
  #NB### Set your correlation method here, default is PCC
  write.table(cor(geneCellMtx), sep = ",",file=flname, append=T, col.names=NA)
  #NB###
})
##########################################################################################
#NEXT STEP: 
# 1. go to https://github.com/VinLau/expressolog-guide
# 2. Parse correlation matrices: clean up the above csv file by using the makeTSVFromCorrMtcs.py file.
# This file will take the correlation matrices and can find the best matching (best expressolog)
# in each gene in each cluster and produce a TSV file.