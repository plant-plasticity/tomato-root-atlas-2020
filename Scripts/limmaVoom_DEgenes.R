## RNA seq analysis using limma-voom
# jrodriguezm at ucdavis dot edu
# github: rodriguezmDNA
# brady lab

# last update 2020.03.03 - LSM

######################## User defined options ########################
######################################################################
# Main directory:
setwd("/mywd/") #Full path of the working directory. 
#It must contain: 
# ** A directory named "Counts" with a delimited file for raw counts. Rows are genes and columns samples.
# ** A directory named "meta" with a metadata file with information about the samples. Ideally the number of rows in the metadata is the same as in the raw counts.
# ** A directory named Scripts with this script and the 'functions.R' script.


## Metadata options
metaFile <- "20180814_atlas-SL-meta_final.csv" #Name of metadata file in csv format.
doFilter <- F #T
#whichFilter <- c("PH.T","XY.TQ","WOX.Ta","XY.I","WOX.Tb","MCO.S")

## Counts file name (with extension) 
countsFile <- "20180814_ATLAS_ITAG3.2_Kallisto_raw_counts.csv" #Name of counts file

shortName <- "191114-Sl_BatchSUBTRACTED" #Short name to append at the end of the filenames.#191107-Sl_BatchCorrect
#If missing it will append the name of the folder where the scripts where run.

## Filter genes with low expression using CPM counts.
filterByCPM <- T #T
CPMcutoff <- 0.5 #2


## pValue (default = 0.05 ) and absolute logFC (default = 2) to color genes on volcano plots
pValCut=0.05  #0.05  
logCut=2 #2

########################
########################

###

library(edgeR) #which loads limma as a dependency
library(reshape)
library(gplots)
library(RColorBrewer)
library(calibrate)
library(Glimma) 
library(viridis)
library(pheatmap)
library(ggplot2)

# ## Output
 outDir = "outDirName/" 
 dir.create(outDir, showWarnings=T)
# 
 sink('sinkName.txt')
# 
 geneListsDir = "geneListsDirName/GeneLists"
 dir.create(geneListsDir, showWarnings=T)
# #
 imgDir = "imgDirName/images/"
 dir.create(imgDir, showWarnings=T)
# ## --

if (is.na(shortName)){
  shortName <- basename(getwd())
}

# Load functions
source("~/Data/functions.R")
######## --- --- --- 


## Start of analysis
####################################################################################
####################################################################################
cat("Reading metadata file \n")

meta <- metaDataProcessing(metaFile,doFilter,whichFilter)
head(meta)

#
cat("Reading counts file:",countsFile,"\n")

GeneCounts <- read.csv(paste0("Counts/",countsFile),row.names = 1)
dim(GeneCounts)

## Check that samples in both counts and metadata are the same.
## Use function filterCounts(counts,meta)
tmp <- filterCounts(GeneCounts,meta)
GeneCounts <- tmp[["counts"]]
meta <- tmp[["meta"]]
rm(tmp)
## --


###### Design matrix
## Convert experimental metadata to factors for the design
experimentFactors <- lapply(apply(meta,2,split,""),unlist)
experimentFactors <- as.data.frame(lapply(experimentFactors,as.factor))

cat ("Create the design with these factors:\n")
print(head(experimentFactors))

###  User modified:
####Simplest design taking into account all possible interactions
Groups <- as.factor(paste0(experimentFactors$Promoter))

###########################################################
# Accounting for batch effect 
###########################################################
cat("accounting for batch effect \n")
Batch <- as.factor(paste0(experimentFactors$Replicate))
print(head(Batch))

design <- model.matrix(~0+Groups) 
# Example of an interaction
#design <- model.matrix(~0+experimentFactors$Sample*experimentFactors$Treatment) #Sample*Treatment interaction

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("Groups","experimentFactors","\\$","\\:","\\-",
                  colnames(experimentFactors)),sep="",collapse = "|")

colnames(design) <- gsub(fixCols,"",colnames(design))
head(design)


####################################################################################
cat("Removing genes with 0 counts on all conditions \n")
cat("Initial number of genes:",nrow(GeneCounts),"\n")
rmIDX <- which(rowSums(GeneCounts) == 0)
cat("Removing",length(rmIDX),"genes \n")
GeneCounts <- GeneCounts[-rmIDX,]
cat("Remaining number of genes:",nrow(GeneCounts),"\n")

cat("Use cpms to uncover lowly expressed genes \n")
dge <- DGEList(counts=GeneCounts,remove.zeros = F)


# Filter genes with low CPMs accross replicates 
cat("Replicates of samples range between:", range(table(Groups)),"\n")

#
if (filterByCPM){
  
  sampleMin <- min(table(Groups))
  cat("Filtering reads with low CPMs ( <",CPMcutoff,") in at least",sampleMin,"replicates \n")
  #
  cpm <- cpm(dge)
  keep.exprs <- rowSums(cpm>CPMcutoff)>=sampleMin
  table(keep.exprs)
  

  cat("Removing",table(keep.exprs)[1],"genes \n")
  cat("Remaining number of genes:",table(keep.exprs)[2],"\n")
  
  #
  
  y <- dge[keep.exprs, , keep.lib.size = FALSE]
} else {
  cat("Not doing CPM filtering")
}

## Normalization with TMM
#y <- calcNormFactors(y, method = "TMM")
# A paper from Gordon K Smyth, who developed limma uses TMM before voom:
# http://dx.doi.org/10.12688/f1000research.9005.1
#####


normalizedExpression <- cpm(y)
colnames(normalizedExpression) <- meta$PromoterLibrary

##########################################################################################################
# Convert to CPM and log2 transformation- FOR VISUALIZATION PURPOSE
##########################################################################################################
cat ("Convert to CPM and log2 transformation- FOR VISUALIZATION PURPOSE \n") 
normalizedlogExpression <- cpm(y, log=TRUE, prior.count=3) # log=TRUE, prior.count=3 based on: https://support.bioconductor.org/p/76837/
#we use a larger prior count to down-weight the contribution of low-abundance genes by squeezing their relative abundances
#toward the mean for each gene. This allows algorithms to focus more on the higher-abundance genes,
#whose measurements have better precision, without explicitly weighting the values.
#IMPORTANT NOTE: Source for differences between MDS plots of Kaisa and mine!?!

colnames(normalizedlogExpression) <- meta$PromoterLibrary

#Save a file with log2 normalized CPM values- FOR VISUALIZATION PURPOSE
tmpSave <- paste(outDir,"normalizedlogExpression_WITH batch","_",shortName,".csv",sep="")
cat("Saving log2normalized data to:", tmpSave, "\n")
save(normalizedlogExpression,file = "normalizedlogCPM_WITH batch.RData")
write.csv(x=normalizedlogExpression,tmpSave,quote = F,row.names = T)

#Save a file with normalized CPM values
tmpSave <- paste(outDir,"normCPM_WITH batch ","_",shortName,".csv",sep="")
cat("Saving normalized data to:", tmpSave, "\n")
save(normalizedExpression,file = "normCPM_WITH batch.RData")
write.csv(x=normalizedExpression,tmpSave,quote = F,row.names = T)

###########################################################
# Remove batch effect
###########################################################
cat("Remove batch effect \n")

logCPM_batchRemoved <- removeBatchEffect(normalizedlogExpression, batch=Batch, design=design) 

colnames(logCPM_batchRemoved) <- meta$PromoterLibrary

tmpSave <- paste(outDir,"logCPM_batchRemoved","_",shortName,".csv",sep="")
cat("Saving batch corrected data to:", tmpSave, "\n")
save(logCPM_batchRemoved,file = "normlogCPM_batchRemoved.RData")
write.csv(x=logCPM_batchRemoved,tmpSave,quote = F,row.names = T)


############################################################################################
# Exploratory analysis BEFORE batch effect correction- FOR VISUALIZATION PURPOSE
############################################################################################
cat("Exploratory analysis BEFORE batch effect correction: MDS & PCA plots, expression correlation, mean center and expression heatmaps \n")

## Easier visualization for MDS plots

cat("Using glimma for MDS plot visualization - Normalized data \n")

glMDSPlot(normalizedlogExpression, labels=meta$PromoterLibrary,
          groups=meta, launch=T, folder = "WithBatch/glimmaPlots_normlogCPM_WITHbatch")

##################################
# PCA

transposed_normalizedlogExpression <- t(normalizedlogExpression)
# perform prcomp on the transposed matrix from which columns (genes) of zero variance have been removed
pca_proc_withoutbatchcorr <- prcomp(transposed_normalizedlogExpression[,apply(transposed_normalizedlogExpression, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE)

cat("summary of pca_proc_withoutbatchcorr \n")
summary(pca_proc_withoutbatchcorr)

#Importance of components:
#                           PC1      PC2      PC3      PC4      PC5      PC6      PC7     PC8      PC9     PC10     PC11     PC12
# Standard deviation     70.1623 32.78613 31.5986 27.66331 26.1211 25.73227 24.74617 24.42381 24.27219 23.4196
# Proportion of Variance  0.2243  0.04899  0.0455  0.03488  0.0311  0.03018  0.02791  0.02719  0.02685  0.0250
# Cumulative Proportion   0.2243  0.27334  0.3189  0.35372  0.3848  0.41500  0.44291  0.47009  0.49694  0.5219

plotData_withoutbatchcorr = meta[,c("Name","Replicate","Promoter","PromoterLibrary")]
plotData_withoutbatchcorr$PC1 <- pca_proc_withoutbatchcorr$x[,1]
plotData_withoutbatchcorr$PC2 <- pca_proc_withoutbatchcorr$x[,2]
plotData_withoutbatchcorr$PC3 <- pca_proc_withoutbatchcorr$x[,3]
plotData_withoutbatchcorr$PC4 <- pca_proc_withoutbatchcorr$x[,4]
plotData_withoutbatchcorr$PC5 <- pca_proc_withoutbatchcorr$x[,5]

# 2D plots of pairs of principal components

pdf("WithBatch/PCA of Sl atlas norm log2CPM.pdf")
qplot(PC1,PC2,data=plotData_withoutbatchcorr,color=Replicate,shape=Promoter) +
  labs(x="PC1 (22% variability)",y="PC2 (5% variability)") +
  ggtitle("PCA without batch effect correction") +
  scale_shape_manual(values=c("EN"=0,"MCO"=1,"MZ"=2,"V"=3,"XY"=4,"EP"=5,"PH"=6,"EXO"=7,"ACT"=8,"X35S"=9,"WOX"=10,"COR"=11)) +
  geom_point(size = 3) 
dev.off()
############################################################################################
# Exploratory analysis AFTER batch effect correction- FOR VISUALIZATION PURPOSE
############################################################################################
logCPM_batchRemoved <- read.csv("./SlAtlas-DiffExprs-BatchCorrect/logCPM_batchRemoved_20191107_visualPerp.csv", header = T, row.names = 1)

cat("Exploratory analysis following batch correction \n")


glMDSPlot(logCPM_batchRemoved, labels=meta$PromoterLibrary,
          groups=meta, launch=T, folder = "WITHOUT_Batch/glimmaPlots_normlogCPM_batchRemoved")

##################################
# PCA
logCPM_batchRemoved <- read.csv("~/Postdoc/Projects/Plasticity/Atlas/Cell-type lists/limma/SlAtlas-DiffExprs-BatchCorrect/logCPM_batchRemoved_20191107_visualPerp.csv", header = T, row.names = 1)
meta <- read.csv("~/Postdoc/Projects/Plasticity/Atlas/Cell-type lists/limma/meta/20180814_atlas-SL-meta_final.csv", header = T)

transposed_logCPM_batchRemoved <- t(logCPM_batchRemoved)
# perform prcomp on the transposed matrix from which columns (genes) of zero variance have been removed
pca_proc_batchRemoved <- prcomp(transposed_logCPM_batchRemoved[,apply(transposed_logCPM_batchRemoved, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE) #,retX=TRUE

cat("summary of pca_proc_withbatchcorr \n")
summary(pca_proc_batchRemoved)
#Importance of components:
#                         PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10     PC11
# Standard deviation     73.5893 36.01638 35.63207 33.6498 29.20187 28.43911 27.7112 26.3310 25.77736 25.656 24.54246
# Proportion of Variance  0.2468  0.05912  0.05786  0.0516  0.03886  0.03686  0.0350  0.0316  0.03028  0.030  0.02745
# Cumulative Proportion   0.2468  0.30592  0.36379  0.4154  0.45426  0.49112  0.5261  0.5577  0.58799  0.618  0.64544

plotData_batchRemoved = meta[,c("Name","Replicate","Promoter","PromoterLibrary")]
plotData_batchRemoved$PC1 <- pca_proc_batchRemoved$x[,1]
plotData_batchRemoved$PC2 <- pca_proc_batchRemoved$x[,2]
plotData_batchRemoved$PC3 <- pca_proc_batchRemoved$x[,3]
plotData_batchRemoved$PC4 <- pca_proc_batchRemoved$x[,4]
plotData_batchRemoved$PC5 <- pca_proc_batchRemoved$x[,5]

# 2D plots of pairs of principal components
library(ggplot2)
library(ellipse)
#https://stackoverflow.com/questions/24268843/ggplot2-stat-ellipse-draw-ellipses-around-multiple-groups-of-points

#Draw ellipses around multiple groups of points
centroids <- aggregate(cbind(PC1,PC2)~Promoter,plotData_batchRemoved,median)

conf.rgn  <- do.call(rbind,lapply(unique(plotData_batchRemoved$Promoter),function(t)
  data.frame(Promoter=as.character(t),
             ellipse(cov(plotData_batchRemoved[plotData_batchRemoved$Promoter==t,5:6]),
                     centre=as.matrix(centroids[t,2:3]),
                     level=0.95),
             stringsAsFactors=FALSE)))

ctCols <- c("COR"="#2E5EAC","EP"="#016E82","EN"="#44C3D0","EXO"="#7D2C91","MCO"="#CD85B9","MZ"="#AA1D3F",
            "PH"="#ABD379","V"="#19B35A","WOX"="#F47752","XY"="#EDE83B","X35S"="#939598","ACT"="#414141")            

pdf("WITHOUT_Batch/PCA of Sl atlas norm log2CPM BatchRemoved w medianEllipse.pdf")
qplot(PC1,PC2,data=plotData_batchRemoved,color=Promoter) + #,color=Replicate,
  labs(x="PC1 (25% variability)",y="PC2 (6% variability)") + #CHANGE!!!!
  #ggtitle("PCA WITH batch effect correction") +
  scale_color_manual(values=ctCols)+
  #scale_shape_manual(values=c("EN"=0,"MCO"=1,"MZ"=2,"V"=3,"XY"=4,"EP"=5,"PH"=6,"EXO"=7,"ACT"=8,"X35S"=9,"WOX"=10,"COR"=11)) +
  geom_path(data=conf.rgn)+ 
  geom_point(size = 3)+
  theme_classic()
dev.off()


pdf("WITHOUT_Batch/PCA of Sl atlas norm log2CPM BatchRemoved colorCode.pdf")
qplot(PC1,PC2,data=plotData_batchRemoved,color=Promoter,shape=Replicate) + #,color=Replicate,
  labs(x="PC1 (25% variability)",y="PC2 (6% variability)") + #CHANGE!!!!
  #ggtitle("PCA WITH batch effect correction") +
  scale_shape_manual(values=c("G_S"=0,"H"=1,"I"=2,"J"=3,"K"=4,"L"=5,"M"=6,"N"=7,"O"=8,"P"=9,"Q"=10,"R"=11,"S"=12,"T"=13)) +
  scale_color_manual(values=ctCols)+
  geom_point(size = 4)+
  theme_classic()
dev.off()



################################################################################################
################################################################################################

#### Start PDF
tmpSave <- paste(imgDir,"DEG_Analysis_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")

### Use voom on the dge object.
cat("Define the model metrix with replicates as batch effect \n")
design2 <- model.matrix(~0+Groups+Batch)

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("Groups", "Batch","experimentFactors","\\$","\\:","\\-",
                   colnames(experimentFactors)),sep="",collapse = "|")

colnames(design2) <- gsub(fixCols,"",colnames(design2))
head(design2)

v <- voom(y, design2, plot = TRUE,normalize.method ="quantile")

#cat("Analyzing",nrow(v),"with",ncol(v),"libraries \n")

######## Visualization and quality control

##################
## Correlation between replicates of samples belonging to same group
##################
colnames(v$E) <- meta$PromoterLibrary
corrSamples <- cor(v$E)
## --
tmpSave <- paste(imgDir,"CorrelationBetweenReplicates_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")#width = 8,height = 6)
#colors <- colorRampPalette(c("darkgoldenrod4","darkgoldenrod1","white","white","steelblue1","steelblue4"))
#pdf("BATCH_CORRECTED_CorrelationBetweenReplicates.pdf")
for (each in (levels(Groups))){
  hmData <- corrSamples[grep(each, rownames(corrSamples)),grep(each, colnames(corrSamples))]
  #hmData <- corrSamples[,grep(each, colnames(corrSamples))]
  hm <- T
  if(!hm){
    cat("Heatmaps with NMF \n")
    NMF::aheatmap(hmData,col=viridis(20),
                  txt = ifelse(hmData<0.8,"<",NA),#Rowv = F,Colv = F,
                  main=paste0("Correlation between samples of group ",each))
    
  } else {
    cat("Heatmaps with heatmap.2 \n")
    heatmap.2(hmData,col=viridis(20), keysize = 0.75,
              cellnote = ifelse(hmData<0.8,"*",NA), notecol = "black",
              
              #margins = c(16,16),
              
              dendrogram = "none", trace = "none",density.info='none',
              cexCol  = 0.8 ,cexRow = 0.8,
              lmat=rbind(c(4, 3, 9),
                         c(2, 1, 6),
                         c(8, 5, 7)),
              lhei=c(0.3, 0.6,0.8),
              lwid=c(0.25, 0.4,0.2),
              main=paste0("Correlation\n",each))
    legend("bottomleft",legend = "* means correlation < 0.8",bty = "n")
  }
  
}
dev.off()
################## --
##################

## Assign colors to each of the experimental factors. 
ColorTable <- assignColorsToFactors(experimentFactors)

## Boxplot of normalized counts ordered by Groups
pdf("SlAtlas-DiffExprs-BatchSUBTRACTED/images/boxplots of quantileNorm Counts.pdf")
boxplot(v$E[,order(Groups)], range=0,col=customColors[Groups[order(Groups)]], 
        ylab="log2[counts]", xlab="sample", main="Quantile normalized Counts",
        cex.axis=0.5,las=2)
dev.off()

#####################
### Do contrasts
#####################
cat("Do contrasts- with design as levels")#If modeling the batch effect use design2!!!!!!!!!!
Groups
cont.matrix= makeContrasts(
    "EXO-EP"=EXO-EP,
    "EXO-COR"=EXO-COR,
    "EXO-MCO"=EXO-MCO,
    "EXO-EN"=EXO-EN,
    "EXO-V"=EXO-V,
    "EXO-XY"=EXO-XY,
    "EXO-PH"=EXO-PH,
    "EXO-WOX"=EXO-WOX,
    "EXO-MZ"=EXO-MZ,
    "EXO-ACT"=EXO-ACT,
    "EXO-35S"=EXO-X35S,
    
    "EP-EXO"=EP-EXO,  
    "EP-COR"=EP-COR,
    "EP-MCO"=EP-MCO,
    "EP-EN"=EP-EN,
    "EP-V"=EP-V,
    "EP-XY"=EP-XY,
    "EP-PH"=EP-PH,
    "EP-WOX"=EP-WOX,
    "EP-MZ"=EP-MZ,
    "EP-ACT"=EP-ACT,
    "EP-35S"=EP-X35S,
    
    "COR-EP"=COR-EP,
    "COR-EXO"=COR-EXO,  
    "COR-MCO"=COR-MCO,
    "COR-EN"=COR-EN,
    "COR-V"=COR-V,
    "COR-XY"=COR-XY,
    "COR-PH"=COR-PH,
    "COR-WOX"=COR-WOX,
    "COR-MZ"=COR-MZ,
    "COR-ACT"=COR-ACT,
    "COR-35S"=COR-X35S,
    
    "MCO-EP"=MCO-EP,
    "MCO-EXO"=MCO-EXO,  
    "MCO-COR"=MCO-COR,
    "MCO-EN"=MCO-EN,
    "MCO-V"=MCO-V,
    "MCO-XY"=MCO-XY,
    "MCO-PH"=MCO-PH,
    "MCO-WOX"=MCO-WOX,
    "MCO-MZ"=MCO-MZ,
    "MCO-ACT"=MCO-ACT,
    "MCO-35S"=MCO-X35S,
    
    "EN-EP"=EN-EP,
    "EN-EXO"=EN-EXO,  
    "EN-COR"=EN-COR,
    "EN-MCO"=EN-MCO,
    "EN-V"=EN-V,
    "EN-XY"=EN-XY,
    "EN-PH"=EN-PH,
    "EN-WOX"=EN-WOX,
    "EN-MZ"=EN-MZ,
    "EN-ACT"=EN-ACT,
    "EN-35S"=EN-X35S,
    
    "V-EP"=V-EP,
    "V-EXO"=V-EXO,  
    "V-COR"=V-COR,
    "V-MCO"=V-MCO,
    "V-EN"=V-EN,
    "V-XY"=V-XY,
    "V-PH"=V-PH,
    "V-WOX"=V-WOX,
    "V-MZ"=V-MZ,
    "V-ACT"=V-ACT,
    "V-35S"=V-X35S,
    
    "XY-EP"=XY-EP,
    "XY-EXO"=XY-EXO,  
    "XY-COR"=XY-COR,
    "XY-MCO"=XY-MCO,
    "XY-EN"=XY-EN,
    "XY-V"=XY-V,
    "XY-PH"=XY-PH,
    "XY-WOX"=XY-WOX,
    "XY-MZ"=XY-MZ,
    "XY-ACT"=XY-ACT,
    "XY-35S"=XY-X35S,

    "PH-EP"=PH-EP,
    "PH-EXO"=PH-EXO,  
    "PH-COR"=PH-COR,
    "PH-MCO"=PH-MCO,
    "PH-EN"=PH-EN,
    "PH-V"=PH-V,
    "PH-XY"=PH-XY,
    "PH-WOX"=PH-WOX,
    "PH-MZ"=PH-MZ,
    "PH-ACT"=PH-ACT,
    "PH-35S"=PH-X35S,
    
    "WOX-EP"=WOX-EP,
    "WOX-EXO"=WOX-EXO,  
    "WOX-COR"=WOX-COR,
    "WOX-MCO"=WOX-MCO,
    "WOX-EN"=WOX-EN,
    "WOX-V"=WOX-V,
    "WOX-XY"=WOX-XY,
    "WOX-PH"=WOX-PH,
    "WOX-MZ"=WOX-MZ,
    "WOX-ACT"=WOX-ACT,
    "WOX-35S"=WOX-X35S,
    
    "MZ-EP"=MZ-EP,
    "MZ-EXO"=MZ-EXO,  
    "MZ-COR"=MZ-COR,
    "MZ-MCO"=MZ-MCO,
    "MZ-EN"=MZ-EN,
    "MZ-V"=MZ-V,
    "MZ-XY"=MZ-XY,
    "MZ-PH"=MZ-PH,
    "MZ-WOX"=MZ-WOX,
    "MZ-ACT"=MZ-ACT,
    "MZ-35S"=MZ-X35S,
    
    "ACT-EP"=ACT-EP,
    "ACT-EXO"=ACT-EXO,  
    "ACT-COR"=ACT-COR,
    "ACT-MCO"=ACT-MCO,
    "ACT-EN"=ACT-EN,
    "ACT-V"=ACT-V,
    "ACT-XY"=ACT-XY,
    "ACT-PH"=ACT-PH,
    "ACT-WOX"=ACT-WOX,
    "ACT-MZ"=ACT-MZ,
    "ACT-35S"=ACT-X35S,
    
    "X35S-EP"=X35S-EP,
    "X35S-EXO"=X35S-EXO,  
    "X35S-COR"=X35S-COR,
    "X35S-MCO"=X35S-MCO,
    "X35S-EN"=X35S-EN,
    "X35S-V"=X35S-V,
    "X35S-XY"=X35S-XY,
    "X35S-PH"=X35S-PH,
    "X35S-WOX"=X35S-WOX,
    "X35S-MZ"=X35S-MZ,
    "X35S-ACT"=X35S-ACT,
    levels=design) #design2
#####################

#### Fit and do Diff Expression
cat("Fit and do Diff Expression with batch included in the model metrix (i.e. design2) \n")
fit <- lmFit(v, design2) #design
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

## -- Summary and Venn diagrams , only good for up to 5 comparisons.
results <- decideTests(fit2)
summary(results)
if (ncol(results) <= 5){
  cat ("Doing Venn Diagrams \n")
  vennDiagram(results,include = c("up","down"), main="DE")
} else {
  cat ("More than 5 comparisons, skipping Venn Diagrams  \n")
}
DESummary <- t(summary(decideTests(fit2)))[,-2]
colnames(DESummary) = c("Downregulated","Upregulated")

# Save as csv
tmpSave <- paste(outDir,"DESummary_pairwise_BatchRemoved_",shortName,".csv",sep="")
write.csv(x=DESummary,tmpSave,quote = F,row.names = T)
write.csv(x=DESummary,"DESummary_pairwise_BATCH_SUBTRACTED.csv",quote = F,row.names = T)

# Write to PDF
plotData <- t(DESummary)
yMax <- max(colSums(plotData))
rownames(plotData) <- c("Down","Up")

pdf("/images/Barplot DEG by contrast BatchSubtracted.pdf")
barplot(plotData,legend.text = rownames(plotData),col=c("orange","steelblue4"),
        xlab = "Contrast", ylab = "Number of genes",
        beside = T,
        ylim = c(0,yMax*1.2), #Comment out? #ylim = c(0,yMax*0.8)
        las=2,
        cex.names = 0.6, border = T, bty="n",
        main="DE genes per contrast Batch Removed")


## Prepare for gene annotation
annotationAvail <- F
if (annotationAvail){
  cat("Reading annotation file \n")
  genealiasfile <- "gene_aliases.txt"
  ID2Symbol <- getGeneSymbols(genealiasfile)
} else{cat("Annotation file unavailable \n")}


## --

DEList <- list()
for (contrast in colnames(cont.matrix)){
  print(contrast)
  ## Sorting by none ensures all contrasts will be in the same order
  tmp <- topTable(fit2, coef=contrast,number = Inf,sort.by = "none")
  #
  pValpassed <- table(tmp$adj.P.Val < 0.05)[2]
  cat ("Number of genes with pVal < 0.05 on ",contrast,":",pValpassed,"\n")
  
  
  ## Write genes that are up or downregulated (logFC > 0; logFC < 0)
  upGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC > 0,]))
  tmpSave <- paste(geneListsDir,"/",contrast,"_up",".csv",sep="")
  write.csv(x=upGenes,tmpSave,quote = F,row.names = T)
  #
  downGenes <- as.data.frame(rownames(tmp[tmp$adj.P.Val < 0.05 & tmp$logFC < 0,]))
  tmpSave <- paste(geneListsDir,"/",contrast,"_down",".csv",sep="")
  write.csv(x=downGenes,tmpSave,quote = F,row.names = T)
  #####
  
  #-- Add gene symbols if available
  tmp[,"Symbol"] <- rownames(tmp)
  if (annotationAvail){
    cat("Adding annotation \n")
    Genes <- rownames(tmp)
    idx <- intersect(names(AGI2Symbol),Genes)
    tmp[idx,"Symbol"] <- AGI2Symbol[idx]
    Genes
  }
  #--
  
  ## Add contrast name to the column names, in case of multiple contrasts.
  colnames(tmp) <- paste(colnames(tmp),contrast,sep = ".")
  
  # Write each contrast to file
  tmpSave <- paste(outDir,contrast,"_",shortName,".csv",sep="")
  write.csv(x=tmp,tmpSave,quote = F,row.names = T)

  # Save result to list
  DEList[[contrast]] <- tmp 
}

tmpSave <- paste(outDir,"DEListBatchSUBTRACTED_",shortName,".RData",sep="")
save(DEList,file = tmpSave)
### ------------


tmpSave <- paste(imgDir,"VolcanoPlots_",shortName,".pdf",sep="")
pdf(tmpSave,paper = "USr")
makeVolcanoPlots(DEList,pValCut=0.01,logCut=2,plotGenes=F) #plotGenes=T to print genes in the plot
dev.off()

### Condense into a single list
DE_All <- condenseListTables(DEList) ## Use a custom function
DE_All <- DE_All[,-grep("t.|B.|P.Value|AveExpr",colnames(DE_All))] #Remove unwanted columns

tmpSave <- paste(outDir,"DEG_All_BatchSUBTRACTED_",shortName,".csv",sep="")
write.csv(x = DE_All,file = tmpSave,quote = F,row.names = T)

## Close main img
dev.off()
sink()

