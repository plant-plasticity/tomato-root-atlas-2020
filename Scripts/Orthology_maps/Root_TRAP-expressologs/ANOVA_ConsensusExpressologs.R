##################################################################################################################
# 
#  This script includes the following steps:
# - Expression filtering and transformation/batch correction for tomato and rice
# - PCA for each species
# - Merge log2 expression data based on the 1.6K positive consensus expressologs
# - Quantile normalization, batch effect correction and PCA
# - ANOVA test to identify genes with a similar expression pattern in one cell type across the three species
# 
##################################################################################################################

setwd("~/mywd/")

library(edgeR) #which loads limma as a dependency
library(reshape2)
library(gplots)
library(RColorBrewer)
library(calibrate)
library(Glimma) 
library(viridis)
library(pheatmap)
library(ggplot2)
library(tidyverse)

##################################################################################################################
#STEP 1: Analyze each species separately. Use TPM/intensity data
#(references for using logs(TPM+x) followed by RemoveBatchEffect: https://www.proteinatlas.org/about/assays+annotation, https://static1.squarespace.com/static/5d93d88fcd5ddb320534ef42/t/5e2b4998afedbe47044dabc2/1579895206870/Trevino+Science+2020.pdf)
##################################################################################################################
#TOMATO
##############
## Metadata file
SlmetaFile <- read.csv("/Data/20180814_atlas-SL-meta_final.csv", header=T) #Name of metadata file in csv format.

## Counts file name (with extension) 
SlTPMfile <- read.csv("/Data/20180811_ATLAS_ITAG3.2_Kallisto_quantile_norm_tpm.csv", header = T,row.names = 1)
dim(SlTPMfile)

SlTPMfileDrop <- SlTPMfile[,-20] #remove sapmle M02 (not present in the meta file)

## Filter genes with low expression using TPM counts.
filterByTPM <- T 
TPMcutoff <- 2 

# Load functions
source("/Data/functions.R")

## Start of analysis
####################################################################################
####################################################################################
## Check that samples in both counts and metadata are the same.
## Use function filterCounts(counts,meta)
tmp <- filterCounts(SlTPMfileDrop,SlmetaFile)
SlGeneTPMs <- tmp[["counts"]]
Slmeta <- tmp[["meta"]]
rm(tmp)
## --

###### Design matrix
## Convert experimental metadata to factors for the design
SlexperimentFactors <- lapply(apply(Slmeta,2,split,""),unlist)
SlexperimentFactors <- as.data.frame(lapply(SlexperimentFactors,as.factor))

#Create the design with these factors
print(head(SlexperimentFactors))

###  User modified:
####Simplest design taking into account all possible interactions
SlGroups <- as.factor(paste0(SlexperimentFactors$Promoter))

#accounting for batch effect
SlBatch <- as.factor(paste0(SlexperimentFactors$Replicate))
print(head(SlBatch))

Sldesign <- model.matrix(~0+SlGroups) 

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("SlGroups","SlexperimentFactors","\\$","\\:","\\-",
                   colnames(SlexperimentFactors)),sep="",collapse = "|")

colnames(Sldesign) <- gsub(fixCols,"",colnames(Sldesign))
head(Sldesign)

####################################################################################
#Removing genes with 0 counts on all conditions
cat("Initial number of genes:",nrow(SlTPMfileDrop),"\n") #Initial number of genes: 35768
SlrmIDX <- which(rowSums(SlTPMfileDrop) == 0)
cat("Removing",length(SlrmIDX),"genes \n") #Removing 8321 genes
SlTPMfileDrop <- SlTPMfileDrop[-SlrmIDX,]
cat("Remaining number of genes:",nrow(SlTPMfileDrop),"\n") #Remaining number of genes: 27447 

# Filter genes with low TPMs accross replicates 
cat("Replicates of samples range between:", range(table(SlGroups)),"\n") #Replicates of samples range between: 4 4 

#
if (filterByTPM){
  
  SlsampleMin <- min(table(SlGroups))
  cat("Filtering reads with low TPMs ( <",TPMcutoff,") in at least",SlsampleMin,"replicates \n")
  #
  #cpm <- cpm(dge)
  Slkeep.exprs <- rowSums(SlTPMfileDrop>TPMcutoff)>=SlsampleMin
  table(Slkeep.exprs)
  
  
  cat("Removing",table(Slkeep.exprs)[1],"genes \n")
  cat("Remaining number of genes:",table(Slkeep.exprs)[2],"\n")
  
  #
  
  Sly <- SlTPMfileDrop[Slkeep.exprs, ,]
} else {
  cat("Not doing TPM filtering")
}

# Filtering reads with low TPMs ( < 2 ) in at least 4 replicates 
# Removing 6630 genes 
# Remaining number of genes: 20817 

colnames(Sly) <- Slmeta$PromoterLibrary

##########################################################################################################
# log2 transformation of TPM+3 & batch effect correction
##########################################################################################################
#log2 transformation of TPMs
SllogTPM <- log2(Sly+3) # for CPM: log=TRUE, prior.count=3 based on: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf & https://support.bioconductor.org/p/76837/
#we use a larger prior count to down-weight the contribution of low-abundance genes by squeezing their relative abundances
#toward the mean for each gene. This allows algorithms to focus more on the higher-abundance genes,
#whose measurements have better precision, without explicitly weighting the values.

# Remove batch effect. 
SllogTPM_batchRemoved <- removeBatchEffect(SllogTPM, batch=SlBatch, design=Sldesign) 

removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}

row.names(SllogTPM_batchRemoved) <- removeDotGene(row.names(SllogTPM_batchRemoved))
dim(SllogTPM_batchRemoved)

#create a table with MEDIAN per promoter
median_SllogTPM_batchRemoved <- data.frame(row.names = rownames(SllogTPM_batchRemoved), 
                       ACT=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="ACT")],1,median),
                       COR=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="COR")],1,median),
                       EN=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="EN")],1,median),
                       EP=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="EP")],1,median),
                       EXO=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="EXO")],1,median),
                       MCO=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="MCO")],1,median),
                       MZ=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="MZ")],1,median),
                       PH=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="PH")],1,median),
                       V=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="V")],1,median),
                       WOX=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="WOX")],1,median),
                       X35S=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="X35S")],1,median),
                       XY=apply(SllogTPM_batchRemoved[,which(Slmeta$Promoter=="XY")],1,median)
)


####################################################################
# PCA- AFTER batch effect correction- FOR VISUALIZATION PURPOSE
####################################################################

transposed_logTPM_batchRemoved <- t(SllogTPM_batchRemoved)
# perform prcomp on the transposed matrix from which columns (genes) of zero variance have been removed
Slpca_proc_batchRemoved <- prcomp(transposed_logTPM_batchRemoved[,apply(transposed_logTPM_batchRemoved, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE) #,retX=TRUE

summary(Slpca_proc_batchRemoved)
#Importance of components:
#                         PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10
# Standard deviation     54.0143 38.05987 37.10649 34.82149 32.05864 31.23735 29.44683 28.74517 27.87447 27.23594
# Proportion of Variance  0.1401  0.06959  0.06614  0.05825  0.04937  0.04687  0.04165  0.03969  0.03732  0.03563
# Cumulative Proportion   0.1401  0.20974  0.27588  0.33413  0.38350  0.43037  0.47203  0.51172  0.54904  0.58468

SlplotData_batchRemoved = Slmeta[,c("Name","Replicate","Promoter","PromoterLibrary")]
SlplotData_batchRemoved$PC1 <- Slpca_proc_batchRemoved$x[,1]
SlplotData_batchRemoved$PC2 <- Slpca_proc_batchRemoved$x[,2]
SlplotData_batchRemoved$PC3 <- Slpca_proc_batchRemoved$x[,3]
SlplotData_batchRemoved$PC4 <- Slpca_proc_batchRemoved$x[,4]
SlplotData_batchRemoved$PC5 <- Slpca_proc_batchRemoved$x[,5]

# 2D plots of pairs of principal components

pdf("Sl_TPM/PCA Sl atlas norm log2TPM BatchRemoved.pdf")
qplot(PC1,PC2,data=SlplotData_batchRemoved,color=Replicate,shape=Promoter) +
  labs(x="PC1 (14% variability)",y="PC2 (7% variability)")
  ggtitle("PCA log2TPM WITH batch effect correction") +
  scale_shape_manual(values=c("EN"=0,"MCO"=1,"MZ"=2,"V"=3,"XY"=4,"EP"=5,"PH"=6,"EXO"=7,"ACT"=8,"X35S"=9,"WOX"=10,"COR"=11)) +
  geom_point(size = 3) 
dev.off()

################################################################
################################################################

##############
#RICE
##############
## Metadata file
SlmetaFile <- read.csv("/Data/20180814_atlas-SL-meta_final.csv", header=T) #Name of metadata file in csv format.


## Metadata file
OsmetaFile <- "Os_Elide_metadata.csv" #Name of metadata file in csv format.
doFilter <- T #F
OswhichFilter <- c("Root_total_1", "Root_total_2", "Root_total_3", "QC_1", "QC_2", "QC_3")

#If there are libraries that need to be filtered out (Avoid removing columns manually from the raw counts file)
OswhichFilter


## Filter genes with low expression using TPM counts.
filterByTPM <- T 
TPMcutoff <- 2 

# Load functions
source("/Data/functions.R")

## Start of analysis
####################################################################################
####################################################################################
#Reading metadata file
#Reading counts file
OsmetaFile <- read.csv(paste0("/Data/",OsmetaFile))
Osmeta <- metaDataProcessing(OsmetaFile,doFilter,OswhichFilter)
head(Osmeta)

#Reading counts file
OsTPMfile <- read.csv("/Data/20200311_Os_Atlas_IRGSP-1.0_TRANSCRIPT_Kallisto_quantile_norm_tpm_genes.csv", header = T,row.names = 1)
dim(OsTPMfile)

## Check that samples in both counts and metadata are the same.
## Use function filterCounts(counts,meta)
Ostmp <- filterCounts(OsTPMfile,Osmeta)
OsGeneTPMs <- Ostmp[["counts"]]
Osmeta <- Ostmp[["meta"]]
rm(Ostmp)
## --

###### Design matrix
## Convert experimental metadata to factors for the design
OsexperimentFactors <- lapply(apply(Osmeta,2,split,""),unlist)
OsexperimentFactors <- as.data.frame(lapply(OsexperimentFactors,as.factor))

#Create the design with these factors
print(head(OsexperimentFactors))

###  User modified:
####Simplest design taking into account all possible interactions
OsGroups <- as.factor(paste0(OsexperimentFactors$Promoter))

Osdesign <- model.matrix(~0+OsGroups) 

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("OsGroups","OsexperimentFactors","\\$","\\:","\\-",
                   colnames(OsexperimentFactors)),sep="",collapse = "|")

colnames(Osdesign) <- gsub(fixCols,"",colnames(Osdesign))
head(Osdesign)

####################################################################################
#Removing genes with 0 counts on all conditions
cat("Initial number of genes:",nrow(OsTPMfile),"\n") #Initial number of genes: 37849
OsrmIDX <- which(rowSums(OsTPMfile) == 0)
cat("Removing",length(OsrmIDX),"genes \n") #Removing 3979 genes
OsTPMfile <- OsTPMfile[-OsrmIDX,]
cat("Remaining number of genes:",nrow(OsTPMfile),"\n") #Remaining number of genes: 33870 

# Filter genes with low TPMs accross replicates 
cat("Replicates of samples range between:", range(table(OsGroups)),"\n") #Replicates of samples range between: 3 4 

#
if (filterByTPM){
  
  OssampleMin <- min(table(OsGroups))
  cat("Filtering reads with low TPMs ( <",TPMcutoff,") in at least",OssampleMin,"replicates \n")
  #
  Oskeep.exprs <- rowSums(OsTPMfile>TPMcutoff)>=OssampleMin
  table(Oskeep.exprs)
  
  
  cat("Removing",table(Oskeep.exprs)[1],"genes \n")
  cat("Remaining number of genes:",table(Oskeep.exprs)[2],"\n")
  
  #
  
  Osy <- OsTPMfile[Oskeep.exprs, ,]
} else {
  cat("Not doing TPM filtering")
}

# Filtering reads with low TPMs ( < 2 ) in at least 3 replicates 
# Removing 11432 genes 
# Remaining number of genes: 22438 

colnames(Osy) <- Osmeta$PromoterReplicate

##########################################################################################################
# log2 transformation of TPM+3 & calculate mean/median TPM
##########################################################################################################
#log2 transformation of TPMs
OslogTPM <- log2(Osy+3) # for CPM: log=TRUE, prior.count=3 based on: https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf & https://support.bioconductor.org/p/76837/

#create a table with mean/ MCOmedian per promoter
mtdt4ct<- read.csv("~/Postdoc/Projects/Plasticity/Atlas/Cell-type lists/CSA/Rice/Limma-TRANSCRIPT/meta/Os_atlas_meta.csv", header=T)
meanmedian_OslogTPM <- data.frame(row.names = rownames(OslogTPM), 
                            MCO=apply(OslogTPM[,which(mtdt4ct$Promoter=="MCO")],1,median),
                            EN=apply(OslogTPM[,which(mtdt4ct$Promoter=="EN")],1,mean),
                            V=apply(OslogTPM[,which(mtdt4ct$Promoter=="V")],1,mean),
                            MZ=apply(OslogTPM[,which(mtdt4ct$Promoter=="MZ")],1,mean),
                            X35S=apply(OslogTPM[,which(mtdt4ct$Promoter=="X35S")],1,mean)
                            
)

#change geneID (from Os0xt to Os0xg)

############################################################################################
# PCA
##################################

transposed_OslogTPM <- t(OslogTPM)
# perform prcomp on the transposed matrix from which columns (genes) of zero variance have been removed
Ospca_proc <- prcomp(transposed_OslogTPM[,apply(transposed_OslogTPM, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE) #,retX=TRUE

summary(Ospca_proc)
#Importance of components:
#                         PC1      PC2      PC3      PC4      PC5      PC6      PC7      PC8      PC9     PC10
# Standard deviation     93.4638 50.2531 39.70272 37.12237 34.27704 34.00471 33.43278 31.27435 30.28298 29.79571 29.66036
# Proportion of Variance  0.3893  0.1125  0.07025  0.06142  0.05236  0.05153  0.04982  0.04359  0.04087  0.03957  0.03921
# Cumulative Proportion   0.3893  0.5019  0.57212  0.63353  0.68590  0.73743  0.78725  0.83084  0.87171  0.91127  0.95048

OsplotData = Osmeta[,c("Name","replicate","Promoter","PromoterReplicate")]
OsplotData$PC1 <- Ospca_proc$x[,1]
OsplotData$PC2 <- Ospca_proc$x[,2]
OsplotData$PC3 <- Ospca_proc$x[,3]
OsplotData$PC4 <- Ospca_proc$x[,4]
OsplotData$PC5 <- Ospca_proc$x[,5]

# 2D plots of pairs of principal components

pdf("../OsTPM_Fig/PCA Os atlas norm log2TPM.pdf")
qplot(PC1,PC2,data=OsplotData,color=replicate,shape=Promoter) +
  labs(x="PC1 (39% variability)",y="PC2 (11% variability)") + #CHANGE!!!!
  ggtitle("PCA log2TPM") +
  scale_shape_manual(values=c("EN"=0,"MCO"=1,"MZ"=2,"V"=3,"X35S"=9)) + #,"QC"=10
  geom_point(size = 3) 
dev.off()

################################################################
################################################################

#############################################
#ARABIDOPSIS- load reMeanlog_plastiDrop
#############################################
Atlog <- read.csv("/Data/reMeanlog_plastiDrop.csv", header = T, row.names = 1)
#Add Species name to colnames
colnames(Atlog) <- c("At_X35S","At_V","At_EN","At_MZ","At_MCO")
Atlog<- Atlog[,c(5,3,2,4,1)]
colnames(Atlog) <- c("At_MCO",	"At_EN",	"At_V",	"At_MZ",	"At_X35S")

##########################################################################################################
#STEP 2: Combine log2 expression data based on the 1.6K positive consensus expressologs map 
##########################################################################################################
consenExpCor <- read.csv("/Data/20200315 1640 consensusExpressologs positive CorrValues_TRANSCRIPT.csv", header = T)

consenExp <- consenExpCor[, -c(1,5:7)]

#Subst Sly atlas for the 5 comparable CT
median_SllogTPM_bRemov5Comp <- median_SllogTPM_batchRemoved[,c(6,3,9,7,11)]
#Add Species name to colnames
colnames(median_SllogTPM_bRemov5Comp) <- c("Sl_MCO","Sl_EN","Sl_V","Sl_MZ","Sl_X35S")

colnames(meanmedian_OslogTPM) <- c("Os_MCO","Os_EN","Os_V","Os_MZ","Os_X35S")

consenExpSl <- merge(consenExp, median_SllogTPM_bRemov5Comp, by.x = "SlGeneID", by.y = "row.names") #1639 genes
consenExpSlAt <- merge(consenExpSl,Atlog, by.x = "AtGeneID", by.y = "row.names") #1557 genes
consenExpSlAtOs <- merge(consenExpSlAt,meanmedian_OslogTPM, by.x = "OsGeneID", by.y = "row.names") #1535 genes

##########################################################################################################
#STEP 3: CSA exploratory data analysis (before running ANOVA) 
##########################################################################################################
# # Include EN, MCO, MZ, V & 35S for normalizeBetweenArrays & removeBatchEffect
# Remove 35S for EDA
####################################################################################
#Exclude geneIDs and keep only AtGeneID as row names
csa_consensusPos5CT <- consenExpSlAtOs[,-c(1:3)]
rownames(csa_consensusPos5CT) <- consenExpSlAtOs$AtGeneID
#Use the function normalizeBetweenArrays to align expression values between samples using quantile normalization
csa_consensusPos5CTNorm <-normalizeBetweenArrays(csa_consensusPos5CT, method="quantile")

#Metadata: 
###### 
meta <- read.csv("/Data/CSA_Metadata.csv", header = T)
head(meta)
meta4CT <- read.csv("/Data/ANOVAmetadata4CT.csv", header = T)
head(meta4CT)
###### 
###### Define experimental factors and design matrix
###### 
## Convert experimental metadata to factors for the design
Groups <- as.factor(meta$Tissue)

Batch <- as.factor(meta$Species)

design <- model.matrix(~0+Groups)

## Ensures column names are optimal for the contrast design
fixCols <- paste(c("Groups","Batch", "meta","\\$","\\:","\\-",
                   colnames(meta)),sep="",collapse = "|") #"Batch",

colnames(design) <- gsub(fixCols,"",colnames(design))
head(design)
####################################      
# Remove batch effect
####################################    
csa_consensusPos5CTNorm_bRemov <- removeBatchEffect(csa_consensusPos5CTNorm, batch=Batch, design=design) 

############################################################################################
# Remove 35S for EDA
####################################################################################
csa_consensusPos4CTNorm_bRemov <- csa_consensusPos5CTNorm_bRemov[,-c(5,10,15)]

##################################
# PCA
##################################
transposed_csa_consensusPos4CTNorm_bRemov <- t(csa_consensusPos4CTNorm_bRemov)
# perform prcomp on the transposed matrix from which columns (genes) of zero variance have been removed
pca_proc_consensusPos4CTNorm_bRemov <- prcomp(transposed_csa_consensusPos4CTNorm_bRemov[,apply(transposed_csa_consensusPos4CTNorm_bRemov, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE) #,retX=TRUE

summary(pca_proc_consensusPos4CTNorm_bRemov)
#Importance of components:
#                          PC1     PC2     PC3     PC4      PC5      PC6     PC7    PC8      PC9     PC10     PC11      PC12
# Standard deviation     21.5266 16.1990 15.3259 13.1102 10.85947 8.57882 8.02520 6.85678 6.76111 5.76110 4.53794 2.646e-14
# Proportion of Variance  0.3019  0.1709  0.1530  0.1120  0.07683 0.04795 0.04196 0.03063 0.02978 0.02162 0.01342 0.000e+00
# Cumulative Proportion   0.3019  0.4728  0.6259  0.7378  0.81465 0.86260 0.90455 0.93518 0.96496 0.98658 1.00000 1.000e+00

plotData_consensusPos4CTNorm_bRemov = meta4CT[,c("Name","Tissue","Species")]
plotData_consensusPos4CTNorm_bRemov$PC1 <- pca_proc_consensusPos4CTNorm_bRemov$x[,1]
plotData_consensusPos4CTNorm_bRemov$PC2 <- pca_proc_consensusPos4CTNorm_bRemov$x[,2]
plotData_consensusPos4CTNorm_bRemov$PC3 <- pca_proc_consensusPos4CTNorm_bRemov$x[,3]
plotData_consensusPos4CTNorm_bRemov$PC4 <- pca_proc_consensusPos4CTNorm_bRemov$x[,4]
plotData_consensusPos4CTNorm_bRemov$PC5 <- pca_proc_consensusPos4CTNorm_bRemov$x[,5]

# 2D plots of pairs of principal components
pdf("PCA qNormlog2TPM 4CT consensusPOSExp bRemoved.pdf")
qplot(PC1,PC2,data=plotData_consensusPos4CTNorm_bRemov,color=Tissue,shape=Species) +
  labs(x="PC1 (30% variability)",y="PC2 (17% variability)") +
  ggtitle("PCA consensusPOSExp 4CT log2TPM/intensity qnorm AFTER batch correction") +
  scale_color_manual(values=c("EN"="#44C3D0","MCO"="#CD85B9","MZ"="#AA1D3F","V"="#19B35A")) +#,"X35S"="#939598"
  scale_shape_manual(values=c("Sly"=15,"Ath"=16,"Osa"=17)) +
  geom_point(size = 5)+
  theme_classic()
dev.off()



##########################################################################################################
#STEP 4: ANOVA
##########################################################################################################

csaMelt_NormbRemov <- melt(data = as.matrix(csa_consensusPos4CTNorm_bRemov), id.vars = rownames(csa_consensusPos4CTNorm_bRemov), 
                       measure.vars = c("Sl_MCO","Sl_EN","Sl_V","Sl_MZ","At_MCO","At_EN","At_V","At_MZ","Os_MCO","Os_EN","Os_V","Os_MZ"))

csaMelt_NormbRemov$Tissue <-  substr(csaMelt_NormbRemov$Var2, 4, 7)
csaMelt_NormbRemov$Species <-  substr(csaMelt_NormbRemov$Var2, 1, 2)
csaANOVA_NormbRemov <- csaMelt_NormbRemov[,-2]

library(tidyverse)
library(agricolae)

#create a character of unique gene names: nb= q-normed & batch corrected
uniqu_gene <- as.character(unique(csaANOVA_NormbRemov$Var1))
gene_tpm_df <- data.frame(row.names=uniqu_gene) # create a dataframe with the unique gene names as row names

###################################################
# HSD.test: alpha = 0.1
###################################################

# create a list of outputs for every step of the analysis

tukeyResults_a01 <- list()

for(i in 1:length(rownames(gene_tpm_df))){
  geneName <- uniqu_gene[i] #uniqe gene names
  lmResults_nb[[geneName]] <- lm(value ~ Tissue, data=csaANOVA_NormbRemov[csaANOVA_NormbRemov$Var1 == uniqu_gene[i], ])
  anovaResults_nb[[geneName]] <- aov(lmResults_nb[[geneName]])
  
  tukeyResults_a01[[geneName]] <- HSD.test(anovaResults_nb[[geneName]],trt = "Tissue",group=T, alpha = 0.1)
  
  boxplotResults_nb[[geneName]] <- csaANOVA_NormbRemov[csaANOVA_NormbRemov$Var1 == uniqu_gene[i], ] %>% 
    ggplot(aes(x=Tissue,y=value)) + 
    geom_boxplot() + ggtitle(geneName) + 
    geom_jitter(aes(shape=Species),size=3) + 
    theme_classic()
  
}

#anovaResults[[1]]$coefficients
#example for extracting Fvalues:
summary(anovaResults_nb$AT1G22860)[[1]]$`F value`[1]

# Extract the F values for each gene
FvaluesResults_nb <- sapply(anovaResults_nb,function(x){
  summary(x)[[1]]$`F value`[1]
})

#create a data frame with gene names & F values
FvaluesResults_nb <- data.frame("Gene"=names(FvaluesResults_nb), "Fvalue"=FvaluesResults_nb,stringsAsFactors = F)
head(FvaluesResults_nb)
# Gene    Fvalue
# AT3G23900 AT3G23900 26.1117384
# AT1G07230 AT1G07230  0.4465772
# AT5G14030 AT5G14030  5.0838739
# AT1G20050 AT1G20050  4.7721490
# AT1G52310 AT1G52310  4.4559665
# AT5G15020 AT5G15020  0.5818210

#Sort the dataframe beased on the order of the F values (highest value at the top of the table)
FvaluesResults_nb_sorted <- FvaluesResults_nb[order(FvaluesResults_nb$Fvalue,decreasing = T),,drop=F]

#Extract Tukey's groups while keeping the tissues in the same order
groups_Tukey_a01 <- lapply(tukeyResults_a01,function(x){
  #x <- tukeyResults_a01[[1]]
  tmp <- x$groups[,2,drop=F]
  data.frame(tmp[sort(rownames(tmp)),,drop=F],stringsAsFactors = F)
})

#Extract mean expression while keeping the tissues in the same order
means_Tukey_a01 <- lapply(tukeyResults_a01,function(x){
  #x <- tukeyResults_a01[[1]]
  tmp <- x$groups[,1,drop=F]
  data.frame(tmp[sort(rownames(tmp)),,drop=F],stringsAsFactors = F)
})

#Create a dataframe from mean expression
tukeyMeansDF_a01 <- t(do.call("cbind",means_Tukey_a01))
rownames(tukeyMeansDF_a01) <- names(means_Tukey_a01)
head(tukeyMeansDF_a01)

#Create a dataframe from Tukey's groups
tukeyGroupsDF_a01 <- t(do.call("cbind",groups_Tukey_a01))
rownames(tukeyGroupsDF_a01) <- names(groups_Tukey_a01)
head(tukeyGroupsDF_a01)

#Merge the two data frames
tukeyMeanGroup_a01 <- merge(tukeyMeansDF_a01, tukeyGroupsDF_a01, by = "row.names")
tukeyMeanGroupFv_a01 <- merge(tukeyMeanGroup_a01, FvaluesResults_nb, by.x="Row.names", by.y = "Gene")

tukeyMeanGroupFv_a01 <-data.frame(tukeyMeanGroupFv_a01[,-1], row.names = tukeyMeanGroup_a01$Row.names)
colnames(tukeyMeanGroupFv_a01) <- c("exp_EN", "exp_MCO", "exp_MZ", "exp_V","grp_EN", "grp_MCO", "grp_MZ", "grp_V", "Fvalue")
write.csv(tukeyMeanGroupFv_a01, "1535 consExp log2 qnorm-bRemoved meanExp, Tukey groups & Fvalue p01.csv", quote = F, row.names = T)

# #Test specific genes
# boxplotResults_nb[["AT1G22860"]]
# tukeyResults_a01[["AT1G22860"]]
# summary(anovaResults_nb[["AT1G22860"]])

summary(FvaluesResults_nb_sorted$Fvalue)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 0.03552   1.42190   2.64548   4.03275   4.70954 163.99641

#Create a table of genes with top 15% Fvalues
top15_nb <- FvaluesResults_nb_sorted[1:floor(nrow(FvaluesResults_nb_sorted) * .15),]

#Create a table of top15% F values with mean expression and Tukey's group
# to select for conserved cell type (enriched) genes
top15expGroup_a01 <- merge(top15_nb, tukeyMeanGroupFv_a01, by= "row.names")
top15expGroup_a01 <- data.frame(top15expGroup_a01[,-c(1,2,12)], row.names = top15expGroup_a01$Row.names)
write.csv(top15expGroup_a01, "Top 15% Fvalues meanExp Tukey groups p01.csv", quote = F, row.names = T)

# Add MapMan annotations to top15expGroup
mapman <- read.csv("/Data/Ath_AGI_LOCUS_TAIR10_Aug2012_drop_CapLock.csv", header = T)

top15mapman_a01 <- merge(top15expGroup_a01,mapman, by.x = "row.names", by.y = "IDENTIFIER", all.x = T)
write.csv(top15mapman_a01, "Top 15% Fvalues meanExp Tukey groups p01 MapMan annot.csv", quote = F, row.names = F)

###########################
# 1. Identify tomato and rice expressologs for the conserved cell type-specific genes
# 2. Test the overlap of these genes with CEGs (constitutively expressed genes) from each species

###############################
#Heatmap of 32 sig. cell type enriched genes based on the ANOVA test
###############################
ctenrichedHeatmap <- read.csv("Data/conservedEnriched_ANOVA_heatmap.csv", header = T)
ctenrichedHeatmap$GeneAnnotation <-paste(ctenrichedHeatmap$GeneName, ctenrichedHeatmap$Annotation,sep="; ")
subctenrichedHeatmap <-data.frame(ctenrichedHeatmap[,c(3:6)], row.names = ctenrichedHeatmap$GeneAnnotation)
subctTukeygroups <-data.frame(ctenrichedHeatmap[,c(7:10)], row.names = ctenrichedHeatmap$GeneAnnotation)

##########################################
#split columns into 4 groups based on the cell type
row_groups <- ctenrichedHeatmap$Tissue
table(row_groups)
# EN MCO  MZ   V 
#  5   2  25   5 

# Data frame with column annotations.
mat_row <- data.frame(group = row_groups)#, order(ctenrichedHeatmap$Tissue))
rownames(mat_row) <- rownames(subctenrichedHeatmap)

# List with colors for each annotation.
mat_colors <- list(group = c("#CD85B9","#44C3D0","#19B35A","#AA1D3F"))
names(mat_colors$group) <- unique(row_groups)

############################################
#heatmap WITH scaling expression by genes
#unlog first and then use the log2 of tissue vs MEDIAN TPM
############################################
#Unlog
subctenrichedHeatmapunlog <- 2^(subctenrichedHeatmap)

#Use the log2 of tissue vs MEDIAN
ctenrichedunlogCenterrelog <- log2(subctenrichedHeatmapunlog/apply(subctenrichedHeatmapunlog,1,median, na.rm = T))
max(ctenrichedunlogCenterrelog) #1.342216
min(ctenrichedunlogCenterrelog) #-1.150866

MedianCenterrelog_colors <- c(min(ctenrichedunlogCenterrelog),seq(-1.15,1.35,by=0.2),max(ctenrichedunlogCenterrelog)) 
medianCenterrelog_palette <- viridis(length(MedianCenterrelog_colors)-2)

# subctenrichedHeatmapCenter <- subctenrichedHeatmap/apply(subctenrichedHeatmap,1,median)
# 
pdf("37 ANOVA CTenriched consensExp-Center qNorm bRemov.pdf")
pheatmap::pheatmap(
  ctenrichedunlogCenterrelog,
  cellwidth = 40,
  cellheight = 12,
  color = medianCenterrelog_palette,
  cluster_rows = F,
  cluster_cols = F,
  fontsize_row = 7,
  show_rownames     = T,
  # display_numbers = subctTukeygroups,
  # fontsize_number = 4,
  fontsize_col = 9,
  angle_col = 45,
  annotation_row    = mat_row,
  annotation_colors = mat_colors,
  border_color      = "gray60")

dev.off()
