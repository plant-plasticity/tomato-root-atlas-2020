####################################################################################
#  This script is aimed to:
#
# 1. Identify reciprocal expressologs high confidence expressologs that have
#    a reciprocal match between Sl&At, Sl&Os and At&Os
#
# 2. Identify consensus expressologs that  have identical expressolog relationships 
#    independent of the reference species and positive expression correlations
####################################################################################
setwd("~/mywd/") 

####################################################################################
# STEP 1.	For each pairwise expressologs (e.g. Sly-Ath_expressologs.tsv & Sly-Osa_expressologs.tsv) 
#        find ortholog pairs with max correlation for each species separately 
#        (i.e. max orthology based on Sl & max orthlogy based on At)

#*curated = remove dupication and kept only one pairwise expression correlation.
#           remove Corr_value = NA, which results from genes that are expressed in only one CT (~17 in Ath&Osa and >500 in Sly&Osa) ) 
####################################################################################

#############################
# Sly-Ath_expressologs.csv
#############################
slat <- read.csv("Sly-Ath_expressologs_curated.csv", header =T , stringsAsFactors = F)

# ARABIDOPSIS
#create a character of unique At genes from Sl-At expressologs
uniqu_Atgene <- unique(slat$SpeciesA)

#create a dataframe with the unique gene names as row names
Atgene_slat <- data.frame(row.names=uniqu_Atgene) #13214 unique genes
Atgene_slat[,c("SpeciesA","SpeciesB", "Corr_value")]<- NA

for (i in 1:length(rownames(Atgene_slat))){
  tmp <- slat[which(slat$SpeciesA==rownames(Atgene_slat[i,])),]
  Atgene_slat[i,] <- lapply(tmp[which.max(tmp$Corr_value),],as.character)
}


# TOMATO *Sly-Ath_expressologs.csv
#create a character of unique Sl genes from Sl-At expressologs
uniqu_Slgene <- unique(slat$SpeciesB) #13199 unique genes

#create a dataframe with the unique gene names as row names
Slgene_slat <- data.frame(row.names=uniqu_Slgene)
#colnames(data)
Slgene_slat[,c("SpeciesA","SpeciesB", "Corr_value")]<- NA

for (i in 1:length(rownames(Slgene_slat))){
  tmp <- slat[which(slat$SpeciesB==rownames(Slgene_slat[i,])),]
  Slgene_slat[i,] <- lapply(tmp[which.max(tmp$Corr_value),],as.character)
}

#############################
# Sly-Osa_expressologs.csv
#############################
slos <- read.csv("20200314_Sly-Osa_expressologs_curated.csv", header =T , stringsAsFactors = F)

# RICE 
#create a character of unique Os genes from Sl-Os expressologs
uniqu_Osgene <- unique(slos$SpeciesA) # Rice TRANSCRIPT=11,794 unique genes (Rice CDS=12051 unique genes)

#create a dataframe with the unique gene names as row names                     
Osgene_slos <- data.frame(row.names=uniqu_Osgene) 
Osgene_slos[,c("SpeciesA","SpeciesB", "Corr_value")]<- NA

for (i in 1:length(rownames(Osgene_slos))){
  tmp <- slos[which(slos$SpeciesA==rownames(Osgene_slos[i,])),]
  Osgene_slos[i,] <- lapply(tmp[which.max(tmp$Corr_value),],as.character)
}


# TOMATO *Sly-Osa_expressologs.csv
#create a character of unique Sl genes from Sl-At expressologs
uniqu_Slgene_slos <- unique(slos$SpeciesB) #Rice TRANSCRIPT=12,567 unique genes (Rice CDS=12682 unique genes)

#create a dataframe with the unique gene names as row names                     
Slgene_slos <- data.frame(row.names=uniqu_Slgene_slos) 
#colnames(data)
Slgene_slos[,c("SpeciesA","SpeciesB", "Corr_value")]<- NA

for (i in 1:length(rownames(Slgene_slos))){
  tmp <- slos[which(slos$SpeciesB==rownames(Slgene_slos[i,])),]
  Slgene_slos[i,] <- lapply(tmp[which.max(tmp$Corr_value),],as.character)
}


#############################
# Ath-Osa_expressologs.csv
#############################
atos <- read.csv("20200314_Ath-Osa_expressologs_curated.csv", header =T , stringsAsFactors = F)

# RICE 
#create a character of unique Osa genes from At-Os expressologs
uniqu_Osgene_atos <- unique(atos$SpeciesB) #Rice TRANSCRIPT=11,675 unique genes (Rice CDS=11924 unique genes)

#create a dataframe with the unique gene names as row names                     
Osgene_atos <- data.frame(row.names=uniqu_Osgene_atos) 
Osgene_atos[,c("SpeciesA","SpeciesB", "Corr_value")]<- NA

for (i in 1:length(rownames(Osgene_atos))){
  tmp <- atos[which(atos$SpeciesB==rownames(Osgene_atos[i,])),]
  Osgene_atos[i,] <- lapply(tmp[which.max(tmp$Corr_value),],as.character)
}


# ARABIDOPSIS *Ath-Osa_expressologs.csv
#create a character of unique Ath genes from Ath-Osa expressologs
uniqu_Atgene_atos <- unique(atos$SpeciesA) #Rice TRANSCRIPT=12,441 unique genes (Rice CDS=12578 unique genes)

#create a dataframe with the unique gene names as row names                     
Atgene_atos <- data.frame(row.names=uniqu_Atgene_atos) 
Atgene_atos[,c("SpeciesA","SpeciesB", "Corr_value")]<- NA

for (i in 1:length(rownames(Atgene_atos))){
  tmp <- atos[which(atos$SpeciesA==rownames(Atgene_atos[i,])),]
  Atgene_atos[i,] <- lapply(tmp[which.max(tmp$Corr_value),],as.character)
}

####################################################################################################
# STEP 2.	Identify identical expressolog pairs- expessologs that have a reciprocal max correlation 
####################################################################################################
#between Sl & At
#####################
slat_reciprocalExpr <- Slgene_slat[Slgene_slat$SpeciesA %in% Atgene_slat$SpeciesA & Slgene_slat$SpeciesB %in% Atgene_slat$SpeciesB,]
#9436 rows
#How many unique genes from each species?
sapply(slat_reciprocalExpr, function(x) length(unique(x)))
# SpeciesA(At)   SpeciesB(Sl)
# 8899            9436

atsl_reciprocalExpr <- Atgene_slat[Atgene_slat$SpeciesA %in% Slgene_slat$SpeciesA & Atgene_slat$SpeciesB %in% Slgene_slat$SpeciesB,]
#9497 rows
sapply(atsl_reciprocalExpr, function(x) length(unique(x)))
# SpeciesA(At)   SpeciesB(Sl)
# 9497            8880

#####################
#between Sl & Os
#####################
slos_reciprocalExpr <- Slgene_slos[Slgene_slos$SpeciesA %in% Osgene_slos$SpeciesA & Slgene_slos$SpeciesB %in% Osgene_slos$SpeciesB,]
#Rice TRANSCRIPT=8,571 unique genes (Rice CDS=8729 rows)
#How many unique genes from each species?
sapply(slos_reciprocalExpr, function(x) length(unique(x)))
# SpeciesA(Os)   SpeciesB(Sl) 
# 8168            8571
ossl_reciprocalExpr <- Osgene_slos[Osgene_slos$SpeciesA %in% Slgene_slos$SpeciesA & Osgene_slos$SpeciesB %in% Slgene_slos$SpeciesB,]
#Rice TRANSCRIPT=8,753 unique genes (Rice CDS=8899 rows)
sapply(ossl_reciprocalExpr, function(x) length(unique(x)))
# SpeciesA(Os)   SpeciesB(Sl)
# 8753            8144

#####################
#between At & Os
#####################
atos_reciprocalExpr <- Atgene_atos[Atgene_atos$SpeciesA %in% Osgene_atos$SpeciesA & Atgene_atos$SpeciesB %in% Osgene_atos$SpeciesB,]
#Rice TRANSCRIPT=8,320 unique genes (Rice CDS=8501 rows)

#How many unique genes from each species?
sapply(atos_reciprocalExpr, function(x) length(unique(x)))
# SpeciesA(At)   SpeciesB(Os) 
# 8320            7963
osat_reciprocalExpr <- Osgene_atos[Osgene_atos$SpeciesA %in% Atgene_atos$SpeciesA & Osgene_atos$SpeciesB %in% Atgene_atos$SpeciesB,]
#Rice TRANSCRIPT=8,554 unique genes (Rice CDS=8708 rows)

sapply(osat_reciprocalExpr, function(x) length(unique(x)))
# SpeciesA(At)   SpeciesB(Os)
# 7975            8554 

####################################################################################################
# STEP 3.	
#Find the union set between Sly-Ath reciprocal expressologs (9436 uniqe genes) 
#and Sly-Osa reciprocal expressologs (8729 uniqe genes)
####################################################################################################
#TOMATO AS A REFERENCE
######
unionSlreciprocalExpr <- merge(slat_reciprocalExpr,slos_reciprocalExpr, by = "SpeciesB") ##Rice TRANSCRIPT= 6293 Sly unique genes (Rice CDS=6403 Sl unique genes)
summary(duplicated(unionSlreciprocalExpr$SpeciesB)) #Sl
#   Mode   FALSE
# logical    6293
summary(duplicated(unionSlreciprocalExpr$SpeciesA.x))#At
#   Mode   FALSE    TRUE
# logical    6173     120
summary(duplicated(unionSlreciprocalExpr$SpeciesA.y))#LOC_Os
#   Mode   FALSE    TRUE
# logical    6155     138

colnames(unionSlreciprocalExpr) <- c("SlGeneID","AtGeneID","SlAt_corrValue","LOC_OsGeneID","SlOs_corrValue")

##########################
#Find the union set between Sly-Ath reciprocal expressologs (atsl_reciprocalExpr:9497 uniqe genes)
#and Ath-Osa reciprocal expressologs (8501 uniqe genes)
#ARABIDOPSIS AS A REFERENCE
#########
unionAtreciprocalExpr <- merge(atsl_reciprocalExpr,atos_reciprocalExpr, by = "SpeciesA") #Rice TRANSCRIPT= 6470 Ath unique genes (Rice CDS=6618 Ath unique genes)
summary(duplicated(unionAtreciprocalExpr$SpeciesA)) #Ath
#   Mode   FALSE 
# logical    6470 
summary(duplicated(unionAtreciprocalExpr$SpeciesB.x))#Sly  
#   Mode   FALSE    TRUE 
# logical    6329     141 
summary(duplicated(unionAtreciprocalExpr$SpeciesB.y))#LOC_Os
#   Mode   FALSE    TRUE 
# logical    6340     130 

colnames(unionAtreciprocalExpr) <- c("AtGeneID","SlGeneID","AtSl_corrValue","LOC_OsGeneID","AtOs_corrValue")


##############################
#Find the union set between Osa-Sly reciprocal expressologs (8899 uniqe genes)
#and Osa-Ath reciprocal expressologs (8708 uniqe genes)
#RICE AS A REFERENCE
###############
unionOsreciprocalExpr <- merge(ossl_reciprocalExpr,osat_reciprocalExpr, by.x = "SpeciesA", by.y = "SpeciesB") #Rice CDS=6648 Osa unique genes
summary(duplicated(unionOsreciprocalExpr$SpeciesA)) #LOC_Osa 
#   Mode   FALSE 
# logical    6516 
summary(duplicated(unionOsreciprocalExpr$SpeciesB))#Sly
#   Mode   FALSE    TRUE 
# logical    6342     174 
summary(duplicated(unionOsreciprocalExpr$SpeciesA.y))#Ath
#   Mode   FALSE    TRUE 
# logical    6347     169 
colnames(unionOsreciprocalExpr) <- c("LOC_OsGeneID","SlGeneID","OsSl_corrValue","AtGeneID","OsAt_corrValue")

##########################################################################################################
#STEP 4. Compile a list of consensus expressologs with positive corr between all pairwise comparisons
# or- how many orthologs have the same match when using Sly, Ath or Osa as reference species?
##########################################################################################################
library(sqldf)

#subset tables for columns with GeneIDs
SubunionSlreciprocalExpr <- unionSlreciprocalExpr[,c(1,2,4)]

SubunionAtreciprocalExpr <- unionAtreciprocalExpr[,c(2,1,4)]

SubunionOsreciprocalExpr <- unionOsreciprocalExpr[,c(2,4,1)]

Sl_At_reciprocalExpr <- sqldf('SELECT * FROM SubunionSlreciprocalExpr INTERSECT SELECT * FROM SubunionAtreciprocalExpr') #4338 Rice CDS=4427 orthologs

Sl_At_OsreciprocalExpr <- sqldf('SELECT * FROM Sl_At_reciprocalExpr INTERSECT SELECT * FROM SubunionOsreciprocalExpr') #4250 Rice CDS=4321 orthologs

#Add RAP names for rice
#Convert RAP to MSU's LOC_Os IDs using the conversion table in https://rapdb.dna.affrc.go.jp/download/irgsp1.html
Osconvert <- read.csv("Data/RAP-MSU_2019-08-29.txt", header = F, sep = "\t",stringsAsFactors = F)

#Remove dots from LOC_Os ID to match gene Id in unionSlreciprocalExpr
removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}

Osconvert$V2 <- removeDotGene(Osconvert$V2)

Sl_At_OsreciprocalExprRAP <- merge(Sl_At_OsreciprocalExpr,Osconvert, by.x = "LOC_OsGeneID", by.y = "V2", all.x = T)
#EXPLORE THE DATA (4296)
summary(duplicated(Sl_At_OsreciprocalExprRAP$SlGeneID)) #Sl, AtGeneID or LOC_OsGeneID
#    Mode   FALSE    TRUE
# logical    4250      46
summary(duplicated(Sl_At_OsreciprocalExprRAP$V1)) #Os
#    Mode   FALSE    TRUE 
#logical    4295       1 

colnames(Sl_At_OsreciprocalExprRAP)[colnames(Sl_At_OsreciprocalExprRAP)=="V1"] <- "OsGeneID"
write.csv(Sl_At_OsreciprocalExprRAP, "4296 duplicated consensus expressologs_TRANSCRIPT.csv", quote = F, row.names = F)

#Converting to RAP IDs resulted in duplication of 46 LOC_Os genes and and 2 Loc_Os without Os ID. remove these genes!



################################################################################################
# Merge corrlation values and filter out negative values.
# Keep only consensus expressologs with positive corr between all pairwise comparisons
################################################################################################
consusExp_Slrecip <- merge(deduprecipExprRAP, unionSlreciprocalExpr, by = "SlGeneID") # 4208
consusExp_Slrecipdrop <- consusExp_Slrecip[,-c(5,7)]
consusExp_SlAtrecip <- merge(consusExp_Slrecipdrop, unionAtreciprocalExpr, by.x = "AtGeneID.x", by.y = "AtGeneID") #4208
consusExp_SlAtrecipdrop <- consusExp_SlAtrecip[, -c(7,9)]
consusExp_SlAtOsrecip <- merge(consusExp_SlAtrecipdrop, unionOsreciprocalExpr, by.x = "LOC_OsGeneID.x", by.y = "LOC_OsGeneID")#4208
consusExp_SlAtOsrecipdrop <- consusExp_SlAtOsrecip[, -c(7,9:12)]

#Filter duplicated rows and keep only rows with positive corr values

library(dplyr)
consusExp_SlAtOsPos <- consusExp_SlAtOsrecipdrop %>% 
  filter(SlAt_corrValue > 0, SlOs_corrValue > 0, AtOs_corrValue > 0) #1640

colnames(consusExp_SlAtOsPos)[colnames(consusExp_SlAtOsPos)=="LOC_OsGeneID.x"] <- "LOC_OsGeneID"
colnames(consusExp_SlAtOsPos)[colnames(consusExp_SlAtOsPos)=="AtGeneID.x"] <- "AtGeneID"
colnames(consusExp_SlAtOsPos)[colnames(consusExp_SlAtOsPos)=="SlGeneID.x"] <- "SlGeneID"

write.csv(consusExp_SlAtOsPos,"20200315 1640 consensusExpressologs positive CorrValues_TRANSCRIPT.csv", row.names = F, quote = F)

####################################################################################
# STEP 5.	Use this list to create an expression matrix of the three species
# and test expression correlation among comparable cell types- based on CPM.
# (IMPORTENT NOTE: expressed genes were filtered based on TPM (>2 in at least one cell type) -> 6403 Sl expressologs
#                 visualization of exp cor in comparable cell type is based on CPM -> 6323 Sl expressologs) 
######################################################################################################################
library(limma)
library(Glimma) 
library(pheatmap)
library(viridis)
library(ggplot2)
library(RColorBrewer)
##################################
#Upload expression data, remove unnecessary cell type and add species initials to each column
##################################
################
###   RICE   ###
################
Os_meancpm <- read.csv("/Data/normlog_mean-medianExp_20200311_TRNASCRIPT.csv", row.names = 1,header = T)

################
###  TOMATO  ###
#################################
#Sl logCPM AFTER batch removal
SllogCPM_batchRemoved <- read.csv("/Data/medianlogCPM_BatchRemoved191107_visualPerp.csv", header = T, row.names = 1)

removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}
row.names(SllogCPM_batchRemoved) <- removeDotGene(row.names(SllogCPM_batchRemoved))

SllogCPM_bRem4CT <- SllogCPM_batchRemoved[,c(6,3,9,7,11)]
colnames(SllogCPM_bRem4CT) <- c("Sl_MCO",	"Sl_EN",	"Sl_V",	"Sl_MZ",	"Sl_X35S")

################
# ARABIDOPSIS  #
################
At <- read.csv("/Data/reMeanlog_plastiDrop.csv", header = T, row.names = 1)
At<- At[,c(5,3,2,4,1)]
colnames(At) <- c("At_MCO",	"At_EN",	"At_V",	"At_MZ",	"At_X35S")


###################################################################################################
# CONSENSUS EXPRESSOLOGS WITH POSITIVE CORRELATION VALUES & CELL TYPE EXPRESSION (LogCPM) VALUES
# Create a dataframe that contains the expression of all tomato, rice and  arabidopsis cell types
###################################################################################################

consensus_slepx <- merge(consusExp_SlAtOsPos, SllogCPM_bRem4CT, by.x = "SlGeneID", by.y = "row.names") #1639
consensus_sl_at_epx <- merge(consensus_slepx, At, by.x = "AtGeneID", by.y = "row.names") #1557
consensus_sl_at_os_epx <- merge(consensus_sl_at_epx, Os_meancpm, by.x = "OsGeneID", by.y = "row.names") #1555
write.csv(consensus_sl_at_os_epx,"1555 consensusExpressologs positive CorrValues& expression data_TRANSCRIPT.csv", row.names = F,quote = F)
####################################################################################
## Start of analysis
####################################################################################
# Include EN, MCO, MZ, V & 35S for normalizeBetweenArrays & removeBatchEffect
# Remove 35S for EDA
####################################################################################
#Exclude geneIDs & corr values, and keep only AtGeneID as row names
csa_consensusPos5CT <- consensus_sl_at_os_epx[,c(8:22)]
rownames(csa_consensusPos5CT) <- consensus_sl_at_os_epx$AtGeneID

#Use the function normalizeBetweenArrays to align expression values between samples using quantile normalization
csa_consensusPos5CTNorm <-normalizeBetweenArrays(csa_consensusPos5CT, method="quantile")

#Metadata: 
###### 
meta <- read.csv("/Data/CSA_Metadata.csv", header = T)
head(meta)
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
###### 
# Remove batch effect
###### 
csa_consensusPos5CTNorm_bRemov <- removeBatchEffect(csa_consensusPos5CTNorm, batch=Batch, design=design) 

############################################################################################
# Remove 35S for EDA
####################################################################################
csa_consensusPos4CTNorm_bRemov <- csa_consensusPos5CTNorm_bRemov[,-c(5,10,15)]

meta4CT <- read.csv("//Data/CSA_Metadata4CT.csv", header = T)
head(meta4CT)

#################################
# PCA

transposed_csa_consensusPos4CTNorm_bRemov <- t(csa_consensusPos4CTNorm_bRemov)
# perform prcomp on the transposed matrix from which columns (genes) of zero variance have been removed
pca_proc_consensusPos4CTNorm_bRemov <- prcomp(transposed_csa_consensusPos4CTNorm_bRemov[,apply(transposed_csa_consensusPos4CTNorm_bRemov, 2, var, na.rm=TRUE) != 0],scale=TRUE,center=TRUE) #,retX=TRUE

summary(pca_proc_consensusPos4CTNorm_bRemov)
#Importance of components:
#                          PC1     PC2     PC3     PC4      PC5      PC6     PC7    PC8      PC9     PC10     PC11      PC12
# Standard deviation     22.0305 16.1432 15.1413 13.1654 10.88011 8.74504 7.95363 6.88476 6.81112 5.90073 4.4437 2.381e-14
# Proportion of Variance  0.3121  0.1676  0.1474  0.1115  0.07613 0.04918 0.04068 0.03048 0.02983 0.02239 0.0127 0.000e+00
# Cumulative Proportion   0.3121  0.4797  0.6271  0.7386  0.81473 0.86391 0.90459 0.93508 0.96491 0.98730 1.0000 1.000e+00

plotData_consensusPos4CTNorm_bRemov = meta4CT[,c("Name","Tissue","Species")]
plotData_consensusPos4CTNorm_bRemov$PC1 <- pca_proc_consensusPos4CTNorm_bRemov$x[,1]
plotData_consensusPos4CTNorm_bRemov$PC2 <- pca_proc_consensusPos4CTNorm_bRemov$x[,2]
plotData_consensusPos4CTNorm_bRemov$PC3 <- pca_proc_consensusPos4CTNorm_bRemov$x[,3]
plotData_consensusPos4CTNorm_bRemov$PC4 <- pca_proc_consensusPos4CTNorm_bRemov$x[,4]
plotData_consensusPos4CTNorm_bRemov$PC5 <- pca_proc_consensusPos4CTNorm_bRemov$x[,5]
plotData_consensusPos4CTNorm_bRemov$PC6 <- pca_proc_consensusPos4CTNorm_bRemov$x[,6]

# 2D plots of pairs of principal components
pdf("20200315 PCA qNormlog2CPM 4CT consensusPOSExp bRemoved.pdf")
qplot(PC1,PC2,data=plotData_consensusPos4CTNorm_bRemov,color=Tissue,shape=Species) +
  labs(x="PC1 (31% variability)",y="PC2 (17% variability)") +
  ggtitle("PCA consensusPOSExp 4CT log2 qnorm bRemoved TRANSCRIPT") +
  scale_color_manual(values=c("EN"="#44C3D0","MCO"="#CD85B9","MZ"="#AA1D3F","V"="#19B35A")) + #,"X35S"="#939598"
  scale_shape_manual(values=c("Sly"=15,"Ath"=16,"Osa"=17)) +
  geom_point(size = 5)+
  theme_classic()
dev.off()