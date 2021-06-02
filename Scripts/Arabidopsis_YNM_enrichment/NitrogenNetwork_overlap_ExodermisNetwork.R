library(tidyverse)
library(UpSetR)
setwd("~/Desktop/NitrogenNetwork_Solyc/Scripts/")
source("~/Documents/PhD/code/omicFunctions/generic_functions.R")
#### Get orthologues from YNM (Gaudinier et al 2018)

YNM <- read.delim("../Data/YNM_S03.txt",sep = "\t",header = T,stringsAsFactors = F)
head(YNM)


## Read orthologue maps
# first, try with expressologs
xpressologs <- read.delim("~/Desktop/Orthology_maps/Expressologs/Sly-Ath_expressologs_curated.csv",
                           sep = ",",header = T,stringsAsFactors = F) #91490 orthologous pairs


## Filter to get only positively correlated xpressologs
xpressologs_positive <- xpressologs %>% filter(Corr_value > 0)

###############################################################
###############################################################
### Q1. Which list to use?
SlyAth_xpressologs <- read.delim(
  "~/Desktop/Orthology_maps/Expressologs/20190927_Sl-At_reciprocal_expressologs-phytozomev12.csv",
  sep = ",",header = T,stringsAsFactors = F) #9436 orthologues

commonXpressAth <-intersect(xpressologs$SpeciesA,SlyAth_xpressologs$SpeciesA)
commonXpressSly <-intersect(xpressologs$SpeciesB,SlyAth_xpressologs$SpeciesB)
length(commonXpressAth)
length(commonXpressSly)

## corr values appear to be the same. Ask lidor and come back to this later
xpressologs[xpressologs$SpeciesB == "Solyc04g009650",]
SlyAth_xpressologs[SlyAth_xpressologs$SpeciesB == "Solyc04g009650",]

xpressologs[xpressologs$SpeciesA == "AT4G12310",]
SlyAth_xpressologs[SlyAth_xpressologs$SpeciesA == "AT4G12310",]

###############################################################
### Q2. keep all possible pairs or only those with the highest correlation?
head(xpressologs_positive)

#head(xpressologs_positive[which(duplicated(xpressologs_positive$SpeciesB)),])

xpressologs_positive[xpressologs_positive$SpeciesB == "Solyc03g121050",]
xpressologs_positive[xpressologs_positive$SpeciesA == "AT2G33050",]

###############################################################
###############################################################


# Find expressologs for TF and promoter
YNM_xpressologs <- YNM %>% left_join(xpressologs_positive %>% select(SpeciesA,SpeciesB),
                  by = c("TF_AGI"="SpeciesA")) %>% rename(.,TF_xpressolog_Solyc=SpeciesB) %>% 
        left_join(xpressologs_positive %>% select(SpeciesA,SpeciesB),
                  by = c("Promoter_AGI"="SpeciesA")) %>% rename(.,Promoter_xpressolog_Solyc=SpeciesB) 

nrow(YNM_xpressologs); nrow(YNM)

length(union(YNM_xpressologs$TF_xpressolog_Solyc,YNM_xpressologs$Promoter_xpressolog_Solyc))

## Remove interactions missing any expressolog (TF or Target)
YNM_xpressologs <- YNM_xpressologs %>% filter(!is.na(TF_xpressolog_Solyc) & !is.na(Promoter_xpressolog_Solyc))

length(union(YNM_xpressologs$TF_xpressolog_Solyc,YNM_xpressologs$Promoter_xpressolog_Solyc))

nrow(YNM_xpressologs)
head(YNM_xpressologs)

SolycTFs <- YNM_xpressologs$TF_xpressolog_Solyc
SolycPromoters <- YNM_xpressologs$Promoter_xpressolog_Solyc

Solyc_YNM <- unique(union(SolycTFs,SolycPromoters))

length(unique(union(YNM$TF_AGI,YNM$Promoter_AGI)))#Genes than in the original dataset
length(Solyc_YNM)  # Ath genes with a Sly expressolog in both TF and target


###### Overlap with CTE genes
SolycCTE <- read.csv("~/Desktop/TRAP_sansCTE/Data/20191127-2TPM_BatchModel_CellTypeEnriched.csv",fill = NA,stringsAsFactors = F,header = T)
head(SolycCTE)

SolycCTE_list <- apply(SolycCTE,2,function(x){x[x!='']})


overlayThis <- c(SolycCTE_list,list("Solyc_YNM"=Solyc_YNM))
upset(fromList(overlayThis),
      nsets = length(SolycCTE_list)+1,nintersects = 70)


###### Overlap with Exo network
ExodermisNetwork <- read.delim("../Data/AtlasNetworks/EXO_1kb_promoter_THS_012120.txt",sep="\t",stringsAsFactors = F,header = T)
head(ExodermisNetwork)

Exo_TF <- ExodermisNetwork$Solyc_ID
Exo_Target <- ExodermisNetwork$Solyc_target

ExoNetworkGenes <- unique(union(Exo_TF,Exo_Target))
length(ExoNetworkGenes)

library(gplots)
venn(list("YNM_Sly"=Solyc_YNM,
          "Sly_ExodermisNetwork"=ExoNetworkGenes,
          "CTE_Exodermis"=SolycCTE_list[["EXO"]]))

pdf("../Overlap_Exodermis_NitroSly.pdf",paper = "a4r")
venn(list("YNM_Sly"=Solyc_YNM,
          "Sly_ExodermisNetwork"=ExoNetworkGenes))
dev.off()
### Is this overlap significant?
#Background: genes with a Ath-Sly expressolog & positively correlated
universe <- unique(xpressologs_positive$SpeciesB)
length(universe)

inExo <- universe[universe %in% ExoNetworkGenes]
notExo <- universe[!universe %in% ExoNetworkGenes]


inExo_inNitro <- inExo[inExo %in% Solyc_YNM]
inExo_notNitro <- inExo[!inExo %in% Solyc_YNM]

notExo_inNitro <- notExo[notExo %in% Solyc_YNM]
notExo_notNitro <- notExo[!notExo %in% Solyc_YNM]

## Should be TRUE
length(inExo_inNitro) + length(inExo_notNitro) + 
  length(notExo_inNitro) + length(notExo_notNitro) == length(universe)

(length(inExo_inNitro) / length(inExo_notNitro)  ) / (length(notExo_inNitro) / length(notExo_notNitro))

fisher.test(matrix(c(length(inExo_inNitro) , length(inExo_notNitro),
                     length(notExo_inNitro) , length(notExo_notNitro)),ncol = 2),alternative = "greater")

matrix(c(length(inExo_inNitro) , length(inExo_notNitro),
         length(notExo_inNitro) , length(notExo_notNitro)),ncol=2)

         
#### Methods

# We identified Solanum lycopersicum expressologs of Arabidopsis genes from a transcriptional network controlling nitrogen-associated processes (Gaudinier et al, 2018; doi:10.1038/s41586-018-0656-3). 
# The Arabidopsis nitrogen network contains a total of 429 genes. Of these, 362 have at least one positively correlated expressolog in S. lycopersicum. A total 301 genes have expressolog in both TF-promoter target pairs. We calculated if the overlap between the Sly Exodermis-predicted network genes () and the nitrogen expressologs in tomato was greater than expected by chance using the fisher.test() function in R with alternative="greater".
# 
# 

YNM_inExo <- YNM_xpressologs[YNM_xpressologs$TF_xpressolog_Solyc %in% inExo_inNitro | YNM_xpressologs$Promoter_xpressolog_Solyc %in% inExo_inNitro,]


library(tidyverse)

TFs_inExo <- YNM_xpressologs[YNM_xpressologs$TF_xpressolog_Solyc %in% inExo_inNitro,]
Promoters_inExo <- YNM_xpressologs[YNM_xpressologs$Promoter_xpressolog_Solyc %in% inExo_inNitro,]

TFs_inExo %>% filter(!duplicated(TF_xpressolog_Solyc)) %>% select(c(TF_AGI,TF_Gene_Name,TF_xpressolog_Solyc,Promoter_AGI,Promoter_Gene_Name,Promoter_xpressolog_Solyc))
Promoters_inExo %>% filter(!duplicated(Promoter_xpressolog_Solyc)) %>% select(c(TF_AGI,TF_Gene_Name,TF_xpressolog_Solyc))

intersect(TFs_inExo$TF_xpressolog_Solyc,Promoters_inExo$Promoter_xpressolog_Solyc)
union(TFs_inExo$TF_xpressolog_Solyc,Promoters_inExo$Promoter_xpressolog_Solyc)


### Read expression data
SolycExpression <- read.csv("~/Desktop/dTHS_analysis/Data/20200225_Sly_atlas_log2_median_TPM_batch_corrected_3priorCount.csv",
                            header = T,row.names = 1,sep=",",stringsAsFactors = F)
SolycExpression <- SolycExpression[,-grep("ACT|35S",colnames(SolycExpression))]
#rownames(SolycExpression) <- removeDotGene(rownames(SolycExpression))
SolycExpression <- SolycExpression[,c("EP","EXO","MCO","COR","EN","PH","XY","V","WOX","MZ")] #Reorde



####
TFs_inExo %>% filter(!duplicated(TF_xpressolog_Solyc)) %>% select(c(TF_AGI,TF_Gene_Name,TF_xpressolog_Solyc,Promoter_AGI,Promoter_Gene_Name,Promoter_xpressolog_Solyc))

Promoters_inExo %>% filter(!duplicated(Promoter_xpressolog_Solyc)) %>% select(c(TF_AGI,TF_Gene_Name,TF_xpressolog_Solyc,Promoter_AGI,Promoter_Gene_Name,Promoter_xpressolog_Solyc))

pheatmap::pheatmap(SolycExpression_log2[union(TFs_inExo$TF_xpressolog_Solyc,Promoters_inExo$Promoter_xpressolog_Solyc),],
                   scale = "row")

intersect(TFs_inExo$TF_xpressolog_Solyc,Promoters_inExo$Promoter_xpressolog_Solyc)

union(TFs_inExo$TF_xpressolog_Solyc,Promoters_inExo$Promoter_xpressolog_Solyc)






