###########################################################################
# Detection of tomato expressologs (homolog genes with the best SCC values)
###########################################################################

#Process Sl expressolog data
#Provart's lab has data for 14,939 of 35,768 tomato genes.
#For every tomato gene, rice and Arabidopsis expressologs (and non-expressologs) were retrieved. 
#The % identity and SCC score for each species gene pair are provided.
#The expressolog is the homolog with the best SCC value (2nd last column).

setwd("/mywd/")

data <- read.csv("Data/Expressologs-Sl_At_Os.csv", header = T)

#Subset the data for Arabidopsis and rice
dataAt <- data[which(data$Genome_B == "Arabidopsis"),]
dataOs <- data[which(data$Genome_B == "Rice"),]

#create a character of unique Sl genes with At/Os orthologs
uniqu_Atgene <- unique(dataAt$Gene_A)
uniqu_Osgene <- unique(dataOs$Gene_A)

# ARABIDOPSIS 
#create a dataframe with the unique gene names as row names                     
Atgene_df <- data.frame(row.names=uniqu_Atgene) 
#colnames(data)
Atgene_df[,c("Gene_A","Gene_B", "Probeset_B", "Genome_B", "SCC_Value", "Seq_Similarity")]<- NA

for (i in 1:length(rownames(Atgene_df))){
  tmp <- dataAt[which(dataAt$Gene_A==rownames(Atgene_df[i,])),]
  Atgene_df[i,] <- lapply(tmp[which.max(tmp$SCC_Value),],as.character)
}

# RICE 
# create a dataframe with the unique gene names as row names                     
Osgene_df <- data.frame(row.names=uniqu_Osgene) 
#colnames(data)
Osgene_df[,c("Gene_A","Gene_B", "Probeset_B", "Genome_B", "SCC_Value", "Seq_Similarity")]<- NA

for (i in 1:length(rownames(Osgene_df))){
  tmp <- dataOs[which(dataOs$Gene_A==rownames(Osgene_df[i,])),]
  Osgene_df[i,] <- lapply(tmp[which.max(tmp$SCC_Value),],as.character)
}

#load rice NetAffx Annotation File with gene name to be matched to probeses id
riceanno <- read.csv("Data/Rice-na36-annot-csv/Rice.na36.annot.csv", skip = 18, header = T)
riceannodrop <- riceanno[,c(1,14,15,31,32,33)]
#Match probeses IDs and creat a table that contains ric gene ID (RAP-DB) and GOs
Osgene_dfAnno <- merge(Osgene_df, riceannodrop, by.x = "Probeset_B", by.y = "Probe.Set.ID")

#Sly genes and their Ath positivCorr expressologs

library(dplyr)

Atgene_dfPos <- Atgene_df %>% select(Gene_A,Gene_B, SCC_Value, Seq_Similarity) %>%
  filter(SCC_Value > 0) #11250 rows

colnames(Atgene_dfPos)[1] <- "SlGeneID"
colnames(Atgene_dfPos)[2] <- "AtGeneID"
colnames(Atgene_dfPos)[3] <- "SlAt_SCC_Value"
colnames(Atgene_dfPos)[4] <- "SlAt_Seq_Similarity"


#Sly genes and their Osa positivCorr expressologs
Osgene_dfAnnoPos <- Osgene_dfAnno %>% select(Gene_A,Gene.Title,Gene_B,SCC_Value, Seq_Similarity) %>%
  filter(SCC_Value > 0) #9383 rows

colnames(Osgene_dfAnnoPos)[1] <- "SlGeneID"
colnames(Osgene_dfAnnoPos)[2] <- "OsGeneID"
colnames(Osgene_dfAnnoPos)[3] <- "LOC_OsGeneID"
colnames(Osgene_dfAnnoPos)[4] <- "SlOs_SCC_Value"
colnames(Osgene_dfAnnoPos)[5] <- "SlOs_Seq_Similarity"

#151 rice genes are not unambiguously defined and 496 genes don't have orthologs. 
#In total 647 genes were excluded from further analysis
#grepl is like grep, but it returns a logical vector (use TRUE/FALSE subsetting instead of numeric). grepl works with one pattern at the time 

Slexpressolog_Os0 <- Osgene_dfAnnoPos[(!grepl(c("///"),Osgene_dfAnnoPos$OsGeneID)),] #9261 expressologs
Osgene_dfAnnoPosCur <- Slexpressolog_Os0[(!grepl(c("---"),Slexpressolog_Os0$OsGeneID)),] #8874 expressologs

#Merge rice and Arabidopsis expressolog tables
At_Osgene_df <- merge(Atgene_dfPos, Osgene_dfAnnoPosdrop, by = "Gene_A") #6968 rows
At_Osgene_dfdrop <- At_Osgene_df[,-c(3,4,7,9)]
colnames(At_Osgene_dfdrop)[1] <- "SlGeneID"
colnames(At_Osgene_dfdrop)[2] <- "AtGeneID"
colnames(At_Osgene_dfdrop)[5] <- "LOC_OsGeneID"
colnames(At_Osgene_dfdrop)[8] <- "OsGeneID"


write.csv(At_Osgene_dfdrop, "Sly genes and their Ath & Osa positive corr expressologs.csv", quote = F, row.names = F)
