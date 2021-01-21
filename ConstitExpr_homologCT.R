############################################################################################################################
# Identify constitutively expressed genes based on TPM values for Sly and Osa & unlogged RMA dataset for Ath
############################################################################################################################
setwd("/mywd")

library(pheatmap)
library(viridis)
library(ggplot2)

bgfold=1.5 # fold change used to define background gene set

### Calculate fold differences
getfold<-function(x){
  x<-x+0.5 # add pseudo count to avoid 0s
  return(max(x)/min(x))
}

######################################
# RICE constitutively expressed genes_TRANSCRIPT
######################################
os_normTPM <- read.csv("~/Postdoc/Projects/Plasticity/Atlas/Cell-type lists/CSA/Rice/KALLISTO_COUNTS_transcript/20200311_Os_Atlas_IRGSP-1.0_TRANSCRIPT_Kallisto_qnorm_MeanMedian_tpm_genes.csv",header = T, row.names = 1)

geneFoldTPM<-apply(os_normTPM,1,getfold)
head(geneFoldTPM)
# Os01g0100100 Os01g0100200 Os01g0100300 Os01g0100400 Os01g0100466 Os01g0100500 
# 1.592720     3.572647     2.153418     1.776758     1.683949     1.532848 

geneMedianTPM<-apply(os_normTPM,1,median)
 head(geneMedianTPM)
 # Os01g0100100 Os01g0100200 Os01g0100300 Os01g0100400 Os01g0100466 Os01g0100500 
 # 9.2228949    0.3799700    0.1457228    3.6396271    0.1160011    8.0581832 
 
 quantile(geneMedianTPM)
 # 0%          25%          50%          75%         100% 
 # 0.000000e+00 8.019725e-03 1.611994e+00 7.513434e+00 7.849941e+05 

##########################################################################################
# constitutive expression (not changing significantly between cell types)
# genes within 1.5 fold change and median expression > 1.611994 (to filter out lowly expressed genes)
##########################################################################################
OsconstitgeneTPM50_1.5<-names(geneFoldTPM[geneFoldTPM<1.5&geneMedianTPM>1.611994]) ##get genes within 1.5 fold change & median of median expression > 1.64268
constitgeneOsTPM50_1.5<-os_normTPM[OsconstitgeneTPM50_1.5,] #constitgeneOsTPM=1523

##################
log2(min(constitgeneOsTPM50_1.5))
log2(max(constitgeneOsTPM50_1.5))

pdf("/Osa_CEGs.pdf")
colors <- c(min(log2(constitgeneOsTPM50_1.5)),seq(0.075,9.4,by=0.5),max(log2(constitgeneOsTPM50_1.5)))
my_palette <- viridis(length(colors)-2,begin = 0.5, end = 1, direction = 1)

pheatmap::pheatmap(log2(constitgeneOsTPM50_1.5),
                   #clustering_distance_cols = "correlation",
                   color = my_palette,
                   cellwidth = 40,
                   #cellheight = 0.5,
                   cluster_rows = T,
                   cluster_cols = F,
                   fontsize = 10,
                   angle_col = 45,
                   fontsize_row = 0.005,
                   breaks = colors)
dev.off()


######################################
# TOMATO constitutively expressed genes. TPM WITHOUT BATCH CORRECTION!!!!!!!!!!!!!!!
######################################

sl_normTPM <- read.csv("~/Postdoc/Projects/Plasticity/Atlas/Resources/Atlas count and TPM/20180920_ATLAS_ITAG3.2_Kallisto_quantile_norm_median_tpm.csv",header = T, row.names = 1)
sl_normTPM <- sl_normTPM[,grep("MCO|MZ|EN|V|X35S",colnames(sl_normTPM))]

 slgeneFoldTPM<-apply(sl_normTPM,1,getfold)
 head(slgeneFoldTPM)
# Solyc00g005000.3.1 Solyc00g005040.3.1 Solyc00g005050.3.1 Solyc00g005060.1.1 Solyc00g005080.2.1 Solyc00g005084.1.1
# 16.469612           1.000000           3.279315           1.000000           7.200374           1.000000

 slgeneMedianTPM<-apply(sl_normTPM,1,median)
 head(slgeneMedianTPM)
# Solyc00g005000.3.1 Solyc00g005040.3.1 Solyc00g005050.3.1 Solyc00g005060.1.1 Solyc00g005080.2.1 Solyc00g005084.1.1
# 16.61511            0.00000           33.16510            0.00000            0.00000            0.00000
 quantile(slgeneMedianTPM)
 # 0%          25%          50%          75%         100%
 # 0.000000e+00 0.000000e+00 1.126064e-01 9.097249e+00 6.479046e+04

 ##########################################################################################
# constitutive expression (not changing significantly between cell types)
# genes within 1.5 fold change and median expressoin > 1.126064. Filter out lowly expressed genes
##########################################################################################
slconstitgeneTPM50_1.5<-names(slgeneFoldTPM[slgeneFoldTPM<1.5&slgeneMedianTPM>1.126064]) #bgmedian, genes within 2 fold change & median expression > median
constitgeneSlTPM50_1.5<-sl_normTPM[slconstitgeneTPM50_1.5,] #308 genes

#Remove genome build
removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}

rownames(constitgeneSlTPM50_1.5) <- removeDotGene(rownames(constitgeneSlTPM50_1.5))

################
constitgeneSlTPM50_1.5 <- data.frame(constitgeneSlTPM50_1.5[,c(3,2,5,4,6)], row.names = constitgeneSlTPM50_1.5$X)
log2(min(constitgeneSlTPM50_1.5)) 
log2(max(constitgeneSlTPM50_1.5)) 

pdf("../Figures/Sly CEGs.pdf")
colors <- c(min(log2(constitgeneSlTPM50_1.5)),seq(0.28,9.56,by=0.5),max(log2(constitgeneSlTPM50_1.5)))
my_palette <- viridis(length(colors)-2,begin = 0.5, end = 1, direction = 1)

pheatmap::pheatmap(log2(constitgeneSlTPM50_1.5),
                   #clustering_distance_cols = "correlation",
                   color = my_palette,
                   cellwidth = 40,
                   #cellheight = 0.5,
                   cluster_rows = T,
                   cluster_cols = F,
                   fontsize = 10,
                   fontsize_row = 0.01,
                   angle_col = 45,
                   breaks = colors)
dev.off()

###############################################
# ARABIDOPSIS constitutively expressed genes
###############################################
at_CT <- read.csv("~/Postdoc/Projects/Plasticity/Atlas/Cell-type lists/CSA/Arabidopsis/Mustroph2009_TRAP/reMustroph/reMeanUNlog_plastiDrop.csv",header = T, row.names = 1)

atgeneFold<-apply(at_CT,1,getfold)

head(atgeneFold)
# AT3G52770 AT3G19780 AT1G51000 AT5G45570 AT1G48720 AT1G23910 
# 1.208352  1.084060  1.013325  1.136026  1.136927  1.117190 

atgeneMedian<-apply(at_CT,1,median)

head(atgeneMedian)
# AT3G52770 AT3G19780 AT1G51000 AT5G45570 AT1G48720 AT1G23910 
#13.60054  13.45748  12.76526  14.25249  12.92506  14.16599 
quantile(atgeneMedian)
# 0%         25%         50%         75%        100% 
# 12.42858    40.42747   132.37255   387.31739 18670.97750

##########################################################################################
# constitutive expression (not changing significantly between cell types)
# genes withing 1.5 fold change and median expressoin > median. Filter out lowly expressed genes
##########################################################################################
atconstitgene50_1.5<-names(atgeneFold[atgeneFold<1.5&atgeneMedian>132.37255]) #genes within 1.5 fold change & median expression > 132.37255
constitgeneAt50_1.5<-at_CT[atconstitgene50_1.5,] #1154

constitgeneAt50_1.5 <- data.frame(constitgeneAt50_1.5[,c(6,4,3,5,2)], row.names = constitgeneAt50_1.5$X)

log2(min(constitgeneAt50_1.5))
log2(max(constitgeneAt50_1.5)) 

pdf("../Figures/Ath_CEGs.pdf")
colors <- c(min(log2(constitgeneAt50_1.5)),seq(6.53,14.23,by=0.5),max(log2(constitgeneAt50_1.5)))
my_palette <- viridis(length(colors)-2,begin = 0.5, end = 1, direction = 1)

pheatmap::pheatmap(log2(constitgeneAt50_1.5),
                   #clustering_distance_cols = "correlation",
                   color = my_palette,
                   cellwidth = 40,
                   #cellheight = 0.5,
                   cluster_rows = T,
                   cluster_cols = F,
                   fontsize = 10,
                   fontsize_row = 0.01,
                   angle_col = 45,
                   breaks = colors)
dev.off()


############################################################################################################
# Test overlaps among orthologs
############################################################################################################
# Convert Sly and Osa to Ath

# Use Pos. corr. in-house expressologs Ath as a reference test how many genes can be detected for Sly and Osa
# 7,295 "20200117 At-Sl reciprocal positivCorr expressologs-phytozomev12"
AtSlExprss <- read.csv("At-Sl reciprocal positivCorr expressologs-phytozomev12.csv", header = T)
AtSlConstitExpr <- merge(constitgeneSlTPM50_1.5, AtSlExprss, by.x = "row.names", by.y = "SlGeneID") #153 out of 308

# 6,059 "20200316 dedup At-Os reciprocal positivCorr expressologs-phytozomev12"
AtOsExprss <- read.csv("dedup At-Os reciprocal positivCorr expressologs-phytozomev12.csv ", header = T)

summary(duplicated(AtOsExprss$OsGeneID)) #324 duplicated Osa genes.
# Mode   FALSE    TRUE 
# logical    6070     324 

AtOsConstitExpr <- merge(constitgeneOsTPM50_1.5, AtOsExprss, by.x = "row.names", by.y = "OsGeneID") #547 out of 1523

summary(duplicated(AtOsConstitExpr$Row.names)) #Osa
#   Mode   FALSE    TRUE 
# logical    523      24  

#Use http://bioinformatics.psb.ugent.be/cgi-bin/liste/Venn/calculate_venn.htpl to creat a Venn diagram and text file with the intersections
library(systemPipeR)

SlyConstAth=AtSlConstitExpr$AtGeneID
AthConstAth=row.names(constitgeneAt50_1.5)
OsaConstAth=AtOsConstitExpr$AtGeneID

namesAth=c("Sly, 153 (out of 308) CEGs","Ath, 1,154 CEGs","Osa 547 (out of 1,523) CEGs")

Athclus_up_each_any=list(as.character(SlyConstAth),as.character(AthConstAth),as.character(OsaConstAth))
names(Athclus_up_each_any)=namesAth

pdf(paste0("Venn CEGs Ath as ref.pdf"))
vennset_up_anyAth <- overLapper(Athclus_up_each_any,type="vennsets")
vennPlot(list(vennset_up_anyAth), mymain="Costitutively expressed genes Ath as a reference", mysub=paste0(""),lcol=c("#b35806","#e08214","#fdb863","#fdb863","#8073ac"),lines=c("#b35806","#e08214","#fdb863","#fdb863","#8073ac"), ccol=c("black"))
dev.off()
