library(tidyverse)
source("roku_functions.R")
dir.create("../roku_out")



geneTPM <- read.csv("../../Data/allCellTypes.q2.noCh0.masked.25mil.counts.per.uTHS.wX35S.csv",
                    row.names = 1,as.is = T,
                    stringsAsFactors = F)
### Fix missing data 
geneTPM[is.na(geneTPM)] <- 0
geneTPM[geneTPM=="NAN"] <- 0
geneTPM[geneTPM==""] <- 0

geneTPM <- geneTPM[,-grep("35S",colnames(geneTPM))]

peaks <- rownames(geneTPM)
geneTPM <- data.frame(apply(geneTPM,2,as.numeric),row.names = peaks)
head(geneTPM)

uniqCTs <- unique(gsub("\\..*$","",colnames(geneTPM)))
geneTPM_median <- matrix(NA,nrow =nrow(geneTPM) ,ncol = length(uniqCTs),dimnames = list(rownames(geneTPM),uniqCTs))

for (CT in uniqCTs){
  print(CT)
  geneTPM_median[,grep(CT,colnames(geneTPM_median))] <- apply(geneTPM[,grep(CT,colnames(geneTPM))],1,median)
  
}
apply(head(geneTPM[,1:4]),1,median)

print(dim(geneTPM_median))
head(geneTPM_median) 
geneTPM <- geneTPM_median

###
### Do not log2 transform, the functions do that.
delta=1 #quantile(rowMeans(geneTPM),0.5)#0.05
lowexp <- 1 # threshold for low expression values
bgfold=2 # fold change used to define background gene set
bgmedian=3 # threshold for median expression level #quantile(geneTPM,0.25)
pvalue=0.05 # p value

##########
ncelltype<-ncol(geneTPM) # get number of tissue/cell types

# a cell type specific genes can only be specific to less than half of the cell types
nmax<-floor(ncelltype/2)
tmp<-geneTPM;tmp[tmp>=lowexp]=1;tmp[tmp<lowexp]=0;tmpsum<-apply(tmp,1,sum)
table(tmpsum)
geneTPMkeep<-geneTPM[tmpsum>0,]

### Entropy
getfold<-function(x)
{
  x<-x+0.5 # add pseudo count to avoid 0s
  return(max(x)/min(x))
}
geneFold<-apply(geneTPM,1,getfold)
geneMedian<-apply(geneTPM,1,median)


head(geneFold)
head(geneMedian)

### AKA constitutive expression (not changing significantly between cell tyoes) - jrm/acp/patito
###########################################################
# get background genes as those withing 2 fold change and median expressoin > 0.5 
# this filter is to not use really lowly expressed genes as inputs.
###########################################################

bggene<-names(geneFold[geneFold<bgfold&geneMedian>bgmedian]) # get background genes as those withing 2 fold change
bggeneTPM<-geneTPM[bggene,]
print(paste('number of genes in the background set:', 
            nrow(bggeneTPM)))

print(nrow(bggeneTPM)/nrow(geneTPM))

# wh <- 1
# pdf("../roku_out/background_log2-TPM.pdf",paper = "A4")
# pheatmap::pheatmap(log2(bggeneTPM),
#                    cellwidth = wh*4,cellheight = wh,border_color = NA,
#                    show_rownames = F,fontsize_col = 4)
# dev.off()


# get the WH for for shaved exprssion
bgWHmat<-weightcenterentropy_log_shave(bggeneTPM, ntoshave=nmax,3)

# get threshold at a pvalue, default pvalue = 0.001
TR_0001<-apply(bgWHmat,2,quantile,prob=pvalue)


###########################################################
# Get tissue specific genes.
###########################################################
geneCellSpecific_P0001<-weightcenterentropy_log_shave_identify(geneTPM,TR_0001,delta)
dim(geneCellSpecific_P0001$index)
table(rowSums(geneCellSpecific_P0001$index))

# -4     -3     -2     -1      0      1      2      3      4 
# 14     87    192    411 106812    523    214     71     11 
# 

###########################################################
# Get tissue specific genes.
###########################################################

geneCellSpecific_P0001<-weightcenterentropy_log_shave_identify(geneTPM,TR_0001,delta)
head(geneCellSpecific_P0001$index)


outputfn <- '../roku_out/Results_roku.RData' 
save(geneCellSpecific_P0001,file=outputfn)
outputfn <- '../roku_out/outputRoku.txt' 
write.table(geneCellSpecific_P0001$index,file=outputfn,sep='\t')

####
outliers <- apply(geneCellSpecific_P0001$index==0,1,all)
outliers <- geneCellSpecific_P0001$index[!outliers,]
uniqueOutliers <- outliers[rowSums(outliers != 0) == 1,]

enrichedCT <- list()
for (CT in colnames(outliers)){
  enrichedCT[[CT]] <- uniqueOutliers[uniqueOutliers[,CT] == 1,]
}

nCTE <- sapply(enrichedCT,nrow)
tmp <- data.frame("CT"=factor(names(nCTE),levels = rev(names(nCTE))),"n"=nCTE)
pdf("../roku_out/roku_uniqueOutliers.pdf")
tmp %>%
  ggplot(aes(x=CT)) +
  geom_bar(stat = "identity",aes(y=n),width = 0.5,fill="skyblue")  +
  geom_text(aes(x=CT,y=n*.70,label=n)) +
  coord_flip() +
  theme_classic() +
  theme(aspect.ratio = 24/8) +
  NULL
dev.off()

allEnrichedRoku <- do.call("rbind",enrichedCT)
write.table(file="../roku_out/allEnrichedRoku.txt",allEnrichedRoku,
            sep="\t",quote=FALSE,col.names=NA)


myColor <- RColorBrewer::brewer.pal(9,"BuPu")
wh <- 4
pdf("../roku_out/roku_enrichedByCT.pdf")
for (CT in names(enrichedCT)) {
  pheatmap::pheatmap(log2(geneTPM_median[rownames(enrichedCT[[CT]]),]+1),
                     border_color = NA,#scale = "row", 
                     cellheight = wh,cellwidth = wh*2,  
                     show_rownames = F,
                     color = c(myColor), main=paste0("roku enriched",CT," (log2)"))
  
  pheatmap::pheatmap( geneTPM_median[rownames(enrichedCT[[CT]]),] ,
                      border_color = NA,#scale = "row", 
                      cellheight = wh,cellwidth = wh*2, 
                      show_rownames = F,
                      color = c(myColor), main=paste0("roku enriched",CT))
  
}

dev.off()
