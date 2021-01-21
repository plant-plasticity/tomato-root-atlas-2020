#################################
# Automating MapMan reports to be use with FunRich report
# In FunRich tool use custom database with "Slyc_ITAG2.3_FunRicAdj_DotCurat.xlsx" file,
# Atlas cell type enriched gene lists are provided in Table S1
#################################

### SET YOUR WD HERE
setwd("~/mywd/")

## INSTALL THESE PACKAGES IF NOT ALREADY INSTALLED
#https://readxl.tidyverse.org/
#NOTE: you will still need to load readxl explicitly, because it is not a core tidyverse package loaded via 
library(tidyverse)
library(readxl)
library(pheatmap)
library(viridis)
library(dplyr)
  

#####################
#    CELL TYPES     #
#####################

#Read excel workbook
path <- "FunRich report.xlsx"  ############# CHANGE CELL TYPE FILE

celltype<- path %>% 
         excel_sheets() %>% 
         set_names() %>% 
         map(read_excel, path = path, skip = 8,col_names = TRUE)


celltype_fdrCorrected <- lapply(celltype,function(x){

  ## Add a column with the cell type name
  
  x$CellType <- "cell type"    ############# CHANGE CELL TYPE NAME
  
  ## Add column with the hierarchy/ level of the terms
  
  x$hierarchy <- colnames(x[1])
  
  ## Fix colname of the first column
  colnames(x)[1] <- "Term"
  
  ## Filter by FE
  idxKeep <- x[,5] > 1
  
  ## Recalculate FDR
  x$CorrFDR_FE <- 1
  x[idxKeep,"CorrFDR_FE"] <- p.adjust(as.numeric(x[idxKeep,6][[1]]),method = "fdr")
  return(x)
})

##Concatenate worksheets into one data frame

CELLTYPE_FDRcor <- do.call("rbind", celltype_fdrCorrected) ############# CHANGE CELL TYPE NAME
write.csv(CELLTYPE_FDRcor, "FILENAME.csv") ############# CHANGE CELL TYPE NAMEX2


#Concatenate dataframe of different cell types into one data frame

ALLcelltypes <-rbind(CELLTYPE1_FDRcor, CELLTYPE2_FDRcor, CELLTYPE3_FDRcor)
write.csv(ALLcelltypes, "FILENAME.csv", quote = F, row.names = T)


############################################################################
#SELECT SIGNIFICANT TERMS FROM ALL FOUR HIERARCHIES & CELL TYPES

sigAllCT <- subset(ALLcelltypes, select = c("Term", "CellType", "CorrFDR_FE"), 
                     CorrFDR_FE <= 0.15) # SET CorrFDR_FE THRESHOLD

write.csv(sigAllCT, "2.FilteredReport/20191220 Sly atlas FE1 FDR0.15 2TPM batchModel ALL terms & CT.csv", quote = F, row.names = F)
