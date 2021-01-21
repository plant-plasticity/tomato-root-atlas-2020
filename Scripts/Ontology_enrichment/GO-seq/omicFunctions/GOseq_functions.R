## To install goseq:
# Make sure a SQL client is installed.
# for mac use `brew install mariadb-connector-c``
### 
##source("https://bioconductor.org/biocLite.R")
# install.packages('RMySQL', type='source')
#biocLite(c("goseq","GenomicRanges","Rsamtools","rtracklayer"))
library(goseq)
library(GenomicRanges)
library(Rsamtools)
library(rtracklayer)

##### About the readGO function #####
### The file format that the readGO function processes is as follows:
# First column are the gene identifiers
# Second column the GO terms associated with the gene, separated by a coma and might contain a space or not

#                 ITAG                        GO
# 1 Solyc00g005000.2.1 0019538, 0009056, 0016787
# 1 Solyc00g005000.2.1 0006508,0004190

############ !!!!!! ############
### If the file you're trying to read has a different format, modify the regex into an appropiate form to match your format. - DON'T PANIC.


##### About the call_goseq function #####
# The removeDots is a logical argument (TRUE/FALSE). It can be turned on if the GFF or GO annotation is from a different genome release than the used in the experiments 
# (ie, when GO annotation is from a previous release -gene id ID#####.2.1- and RNASeq was done in a release v3.1 -gene id ID#####.3.1-)
# The parameter is intended to 'normalize' the annotations to avoid resulting in few overlaps 
####


######################## Here be functions ########################

######## Read GFF and get gene lengths ########

## This function takes in the name of a GFF file to read and process the lengths of the genes to be used by goseq 
# the name must be quoted and must be inside the GOdata/ directory 
# GFFfile <- "Name_of_file.gff"

getGeneLengthfromGFF <- function(GFFPath) {
  ## Requires rtracklayer
  require(rtracklayer)
  ###
  GFFPath <- #paste0("GOdata/",gffFileName)
  GFF <- import.gff(GFFPath,version="3",feature.type="gene")
  grl <- GenomicRanges::reduce(split(GFF, mcols(GFF)$Name)) #Conflict with reduce between GRanges and tidyverse
  reducedGTF <- unlist(grl, use.names=T)
  mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
  reducedGTF
  AllLengths <- width(reducedGTF)
  names(AllLengths) <- mcols(reducedGTF)$Name
  return(AllLengths)
}



######## Read GO file and process it ########
## This function takes in the name of a GO file to read and process it to a format that goseq likes
# the name must be quoted and must be inside the "GOdata/" directory 
# goFileName <- "Name_of_file.gff"

readGO <- function(goFileName){
  
  GOitag <- read.table( goFileName ,
                        stringsAsFactors = F,header = T)
  
  GOitagSplit <- split(GOitag,GOitag[,1])
  tmp <- GOitagSplit[[1]]
  tst <- lapply(GOitagSplit,function(tmp){cbind( 
    tmp[,1],unlist(     # Split by regex:
      strsplit(tmp[,2],",[[:space:]]{0,}"))) #Split by a coma and any or no space after. 
    # This makes it work for both files I got from Ted's data, but just in case any other GO file might be using tabs or any type of space.
    # If
  })
  GOitagNew <- data.frame(do.call(rbind,tst),stringsAsFactors = F)
  colnames(GOitagNew) <- c("ITAG","GOID")
  head(GOitagNew)
  GOitagNew[,2] <- paste0("GO:",GOitagNew[,2])
  cat ("Number of GO terms in data:",nrow(GOitagNew),"\n")
  return(GOitagNew)
}




######## Call goseq functions and perform enrichment testing ########

## Calls go seq functions
call_goseq <- function(genesToAnalyze,assayed.genes,AllLengths,go.goseq,removeDots=TRUE,plotFit=FALSE){
  ## Intended to normalize names if using info from different builds. Doesn't hurt if the genes names are the same across different data sets
  if (removeDots){
    cat ("--- \n Removing genome version info from gene names (GeneId1234.genomeVersion3 --> GeneId1234)\n---\n")
    
    assayed.genes <- removeDotGene(assayed.genes)
    genesToAnalyze <- removeDotGene(genesToAnalyze)
    names(AllLengths) <- removeDotGene ( names(AllLengths) )
    go.goseq$ITAG <- removeDotGene (go.goseq$ITAG)
  } else cat ("Leaving gene ID's as is")
  
  
  cat("Total number of genes analyzed (universe):", length(assayed.genes),"\n")
  cat("Number of genes used for testing", length(genesToAnalyze),"\n")
  
  #
  gene.vector=as.integer(assayed.genes%in%genesToAnalyze)
  names(gene.vector)=assayed.genes
  #head(gene.vector)
  #table(gene.vector)
  
  ### Filter length info from all the genes in the genome to those analyzed (present in the universe)
  GeneLengths <- AllLengths[names(gene.vector)]
  
  ### Calls go seq function
  pwf=nullp(gene.vector, bias.data=GeneLengths,plot.fit = plotFit) #By default doens't plot the fit 
  go.all <- goseq(pwf, gene2cat=go.goseq)
  
  # Return results
  return(go.all)
}



########## Functions
############################################################

##### Extract significant GO terms from a list of GO enrichment tests // output from GOseq
#########################
get_significant_GOterms <- function(GOresults,pValColumn="over_represented_pvalue",pThreshold = 0.05) {
  lapply(GOresults, function(x){  
    pIdx <- which(colnames(x) %in% pValColumn)
    filterP <- x[,pIdx] < pThreshold
    return(x[filterP,])})
}
####################################################################################################

##### Extract a matrix of p values (or FDR,etc) to prepare as input for a heatmap/tileplot
#########################
get_pVals_matrix <- function(GOresults,pValColumn="over_represented_pvalue") {
  ## Function to return a matrix with all p values from the GO enrichment testing
  # The output can be processed for p value adjustment and further filtered and passed to a heatmap/tile plot function for vizualization
  pValsList <- lapply(GOresults, function(y){  
    
    byOnto <- split(y,y$ontology)
    
    pVals_byOnto <- lapply(byOnto,function(x) {
      pIdx <- which(colnames(x) %in% pValColumn)
      tmp <- x[,pIdx,drop=F]
      rownames(tmp) <- ifelse(is.na(x$term),x$category,x$term)
      return(tmp)
    })
    return(pVals_byOnto)
  })
  
  pValsMatrix_BP <- aggregate_list_tables(lapply(pValsList,"[[","BP"))
  pValsMatrix_MF <- aggregate_list_tables(lapply(pValsList,"[[","MF"))
  pValsMatrix_CC <- aggregate_list_tables(lapply(pValsList,"[[","CC"))
  return(
    list("BP"=pValsMatrix_BP,
         "MF"=pValsMatrix_MF,
         "CC"=pValsMatrix_CC)
  )
}
####################################################################################################

##### After testing multiple gene sets for enriched GO categories, see the most frequent terms by ontology
#########################
get_significanGO_frequencies <- function(GOsignificantList,GOtermColumn="term",GOontologyColumn="ontology"){
  # Check libs
  require(tidyverse)
  ## Aggregate all into
  aggregateAll <- do.call("rbind",GOsignificantList)
  nrow(aggregateAll)==sum(sapply(GOsignificantList, nrow)) # Should be TRUE
  ## Get frequency of terms
  freqTerms <- aggregateAll %>% filter(!is.na(.[GOtermColumn])) %>% 
    select(GOtermColumn,GOontologyColumn) %>% 
    group_by_all() %>% 
    summarise("Count"=n()) %>% arrange(.[[2]],desc(Count))
  return(freqTerms)
}
####################################################################################################



######## Find which genes were tested on each category
## signGOList is a list of all the GO terms that were significant when testing genes. 
## listOfGenes is a table where each column is a list of genes of interest that were used for the GO enrichment testing
#########################
get_Genes_Per_Category <- function(signGOList,listOfGenes,GOdata,GOColumn="category") {
  require(tidyverse)
  
  ##### Cats by Term
  # Dictionary of terms to GOID
  test <- lapply(signGOList,"[",c("category","term"))
  catBYterm <- lapply(test,function(x){
    split(x$category,x$term)  
  })
  uniqueTerms <- Reduce(union,lapply(catBYterm,names))
  uniqueGOID <- unlist(Reduce(union,lapply(catBYterm,"[")))
  names(uniqueTerms) <- uniqueGOID
  #########################
  getGenesPerCategory <- lapply(seq_along(signGOList),function(eachSet){
    
    cat("Genes per category in pattern", names(signGOList)[eachSet],"\n")
    ## Get DE genes
    de.genes <- listOfGenes[listOfGenes[,eachSet]!="",eachSet]
    
    ### Get categories
    ############
    if (nrow(signGOList[[eachSet]]) !=0){
      
      categories <- signGOList[[eachSet]][,GOColumn]
      categories <- categories[!is.na(categories)]
      
      ## From the full GO seq list, get the genes that overlap with the DE
      DEgenesInCategory <- GOdata %>%
        ## GOColumn="category" is category and not term because the column in the GO table is GOID
        filter(.$GOID %in% categories) %>% group_by(GOID) %>%  
        filter(removeDotGene(ITAG) %in% removeDotGene(de.genes)) %>% nest()
      
      #### Extract the lists into a matrix where the column names are the categories and contains the genes used for each test
      names(DEgenesInCategory$data) <- DEgenesInCategory$GOID
      GenesPerCategory <- condenseGeneList_toMatrix(lapply(DEgenesInCategory$data,as.matrix))
      
      
      colnames(GenesPerCategory) <- uniqueTerms[colnames(GenesPerCategory)]
      head(GenesPerCategory)
      return(GenesPerCategory)
    }
  })
  names(getGenesPerCategory) <- names(signGOList)
  return(getGenesPerCategory)
}
####################################################################################################

make_GOfrequency_plot <- function(GOfrequencies,plotTitle="Frequency of Categories",
                                  minFreq=2,
                                  xtraName="",
                                  filterByKeyWord=F,
                                  filterbyOnto="BP",
                                  keywords="*"){
  require(ggplot2)
  require(tidyverse)
  ### Adjust plot name
  if (xtraName != "") plotTitle <- paste(plotTitle,xtraName,sep=" - ")
  
  ######## GO term frequency
  ## Filter by a minimum frequency (default: 2) - ie, gets rid of categories that are unique to a single tested set.
  tmpCounts <- GOfrequencies
  if (max(GOfrequencies$Count)>=minFreq) tmpCounts <- GOfrequencies %>% filter(Count >= minFreq) 
  
  ###### Filter by Keywords
  if (filterByKeyWord) {
    searchThis <- paste(keywords,collapse = "|") 
    plotTitle <- paste(plotTitle, paste("keywords: ", paste(keywords,collapse=",") ,sep = ""), sep="\n")
    tmpCounts <- tmpCounts %>%  filter(grepl(searchThis,term,ignore.case = T))
  } else {
    cat ("Not filtering by keywords \n")
  }
  
  ###### Filter by ontology
  tmpCounts <- tmpCounts %>% filter(ontology %in% filterbyOnto)
  
  ## Make plot
  plotGOFreq <- tmpCounts %>% 
    ggplot(.,aes(x=reorder(term,Count), # Order X axis labels by Count
                 y=Count, #Counts go in the Y axis
                 fill=ontology)) + #Color by ontology
    geom_bar(stat="identity") + coord_flip() + # Do a bar plot 
    xlab("GO category") + xlab("Frequency") + labs(title=plotTitle)  + 
    facet_grid(.~ontology,scales="free_y") #Separate by ontology
  return(plotGOFreq) 
}



make_pHeatmap <- function(matrixData,plotTitle="Heatmap Of Categories",col.pal=RColorBrewer::brewer.pal(9,"Blues"),
                          filterByKeyWord=F, xtraName="", keywords="*",
                          filterbyOnto=c("BP","MF","CC")){
  require(ggplot2)
  require(tidyverse)
  require(pheatmap)
  
  matrixData <- matrixData[names(matrixData) %in% filterbyOnto]
  cat("pheatmap on p values of:",names(matrixData),"\n")
  
  lapply(seq_along(matrixData), function(each){
    onto <- names(matrixData)[each]
    tmpData <- matrixData[[each]]
    
    ### Set title
    if (xtraName != "") {plotTitle <- paste(plotTitle,xtraName,onto,sep=" - ")
    } else plotTitle <- paste(plotTitle,onto,sep=" - ")
    
    
    ###### Filter by Keywords
    if (filterByKeyWord) {
      cat ("Filtering by keywords: ",keywords,"\n")
      searchThis <- paste(keywords,collapse = "|") 
      plotTitle <- paste(plotTitle, paste("keywords: ", paste(keywords,collapse=",") ,sep = ""), sep="\n") #Update plot title
      tmpData <- tmpData[grepl(searchThis,rownames(tmpData),ignore.case = T),] #Filter by keyword
    } else  { cat ("Not filtering by keywords \n")
    }
    
    ## Check that the filtered data is not empty (ie, that at least one word matched)
    if (nrow(tmpData)==0){
      cat ("None of the keywords were found \n")
      cat ("Not filtering by keywords \n")
      matrixData[[each]]
      tmpData <- matrixData[[each]] #Return data to the original, unfiltered
    }
    
    
    ### Filter for terms that are significant in at least 1
    cat ("Showing genes tha are significant in at least one term \n")
    tmpData <- tmpData[rowSums(tmpData <= 0.05)>=1,]
    
    ### Significant?
    signifData <- tmpData
    signifData[signifData<=0.05] <- "*"
    signifData[signifData>0.05] <- ""
    
    tmpData <- -log10(tmpData)
    size = 10
    pheatmap::pheatmap(tmpData,
                       display_numbers = signifData,number_color = "black",
                       fontsize = size/2.5,
                       cellwidth = size/2, cellheight = size/2, 
                       scale = "none",main = plotTitle)
    
    
  })
}



######## Graveyard of old functions:

######## Find which genes from each ########
## Calls the findOrtho_Sly2AGI from an accompanying script
# Need to add > source(orthologueFindingFunctions.R) to the script
findOrthologues_GOgenes <- function(significantGOtibble){
  significantGOtibble <- significantGOtibble %>%  # To this table
    left_join( # Join the following transformation:
      significantGOtibble %>% select(uniqueID,data) %>%
        # Group by category (maybe this isn't necessary)
        mutate(                  # Create a new column that is
          "Sly2AGI_orthologues" = map(data, .f= findOrtho_Sly2AGI)) %>% 
        select(uniqueID,3) 
      , by = c("uniqueID" = "uniqueID")
    ) %>% rename("DE_ITAGsinCats"="data")
  return(significantGOtibble)
}


gg_makeTileMap <- function(significantGOtibble,extraName="",keyword=NULL,filterbyOnto=c("MF","BP","CC")){
  require(gtools)
  plotTitle <- "Heatmap"
  if (extraName != "") plotTitle <- paste(plotTitle,extraName,sep=" - ")
  if (!is.null(keyword)) {
    searchThis <- paste(keyword,collapse = "|") 
    plotTitle <-paste(plotTitle, paste("keywords: ", paste(keyword,collapse=",") ,sep = ""), sep="\n")
  } else searchThis <- "*"
  ## Heatmap
  plotHeat <- significantGOtibble %>%  filter(grepl(searchThis,term,ignore.case = T)) %>%
    ## Do some filtering
    select(Set,term, category,ontology,2) %>% 
    mutate(Set=factor(Set,levels = mixedsort(unique(Set),decreasing = F))) %>%
    filter(ontology %in% filterbyOnto) %>% #If only one ontology term is wanted
    #### Dealing with NA's
    #mutate("term" = ifelse(is.na(term),category,term)) %>% # Either assign the GOcategory to a term that is NA,
    filter(!is.na(term)) %>%# or filter them out
    
    ## Group
    group_by(term,ontology,Set) %>%
    #gather(., key=contrast,value = pVal,over_represented_pvalue)
    ## In case we want significance points
    mutate("signif"=ifelse(over_represented_pvalue < 0.05,"*","")) %>% 
    ## Transform 
    mutate("log10pVal"=-log10(over_represented_pvalue)) %>%
    ### Make the heatmap
    ggplot(., aes(Set, term)) + 
    geom_tile(aes(fill=log10pVal),colour="white",width=1, height=1) + 
    #scale_fill_distiller(palette = "BuGn",direction = 1) +
    scale_fill_distiller(palette = "Blues",direction = 1,na.value="white",name="-log10(P)")+
    facet_grid(ontology~.,scales="free") +
    
    ### Add text
    #geom_text(aes(label=signif),size=2.5, hjust = 0.5) + 
    ## Names of axes and other graph aesthetics
    ggtitle(plotTitle) + xlab("Set") + ylab("GO terms") + 
    theme_gray() + 
    theme(#axis.title.x = element_blank(),
      axis.ticks = element_blank(), 
      panel.background=element_rect(fill="white", colour="lightgray"),
      panel.grid.major.x =  element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color="black", size=8),
      axis.text.x = element_text(angle = 15, hjust = 1),
      legend.position="bottom")
    #scale_alpha(guide = 'none') + 
    #guides(fill=guide_legend(title="-log10(pVal)"))
  return(plotHeat)
}

# ######## Make a wrapper to deal with the results ########
# wrapSignificantGOterms <- function(GOresults,GeneList,saveToFile=F,resultsPath="GenesPerCategory",extraName="",keyword=NULL,filterbyOnto=c("MF","BP","CC"),minFreq=1,findOrtho=F,go.goseq=go.goseq) {
#   require(purrrlyr)
#   ##### Filter by p-value of overrepresented (enriched) categories
#   listSignificantGOterms <- lapply(GOresults, function(x){  x[x$over_represented_pvalue < 0.05,] })
#   
#   
#   ### Calls a function to get the genes that were present in each category
#   getGenesPerCategory <- get_Genes_Per_Category(listSignificantGOterms,GeneList,go.goseq)
#   
#   # Make a table out of the lists
#   significantGOtibble <- bind_rows(getGenesPerCategory) %>%  # Bind lists into a single table
#     mutate(uniqueID = row_number())
#   
#   ### Get the orthologues of the genes
#   if (findOrtho) {
#   ### Find the Ath orthologues for each Sly gene
#   significantGOtibble <- findOrthologues_GOgenes(significantGOtibble) 
#   }
#   significantGOtibble
#   
#   
#   ####### Make pretty graphs ########
#   plotTitle <- "Frequency of Categories"
#   if (extraName != "") plotTitle <- paste(plotTitle,extraName,sep=" - ")
#   ## GO term frequency
#   
#   tmpCounts <- significantGOtibble %>% filter(!is.na(term)) %>% 
#     filter(ontology %in% filterbyOnto) %>% # If only selected ontology term are wanted
#     group_by(ontology,term) %>% summarise("Count"=n()) %>% arrange(ontology,desc(Count)) 
#   ## 
#   if (max(tmpCounts$Count)>=minFreq) tmpCounts <- tmpCounts %>% filter(Count >= minFreq) 
#   ## Make plot
#   plotByOnthology <- tmpCounts %>% 
#     ggplot(.,aes(x=reorder(term,Count), # Order X axis labels by Count
#                  y=Count, #Counts go in the Y axis
#                  fill=ontology)) + #Color by ontology
#     geom_bar(stat="identity") + coord_flip() + # Do a bar plot 
#     xlab("GO category") + xlab("Frequency") + labs(title=plotTitle)  + 
#     facet_grid(.~ontology,scales="free") #Separate by ontology
#   
#   plotTitle <- "Number of enriched GO terms per set"
#   if (extraName != "") plotTitle <- paste(plotTitle,extraName,sep="\n")
#   ## Terms per Set
#   plotOfCatNumber <- significantGOtibble %>% select(Set) %>% 
#     mutate(Set=factor(Set,levels = mixedsort(unique(Set),decreasing = F))) %>%
#     group_by(Set) %>% count() %>%
#     ggplot(.,aes(x=as.factor(Set),y=n)) + labs(title=plotTitle)  + 
#     geom_bar(stat="identity") + #,width = 10/length(unique(significantGOtibble$Set)) ) + # Do a bar plot 
#     xlab("Set") + ylab("# of significant GO terms") + theme_gray()
#   
#   
#   ## Heatmap
#   plotHeat <- gg_makeTileMap(significantGOtibble,extraName,keyword,filterbyOnto)
#   
#   ######## Save to files (or not) ########
#   if (saveToFile) { 
#     cat("---\nSaving to file\n---\n")
#     ##### Create directories to save outputs
#     if (extraName != "") resultsPath <- paste(resultsPath,extraName,sep = "-")
#     
#     itagPath <- paste0(resultsPath,"/ITAG_byCategory/")
#     orthoPath <- paste0(resultsPath,"/AthOrthologues_byCategory/")
#     dir.create(itagPath,recursive = T) 
#     dir.create(orthoPath,recursive = T)
#     
#     unite_(significantGOtibble,"tmpName",c("term","ontology","Set"),sep="_",remove=F) %>%  mutate("tmpName"=gsub("[\\>'()\\[\\\\/]","",tmpName)) %>%
#       ## Create new names for the files
#       mutate("itagFile"=paste0(itagPath,"/",gsub("/|[[:space:]]","-",tmpName),".txt")) %>%
#       by_row(~write.table(.$DE_ITAGsinCats, file = .$itagFile,sep = "\t",row.names = F,quote = F)) %>%
#       mutate("orthoFile"=paste0(orthoPath,"/",gsub("/|[[:space:]]","-",tmpName),".txt")) %>%
#       by_row(~write.table(.$Sly2AGI_orthologues, file = .$orthoFile,sep = "\t",row.names = F,quote = F))
#   } else cat("---\nNot saving to file\n---\n")
#   
#   return(list("significantGOtibble"=significantGOtibble,
#               "plotByOnthology"=plotByOnthology,
#               "plotOfCatNumber"=plotOfCatNumber,
#               "plotHeatmap"=plotHeat))
# }
# 
