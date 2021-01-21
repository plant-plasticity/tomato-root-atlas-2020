##### Generic Genomic functions


### From a tibble of gene identifieres remove everything after a .
# Use this to remove information from genomic build (ie 
# SolycXXXX.3.2,SolycYYY.2.5 -> SolycXXXX,SolycYYYY )

removeDotGene <- function(dotGeneName,ColumnName="",verbose=F){
  if (verbose) cat ("removeDotGene function called\n")
  tmp <- gsub("\\.[0-9].*$","",unlist(dotGeneName))
  return(tmp)
}


aggregate_list_tables <- function(listDFs) {
  cat ("-- aggregate list of tables function called \n")
  # First make an empty table
  uniqRows <- Reduce(union,lapply(listDFs,rownames))
  uniqCols <- names(listDFs)
  naTable <- as.data.frame(matrix(NA,
                                  nrow = length(uniqRows),
                                  length(uniqCols)),
                           row.names =uniqRows)
  colnames(naTable) = uniqCols
  # Then Fill it
  for (each in names(listDFs)){
    cat ("Filling binary table:", each,"\n")
    naTable[rownames(listDFs[[each]]),
            each] <- listDFs[[each]]
  }
  #zeroTable <- apply(zeroTable,c(1,2),as.numeric)
  #zeroTable[is.na(zeroTable)] <- 0
  return(naTable)
  cat ("-- aggregate table done \n")
  cat ("\n")
}
