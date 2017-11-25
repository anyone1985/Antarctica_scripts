# antarctica_goseqfdr_corr.R
# The antarctica_goseq.R file originally outputs terms on their own if their FDR < 0.01. This script was intended to read in the results files from this first script and used it to create a key:value output where each output term is associated with its FDR value to enable the signature_terms script to utilise FDR values as part of its process.

library(plyr)  # Can't remember if this was necessary or not...
# Go to directory

# Specify sample ordering
samples <- c("AR", "DW", "GH1", "LS", "OF", "STP1", "STP11A", "STP12B", "STP13A", "STP15C", "STP2", "WS")

##### Loop through each pairwise combination
N = 12
for(i in 1:(N-1)) {
  for(x in (i+1):N) {
    word <- 'DOWN'
    setwd("/directory/csv_results")
    matrixname <- paste(samples[i],".vs.",samples[x],".",word,".goSeqResults.csv",sep="")  # Modify the ".goSeqResults.csv" section to be applicable to the data type e.g., FOAM, Pfam, or Resfams.
    #### Load in vector matrix
    matrixtable <- read.csv(matrixname, as.is=T, header = T)
    ### Process GO output
    setwd("/directory/fdr_corrected")
    ### FDR correction
    GO_terms = matrixtable$category[p.adjust(matrixtable$over_represented_pvalue,method = "BH")<.01]
    FDR = p.adjust(matrixtable$over_represented_pvalue,method = "BH")
    tempList = list(category = GO_terms, over_represented_FDR = FDR[1:length(GO_terms)])
    outTable = data.frame(tempList)
    write.table(outTable, file=paste(strsplit(matrixname, '.goSeqResults.csv')[[1]][1],".fdrCorr.txt", sep=""), row.names=FALSE, sep = "\t", quote=FALSE) # As above, change file name
  }
}
