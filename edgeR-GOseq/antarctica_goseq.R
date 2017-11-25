############## START

# Goseq workaround [I had issues running it on the HPC cluster because I lack permissions for installing libraries permanently]
library(goseq, lib.loc="/tmp/")
library(GO.db, lib.loc="/tmp/")

# Go to directory
setwd("/directory/mapping_files")  # Make this wherever the file indicated by filename below is located

# Specify sample ordering
samples <- c("AR", "DW", "GH1", "LS", "OF", "STP1", "STP11A", "STP12B", "STP13A", "STP15C", "STP2", "WS")  # This is used for looping through the files. These must be the site name identifies used as part of the genematrix file naming scheme, and must be ordered correctly such that files that start with "AR.vs." occurs 11 times for UP and 11 times for DOWN, "DW.vs." occurs 10 times for UP/DOWN, etc...

# Load in annotations and gene lengths
  #>GO annotations
goannotfile = "antarctica_go_revised_mapping.out"   # Change file names as appropriate... this file is a key:value\n formatted text file where a value of 0 means no GO annotations associated with the key (i.e., sequence ID)
goannot=read.table(goannotfile, as.is=T, header=T, sep='\t')
goannot[goannot == 0] <- NA
  #>Generate GO mapping with gene2cat
splitgoannot = strsplit(goannot[,2], split='; ')
name.list = as.vector(goannot[,1])
names(splitgoannot) = name.list
rm(goannot) # tidy up

  #>PFAM annotations
pfamannotfile = "antarctica_pfam_revised_mapping.out"  # Change file names as appropriate... this file is a key:value\n formatted text file where a value of 0 means no pfam annotations associated with the key (i.e., sequence ID) [you get the picture, I won't paste this next to other files]
pfamannot=read.table(pfamannotfile, as.is=T, header=T, sep='\t')
pfamannot[pfamannot == 0] <- NA
  #>Generate PFAM mapping with gene2cat
splitpfamannot = strsplit(pfamannot[,2], split='; ')
name.list = as.vector(pfamannot[,1])
names(splitpfamannot) = name.list
rm(pfamannot) # tidy up

  #>FOAM annotations
foamannotfile = "antarctica_foam_revised_mapping.out"
foamannot=read.table(foamannotfile, as.is=T, header=T, sep='\t')
foamannot[foamannot == 0] <- NA
  #>Generate foam mapping with gene2cat
splitfoamannot = strsplit(foamannot[,2], split='; ')
name.list = as.vector(foamannot[,1])
names(splitfoamannot) = name.list
rm(foamannot) # tidy up

  #>Resfams annotations
resfamsannotfile = "antarctica_resfams_revised_mapping.out"
resfamsannot=read.table(resfamsannotfile, as.is=T, header=T, sep='\t')
resfamsannot[resfamsannot == 0] <- NA
  #>Generate Resfams mapping with gene2cat
splitresfamsannot = strsplit(resfamsannot[,2], split='; ')
name.list = as.vector(resfamsannot[,1])
names(splitresfamsannot) = name.list
rm(resfamsannot) # tidy up

  #>Lengths
lengthfile = "antarctica.revised.gene.lengths"      # The gene lengths file is a key:value\n formatted text file where the value == the length of the gene as a nucleotide, and the key == sequence ID
genelengths=read.table(lengthfile, as.is=T, header=T)
prevector = genelengths[2]
rownames(prevector) = genelengths[,1] 

length.vector = as.integer(prevector[,1])
names(length.vector) = rownames(prevector)
rm(prevector) # tidy up
rm(genelengths) # tidy up

##### Loop through each pairwise combination
N = 12
for(i in 1:(N-1)) {   # Note that goseq can be slow for some pairwise comparisons, so you can modify the looping system to create separate scripts to run in parallel
  for(x in (i+1):N) {
    word <- 'DOWN'
    setwd("/directory/edgeR") # Change wd to matrix file directory
    matrixname <- paste(samples[i],".vs.",samples[x],".",word,".genematrix.txt",sep="")
    #### Load in vector matrix
    matrixtable <- read.table(matrixname, as.is=T, header = F)
    prevector <- matrixtable[2]
    rownames(prevector) <- matrixtable[,1]
    ### Tidy up vector matrix
    gene.vector <- as.integer(prevector[,1])
    names(gene.vector) <- rownames(prevector)
    rm(prevector) # tidy up
    #### Run PWF for length bias
    pwf <- nullp(DEgenes=gene.vector, bias.data=length.vector)
    #### Run GOseq for pairwise up/down category
    go_output <- goseq(pwf, gene2cat = splitgoannot)
    pfam_output <- goseq(pwf, gene2cat = splitpfamannot)
    resfams_output <- goseq(pwf, gene2cat = splitresfamsannot)
    foam_output <- goseq(pwf, gene2cat = splitfoamannot)
    
    ### Process GO output
    setwd("/directory/go") # Change wd to output directory
    write.csv(go_output, file=paste(strsplit(matrixname, '.genematrix')[[1]][1],".goSeqResults.csv", sep=""))
    ### FDR correction
    enriched.GO = go_output$category[p.adjust(go_output$over_represented_pvalue,method = "BH")<.01]
    ### Get GO term descriptions recursively and produce output file
    sink(paste(strsplit(matrixname, '.genematrix')[[1]][1],".goSeqDescriptions.txt", sep=""))
    for(go in enriched.GO) {
      goterm <- go
      print(GOTERM[[goterm]])
      cat("-------------------------------\n")
    }
    sink()
    ### Process PFAM output
    setwd("/directory/pfam") # Change wd to output directory
    write.csv(pfam_output, file=paste(strsplit(matrixname, '.genematrix')[[1]][1],".pfamSeqResults.csv", sep=""))
    ### FDR correction
    enriched.GO = pfam_output$category[p.adjust(pfam_output$over_represented_pvalue,method = "BH")<.01]
    ### Produce output file of FDR corrected terms
    sink(paste(strsplit(matrixname, '.genematrix')[[1]][1],".pfamSeqFDRCorr.txt", sep=""))
    for(go in enriched.GO) {
      print(go)
    }
    sink()
    
    ### Process FOAM output
    setwd("/directory/foam") # Change wd to output directory
    write.csv(foam_output, file=paste(strsplit(matrixname, '.genematrix')[[1]][1],".foamSeqResults.csv", sep=""))
    ### FDR correction
    enriched.GO = foam_output$category[p.adjust(foam_output$over_represented_pvalue,method = "BH")<.01]
    ### Produce output file of FDR corrected terms
    sink(paste(strsplit(matrixname, '.genematrix')[[1]][1],".foamSeqFDRCorr.txt", sep=""))
    for(go in enriched.GO) {
      print(go)
    }
    sink()
    
    ### Process Resfams output
    setwd("/directory/resfams") # Change wd to output directory
    write.csv(resfams_output, file=paste(strsplit(matrixname, '.genematrix')[[1]][1],".resfamsResults.csv", sep=""))
    ### FDR correction
    enriched.GO = resfams_output$category[p.adjust(resfams_output$over_represented_pvalue,method = "BH")<.01]
    ### Produce output file of FDR corrected terms
    sink(paste(strsplit(matrixname, '.genematrix')[[1]][1],".resfamsFDRCorr.txt", sep=""))
    for(go in enriched.GO) {
      print(go)
    }
    sink()
  }
}