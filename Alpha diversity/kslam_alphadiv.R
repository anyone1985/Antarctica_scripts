# kslam_alphadiv.R
# This script will run various alpha diversity metrics from provided table
library(vegan)

# Read data
setwd("/directory/alpha-div")   # Make this wherever the file indicated by filename below is located
filename = "kslam.alphadiv.table"
kslam.dat=read.delim(filename, as.is=T, header=T, row.names = 1)
rownames(kslam.dat) = c('AR', 'DW', 'GH1', 'LS', 'OF', 'STP1', 'STP11A', 'STP12B', 'STP13A', 'STP15C', 'STP2', 'WS')
summary(kslam.dat)

# Convert to matrix
kslam.matrix = as.matrix(kslam.dat)

# Compute simpsons diversity (1/D)
simpsondiv <- diversity(kslam.matrix, index = "invsimpson")