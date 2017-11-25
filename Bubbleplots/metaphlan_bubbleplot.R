# metaphlan_bubbleplot.R
# R script to generate a bubbleplot using a NMDS formatted table as input

library(ggplot2)
library(reshape)

# Read data
setwd("C:/Users/lythl/Desktop/Honours 2016/Project_work/functional_analysis/metaphlan2/bubbleplot")
filename = "metaphlan2.class.nmds"
metaphlan.dat=read.delim(filename, as.is=T, header=T, row.names=1)
rownames(metaphlan.dat) = c('ARL', 'DW', 'GH1', 'LAR1', 'OF', 'STP1', 'STP11A', 'STP12B', 'STP13A', 'STP15C', 'STP2', 'WEDSI')
summary(metaphlan.dat)

# Order alphabetically
metaphlan.dat <- metaphlan.dat[,order(names(metaphlan.dat))]
metaphlan.dat <- metaphlan.dat[rev(names(metaphlan.dat))]

# Remove low count features
metaphlan.dat <- metaphlan.dat[,colSums(metaphlan.dat) > 0.1]  # Number is found by trial and error until only 12 results remain

# Insert descriptive column for melt
metaphlan.dat <- data.frame(sites = c('ARL', 'DW', 'GH1', 'LAR1', 'OF', 'STP1', 'STP11A', 'STP12B', 'STP13A', 'STP15C', 'STP2', 'WEDSI'), metaphlan.dat)

# Rearrange rows to match cluster dendro
target <- c('OF', 'ARL', 'LAR1', 'WEDSI', 'DW', 'STP11A', 'STP13A', 'STP12B', 'STP15C', 'STP2',  'GH1', 'STP1')  # This is entered manually based upon the order from left to right that the samples are grouped by the cluster dendrogram

# Melt into data frame
metaphlan.melt <- melt(metaphlan.dat)
metaphlan.melt$radius <- sqrt( metaphlan.melt$value / pi )

# Plot
graph2 <-ggplot(metaphlan.melt, aes(x = sites, y = variable, label =value ))+
  geom_point(aes(size = radius*40,color=metaphlan.melt$value,fill=metaphlan.melt$value),shape=21)+
  scale_x_discrete(limits=target)+ # Reorder the values
  scale_size_identity()+
  xlab("Sites") + ylab("") + labs(title="MetaPhlAn2 abundance predictions to class or nearest level")+
  scale_fill_gradient(low = "grey", high="red", limits = c(0, 1))+
  scale_color_gradient(low = "grey",  high="red", limits = c(0, 1))+
  theme_bw()+
  theme(plot.title = element_text(size=17), axis.text = element_text(size=13, colour="black"), axis.text.x = element_text(angle=45, hjust = 1), axis.title = element_text(size=15, colour="black"))
graph2
