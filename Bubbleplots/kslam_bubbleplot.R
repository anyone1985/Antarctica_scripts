# kslam_bubbleplot.R
# R script to generate a bubbleplot using a NMDS formatted table as input

library(ggplot2)
library(reshape)

# Read data
setwd("/directory/nmds")   # Make this wherever the file indicated by filename below is located
filename = "kslam.class.nmds"
kslam.dat=read.delim(filename, as.is=T, header=T, row.names=1)
rownames(kslam.dat) = c('ARL', 'DW', 'GH1', 'LAR', 'OF', 'STP1', 'STP11A', 'STP12B', 'STP13A', 'STP15C', 'STP2', 'WEDSI')
summary(kslam.dat)

# Order alphabetically
kslam.dat <- kslam.dat[,order(names(kslam.dat))]
kslam.dat <- kslam.dat[rev(names(kslam.dat))]

# Remove low count features
kslam.dat <- kslam.dat[,colSums(kslam.dat) > 0.15]  # Number is found by trial and error until only 12 results remain

# Insert descriptive column for melt
kslam.dat <- data.frame(sites = c('ARL', 'DW', 'GH1', 'LAR', 'OF', 'STP1', 'STP11A', 'STP12B', 'STP13A', 'STP15C', 'STP2', 'WEDSI'), kslam.dat)    # Site names can be manually specified (like here) or could theoretically just be the row names

# Rearrange rows to match cluster dendro
target <- c('ARL', 'LAR', 'WEDSI', 'STP2', 'STP11A', 'STP13A', 'GH1', 'STP1', 'STP12B', 'STP15C', 'DW', 'OF') # This is entered manually based upon the order from left to right that the samples are grouped by the cluster dendrogram

# Melt into data frame
kslam.melt <- melt(kslam.dat)
kslam.melt$radius <- sqrt( kslam.melt$value / pi )
#kslam.matrix = as.matrix(kslam.dat)

# Plot
graph2 <-ggplot(kslam.melt, aes(x = sites, y = variable, label =value ))+
  geom_point(aes(size = radius*40,color=kslam.melt$value,fill=kslam.melt$value),shape=21)+
  scale_x_discrete(limits=target)+ # Reorder the values
  scale_size_identity()+
  xlab("Sites") + ylab("") + labs(title="K-SLAM abundance predictions to class or nearest level")+
  scale_fill_gradient(low = "grey", high="red", limits = c(0, 0.8))+
  scale_color_gradient(low = "grey",  high="red", limits = c(0, 0.8))+
  theme_bw()+
  theme(plot.title = element_text(size=17), axis.text = element_text(size=13, colour="black"), axis.text.x = element_text(angle=45, hjust = 1), axis.title = element_text(size=15, colour="black"))
graph2
