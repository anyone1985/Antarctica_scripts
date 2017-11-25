# qiime_bubbleplot.R
# R script to generate a bubbleplot using a NMDS formatted table as input

library(ggplot2)
library(reshape)

# Read data
setwd("/directory/nmds")   # Make this wherever the file indicated by filename below is located
filename = "qiime.class.nmds"
qiime.dat=read.delim(filename, as.is=T, header=T, row.names=1)
rownames(qiime.dat) = c('ARL', 'DW', 'GH1', 'LAR', 'OF', 'STP1', 'STP11A', 'STP12B', 'STP13A', 'STP15C', 'STP2', 'WEDSI')
summary(qiime.dat)

# Order alphabetically
qiime.dat <- qiime.dat[,order(names(qiime.dat))]
qiime.dat <- qiime.dat[rev(names(qiime.dat))]

# Remove low count features
qiime.dat <- qiime.dat[,colSums(qiime.dat) > 0.32]  # Number is found by trial and error until only 12 results remain

# Insert descriptive column for melt
names <- toupper(row.names(qiime.dat))
qiime.dat <- data.frame(sites = names, qiime.dat)

# Rearrange rows to match cluster dendro
target <- c('ARL','AZOL','LAR1','SKUARL','WEDSILS','LAR2','WEDRL','ABI','AGI','SKUAZOL','DW','OF1','GH1','STP15B','STP15C','STP0','STP10B','STP6','STP14B','STP10A','STP12B','STP16','STP10C','STP11B','STP13B','STP14A','STP4','STP10D','STP12A','STP17','STP1','STP2','STP11A','STP14','STP13A','STP5','STP7','STP9','STP15A') # This is entered manually based upon the order from left to right that the samples are grouped by the cluster dendrogram

# Melt into data frame
qiime.melt <- melt(qiime.dat)
qiime.melt$radius <- sqrt( qiime.melt$value / pi )

# Plot
graph2 <-ggplot(qiime.melt, aes(x = sites, y = variable, label =value ))+
  geom_point(aes(size = radius*40,color=qiime.melt$value,fill=qiime.melt$value),shape=21)+
  scale_x_discrete(limits=target)+ # Reorder the values
  scale_size_identity()+
  xlab("Sites") + ylab("") + labs(title="QIIME abundance predictions to class or nearest level")+
  scale_fill_gradient(low = "grey", high="red")+
  scale_color_gradient(low = "grey",  high="red")+
  theme_bw()+
  theme(plot.title = element_text(size=17), axis.text = element_text(size=13, colour="black"), axis.text.x = element_text(angle=45, hjust = 1), axis.title = element_text(size=15, colour="black"))
graph2
