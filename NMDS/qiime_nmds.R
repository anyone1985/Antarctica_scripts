# qiime_nmds.R
# Creates NMDS plots based upon NMDS formatted table file
# The outputs are intended to be created as PDFs to be edited in Illustrator manually to organise the labels optimally

library(vegan)
library(ggplot2)
library(ggrepel)

# Read data
setwd("/directory/nmds")
filename = "qiime.class.nmds"
qiime.dat=read.delim(filename, as.is=T, header=T, row.names = 1)
summary(qiime.dat)

# Convert to matrix
qiime.matrix = as.matrix(qiime.dat)

# Define parameters
qiime_NMDS=metaMDS(qiime.matrix, distance = 'bray', k = 2, trymax = 100) # k = number of dimensions, trymax = number of iterations to maximise score
stressplot(qiime_NMDS)  # Stress plot looks good
treat=c('Animal', 'Animal', 'Animal', 'Animal', 'Sewage', 'Sediment', 'Animal', 'Animal', 'Sediment', 'Animal', 'Animal', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Animal', 'Animal')   # Manually entered groupings for the sample types, corresponds to the rownames of the file

data.scores <- as.data.frame(scores(qiime_NMDS))
data.scores$site <- rownames(data.scores)
data.scores$grp <- treat
head(data.scores) 

species.scores <- as.data.frame(scores(qiime_NMDS, "species"))  
species.scores$species <- rownames(species.scores)
head(species.scores)

#### Create detailed plots

#### Plot 4: My adaption of plot3 subplot 2
grp.a <- data.scores[data.scores$grp == "Animal", ][chull(data.scores[data.scores$grp == "Animal", c("NMDS1", "NMDS2")]), ]
grp.b <- data.scores[data.scores$grp == "Sediment", ][chull(data.scores[data.scores$grp == "Sediment", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scores[data.scores$grp == "Sewage", ][chull(data.scores[data.scores$grp == "Sewage", c("NMDS1", "NMDS2")]), ]
hull.data <- rbind(grp.a, grp.b, grp.c)  #combine grp.a and grp.b
hull.data

# Pull out groups above a specified abundance level
abundantspecies <- qiime.dat[,colSums(qiime.dat) > 0.32]    # Arbitrary value: this presented the most abundant species without cluttering the plot
ncol(abundantspecies)
abundantspecies = colnames(abundantspecies)
abun.species.nmds.scores <- subset(species.scores, species.scores[,3] %in% abundantspecies)

graph <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp, colour=grp),alpha=0.30,size=0.5,linetype=1,col='black') + # add the convex hulls
  geom_text_repel(data=abun.species.nmds.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=1, size=4, segment.size = 1, box.padding = unit(0.5, 'lines'), point.padding = NA, nudge_y = ifelse(abun.species.nmds.scores$species == 'Alphaproteobacteria', 0.01, 0), nudge_x = ifelse(abun.species.nmds.scores$species == 'Coriobacteriia', 0.25, 0)) +  # add the species labels [if above 1% in any sample]
  geom_text_repel(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site, colour=grp),size=5, fontface="bold", segment.size = 1, segment.colour = "black", box.padding = unit(0.5, 'lines'), point.padding = NA, nudge_y = ifelse(data.scores$site == 'LS', 0.05, 0)) +  # add the site labels
  labs(y = 'NMDS2', x = 'NMDS1', title = "Bray-curtis dissimilarity NMDS from QIIME results") + # change labels
  scale_colour_manual(values=c("Animal" = "red2", "Sewage" = "sienna2", "Sediment" = "darkgreen"), guide = FALSE) +
  scale_fill_manual(values=c("Animal" = "red2", "Sewage" = "sienna2", "Sediment" = "darkgreen"), name = "Sample origin") +
  coord_equal() +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 15, colour = 'black'), # alter label aesthetics
        axis.title.x = element_text(size = 15, colour = 'black'), 
        axis.text.y = element_text(size = 15, colour = 'black'), 
        axis.title.y = element_text(size = 15, colour = 'black'),
        plot.title = element_text(size = 15, colour = 'black'),
        legend.text = element_text(size = 15, colour = 'black'),
        legend.title = element_text(size = 15, colour = 'black')
  )
graph

# cluster
clus <- hclust(vegdist(qiime.dat))
plot(clus)
rect.hclust(clus, h = 0.6)
