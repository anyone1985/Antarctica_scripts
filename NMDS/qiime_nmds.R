# Actual script

library(vegan)
library(ggplot2)
library(ggrepel)
#library(TeachingDemos)

# Read data
setwd("C:/Users/lythl/Desktop/Honours 2016/Project_work/functional_analysis/k-slam/kslam2krona/nmds/qiime_genus")
filename = "qiime.class.nmds"
kslam.dat=read.delim(filename, as.is=T, header=T, row.names = 1)
#rownames(kslam.dat) = c('AR', 'DW', 'GH1', 'LS', 'OF', 'STP1', 'STP2', 'STP11A', 'STP12B', 'STP13A', 'STP15C', 'WS')
summary(kslam.dat)

# Convert to matrix
kslam.matrix = as.matrix(kslam.dat)
#kslam.matrix = as.matrix(sapply(kslam.dat, as.numeric))
#rownames(kslam.matrix) <- kslam.dat$X

# Define parameters
kslam_NMDS=metaMDS(kslam.matrix, distance = 'bray', k = 2, trymax = 100) # k = number of dimensions, trymax = number of iterations to maximise score
stressplot(kslam_NMDS)  # Stress plot looks good
treat=c('Animal', 'Animal', 'Animal', 'Animal', 'Sewage', 'Sediment', 'Animal', 'Animal', 'Sediment', 'Animal', 'Animal', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Sediment', 'Animal', 'Animal')

data.scores <- as.data.frame(scores(kslam_NMDS))
data.scores$site <- rownames(data.scores)
data.scores$grp <- treat
head(data.scores) 

species.scores <- as.data.frame(scores(kslam_NMDS, "species"))  
species.scores$species <- rownames(species.scores)
head(species.scores)

#### Create detailed plots

#### Plot 4: My adaption of plot3 subplot 2
grp.a <- data.scores[data.scores$grp == "Animal", ][chull(data.scores[data.scores$grp == "Animal", c("NMDS1", "NMDS2")]), ]
grp.b <- data.scores[data.scores$grp == "Sediment", ][chull(data.scores[data.scores$grp == "Sediment", c("NMDS1", "NMDS2")]), ]
grp.c <- data.scores[data.scores$grp == "Sewage", ][chull(data.scores[data.scores$grp == "Sewage", c("NMDS1", "NMDS2")]), ]
#hull.data <- rbind(grp.a, grp.b)  #combine grp.a and grp.b
hull.data <- rbind(grp.a, grp.b, grp.c)  #combine grp.a and grp.b
hull.data

# Pull out groups above a specified abundance level
abundantspecies <- kslam.dat[,colSums(kslam.dat) > 0.32]
ncol(abundantspecies)
abundantspecies = colnames(abundantspecies)
abun.species.nmds.scores <- subset(species.scores, species.scores[,3] %in% abundantspecies)

graph <- ggplot() + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=grp,group=grp, colour=grp),alpha=0.30,size=0.5,linetype=1,col='black') + # add the convex hulls
  geom_text_repel(data=abun.species.nmds.scores,aes(x=NMDS1,y=NMDS2,label=species),alpha=1, size=4, segment.size = 1, box.padding = unit(0.5, 'lines'), point.padding = NA, nudge_y = ifelse(abun.species.nmds.scores$species == 'Alphaproteobacteria', 0.01, 0), nudge_x = ifelse(abun.species.nmds.scores$species == 'Coriobacteriia', 0.25, 0)) +  # add the species labels [if above 1% in any sample]
  #geom_text(data=data.scores,aes(x=NMDS1,y=NMDS2,label=site, colour=grp),size=6,vjust=0, fontface="bold") +  # add the site labels
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
clus <- hclust(vegdist(kslam.dat))
plot(clus)
rect.hclust(clus, h = 0.6)
