# anova_tukeys_stats.R
# This script will compute ANOVA and Tukey's HSD for input files

# Read data
setwd("/directory/alpha-div")   # Make this wherever the file indicated by filename(s) below are located

# ANOVA
filename = "kslam_alpha_anova.csv"
kslam.corr=read.csv(filename, as.is=T, header=T, row.names = 1)
fit <- aov(X1.D ~ type, data=kslam.corr)
summary(fit)  # P = 0.0517, not significant enough to perform Tukey's
#TukeyHSD(fit)

filename = "kslam_bowtie_anova.csv"
kslam.corr=read.csv(filename, as.is=T, header=T, row.names = 1)
fit <- aov(Percent~ type, data=kslam.corr)
summary(fit) # P = 0.00137, significant
TukeyHSD(fit)

filename = "kslam_swiss_anova.csv"
kslam.corr=read.csv(filename, as.is=T, header=T, row.names = 1)
fit <- aov(Percent~ type, data=kslam.corr)
summary(fit) # P = 0.000153, significant
TukeyHSD(fit)

filename = "kslam_kslam_anova.csv"
kslam.corr=read.csv(filename, as.is=T, header=T, row.names = 1)
fit <- aov(percent ~ type, data=kslam.corr)
summary(fit)  # P = 0.00154, significant
TukeyHSD(fit)