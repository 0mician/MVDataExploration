############################
#### Libraries and setup ###
############################

library(ggplot2)
library(corrplot)
library(reshape2)
library(plyr)
setwd(dir = "/home/sid/Dev/MVDataExploration/data/")
options(scipen=1000000)

# Dataset loading and basic info
asm.qc <- read.csv(file = "quast_all_metrics_no_na.csv", header=TRUE)
str(asm.qc)
asm.qc.nolabels <- asm.qc[,4:36] # subsetting no labels
asm.corr <- cor(asm.qc.nolabels)
corrplot(asm.corr, type="upper", order="hclust", tl.col="black", tl.srt=90, tl.cex=0.7)

# outliers
md <- mahalanobis(asm.qc.nolabels, colMeans(asm.qc.nolabels), cov(asm.qc.nolabels))
plot(md)
outliers <- md<100
asm.qc.filters <- asm.qc.nolabels[outliers,]
md1 <- mahalanobis(asm.qc.filters, colMeans(asm.qc.filters), cov(asm.qc.filters))
plot(md1)

# PCA
asm.pca <- princomp(asm.qc.filters, cor = T)
screeplot(asm.pca)
asm.pca$loadings
summary(asm.pca)
biplot(asm.pca)

library(rgl)
plot3d(asm.pca$scores[,1:3]) # , col=asm.qc.filters$Hybrid)
