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

# outliers attempt1 mahalanobis
md <- mahalanobis(asm.qc.nolabels, colMeans(asm.qc.nolabels), cov(asm.qc.nolabels))
plot(md)
outliers <- md<100
asm.qc.filters <- asm.qc.nolabels[outliers,]
md1 <- mahalanobis(asm.qc.filters, colMeans(asm.qc.filters), cov(asm.qc.filters))
plot(md1)

# outliers PCA 

# PCA
asm.pca <- princomp(asm.qc.filters, cor = T)
screeplot(asm.pca, type = "line")
asm.pca$loadings
summary(asm.pca)
biplot(asm.pca)

library(rgl)
plot3d(asm.pca$scores[,1:3]) # , col=asm.qc.filters$Hybrid)

# clustering
asm.pcac <- princomp(asm.qc[,3:36], cor=T)
asm.pkcl <- kmeans(asm.pcac$loadings, 2, 20)
plot3d(asm.pcac$scores[,1:3], col=asm.pkcl$cluster)

# parcoord
asm.prof <- as.matrix(asm.qc[,4:36])
dim(asm.prof)
apply(asm.prof,2,max)
plot(c(0,33), c(0,7500000), type="n")
for(k in (1:3012)){points(1:33, asm.prof[k,],type="l")}

asm.prof.s <- as.matrix(asm.qc[,4:36])
asm.s <- scale(asm.prof.s, center=T, scale=T)
apply(asm.s,2,max)
plot(c(0,33), c(0,33), type="n")
for(k in (1:3012)){points(1:33, asm.s[k,],type="l")}

#scaling down
train <- sample(3012, 150)
asm.s2 <- asm.s[train,]
apply(asm.s2,2,max)
plot(c(0,33), c(0,9), type="n")
for(k in (1:150)){points(1:33, asm.s2[k,],type="l")}

# removing VRFPA04 & contigs > 1000
asm.filter <- subset(asm.qc, Assembly != "VRFPA04" & Ncontigs<300)
asm.prof <- as.matrix(asm.filter[,4:36])
dim(asm.prof)

#scaling down
train <- sample(2814, 150)
asm.s <- scale(asm.prof, center=T, scale=T)
asm.s2 <- asm.s[train,]
apply(asm.s2,2,max)
plot(c(0,33), c(0,6), type="n")
for(k in (1:150)){points(1:33, asm.s2[k,],type="l")}
