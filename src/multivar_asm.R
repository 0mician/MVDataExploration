#######################
# Libraries and setup #
#######################

library(ggplot2)
library(corrplot)
library(pastecs)
library(plotrix)
library(rpart)
library(rgl)

setwd(dir = "/home/sid/Dev/MVDataExploration/data/")
asm <- read.csv(file = "quast_all_metrics_reduced.csv", header=TRUE)
str(asm); dim(asm)

# basic dataset statistics, correlations
stat.desc(asm[,4:36], basic=TRUE, desc=TRUE)
asm.corr <- cor(asm[,4:36])
corrplot(asm.corr, type="upper", order="hclust", tl.col="black", tl.srt=90, tl.cex=0.7)

# Profile plot on subset of data
asm.matrix <- as.matrix(asm[,4:36])
asm.matrix.std <- scale(asm.matrix, center=T, scale=T)
asm.melt <- melt(asm.matrix.std)
colnames(asm.melt) <- c("RowID", "Variable", "Value")

ggplot(asm.melt,aes(x=Variable,y=Value,group=RowID)) +
  geom_line(colour=I("blue"), alpha=0.1) +
  xlab("Variable") + ylab("Std Value") +
  theme(axis.text.x = element_text(angle=270, hjust=0, vjust=0.5), 
        axis.text=element_text(size=10), 
        text = element_text(size=14),
        plot.title = element_text(size=16))
ggsave("parcoord.pdf")

##########################################
# Focus on outliers detected in parcoord #
##########################################
ggplot(asm, aes(GenomeFraction, Nmisassemblies)) +
  geom_point(data=subset(asm, Assembly != "PA7" & Assembly != "VRFPA04"), colour=I("blue"),alpha=0.5) +
  geom_point(data=subset(asm, Assembly == "PA7"), colour=I("sienna1"),alpha=0.5) +
  geom_point(data=subset(asm, Assembly == "VRFPA04"), colour=I("cyan4"),alpha=0.5) +
  annotate("text", x = 90, y = 410, label = "VRFPA04", color="cyan4", size=5) +
  annotate("text", x = 25, y = 42, label = "PA7", color="sienna1", size=5) +
  ggtitle("Genome fraction vs. # misassemblies") +
  ylab("# Misassemblies") +
  xlab("Genome fraction (%)") +
  theme(axis.text=element_text(size=10),
        text = element_text(size=14),
        plot.title = element_text(size=16)) 
ggsave("scatterplot_gfvsmis.pdf")

# subsetting (removal of the 2 pipelines PA7 and VRFPA04, and Ncontigs above 500)
asm.filter <- subset(asm, Assembly != "PA7" & Assembly != "VRFPA04" & Ncontigs < 400)

ggplot(asm.filter, aes(NmismatchesPer100kbp, NindelsPer100kbp)) +
  geom_point(data=subset(asm.filter, StrainID != "9122" & StrainID != "9132"), colour=I("blue"),alpha=0.5) +
  geom_point(data=subset(asm.filter, StrainID == "9122"), colour=I("sienna1"),alpha=0.5) +
  geom_point(data=subset(asm.filter, StrainID == "9132"), colour=I("cyan4"),alpha=0.5) +
  annotate("text", x = 1900, y = 55, label = "Strain 9122", color="cyan4", size=5) +
  annotate("text", x = 1900, y = 58, label = "Strain 9132", color="sienna1", size=5) +
  ggtitle("Nindels vs. #mismatches per 100kbp ") +
  ylab("# mismatchesPer100kbp") +
  xlab("NindelsPer100kbp") +
  theme(axis.text=element_text(size=10),
        text = element_text(size=14),
        plot.title = element_text(size=16)) 
ggsave("scatterplot_indels.pdf")

# further outlier detection and subsetting using mahalanobis
asm.filter <- subset(asm.filter, StrainID != "9122" & StrainID != "9132" & Ncontigs < 400)
md <- mahalanobis(asm.filter[,4:36], colMeans(asm.filter[,4:36]), cov(asm.filter[,4:36]))
plot(md, main="Mahalanobis outlier detection")
identify(md)
asm.filter <- asm.filter[-c(865,254,871,1329,2080,2381,888,1303,866,864),]

################
# PCA analysis #
################
asm.pca <- princomp(asm.filter[,4:36], cor = T)
screeplot(asm.pca, main="scree plot of PCA")
asm.pca$loadings
summary(asm.pca)
plot(asm.pca$scores[,1:2], type="p", pch=19, cex=0.7, col=asm.filter$StrainID)
plot3d(asm.pca$scores[,1:3], col=asm.filter$StrainID)

# biplot
########
asm.matrix <- as.matrix(asm.filter[,4:36])
xm<-apply(asm.matrix,2,mean)
y<-sweep(asm.matrix,2,xm)
ss<-(t(y)%*%y)
s<-ss/(nrow(x)-1)
d<-(diag(ss))^(-1/2)
e<-diag(d,nrow=ncol(asm.matrix),ncol=ncol(asm.matrix))
z<-y%*%e
r<-t(z)%*%z
q<-svd(z)
gfd<-((q$d[1])+(q$d[2]))/sum(q$d)
gfz<-(((q$d[1])^2)+((q$d[2])^2))/sum((q$d)^2)
gfr<-(((q$d[1])^4)+((q$d[2])^4))/sum((q$d)^4)
l<-diag(q$d,nrow=ncol(asm.matrix),ncol=ncol(asm.matrix))
R.B<-q$u        #scores matrix
C.B<-q$v%*%l    #loadings
#possibility to stretch scores by a scale factor
scalefactor<-10
R.B<-q$u *scalefactor

par(mar=c(4,4,4,4),pty='s',oma=c(5,0,0,0),font=2)
plot(R.B[ ,1],R.B[ ,2],axes=F,xlim=c(-1,1),ylim=c(-1,1),xlab=' ',ylab=' ',cex=2.8, pch=".", col="grey")
mtext('First component',side=1,line=3,cex=.8)
mtext('Second component',side=2,line=3,cex=.8)
axis(1,at=c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1),cex=.8)
axis(2,at=c(-1,-.8,-.6,-.4,-.2,0,.2,.4,.6,.8,1),cex=.8)
box( )

text(C.B[,1]-.05,C.B[,2]+.05,as.character(dimnames(asm.matrix)[[2]]),cex=0.7,col="green4")
for (i in seq(1,nrow(C.B),by=1))
  arrows(0,0,C.B[i,1],C.B[i,2],col="green4")

#Draw circle unit
draw.circle(0,0,1,border='black')

results<-list('correlation matrix'=r,'column effects'=C.B,'row effects'=R.B)
cat('The goodness of fit for the correlation matrix is',gfr,'for the centered, standardized design matrix',gfz,'and for the Mahalanobis distances is',gfd,' ') 
results

par(mfrow = c(3,3))
plot(asm.pca$scores[,1]~asm.pca$scores[,2], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,2]~asm.pca$scores[,3], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,3]~asm.pca$scores[,4], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,4]~asm.pca$scores[,5], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,5]~asm.pca$scores[,6], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,6]~asm.pca$scores[,7], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,7]~asm.pca$scores[,8], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,8]~asm.pca$scores[,9], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
plot(asm.pca$scores[,9]~asm.pca$scores[,10], xlim = c(-10,+10), ylim = c(-10,+10) ,cex=0.5)
par(mfrow =c(1,1))

###################
# factor analysis #
###################
asm.cor <- cor(asm.filter[,4:36])
asm.cor.eigen <- eigen(asm.cor)

p <- asm.cor.eigen$vectors[,1:8]
d <- diag(sqrt(asm.cor.eigen$values[1:8]))
b <- p%*%d
rownames(b) <- names(asm.filter[0,4:36])
colnames(b) <- c("Factor1", "Factor2","Factor3", "Factor4", "Factor5", "Factor6", "Factor7", "Factor8")
b

# Variance explained
asm.cor.eigen$values[1:8]

# Final communality estimates
bbt <- b%*%t(b)
diag(bbt)
sum(diag(bbt))

asm.res.cor <- asm.cor - bbt
psi <- diag(diag(asm.res.cor))
res <- (asm.res.cor - psi)^2
overall.rms <- sqrt(sum(res)/(ncol(res)*(ncol(res)-1)))
overall.rms

###################
# Clustering tree #
###################
asm.sub <- sample(2751, 200)
asm.matrix <- as.matrix(asm.filter[asm.sub,4:36])
asm.clust <- hclust(dist(asm.matrix), method="average")
plclust(asm.clust)
asm.gp <- cutree(asm.clust, k=10)
clusplot(asm.matrix, asm.gp, stand=TRUE, labels=10)

#######
# MDS #
#######
library("vegan")
asm.sub <- sample(2751, 200)
asm.matrix <- as.matrix(asm.filter[asm.sub,4:36])
tree.asm <- metaMDS(asm.matrix)
plot(tree.asm, type="t")
