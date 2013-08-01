setwd("~/Documents/Neuroscience/jefferis_lab/LHNAnalysis")
require(RColorBrewer)
require(NMF)
require(gphys)
require(gplots)
require(glmnet)
require(MASS)
require(penalized)
require(scatterplot3d)
library(pvclust)
require(mclust)

# load AnalysisSuite stuff
source("~/projects/AnalysisSuite/R/Code/Startup.R")

# load local functions
# for dealing with ephys data
source("Avg_IgorImport.R")
source("PSTH_FUNC.R")
source("ClusteringFunctions.R")

# for analysis
source("infotheory.R")
source("decoding.R")
source("modeling.R")
source("encoding.R")
source("preprocessing.R")



# load all the data needed for subsequent analysis
source("loaddata.R")



# clustering
## Try clustering principal components

#pcs <- prcomp(lhn.mat.noblank)


# function to color leaves of odor dendrogram by odor type

#labels <- cutree(clust,5)
#heatmap.2(pcs$x[,1:3], scale='none', Colv=NA, trace='none', col=jet.colors, distfun=function(x) dist(x))
#scatterplot3d(pcs$x[,1:3], color=cut(labels,4, labels=c("cyan", "magenta", "black","red")))
# doesn't really work that well...

## Assign weights using LDA
##
## TODO Try to allow only positive weights
#x <- glmnet(pn, factor(c(rep(0, 55),rep(1,55))), family="binomial")
# split into two separate groups
#set0 <- c(1, rep(0, 109))
#set1 <- c(rep(0,55), rep(1,55))
#set2 <- c(rep(0,30), rep(1,30), rep(0, 50))
#set3 <- c(rep(1,40), rep(0,40), rep(1, 30))
set4 <- c(rep(1, 10), rep(0,100))
#g1 <- lda(t(pn), factor(set1))
#g2 <- lda(t(pn), factor(set2))
#g3 <- lda(t(pn), factor(set3))

#plot.lda.response.hist(make.model.responses(pn, set0, 100, TRUE), 30)
#plot.lda.response.hist(make.model.responses(pn, set1, 100, method="nonneg"), 10) # each cell can't discriminate...
#plot.lda.response.hist(make.model.responses(pn, set2, 100, method="nonneg"), 10)
#plot.lda.response.hist(make.model.responses(pn, set3, 100, method="nonneg"), 10)
#plot.lda.response.hist(make.model.responses(pn, rbinom(110, 1, 0.05), 100, method="lda"), 10, x.min=-0.5, x.max=1.5)

# find subsets of highly correlated odors
# show PN 
pn.factors <- prcomp(t(pn))$x
#pn.factor.clust <- hclust(dist(pn.factors[,1:10]), method="ward")
pn.cor.clust <- hclust(as.dist(1-cor(pn)))
num.pn.clusts <- 12
pn.factor.labels <- cutree(pn.cor.clust, num.pn.clusts)
par(mfrow=c(3,4))
for (i in 1:num.pn.clusts) {
    curr <- ifelse(pn.factor.labels == i, 1, 0)
#    g <- lda(t(pn), factor(curr))
    plot.model.response.hist(make.model.responses(pn, curr, 100, TRUE, method="lda"), 30, title=paste("Cluster", i))
}

### Do just for different odor categories
terpenes <- ifelse(rownames(pn) %in% all.odor.type[all.odor.type$Type == "terpene",1], 1, 0)
clust.acids <- c("Propionic acid","Isopentanoic acid", "Pentanoic acid","Butyric acid","Isobutyric acid","Acetic acid","Hexanoic acid", "Methanoic acid","Lactic acid")

acids <- ifelse(rownames(pn) %in% clust.acids, 1, 0)
# TODO Make this work...
#plot.model.response.hist(make.model.responses(t(pn), acids, 100, TRUE, method="lda"), 100)

#plot(prcomp(lhn.common)$x[,1:2],xlim=c(-80,80),ylim=c(-80,80),col='magenta')
#par(new=T)
#plot(prcomp(pn.common)$x[,1:2],xlim=c(-80,80),ylim=c(-80,80),col='cyan')

heatmap.2(orn.lhn.corr, col=jet.colors,scale='none', trace='none')
heatmap.2(pn.lhn.corr, col=jet.colors,scale='none', trace='none')

# LHN heatmap
lhn.cor <- cor(rates.nona.mat, use='complete.obs')
lhn.odors.cor <- cor(t(rates.nona.mat), use='complete.obs')
orn.common.cor <- cor(t(orn.common), use='complete.obs')
pn.common.cor <- cor(t(pn.common), use='complete.obs')
lhn.common.cor <- cor(t(lhn.common), use='complete.obs')
heatmap.2(lhn.cor, hclustfun=function(x) hclust(x, method='complete'), distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')
hist(lhn.cor, main="LHN Cross-correlation", xlab="Correlation Coefficient")
heatmap.2(lhn.odors.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')

# just heatmaps of common odors
heatmap.2(orn.common.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, Rowv=NA, Colv=NA, col=jet.colors, trace='none')
heatmap.2(pn.common.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, Rowv=NA, Colv=NA, col=jet.colors, trace='none')
heatmap.2(lhn.common.cor, distfun=function(x) as.dist(1-x), scale='none',Rowv=NA, Colv=NA, symm=T, col=jet.colors, trace='none')

# TODO factor analysis...



