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

setwd("~/Documents/Neuroscience/jefferis_lab/LHNAnalysis")
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
source("util.R")

# load all the data needed for subsequent analysis
source("loaddata.R")
source("infotheory.R")
source("modeling.R")

source("willanatomy.R")
source("figscripts/clustering.R")
source("figscripts/modeling.R")
source("figscripts/prediction.R")
source("figscripts/responses.R")
source("figscripts/valence.R")

# clustering
## Try clustering principal components

#pcs <- prcomp(lhn.mat.noblank)


# function to color leaves of odor dendrogram by odor type

#labels <- cutree(clust,5)
#heatmap.2(pcs$x[,1:3], scale='none', Colv=NA, trace='none', col=jet.colors, distfun=function(x) dist(x))
#scatterplot3d(pcs$x[,1:3], color=cut(labels,4, labels=c("cyan", "magenta", "black","red")))
# doesn't really work that well...

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



