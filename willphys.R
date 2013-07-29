#setwd("~/Documents/Neuroscience/jefferis_lab/shahar_data/RandIgor/LHNAnalysis")
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

#require(e1071)
#require(nnet)
#require(gptk)

source("Avg_IgorImport.R")
source("PSTH_FUNC.R")
source("ClusteringFunctions.R")
source("infotheory.R")
source("decoding.R")
source("modeling.R")
source("encoding.R")
source("preprocessing.R")


################################################################################
## General setup stuff -- define colors
jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))(256)
bin.colors <- colorRampPalette(c('white', 'black'))(100)
bin.rev.colors <- colorRampPalette(c('black', 'white'))(100)

################################################################################
## LOAD LHN SPIKE DATA
#badcells=c("nm20130329c0","nm20130206c1", "nm20130606c0")
#SpikesT <- list()
#for (x in names(Spikes)) {
#  if (x %in% badcells) {
#  } else {
#    SpikesT[[x]] <- as.repeatedTrain(Spikes[[x]][[1]])
#  }
#}
#save(SpikesT, file="spikeT_summer.rda")

load("spikeT_summer.rda")

## Analysis of LHN data
# compute average firing rates for each recorded LHN
mean.rates <- lapply(SpikesT, cell.mean.rates)

# get odors that overlap with Hallem and Carlson odors
allodors <- unique(unlist(sapply(mean.rates,names)))
goododors <- allodors[1:36]
common.odors <- list(E2Hex="E2Hexenal", GerAc="Geranylacetat", Prpyl="Propylacetate", IPenA="isopentyl acetate", Et3HB="Ethyl-3-hydroxybutyrate", MetSl="Methylsalicylat", PeEtA="Phenylethyl alcohol", EtHex="Ethylhexanoate", BeZal="benzaldehyde", bCitr="b-citronellol", `1HxOl`="1-Hexanol", Cdvrn="cadaverine", MtAct="Methylacetat", AcAcd="Acetic acid", PrpnA="Propionic acid", BtrAc="Butyric acid", Amnia="ammonium hydroxide", PAcHd="Phenylacetaldehyd")

# types are: none, terpenes, alcohols, ester, ketones, lactones, acids, aromatics, aldehydes, sulfur compounds, amines

lhn.odors.types <- data.frame(Odorant=c("OilBl", "E2Hex", "GerAc", "Prpyl", "IPenA", "Et3HB", "Nonnl", 
                                       "CiVAc", "MetSl", "HexAc", "PeEtA", "AceAc", "EtHex", "2PnAc", 
                                       "5OdMx", "BeZal", "bCitr", "1HxOl", "Frnsl", "WatBl", "Cdvrn", 
                                       "Sprmn", "Acoin", "MtAct", "AcAcd", "PrpnA", "BtrAc", "Amnia", 
                                       "Pyrdn", "PAcHd", "HCL36", "PAcAc", "Vingr", "Geosn", "VinGe", 
                                       "PEtAm"),
                              Type=c("none", "aldehyde", "ester", "ester", "ester", "ester", "aldehyde",
                                "ester", "aromatic", "ester", "alcohol", "ketone", "ester", "ester",
                                "ester", "aromatic", "terpene", "alcohol", "terpene", "none", "amine",
                                "amine", "ketone", "ester", "acid", "acid", "acid", "amine", "amine",
                                "aromatic", "acid", "acid", "acid", "alcohol", "acid", "amine"))
rownames(lhn.odors.types) <- lhn.odors.types$Odorant

# figure out toal number of usable trials
total.trials <- sum(unlist(sapply(SpikesT, function(x) {sapply(x, length)})))
total.usable <- length(SpikesT)*36*4

# make matrix of all usable trials
lhn.mat <- make.total.data.matrix(SpikesT, length(mean.rates), 36*4)
rownames(lhn.mat) <- names(SpikesT)

# get labels of each odor
lhn.labels <- factor(colnames(lhn.mat))
# make matrix for each cell if 
lhn.bin <- binarize.mat.per.cell(lhn.mat, lhn.labels)

# make per cell average firing rates
rates.mat <- make.rate.per.cell(mean.rates, goododors)

# make matrix of rates without NAs
rates.nona.mat <- rates.mat[,apply(rates.mat, 2, function(y) !any(is.na(y)))]

# finally, lhn is all the average firing rates without any NAs for just the odors that are common to all trials
lhn <- rates.nona.mat


################################################################################
## LOAD HALLEM AND CARLSON DATA
x.names <- c("X2a", "X7a", "X9a", "X10a", "X19a", "X22a", "X23a", "X33b", "X35a", "X43a", "X43b", "X47a", "X47b", "X49b", "X59b", "X65a", "X67a", "X67c", "X82a", "X85a", "X85b", "X85f", "X88a", "X98a")
orn.names <- c("Or2a", "Or7a", "Or9a", "Or10a", "Or19a", "Or22a", "Or23a", "Or33b", "Or35a", "Or43a", "Or43b", "Or47a", "Or47b", "Or49b", "Or59b", "Or65a", "Or67a", "Or67c", "Or82a", "Or85a", "Or85b", "Or85f", "Or88a", "Or98a")
spontaneous.rates <- c(8, 17, 3, 14, 29, 4, 9, 25, 17, 21, 2, 1, 47, 8, 2, 18, 11, 6, 16, 14, 13, 7, 26, 12) # ordered by OR; from table in supplement of paper

# fine common odors

val.and.ors <- read.csv("val_orn.csv", header=T)
rownames(val.and.ors) <- val.and.ors$Odorant
val <- val.and.ors[,c("Odorant", "Attraction.index")]
all.odor.type <- val.and.odors[, c("Odorant", "Type")]
rownames(all.odor.type) <- all.odor.type$Odorant
orn <- as.matrix(val.and.ors[,x.names])
colnames(orn) <- orn.names
rownames(orn) <- val.and.ors$Odorant
# add baseline firing rates and make all negative rates 0
orn <- apply(orn, 1, function(x) {x + spontaneous.rates})
orn[orn < 0] <- 0
orn <- t(orn)

# make pn rates with and without lateral inhibition (by model)
nl <- make.pn.rates(orn,m=0)
pn <- make.pn.rates(orn)

val.common <- val[val$Odorant %in% unlist(common.odors),2]
orn.common <- orn[rownames(orn) %in% unlist(common.odors),]
pn.common <- pn[rownames(pn) %in% unlist(common.odors),]
lhn.common <- rates.nona.mat[rownames(lhn) %in% unlist(names(common.odors)),]
rownames(lhn.common) <- unlist(common.odors)



################################################################################
# Plot general response properties of each cell type

# Using all odors
pdf("orn_odor_response.pdf")
heatmap.2(t(orn), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors,scale="none", trace="none", cexCol=0.3, sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN All Odors')
dev.off()
pdf("pn_odor_response.pdf")
heatmap.2(t(pn), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", cexCol=0.3, sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN All Odors')
dev.off()
pdf("lhn_odor_response.pdf")
heatmap.2(t(lhn), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", cexRow=0.3,sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN All Odors')
dev.off()


# using common odors
pdf("orn_odor_response_common.pdf")
heatmap.2(t(orn.common), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN Common Odors', Rowv=NA)
dev.off()
pdf("pn_odor_response_common.pdf")
heatmap.2(t(pn.common), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN Common Odors', Rowv=NA)
dev.off()
pdf("lhn_odor_response_common.pdf")
heatmap.2(t(lhn.common), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ),main='LHN Common Odors', Rowv=NA)
dev.off()

# plot all response
pdf("all_response_hist.pdf", height=2.5, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 3, 0))
hist(orn, 20, xlab='Firing Rate',  freq=F, main='', col='darkgreen', border=F)
hist(pn, 20, xlab='Firing Rate',  freq=F, main='', col='blue', border=F)
hist(lhn, 20, xlab='Firing Rate',  freq=F, main='', col='purple', border=F)
mtext("Responses to All Odors", outer = TRUE, cex = 1)
dev.off()


################################################################################
# Plot correlations and changes in correlations between receptors/cells over first three synapses

# using all odors
# cell by cell
pdf("orn_cell_corr.pdf")
heatmap.2(cor(orn), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN - OR x OR')
dev.off()
pdf("pn_cell_corr.pdf")
heatmap.2(cor(pn), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN - OR x OR', cex.main=3)
dev.off()
pdf("lhn_cell_corr.pdf")
heatmap.2(cor(lhn), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN - Cell X Cell', cex.main=3)
dev.off()

# histograms of correlations
pdf("cell_corr_hist.pdf", height=2, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0))

hist(cor(orn), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='darkgreen', border=F)
hist(cor(pn), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='blue', border=F)
hist(cor(lhn), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='purple', border=F)
dev.off()

# odor by odor
pdf("orn_odor_corr.pdf")
heatmap.2(cor(t(orn)), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN - Odor x Odor')
dev.off()
pdf("pn_odor_corr.pdf")
heatmap.2(cor(t(pn)), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN - Odor x Odor')
dev.off()
pdf("lhn_odor_corr.pdf")
heatmap.2(cor(t(lhn)), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN - Odor x Odor')
dev.off()

pdf("odor_corr_hist.pdf", height=2, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0))
hist(cor(t(orn)), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='darkgreen', border=F)
hist(cor(t(pn)), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='blue', border=F)
hist(cor(t(lhn)), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='purple', border=F)
dev.off()

# change in correlations over synapses
orn.pn.corr <- cor(orn.common, pn.common)
orn.lhn.corr <- cor(orn.common, lhn.common)
pn.lhn.corr <- cor(pn.common, lhn.common)

#heatmap.2(cor(orn.common), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN')

pdf("orn_pn_corr.pdf")
heatmap.2(orn.pn.corr, distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=F, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN -> PN')
dev.off()

#heatmap.2(cor(pn.common), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN')

pdf("pn_lhn_corr.pdf")
heatmap.2(pn.lhn.corr,  hclustfun=function(x) hclust(x, method="ward"), symm=F, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN -> LHN')
dev.off()

pdf("orn_lhn_corr.pdf")
heatmap.2(orn.lhn.corr,  hclustfun=function(x) hclust(x, method="ward"), symm=F, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN -> LHN')
dev.off()

#heatmap.2(cor(lhn.common), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN')
pdf("trans_corr.pdf", height=2.5, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0))
hist(orn.pn.corr, 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='gray', border=F)
hist(orn.lhn.corr, 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='gray', border=F)
hist(pn.lhn.corr, 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='gray', border=F)
dev.off()

## higher order correlations
heatmap.2(cor(pn.lhn.corr),  hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN -> LHN')

################################################################################
# Changes in normalization/decorrleation
pdf("norm_pca_plot.pdf", height=4.5, width=8.5, pointsize=10)
par(mfrow=c(2,3), mar=c(4,4,3,4))
s.orn <- apply(orn,1,sum)
s.nl <- apply(nl,1,sum)
s.pn <- apply(pn,1,sum)
s.lhn <- apply(rates.nona.mat[!(rownames(rates.nona.mat) %in% c("WatBl", "OilBl")),], 1, sum)

# Equalization histogram
barplot(s.orn/mean(s.orn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn), main="ORN", col="darkgreen", border=NA)
#barplot(s.nl/mean(s.nl), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn), main="NL",)
barplot(s.pn/mean(s.pn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(pn),main="PN", col="blue", border=NA)
barplot(s.lhn/mean(s.lhn), xlab="Odor", ylab="firing rates/avg", names.arg=1:length(s.lhn), main="LHN", col="purple", border=NA)

# Variance explained plot
pc.orn <- princomp(orn)$sdev
pc.nl <- princomp(nl)$sdev
pc.pn <- princomp(pn)$sdev
pc.lhn <- princomp(t(rates.nona.mat))$sdev
#par(mfrow=c(1,3))
plot(100*pc.orn^2/sum(pc.orn^2), xlab="PC component", ylab="% variance", ylim=c(0,75), col="darkgreen")
#plot(100*pc.nl^2/sum(pc.nl^2), xlab="PC component", ylab="% variance", ylim=c(0,75))
plot(100*pc.pn^2/sum(pc.pn^2), xlab="PC component", ylab="% variance",ylim=c(0,75), col="blue")
plot(100* pc.lhn ^2/sum(pc.lhn ^2), xlab="PC component", ylab="% variance", ylim=c(0,75), col="purple")
dev.off()

# clustering
## Try clustering principal components
orn.pcs <- prcomp(orn)
pn.pcs <- prcomp(pn)
lhn.pcs <- prcomp(lhn[!(rownames(lhn) %in% c("OilBl", "WatBl")),])
#pcs <- prcomp(lhn.mat.noblank)
type.colors <- c(amine="gold", lactone="darkblue", acid="pink", sulfur="black", terpene="lightgreen", aldehyde="gray", ketone="yellow", aromatic="lightblue", alcohol="red", ester="darkgreen")

# function to color leaves of odor dendrogram by odor type

color.leaves.all <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = type.colors[[as.character(all.odor.type[all.odor.type$Odor == a$label,]$Type)]], lab.cex=0.6))
  }
  n
}

color.leaves.lhn <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
   
    col <- type.colors[[as.character(lhn.odors.types[lhn.odors.types$Odorant == a$label,]$Type)]]
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=1))
  }
  n
}
#orn.clust <- dendrapply(as.dendrogram(hclust(dist(orn.pcs$x[,1:10]), method='ward')), color.leaves.all)
#pn.clust <- dendrapply(as.dendrogram(hclust(dist(pn.pcs$x[,1:10]), method='ward')), color.leaves.all)
#lhn.clust <- dendrapply(as.dendrogram(hclust(dist(lhn.pcs$x[,1:10]), method='ward')), color.leaves.lhn)

#orn.clust <- hclust(as.dist(1-cor(t(orn))), method='ward')
#pn.clust <- hclust(as.dist(1-cor(t(pn))), method='ward')
#lhn.clust <- hclust(as.dist(1-cor(t(lhn[!(rownames(rates.nona.mat) %in% c("OilBl", "WatBl")),]))), method='ward')

orn.clust <- hclust(dist(orn), method='ward')
pn.clust <- hclust(dist(pn), method='ward')
lhn.clust <- hclust(dist(lhn[!(rownames(rates.nona.mat) %in% c("OilBl", "WatBl")),]), method='ward')

orn.den <- dendrapply(as.dendrogram(orn.clust), color.leaves.all)
pn.den <- dendrapply(as.dendrogram(pn.clust), color.leaves.all)
lhn.den <- dendrapply(as.dendrogram(lhn.clust), color.leaves.lhn)

pdf("cluster_res.pdf", height=6, width=8.5)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0), mar=c(3,1,1,6))
plot(orn.den, horiz=T, main="ORN", center=T)
legend("topleft", legend=names(type.colors), fill=unlist(type.colors), border=unlist(type.colors), bty='n', cex=0.8)
plot(pn.den, horiz=T, main="PN", center=T)
plot(lhn.den, horiz=T, main="LHN", center=T)
dev.off()

num.clusts <- 10
orn.clust.labels <- cutree(orn.clust, num.clusts)
pn.clust.labels <- cutree(pn.clust, num.clusts)
lhn.clust.labels <- cutree(lhn.clust, 8)

orn.types.num <- as.numeric(factor(all.odor.type[names(orn.clust.labels),]$Type))
pn.types.num <- as.numeric(factor(all.odor.type[names(pn.clust.labels),]$Type))
lhn.types.num <-  as.numeric(factor(lhn.odors.types[names(lhn.clust.labels),]$Type))

per.class.purity <- function(labels, true.cats) {
  p <- c()
  for (i in unique(true.cats)) {
    curr <- labels[true.cats == i]
    most.freq <- as.integer(names(which.max(table(curr))))
    p <- c(p, sum(curr == most.freq)/length(curr))
    # fraction of pairs that are in the same cluster vs different clusters
  }
  return(p)
}

cluster.purity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

# numeric labels for odor class
#orn.ri <- adjustedRandIndex(orn.clust.labels, orn.types.num)
#pn.ri <-  adjustedRandIndex(pn.clust.labels, pn.types.num)
#lhn.ri <-  adjustedRandIndex(lhn.clust.labels, lhn.types.num)

orn.purity <- per.class.purity(orn.clust.labels,orn.types.num)*100
names(orn.purity) <- unique(all.odor.type$Type)
pn.purity <- per.class.purity(pn.clust.labels, pn.types.num)*100
names(pn.purity) <- unique(all.odor.type$Type)
lhn.purity <- per.class.purity(lhn.clust.labels, lhn.types.num)*100
names(lhn.purity) <- unique(factor(lhn.odors.types[names(lhn.clust.labels),]$Type))

pdf("cluster_purity.pdf", height=2.5, width=8.5)
par(mfrow=c(1,3))
barplot(orn.purity, col=sapply(names(orn.purity), function(x) type.colors[[x]]), las=3, ylab="Purity (%)")
barplot(pn.purity, col=sapply(names(pn.purity), function(x) type.colors[[x]]), las=3, ylab="Purity (%)")
barplot(lhn.purity, col=sapply(names(lhn.purity), function(x) type.colors[[x]]), las=3, ylab="Purity (%)")
dev.off()

#######################################################
# make per odor class mean counts

# computes t-statistic
compute.odor.type.score <- function(mat, types) {
  m <- matrix(0, nrow=length(unique(types$Type)), ncol=ncol(mat))
  rownames(m) <- unique(types$Type)
  colnames(m) <- colnames(mat)
  group.means <-  sapply(unique(types$Type), function(x) mean(mat[types$Type == x,]))
  group.sd <-  sapply(unique(types$Type), function(x) sd(mat[types$Type == x,])/length(x))
  for (i in 1:ncol(mat)) {
   m[,i] <- sapply(unique(types$Type), function(x) mean(mat[types$Type == x,i]))
#    m[,i] <- sapply(unique(types$Type), function(x) t.test(mat[types$Type == x,i], mat[,i])$p.value)
#    m[,i] <- (m[,i] - mean(mat[,i]))/(sd(mat[,i])/sqrt(length(mat[,i])))
    m[,i] <- (m[,i] - group.means)/(group.sd)
  }
  return(m)
}


orn.type.mean <- compute.odor.type.score(orn, all.odor.type)
pn.type.mean <-compute.odor.type.score(orn, all.odor.type)
lhn.type.mean <- compute.odor.type.score(lhn[!(rownames(lhn) %in% c("OilBl", "WatBl")),], lhn.odors.types[!(lhn.odors.types$Odorant %in% c("OilBl", "WatBl")),])


orn.cells.types <- apply(orn.type.mean, 2, function(x) all.odor.type$Type[which.max(x)])
pn.cells.types <- apply(pn.type.mean, 2, function(x) all.odor.type$Type[which.max(x)])
lhn.cells.types <- apply(lhn.type.mean, 2, function(x) lhn.odors.types[!(lhn.odors.types$Odorant %in% c("OilBl", "WatBl")),]$Type[which.max(x)])

orn.cells.clust <- hclust(as.dist(1-cor(orn)), method='ward')
pn.cells.clust <- hclust(as.dist(1-cor(pn)), method='ward')
lhn.cells.clust <- hclust(as.dist(1-cor(lhn[!(rownames(rates.nona.mat) %in% c("OilBl", "WatBl")),])), method='ward')

color.leaves.cell.orn <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    x <- orn.cells.types[[a$label]]
    col <- type.colors[[x]] 
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=0.5))
  }
  n
}
color.leaves.cell.pn <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    x <- pn.cells.types[[a$label]]
    col <- type.colors[[x]] 
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=0.5))
  }
  n
}
color.leaves.lhn.cell <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    x <- lhn.cells.types[[a$label]]
    col <- type.colors[[x]] 
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=0.5))
  }
  n
}

orn.cell.den <- dendrapply(as.dendrogram(orn.cells.clust), color.leaves.cell.orn)
pn.cell.den <- dendrapply(as.dendrogram(pn.cells.clust), color.leaves.cell.pn)
lhn.cell.den <- dendrapply(as.dendrogram(lhn.cells.clust), color.leaves.lhn.cell)

lhn.score.pcs <- prcomp(t(lhn.type.mean))$x
pn.score.pcs <- prcomp(t(pn.type.mean))$x
#par(mfrow=c(1,2))
pdf("cell_cluster_res.pdf", height=4, width=8.5)
par(oma = c(0, 0, 0, 0))
plot(lhn.cell.den, horiz=F, main="LHN", center=T)
legend("topleft", legend=names(type.colors), fill=unlist(type.colors), border=unlist(type.colors), bty='n', cex=0.8)
dev.off()

pdf("orn_odor_score.pdf")
heatmap.2(t(orn.type.mean), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN')
dev.off()

pdf("pn_odor_score.pdf")
heatmap.2(t(pn.type.mean), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN')
dev.off()

pdf("lhn_odor_score.pdf")
heatmap.2(t(lhn.type.mean), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN')
dev.off()

#labels <- cutree(clust,5)
#heatmap.2(pcs$x[,1:3], scale='none', Colv=NA, trace='none', col=jet.colors, distfun=function(x) dist(x))
#scatterplot3d(pcs$x[,1:3], color=cut(labels,4, labels=c("cyan", "magenta", "black","red")))
# doesn't really work that well...

## Assign weights using LDA
##
## TODO Try to allow only positive weights
#x <- glmnet(pn, factor(c(rep(0, 55),rep(1,55))), family="binomial")
# split into two separate groups
set0 <- c(1, rep(0, 109))
set1 <- c(rep(0,55), rep(1,55))
set2 <- c(rep(0,30), rep(1,30), rep(0, 50))
set3 <- c(rep(1,40), rep(0,40), rep(1, 30))
set4 <- c(rep(1, 10), rep(0,100))
g1 <- lda(t(pn), factor(set1))
g2 <- lda(t(pn), factor(set2))
g3 <- lda(t(pn), factor(set3))

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


###################
## Data Analysis ##
###################
  
lhn.trial.bin <- binarize.trials(lhn.mat, lhn.labels)
lhn.probs <- make.trial.probs(lhn.trial.bin)
lhn.mat.noblank <- lhn.mat[, !(colnames(lhn.mat) %in% c("OilBl", "WatBl"))]

# plot relative entropy distance matrix
D <- make.D(lhn.probs)
pdf("klmat.pdf")
heatmap.2(D, hclustfun = function(x) { hclust(as.dist(x), method="complete") }, scale='none', symm=T, trace='none', col=bin.rev.colors, main='Relative Entropy')
dev.off()
labels = cutree(hclust(as.dist(D), method="complete"),4)

# Show how clustering affects classification accuracy
success.rates <- sapply(1:nrow(lhn.trial.bin), function(x) mean(cross.validate(lhn.trial.bin, num.classes=x)))

pdf("classification_acc_and_mi.pdf", width=8.5, height=4)
par(mfrow=c(1,2))
plot(success.rates, xlab="Number of Clusters", ylab="Mean % Odors Correctly Identified", type='l', main='Classification Accuracy')

# normalized mutual information scores
#MI.scores <- sapply(1:nrow(lhn.trial.bin), function(x) cluster.mutual.info(lhn.trial.bin,x))/mean(mutual.info(lhn.probs))
plot(MI.scores, ylab="Normalized Information About Stimulus", xlab="Number of Clusters", type='l', main='Mutual Information')
dev.off()

# plot confusion matrices
#confusion <- single.neuron.confusion(lhn.trial.bin)
# TODO Calculate mutual information and show redundancy

#leave.out <- seq(1,ncol(lhn.mat.noblank),4)
#train.lhn.rates <- make.trial.rates(lhn.mat.noblank[, !(1:ncol(lhn.mat.noblank) %in% leave.out)])
#test.lhn.rates <- lhn.mat.noblank[,leave.out]
#train.lhn.bin <- make.trial.probs(lhn.trial.bin[,!(1:ncol(lhn.trial.bin) %in% leave.out)], 3)
#test.lhn.bin <- lhn.trial.bin[,leave.out]


# use single-linkage clustering to get minimum info dists

## Look for combinations of PN factors that predict LHN
# use first 5 factors of odor responses and try to predict?
#pn.factors <- nmf(pn,5)
#odor.factors <- fit(pn.factors)@H

###
# valence stuff

#orn.val <- sapply(summary(cor(orn.common ~ val.common)), function(x) x$r.squared)
#pn.val <- sapply(summary(lm(pn.common ~ val.common)), function(x) x$r.squared)
#lhn.val <- sapply(summary(lm(lhn.common ~ val.common)), function(x) x$r.squared)
orn.val.all <- cor(orn, val$Attraction.index)
pn.val.all <- cor(pn, val$Attraction.index)
orn.val <- cor(orn.common, val.common)
pn.val <- cor(pn.common, val.common)
lhn.val <- cor(lhn.common, val.common)

# Or65a is highly correlated with valence (.321 in ORN, .45 in PN)
# also Or33b
# Or65a is thought to mediate aggression
par(mfrow=c(1,3))
hist(orn.val, main="ORN vs valence", xlab="Correlation Coefficient", xlim=c(-1,1))
hist(pn.val, main="PN vs valence", xlab="Correlation Coefficient", xlim=c(-1,1))
hist(lhn.val, main="LHN vs valence", xlab="Correlation Coefficient", xlim=c(-1,1))

