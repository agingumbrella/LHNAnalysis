################################################################################
# Plot general response properties of each cell type

# Using all odors
pdf("figs/orn_odor_response.pdf")
heatmap.2(t(orn), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors,scale="none", trace="none", cexCol=0.3, sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN All Odors')
dev.off()
pdf("figs/pn_odor_response.pdf")
heatmap.2(t(pn), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", cexCol=0.3, sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN All Odors')
dev.off()
pdf("figs/lhn_odor_response.pdf")
heatmap.2(t(lhn), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", cexRow=0.3,sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN All Odors')
#pheatmap(lhn, col=jet.colors, clustering_method="ward", fontsize_row=5, fontsize_col=6, border_col=NA, scale='none', annotation_colors =  sapply(lhn.odors.types$Type, function(x) type.colors[[x]]))
dev.off()


# using common odors
pdf("figs/orn_odor_response_common.pdf")
heatmap.2(t(orn.common), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN Common Odors', Rowv=NA)
dev.off()

pdf("figs/pn_odor_response_common.pdf")
heatmap.2(t(pn.common), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none',lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN Common Odors', Rowv=NA)
dev.off()

pdf("figs/lhn_odor_response_common.pdf")
heatmap.2(t(lhn.common), hclustfun=function(x) hclust(x, method="ward"), col=jet.colors, scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ),main='LHN Common Odors', Rowv=NA)
dev.off()

# plot all response
pdf("figs/all_response_hist.pdf", height=2.5, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 3, 0))
hist(orn, 20, xlab='Firing Rate',  freq=F, main='', col='darkgreen', border=F)
hist(pn, 20, xlab='Firing Rate',  freq=F, main='', col='blue', border=F)
hist(lhn, 20, xlab='Firing Rate',  freq=F, main='', col='purple', border=F)
mtext("Responses to All Odors", outer = TRUE, cex = 1)
dev.off()

lhn.delta.positive <- lhn.delta
lhn.delta.positive[lhn.delta < 0] <- 0
lifetime.sparseness(lhn.delta.positive)
pdf("figs/sparseness_hist.pdf", height=2.5, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 3, 0))
hist(lifetime.sparseness(orn), 10, xlab='Lifetime Sparseness', xlim=c(0,1), freq=F, main='', col='darkgreen', border=F)
hist(lifetime.sparseness(pn), 10, xlab='Lifetime Sparseness',  xlim=c(0,1), freq=F, main='', col='blue', border=F)
hist(lifetime.sparseness(lhn.delta.positive), 10, xlim=c(0,1), xlab='Lifetime Sparseness',  freq=F, main='', col='purple', border=F)
mtext("Sparseness of All Cells", outer = TRUE, cex = 1) 
dev.off()

################################################################################
# Plot correlations and changes in correlations between receptors/cells over first three synapses

# using all odors
# cell by cell
pdf("figs/orn_cell_corr.pdf")
heatmap.2(cor(orn), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN - OR x OR')
dev.off()

pdf("figs/pn_cell_corr.pdf")
heatmap.2(cor(pn), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN - OR x OR', cex.main=3)
dev.off()

pdf("figs/lhn_cell_corr.pdf")
heatmap.2(cor(lhn), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN - Cell X Cell', cex.main=3)
dev.off()

# histograms of cell type correlations
pdf("figs/cell_corr_hist.pdf", height=2, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0))
hist(cor(orn), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='darkgreen', border=F)
hist(cor(pn), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='blue', border=F)
hist(cor(lhn), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='purple', border=F)
dev.off()

# odor by odor
pdf("figs/orn_odor_corr.pdf")
heatmap.2(cor(t(orn)), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN - Odor x Odor')
dev.off()
pdf("figs/pn_odor_corr.pdf")
heatmap.2(cor(t(pn)), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN - Odor x Odor')
dev.off()
pdf("figs/lhn_odor_corr.pdf")
heatmap.2(cor(t(lhn)), distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN - Odor x Odor')
dev.off()

# histograms of odor correlations
pdf("figs/odor_corr_hist.pdf", height=2, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0))
hist(cor(t(orn)), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='darkgreen', border=F)
hist(cor(t(pn)), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='blue', border=F)
hist(cor(t(lhn)), 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='purple', border=F)
dev.off()

# change in correlations over synapses
orn.pn.corr <- cor(orn.common, pn.common)
orn.lhn.corr <- cor(orn.common, lhn.common)
pn.lhn.corr <- cor(pn.common, lhn.common)

pdf("figs/orn_pn_corr.pdf")
heatmap.2(orn.pn.corr, distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=F, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', Colv=NA, Rowv=NA, dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN -> PN')
dev.off()

pdf("figs/pn_lhn_corr.pdf")
heatmap.2(pn.lhn.corr,  hclustfun=function(x) hclust(x, method="ward"), symm=F, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN -> LHN')
dev.off()

pdf("figs/orn_lhn_corr.pdf")
heatmap.2(orn.lhn.corr,  hclustfun=function(x) hclust(x, method="ward"), symm=F, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN -> LHN')
dev.off()

pdf("figs/trans_corr.pdf", height=2.5, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0))
hist(orn.pn.corr, 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='gray', border=F)
hist(orn.lhn.corr, 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='gray', border=F)
hist(pn.lhn.corr, 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='gray', border=F)
dev.off()

## higher order correlations
#heatmap.2(cor(pn.lhn.corr),  hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN -> LHN')

################################################################################
# Plot changes in normalization/decorrelation and PCA analysis

# compute sums
s.orn <- apply(orn,1,sum)
s.nl <- apply(nl,1,sum)
s.pn <- apply(pn,1,sum)
s.lhn <- apply(lhn[!(rownames(lhn) %in% c("WatBl", "OilBl")),], 1, sum)

# compute eigenvalues
pc.orn <- prcomp(orn)$sdev
pc.nl <- prcomp(nl)$sdev
pc.pn <- prcomp(pn)$sdev
pc.lhn <- prcomp(lhn)$sdev

## set up figure
pdf("figs/norm_pca_plot.pdf", height=4.5, width=8.5, pointsize=10)
par(mfrow=c(2,3), mar=c(4,4,3,4))

# Equalization histogram
barplot(s.orn/mean(s.orn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn), main="ORN", col="darkgreen", border=NA)
#barplot(s.nl/mean(s.nl), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn), main="NL",)
barplot(s.pn/mean(s.pn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(pn),main="PN", col="blue", border=NA)
barplot(s.lhn/mean(s.lhn), xlab="Odor", ylab="firing rates/avg", names.arg=1:length(s.lhn), main="LHN", col="purple", border=NA)


# Variance explained plot by PCA
plot(100*pc.orn^2/sum(pc.orn^2), xlab="PC component", ylab="% variance", ylim=c(0,75), col="darkgreen", pch=20)
plot(100*pc.pn^2/sum(pc.pn^2), xlab="PC component", ylab="% variance",ylim=c(0,75), col="blue", pch=20)
plot(100* pc.lhn^2/sum(pc.lhn ^2), xlab="PC component", ylab="% variance", ylim=c(0,75), col="purple", pch=20)
dev.off()

pdf("figs/pca_project.pdf", width=8.5, height=2.5, pointsize=10)
par(mfrow=c(1,3), mar=c(4,4,3,4))
pc.common.orn <- prcomp(orn)
pc.common.pn <- prcomp(pn)
pc.common.lhn <- prcomp(lhn)

plot(pc.common.orn$x[,1:2], col='darkgreen', pch=20)
plot(pc.common.pn$x[,1:2], col='blue', pch=20)
plot(pc.common.lhn$x[,1:2], col='purple', pch=20)
dev.off()

