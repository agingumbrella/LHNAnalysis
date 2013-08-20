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
hist(orn, 20, xlab='Response Magnitude (spikes/s)',  freq=F, main='', col='darkgreen', border=F)
hist(pn, 20, xlab='Response Magnitude (spikes/s)',  freq=F, main='', col='blue', border=F)
hist(lhn, 20, xlab='Response Magnitude (spikes/s)',  freq=F, main='', col='purple', border=F)
mtext("Responses to All Odors", outer = TRUE, cex = 1)
dev.off()

pdf("figs/eq_hist.pdf", height=2.25, width=8.5, pointsize=10)
par(mfrow=c(1,3), mar=c(4,4,3,4))

# Equalization histogram
barplot(s.orn/mean(s.orn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn), main="ORN", col="darkgreen", border=NA)
#barplot(s.nl/mean(s.nl), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn), main="NL",)
barplot(s.pn/mean(s.pn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(pn),main="PN", col="blue", border=NA)
barplot(s.lhn/mean(s.lhn), xlab="Odor", ylab="firing rates/avg", names.arg=1:length(s.lhn), main="LHN", col="purple", border=NA)
dev.off()


################################################################################
# Plot correlations and changes in correlations between receptors/cells over first three synapses

# using all odors
# cell by cell
pdf("figs/orn_cell_corr.pdf")
cor.orn <- cor(orn)
heatmap.2(cor.orn, distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='ORN - OR x OR')
dev.off()

pdf("figs/pn_cell_corr.pdf")
cor.pn <- cor(pn)
heatmap.2(cor.pn, distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='PN - OR x OR', cex.main=3)
dev.off()

pdf("figs/lhn_cell_corr.pdf")
cor.lhn <- cor(lhn)
heatmap.2(cor.lhn, distfun=function(x) as.dist(1-x), hclustfun=function(x) hclust(x, method="ward"), symm=T, col=jet.colors,scale="none", trace="none", sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ), main='LHN - Cell X Cell', cex.main=3)
dev.off()

# histograms of cell type correlations
pdf("figs/cell_corr_hist.pdf", height=2, width=8.5, pointsize=10)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0))
hist(cor.orn[lower.tri(cor.orn)], 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='darkgreen', border=F)
hist(cor.pn[lower.tri(cor.pn)], 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='blue', border=F)
hist(cor.lhn[lower.tri(cor.lhn)], 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='purple', border=F)
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
hist(cor(t(orn))[lower.tri(cor(t(orn)))], 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='darkgreen', border=F)
hist(cor(t(pn))[lower.tri(cor(t(pn)))], 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='blue', border=F)
hist(cor(t(lhn))[lower.tri(cor(t(lhn)))], 20, xlab='Correlation Coefficient', xlim=c(-1, 1), freq=F, main='', col='purple', border=F)
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
## Plot per-cell and per-odor sparseness scores
lhn.delta.positive <- lhn.delta.rates
lhn.delta.positive[lhn.delta.rates < 0] <- 0

# plot lifetime sparseness and coefficient of variation
pdf("figs/sparseness_hists.pdf", height=4.5, width=8.5, pointsize=10)
par(mfrow=c(1,2), oma = c(0, 0, 3, 0), marc=c(4,4,3,4))

orn.col.alpha <- rgb(0, 100/255, 0, 0.7)
pn.col.alpha <- rgb(0, 0, 1, 0.7)
lhn.col.alpha <- rgb(160/255,32/255,240/255,0.7)

# lifetime sparseness
cell.orn.s <- lifetime.sparseness(orn)
cell.pn.s <- lifetime.sparseness(pn)
cell.lhn.s <- lifetime.sparseness(lhn)

pop.orn.s <- lifetime.sparseness(t(orn))
pop.pn.s <- lifetime.sparseness(t(pn))
pop.lhn.s <- lifetime.sparseness(t(lhn))

d.orn <- density(lifetime.sparseness(orn), bw="SJ")
d.pn <- density(lifetime.sparseness(pn), bw="SJ")
d.lhn <- density(lifetime.sparseness(lhn), bw="SJ")

plot(d.orn, col=orn.col.alpha, xlab='Sparseness Score', xlim=c(0,1), ylim=c(0,3.5), main='Lifetime Sparseness')
polygon(d.orn,  col=orn.col.alpha)
par(new=T)
plot(d.pn, col=pn.col.alpha, xlim=c(0,1), ylim=c(0,3.5), main='', ann=F, xaxt='n', yaxt='n')
polygon(d.pn, col=pn.col.alpha)
par(new=T)
plot(d.lhn, col=lhn.col.alpha, xlim=c(0,1), ylim=c(0,3.5), main='', ann=F, xaxt='n', yaxt='n')
polygon(d.lhn, col=lhn.col.alpha)

#mtext("Sparseness of All Cells", outer = TRUE, cex = 1)
legend('topright', legend=c("ORN", "PN", "LHN"), border="white", fill=c(orn.col.alpha, pn.col.alpha, lhn.col.alpha), col=c(orn.col.alpha, pn.col.alpha, lhn.col.alpha), box.col='white')
       
# coefficient of variations
#cv.orn <- apply(orn, 1, function(x) sd(x)/mean(x))
#cv.pn <- apply(pn, 1, function(x) sd(x)/mean(x))
#cv.lhn <- apply(lhn, 1, function(x) sd(x)/mean(x))

#d.orn <- density(cv.orn, bw="SJ")
#d.pn <- density(cv.pn, bw="SJ")
#d.lhn <- density(cv.lhn, bw="SJ")

d.orn <- density(lifetime.sparseness(t(orn)), bw="SJ")
d.pn <- density(lifetime.sparseness(t(pn)), bw="SJ")
d.lhn <- density(lifetime.sparseness(t(lhn)), bw="SJ")

plot(d.orn, col=orn.col.alpha, xlab='Sparseness Score', xlim=c(0,1), ylim=c(0,5.5), main='Population Sparseness')
polygon(d.orn,  col=orn.col.alpha)
par(new=T)
plot(d.pn, col=pn.col.alpha, xlim=c(0,1), ylim=c(0,5.5), main='', xaxt='n', yaxt='n', ann=F)
polygon(d.pn, col=pn.col.alpha)
par(new=T)
plot(d.lhn, col=lhn.col.alpha, xlim=c(0,1), ylim=c(0,5.5), main='', xaxt='n', yaxt='n', ann=F)
polygon(d.lhn, col=lhn.col.alpha)

#mtext("Sparseness of All Odors", outer = TRUE, cex = 1)
legend('topright', legend=c("ORN", "PN", "LHN"), border="white", fill=c(orn.col.alpha, pn.col.alpha, lhn.col.alpha), col=c(orn.col.alpha, pn.col.alpha, lhn.col.alpha), box.col='white')

dev.off()


################################################################################
## Plot changes in normalization/decorrelation and PCA analysis

# compute sums
s.orn <- apply(orn,1,sum)
s.nl <- apply(nl,1,sum)
s.pn <- apply(pn,1,sum)
s.lhn <- apply(lhn[!(rownames(lhn) %in% c("WatBl", "OilBl")),], 1, sum)

# compute eigenvalues by cell
pc.orn <- prcomp(orn)$sdev
pc.pn <- prcomp(pn)$sdev
pc.lhn <- prcomp(lhn)$sdev

# compute eigenvalues by odor
pc.odor.orn <- prcomp(t(orn))$sdev
pc.odor.pn <- prcomp(t(pn))$sdev
pc.odor.lhn <- prcomp(t(lhn))$sdev

## set up figure
pdf("figs/pca_plot.pdf", height=5.5, width=6.5, pointsize=10)
par(mfrow=c(3,3), mar=c(4,4,3,4))

# Variance explained plot by PCA
plot(100*pc.orn^2/sum(pc.orn^2), xlab="PC component", xlim=c(0,20), ylab="% variance", ylim=c(0,75), col="darkgreen", pch=20, main='ORN Cell')
plot(100*pc.pn^2/sum(pc.pn^2), xlab="PC component", xlim=c(0,20), ylab="% variance",ylim=c(0,75), col="blue", pch=20, main='PN Cell')
plot(100* pc.lhn^2/sum(pc.lhn ^2), xlab="PC component", xlim=c(0,20), ylab="% variance", ylim=c(0,75), col="purple", pch=20, main='LHN Cell')

plot(100*pc.odor.orn^2/sum(pc.odor.orn^2), xlab="PC component", xlim=c(0,20), ylab="% variance", ylim=c(0,75), col="darkgreen", pch=20, main='ORN Odor')
plot(100*pc.odor.pn^2/sum(pc.odor.pn^2), xlab="PC component", xlim=c(0,20), ylab="% variance",ylim=c(0,75), col="blue", pch=20, main='PN Odor')
plot(100* pc.odor.lhn^2/sum(pc.odor.lhn ^2), xlab="PC component", xlim=c(0,20), ylab="% variance", ylim=c(0,75), col="purple", pch=20, main='LHN Odor')


#pdf("figs/pca_project.pdf", width=8.5, height=2.5, pointsize=10)
#par(mfrow=c(1,3), mar=c(4,4,3,4))
pc.common.orn <- prcomp(orn)
pc.common.pn <- prcomp(pn)
pc.common.lhn <- prcomp(lhn)

plot(pc.common.orn$x[,1:2], col='darkgreen', pch=20, main='ORN Odor')
plot(pc.common.pn$x[,1:2], col='blue', pch=20, main='PN Odor')
plot(pc.common.lhn$x[,1:2], col='purple', pch=20, main='LHN Odor')
dev.off()

