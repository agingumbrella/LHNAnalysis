
################################################################################
# utility functions for plotting dendrogram (modified from clustering)
color.leaves.lhn.anatomy <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    col <- anatomy.lines.colors[[a$label]]
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=1))
  }
  n
}

# utility function to plot each cluster of neuron structures -- from output of pvpick
# structs is a neuronlist
plot.structure.clusters <- function(structs, labels, clusts) {
  bg3d('white')
  colors <- list(c("darkgreen","green"), c("red","magenta"), c("gray","black"))
  for (x in 1:length(clusts)) {
    clear3d()
    curr.colors <- colorRampPalette(colors[[x]])(sum(labels %in% clusts[[x]]))
    plot3d(neurons[labels %in% clusts[[x]]], col=curr.colors, WithLine=T, WithNodes=F)
    view3d(zoom=0.75)
    par3d(userMatrix=matrix(c(1,rep(0,4), -1, rep(0,4), -1, rep(0,4), 1),ncol=4))
    snapshot3d(paste("figs/cluster_",x,".png",sep=""))
  }
}

# counts the fraction of times each line was in the same cluster in the two clusterings
per.class.overlap <- function(x,y) {
  ret <- list()
  for (l in unique(names(x))) {
    curr.x <- x[names(x) %in% l]
    curr.y <- y[names(y) %in% l]
    ret[[l]] <- purity(factor(curr.x), factor(curr.y))
  }
  return(ret)
}

anatomy.lines <- unique(colnames(lhn.anatomy.dists))
anatomy.lines.colors <- colorRampPalette(c("red","green","yellow","blue", "purple"))(length(anatomy.lines))
names(anatomy.lines.colors) <- anatomy.lines
#heatmap.2(dists, symm=T, trace='none', scale='none', RowSideColors=make.side.colors(anatomy.cross.ids, cross.colors), ColSideColors=make.side.colors(anatomy.cross.ids, cross.colors), col=bin.colors)
trace.cluster.hclust <- hclust(as.dist(1-cor(lhn.anatomy.dists)), method='ward')
trace.cluster <- pvclust(lhn.anatomy.dists)
#trace.cluster.labels <- pvpick(trace.cluster, alpha=0.9)$cluster
trace.cluster.labels <- cutree(trace.cluster.hclust, 3)

pdf("figs/anatomy_cluster.pdf", height=8.5, width=6)
par(mfrow=c(1,2), oma = c(0, 0, 0, 0), mar=c(2,2,2,7))
plot(dendrapply(as.dendrogram(trace.cluster.hclust),color.leaves.lhn.anatomy), horiz=T, main='Anatomical', cex.main=1)
#par(oma = c(0, 0, 0, 0), mar=c(3,1,1,6))
legend("topleft", legend=anatomy.lines, fill=unlist(anatomy.lines.colors), border=unlist(anatomy.lines.colors), bty='n', cex=0.5)
plot(dendrapply(as.dendrogram(hclust(as.dist(1-cor(lhn.cell.dists)), method='ward')),color.leaves.lhn.anatomy), horiz=T, main='Functional',cex.main=1)
dev.off()

# plot structures
open3d();
plot.structure.clusters(neurons, colnames(lhn.anatomy.dists), trace.cluster.labels)

# Show overlap between clustering by anatomy and by physiology -- correlations and overlapping clusters
lhn.with.anatomy <- lhn[,which(lhn.cell.names %in% traced.info$Cell)]
lhn.cell.dists <- as.matrix(dist(t(lhn.with.anatomy)))
phys.cluster.hclust <- hclust(as.dist(1-cor(lhn.cell.dists)), method='ward')
phys.cluster.labels <- cutree(phys.cluster.hclust, 3)


heatmap(cor(lhn.cell.dists, lhn.anatomy.dists), symm=F, trace='none', scale='none')

rownames(lhn.cell.dists) <- 1:nrow(lhn.cell.dists)
rownames(lhn.anatomy.dists) <- 1:nrow(lhn.anatomy.dists)
mantel.rtest(as.dist(lhn.cell.dists), as.dist(lhn.anatomy.dists))
require(ade4)
mantel.rtest(as.dist(lhn.cell.dists), as.dist(lhn.anatomy.dists))

require(dendroextras)
phys.dendro <- as.dendrogram(phys.cluster.hclust)
anatomy.dendro <- as.dendrogram(trace.cluster.hclust)

pdf("figs/phys_heatmap.pdf",width=4,height=4)
heatmap(lhn.with.anatomy[,phys.cluster.hclust$order], Colv=NA, Rowv=NA, col=jet.colors,  scale='n', cexRow=0.6, cexCol=0.8)
dev.off()
pdf("figs/anatomy_heatmap.pdf",width=4,height=4)
heatmap(lhn.with.anatomy[,trace.cluster.hclust$order], Colv=NA, Rowv=NA, col=jet.colors,  scale='n', cexRow=0.6, cexCol=0.8)
dev.off()

## compute purity scores for overlap between physiological and anatomical clustering
anatomy.purity <- purity(factor(trace.cluster.labels), factor(phys.cluster.labels))

per.class.purity(factor(phys.cluster.labels), factor(trace.cluster.labels))
#per.class <- per.class.purity(trace.cluster.labels, phys.cluster.labels)
trace.labels.over <- cutree(trace.cluster.hclust, length(anatomy.lines))
phys.labels.over <- cutree(phys.cluster.hclust, length(anatomy.lines))

phys.purity <- per.class.purity(phys.labels.over, factor(names(phys.labels.over)))
# Do joint clustering with anatomy and physiology?
#require(iCluster)
#joint.cluster <- iCluster(list(cells=t(lhn), anatomy=lhn.anatomy.dists), 10, lambda=c(0.2,0.2))

# try to predict class labels for left out cells

require(e1071)
preds <- lapply(3:15, function(x) {
  preds <- c()
  for (i in 1:10) {
    idx <- sample(1:ncol(lhn.with.anatomy), x)
    not.idx <- !(1:ncol(lhn.with.anatomy) %in% idx)
    tryCatch(model <- lda(t(lhn.with.anatomy[, idx]), factor(trace.cluster.labels[idx])), error=function(e) NA, finally=preds <- c(preds, sum(predict(model, t(lhn.with.anatomy[, not.idx]))$class == trace.cluster.labels[not.idx])/length(trace.cluster.labels[not.idx])))
  }
  preds
})


pdf('figs/anatomy_pred.pdf', width=4,height=4)
par(oma=c(1,1,1,1))
plot(3:15,sapply(preds,function(x) mean(x[!is.na(x)])), type='l', ylim=c(0,1), xlab='Number Labeled Neurons', ylab='Average Prediction Accuracy (%)', main='Prediction of Class Labels')
dev.off()

