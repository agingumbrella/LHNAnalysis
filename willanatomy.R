
################################################################################
# utility functions for plotting dendrogram (modified from clustering)
color.leaves.lhn.anatomy <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    col <- cross.colors[[a$label]]
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=1))
  }
  n
}

# utility function to plot each cluster of neuron structures -- from output of pvpick
# structs is a neuronlist
plot.structure.clusters <- function(structs, labels, clusts) {
  bg3d('white');
#  par3d(windowRect=c(150,150,150,150));
  colors <- colorRampPalette(c("darkgreen", "magenta", "black"))(length(clusts))
  for (x in 1:length(clusts)) {
    clear3d()
    plot3d(neurons[labels %in% clusts[[x]]], col=colors[x], WithLine=T, WithNodes=F)
   	view3d(zoom=0.75);
    par3d(userMatrix=matrix(c(1,rep(0,4), -1, rep(0,4), -1, rep(0,4), 1),ncol=4))
    snapshot3d(paste("figs/cluster_",x,".png",sep=""))
  }
}



#heatmap.2(dists, symm=T, trace='none', scale='none', RowSideColors=make.side.colors(anatomy.cross.ids, cross.colors), ColSideColors=make.side.colors(anatomy.cross.ids, cross.colors), col=bin.colors)
trace.cluster.hclust <- hclust(as.dist(1-cor(lhn.anatomy.dists)), method='ward')
trace.cluster <- pvclust(lhn.anatomy.dists)
trace.cluster.labels <- pvpick(trace.cluster, alpha=0.9)$cluster
#trace.cluster.labels <- cutree(trace.cluster, num.anatomy.cross.ids)

pdf("figs/anatomy_cluster.pdf", height=8.5, width=3)
par(oma = c(0, 0, 0, 0), mar=c(3,1,1,6))
plot(dendrapply(as.dendrogram(trace.cluster.hclust),color.leaves.lhn.anatomy), horiz=T)
dev.off()

# plot structures
open3d();
plot.structure.clusters(neurons, colnames(lhn.anatomy.dists), trace.cluster.labels)

# Show overlap between clustering by anatomy and by physiology -- correlations and overlapping clusters
lhn.with.anatomy <- lhn[,which(lhn.cell.names %in% traced.info$Cell)]
lhn.cell.dists <- as.matrix(dist(t()))


heatmap.2(cor(lhn.cell.dists, lhn.anatomy.dists), symm=T, trace='none', scale='none')


# Do joint clustering with anatomy and physiology?
#require(iCluster)
#joint.cluster <- iCluster(list(cells=t(lhn), anatomy=lhn.anatomy.dists), 10, lambda=c(0.2,0.2))

