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



# had to do this by hand because the names of the files included date strings and didn't exactly match the cross IDs in the ephys data
anatomy.cross.ids <- c("SF274-JK1742", "SF274-JK1742", "SF274-JK304", "SF274-JK1742", "SF274-JK1742", "SF274-JK1742", "SF274-JK294", "JK843-SF341", "JK843-SF401PC", "SF274-JK1742", "JK843-SF401PC", "JK843-SF401PC", "JK2204Cha", "JK801Cha", "JK801Cha")
rownames(dists) <- colnames(dists) <- anatomy.cross.ids

num.anatomy.cross.ids <- length(unique(rownames(dists)))

#heatmap.2(dists, symm=T, trace='none', scale='none', RowSideColors=make.side.colors(anatomy.cross.ids, cross.colors), ColSideColors=make.side.colors(anatomy.cross.ids, cross.colors), col=bin.colors)
trace.cluster <- hclust(dist(dists), method='ward')
trace.cluster.labels <- cutree(trace.cluster, num.anatomy.cross.ids)

pdf("cell_cluster_res.pdf", height=8.5, width=3)
par(oma = c(0, 0, 0, 0))
plot(dendrapply(as.dendrogram(trace.cluster),color.leaves.lhn.anatomy), horiz=T)
dev.off()

# Show overlap between clustering by anatomy and by physiology

# Do joint clustering with anatomy and physiology
