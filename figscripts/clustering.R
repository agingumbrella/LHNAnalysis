################################################################################
## Utility functions used for dendrogram plots

# there's probably a better way of doing this, but the only way I could figure out to add the right colors to the labels on the leaves of the dendrograms was to write separate functions for each that can then be dendrapply'd

# for odor clustering of ORN, PN
color.leaves.all <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = type.colors[[as.character(all.odor.type[all.odor.type$Odor == a$label,]$Type)]], lab.cex=0.6))
  }
  n
}

# for odor clustering of LHNs
color.leaves.lhn <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
   
    col <- type.colors[[as.character(lhn.odors.types[lhn.odors.types$Odorant == a$label,]$Type)]]
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=1))
  }
  n
}

# for cell clustering of LHNs
color.leaves.lhn.cell <- function(n) {
  if (is.leaf(n)) {
    a <- attributes(n)
    col <- cross.colors[[a$label]]
    attr(n, "nodePar") <- c(a$nodePar, list(lab.col = col, lab.cex=0.4))
  }
  n
}

# compute per class purity score
per.class.purity <- function(labels, true.cats) {
  p <- c()
  for (i in unique(true.cats)) {
    curr <- labels[true.cats == i]
    most.freq <- as.integer(names(which.max(table(curr))))
    # fraction of pairs that are in the same cluster vs different clusters
    p <- c(p, sum(curr == most.freq)/length(curr))
  }
  return(p)
}

# compute cluster purity for all clusters known classes
cluster.purity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

################################################################################
## Dendrograms of cells and odors 
require(dendroextras)
# multiple ways of clustering -- uncomment to try, but they all give similar results

# using PCA
#orn.pcs <- prcomp(orn)
#pn.pcs <- prcomp(pn)
#lhn.pcs <- prcomp(lhn[!(rownames(lhn) %in% c("OilBl", "WatBl")),])
#orn.clust <- dendrapply(as.dendrogram(hclust(dist(orn.pcs$x[,1:10]), method='ward')), color.leaves.all)
#pn.clust <- dendrapply(as.dendrogram(hclust(dist(pn.pcs$x[,1:10]), method='ward')), color.leaves.all)
#lhn.clust <- dendrapply(as.dendrogram(hclust(dist(lhn.pcs$x[,1:10]), method='ward')), color.leaves.lhn)

# using euclidean distances
#orn.clust <- hclust(dist(orn), method='ward')
#pn.clust <- hclust(dist(pn), method='ward')
#lhn.clust <- hclust(dist(lhn[!(rownames(rates.nona.mat) %in% c("OilBl", "WatBl")),]), method='ward')

# using correlation scores -- same as matrix in other figures
orn.clust <- hclust(as.dist(1-cor(t(orn))), method='ward')
pn.clust <- hclust(as.dist(1-cor(t(pn))), method='ward')
lhn.clust <- hclust(as.dist(1-cor(t(lhn[!(rownames(lhn) %in% c("OilBl", "WatBl")),]))), method='ward')

# make and color the dendrograms for each cell class
orn.den <- dendrapply(as.dendrogram(orn.clust), color.leaves.all)
pn.den <- dendrapply(as.dendrogram(pn.clust), color.leaves.all)
lhn.den <- dendrapply(as.dendrogram(lhn.clust), color.leaves.lhn)


num.clusts <- 10
orn.clust.labels <- cutree(orn.clust, num.clusts)
pn.clust.labels <- cutree(pn.clust, num.clusts)
lhn.clust.labels <- cutree(lhn.clust, 8)

orn.types.num <- as.numeric(factor(all.odor.type[names(orn.clust.labels),]$Type))
pn.types.num <- as.numeric(factor(all.odor.type[names(pn.clust.labels),]$Type))
lhn.types.num <-  as.numeric(factor(lhn.odors.types[names(lhn.clust.labels),]$Type))

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

# plow each of the dendrograms with legend
pdf("cluster_res.pdf", height=6, width=8.5)
par(mfrow=c(1,3), oma = c(0, 0, 0, 0), mar=c(3,1,1,6))
plot(orn.den, horiz=T, main="ORN", center=T)
legend("topleft", legend=names(type.colors), fill=unlist(type.colors), border=unlist(type.colors), bty='n', cex=0.8)
plot(pn.den, horiz=T, main="PN", center=T)
plot(lhn.den, horiz=T, main="LHN", center=T)
dev.off()

pdf("cluster_purity.pdf", height=2.5, width=8.5)
par(mfrow=c(1,3))
barplot(orn.purity, col=sapply(names(orn.purity), function(x) type.colors[[x]]), las=3, ylab="Purity (%)")
barplot(pn.purity, col=sapply(names(pn.purity), function(x) type.colors[[x]]), las=3, ylab="Purity (%)")
barplot(lhn.purity, col=sapply(names(lhn.purity), function(x) type.colors[[x]]), las=3, ylab="Purity (%)")
dev.off()

#######################################################
## cluster LHNs and look at similarity with split-Gal4 lines

# cluster LHNs by correlation score
colnames(lhn.delta.rates)[colnames(lhn.delta.rates) %in% "JKSF274-JK671"] <- "SF274-JK671"
lhn.cells.clust <- hclust(dist(t(lhn.delta.rates)), method='ward')

# separate LHNs into clusters (number of different lines)
lhn.cell.labels <- cutree(lhn.cells.clust, length(unique(colnames(lhn.delta.rates))))

# make dendrogram for LHNs
cross.colors <- rainbow(length(unique(colnames(lhn.delta.rates))))
names(cross.colors) <- unique(colnames(lhn.delta.rates))
lhn.cell.den <- dendrapply(as.dendrogram(lhn.cells.clust), color.leaves.lhn.cell)

# plot dendrogram clustering of cells, with cross ID indicated by color
pdf("figs/cell_cluster_res.pdf", height=8.5, width=3, pointsize=10)
par(oma = c(0, 0, 0, 0), mar=c(3,1,1,6))
plot(lhn.cell.den, horiz=T, main="LHN", center=T, xlab='Distance')
legend("topleft", legend=names(cross.colors), fill=unlist(cross.colors), border=unlist(cross.colors), bty='n', cex=0.5)
dev.off()

library(mclust)
apply.collapsed <- function(rates, clusters, odors, func) {
  out <- matrix(0, nrow=length(unique(clusters)), ncol=length(unique(odors)))
  rownames(out) <- unique(clusters)
  colnames(out) <- unique(odors)
  for (i in unique(clusters)) {
    for (j in unique(odors)) {
      curr <- rates[clusters %in% i, odors %in% j]
      out[i,j] <- func(curr)
    }
  }
  return(out)
}




#km <- kmeans(t(lhn.delta.rates), length(unique(colnames(lhn.delta.rates))))
km <- Mclust(t(lhn.delta.rates), G=10)
x = apply.collapsed(t(lhn.delta.rates), km$classification, factor(lhn.odors.types$Type[lhn.odors.types$Type != "none"]), function(x) mean(x))
pdf("figs/celltyperesponse.pdf")
heatmap.2(x, scale='row', col=jet.colors, trace='n', Colv=NA, Rowv=NA, main='Mean Per Cluster, Per Odor Type Responses', sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(0.5, 3.), lhei=c(1, 4, 1), xlab='Odorant Type', ylab='Cluster ID')
dev.off()

plot(km.pc$x[,1:2])
# compute per cross purity score
lhn.cell.purity <- per.class.purity(km$classification, colnames(lhn.delta.rates))
names(lhn.cell.purity) <- unique(colnames(lhn.delta.rates))

pdf("figs/lhn_cell_purity.pdf", width=4, height=4, pointsize=10)
par(oma=c(0,0,0,0), mar=c(8,4,3,3))
barplot(lhn.cell.purity*100, col=cross.colors, las=3, ylab="Purity (%)", main='LHN Line Type Purity Scores')
dev.off()
