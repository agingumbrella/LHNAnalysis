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
#par(mfrow=c(3,4))
#for (i in 1:num.pn.clusts) {
#    curr <- ifelse(pn.factor.labels == i, 1, 0)
#    g <- lda(t(pn), factor(curr))
#    plot.model.response.hist(make.model.responses(pn, curr, 100, TRUE, method="lda"), 30, title=paste("Cluster", i))
#}

### Do just for different odor categories
make.odor.type.label <- function(pn, types, curr.type) {
  return(ifelse(rownames(pn) %in% all.odor.type[all.odor.type$Type == curr.type,1], 1, 0))
}


  odor.types <- unique(all.odor.type$Type)
  
  terpenes <- make.odor.type.label(pn, odor.types, "terpene")
  aldehydes <- make.odor.type.label(pn, odor.types, "amine")
  acid <- make.odor.type.label(pn, odor.types, "acid")
  alcohol <- make.odor.type.label(pn, odor.types, "alcohol")
  ester <- make.odor.type.label(pn, odor.types, "ester")
  aromatic <- make.odor.type.label(pn, odor.types, "aromatic")
  ketone <- make.odor.type.label(pn, odor.types, "ketone")
  lactone <- make.odor.type.label(pn, odor.types, "lactone")
  sulfur <- make.odor.type.label(pn, odor.types, "sulfur")

# TODO Make this work...
pdf("figs/modeled_responses.pdf")
par(mfrow=c(3,3))
plot.model.response.hist(make.model.responses(pn, aldehydes, 100), x.min=-0.005, x.max=0.01, ymax=1000, title='Aldehyde')
plot.model.response.hist(make.model.responses(pn, acid, 100), x.min=-0.01, x.max=0.005, ymax=800, title='Acid')
legend("topright", legend=c("Targeted Odors", "Other Odors"), fill=c("red","blue"), border=c("red","blue"), bty='n', cex=0.5, title="Total Input to Model LHN")
plot.model.response.hist(make.model.responses(pn, alcohol, 100), x.min=-0.01, x.max=0.01, ymax=1000, title='Alcohol')
plot.model.response.hist(make.model.responses(pn, ester, 100), x.min=-0.005, x.max=0.001, ymax=1200, title='Ester')
plot.model.response.hist(make.model.responses(pn, aromatic, 100), x.min=-0.01, x.max=0.005, ymax=800, title='Aromatic')
plot.model.response.hist(make.model.responses(pn, ketone, 100), x.min=-0.0, x.max=0.01, ymax=1000, title='Ketone')
plot.model.response.hist(make.model.responses(pn, lactone, 100), x.min=-0.005, x.max=0.006, ymax=1000, title='Lactone')
plot.model.response.hist(make.model.responses(pn, sulfur, 100), x.min=-0.005, x.max=0.006, ymax=1000, title='Sulfur')
dev.off()

pdf("figs/modeled_terpene.pdf",width=4,height=4)
plot.model.response.hist(make.model.responses(pn, terpenes, 1000), x.min=-0.1, x.max=0.4, ymax=15, title='Terpene')
legend("topright", legend=c("Terpenes", "Other Odorants"), border='white', fill=c('red', 'blue'), col=c('red','blue'), box.col='white')
#plot(prcomp(lhn.common)$x[,1:2],xlim=c(-80,80),ylim=c(-80,80),col='magenta')
#par(new=T)
#plot(prcomp(pn.common)$x[,1:2],xlim=c(-80,80),ylim=c(-80,80),col='cyan')

groups <- list(terpenes, aldehydes, acid, alcohol, ester, aromatic, ketone, lactone, sulfur)
ldas <- lapply(groups, function(x) fisher.lda(pn, x))
weights <- lapply(ldas, function(x) x$w)
names(weights) <- c("terpene", "amine", "acid", "alcohol", "ester", "aromatic" ,"ketone", "lactone", "sulfur")
thresh <- lapply(ldas, function(x) as.vector(t(x$w) %*% (x$m1 +x$m2)/2))
names(thresh) <- names(weights)
classifications = classify.all(pn, groups, weights, thresh)
rownames(classifications) <- names(weights)
colnames(classifications) <- names(weights)
pdf("figs/confusion.pdf", width=4, height=4)
heatmap.2(classifications, Colv=NA, Rowv=NA, scale='n', col=jet.colors, trace='n', main='Response Probability', sepwidth=c(0,0), density.info='none', dendrogram='none', lmat=rbind( c(0, 3), c(2,1), c(0,4) ), lwid=c(1, 4), lhei=c(1, 4, 1 ))
dev.off()
barplot(diag(classifications), main='Prediction Accuracy')
