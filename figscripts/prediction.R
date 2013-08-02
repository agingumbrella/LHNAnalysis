  
lhn.trial.bin <- binarize.trials(lhn.mat, lhn.labels)
lhn.probs.binom <- make.trial.probs(lhn.trial.bin)
#lhn.mat.noblank <- lhn.mat[, !(colnames(lhn.mat) %in% c("OilBl", "WatBl"))]

# plot relative entropy distance matrix
D.binom <- make.D(lhn.probs.binom, binom=T)
pdf("figs/klmat.pdf")
heatmap.2(D.binom, hclustfun = function(x) { hclust(as.dist(x), method="complete") }, scale='none', symm=T, trace='none', col=bin.rev.colors, main='Relative Entropy')
dev.off()
labels = cutree(hclust(as.dist(D.binom), method="complete"),4)

# Show how clustering affects classification accuracy
success.rates <- sapply(1:nrow(lhn.trial.bin), function(x) mean(cross.validate(lhn.trial.bin, num.classes=x)))

# normalized mutual information scores
MI.scores <- sapply(1:nrow(lhn.trial.bin), function(x) cluster.mutual.info(lhn.trial.bin,x))/mean(mutual.info(lhn.probs.binom))

pdf("figs/classification_acc_and_mi.pdf", width=8.5, height=4)
par(mfrow=c(1,2))
plot(success.rates, xlab="Number of Clusters", ylab="Mean % Odors Correctly Identified", type='l', main='Classification Accuracy')
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
