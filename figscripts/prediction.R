# TODO
# [ ] Compute mutual info
# [ ] Do with max MI per cluster
# [ ] Compute category instead of identity
# [ ] Do comparison with LDA for each

################################################################################
## Encode LHN data and compute MI

# note that all matrices in this section are transponsed from what they are in the rest of the program
# i.e. rows are cells and columns are odors

# call significant responses and compute probabilities
lhn.trial.bin <- binarize.trials(lhn.mat)
lhn.probs.binom <- make.trial.probs(bin2data.frame(lhn.trial.bin))

lhn.odor.types.notnone <- factor(subset(lhn.odors.types, Type != "none")$Type)
# testing...

#lhn.probs.odor.types <- make.trial.probs(lhn.trial.bin, odor.class.labels=lhn.odor.types.notnone)
#lhn.probs.cell.types <- make.trial.probs(lhn.trial.bin, cell.class.labels=factor(rownames(lhn.mat)))
#lhn.probs.cell.odor.types <- make.trial.probs(lhn.trial.bin, cell.class.labels=x, odor.class.labels=lhn.odor.types.notnone)

# Compute mutual information for all LHN
lhn.mi <- mutual.info(t(lhn.probs.binom))

plot(cumsum(sort(lhn.mi)), xlab='Number of Cells', ylab='Total Mutual Information (bits)', type='l')
################################################################################
## Cluster cells by relative entropy

# plot relative entropy distance matrix
D.binom <- make.D(lhn.probs.binom, binom=T)
pdf("figs/klmat.pdf")
heatmap.2(D.binom, hclustfun = function(x) { hclust(as.dist(x), method="complete") }, scale='none', symm=T, trace='none', col=bin.rev.colors, main='Average Relative Entropy')
dev.off()

################################################################################

lhn.mat.noblanks <- lhn.mat[,!(colnames(lhn.mat) %in% c("OilBl", "WatBl"))]
lhn.mat.delta <- lhn.mat.noblanks - rep.col(lhn.spontaneous.rates, ncol(lhn.mat.noblanks))
lhn.delta <- avg.by.odor(lhn.mat.delta, colnames(lhn.mat.delta))
indep.gauss <- indep.gauss.model.all(lhn.mat.delta)
dep.gauss <- dep.gauss.model.all(lhn.mat.delta)

lhn.clust.cell <-  hclust(as.dist(1-cor(lhn)), method='ward')
lhn.types <- subset(lhn.odors.types, Type != "none")

################################################################################
## Look at encoding
sim.indep <- c()
for (i in 1:length(indep.gauss)) {
  sim.indep <- rbind(sim.indep, rmvnorm(1, mean=indep.gauss[[i]]$mu, sigma=indep.gauss[[i]]$sigma))
}

sim.dep <- c()
for (i in 1:length(indep.gauss)) {
  sim.dep <- rbind(sim.dep, rmvnorm(1, mean=dep.gauss[[i]]$mu, sigma=dep.gauss[[i]]$sigma))
}

pdf("figs/indep_heatmap.pdf", width=4, height=4)
heatmap(cor(t(sim.indep)), col=jet.colors, scale='n', main='Independent')
dev.off()

pdf("figs/indep_heatmap.pdf", width=4, height=4)
heatmap(cor(t(sim.dep)), col=jet.colors, scale='n', main='Correlated')
dev.off()


################################################################################
## Do predictions with Gaussian
cv.gauss.indep <- c()
for (i in 2:nrow(lhn.mat.delta)) {
  print(i)
  cv.gauss.indep <- c(cv.gauss.indep, cross.validate.gauss(lhn.mat.delta, lhn.clust.cell, i, indep=T))
}
  
cv.gauss.dep <- c()
for (i in 2:nrow(lhn.mat.delta)) {
  print(i)
  cv.gauss.dep <- c(cv.gauss.dep, cross.validate.gauss(lhn.mat.delta, lhn.clust.cell, i, indep=F))
}

cv.gauss.type.indep <- c()
for (i in 2:nrow(lhn.mat.delta)) {
  print(i)
  cv.gauss.type.indep <- c(cv.gauss.type.indep, cross.validate.gauss(lhn.mat.delta, lhn.clust.cell, i, indep=T, do.type=T))
}

cv.gauss.type.dep <- c()
for (i in 2:nrow(lhn.mat.delta)) {
  print(i)
  cv.gauss.type.dep <- c(cv.gauss.type.dep, cross.validate.gauss(lhn.mat.delta, lhn.clust.cell, i, indep=F, do.type=T))
}

save(cv.gauss.type.dep, cv.gauss.type.indep, file="data/gausspredicttype.rda")
#save(cv.gauss.indep, cv.gauss.dep, file="data/gausspredict.rda")
load("data/gausspredict.rda")

x = cross.validate.single.neurons(lhn.mat.delta, lhn.clust.cell, 128)
barplot(100*x, ylim=c(0,20), xlab='Cell Number', ylab='Average Prediction Accuracy (%)', main='Single Cell Predictions')

cv.single <- c()
for (i in 2:nrow(lhn.mat.delta)) {
  print(i)
  cv.single <- c(cv.single, mean(cross.validate.single.neurons(lhn.mat.delta, lhn.clust.cell, i)))
}

cv.single.type <- c()
for (i in 2:nrow(lhn.mat.delta)) {
  print(i)
  cv.single.type <- c(cv.single.type, mean(cross.validate.single.neurons(lhn.mat.delta, lhn.clust.cell, i, do.type=T)))
}

pdf("figs/single_pred.pdf")
plot(100*cv.single, ylim=c(0,30), ylab='Average Prediction Accurary (%)', xlab='Number of Subtypes', col='blue', type='l', main='Single Cell Predictions')
par(new=T)
plot(100*cv.single.type, ylim=c(0,30), ylab='', xlab='', xaxt='n', yaxt='n', ann=F, col='red', type='l')
abline(h=2.9, col='blue', lty=2)
abline(h=11.1, col='red', lty=2)
legend("topleft", legend=c("Identity", "Type"), fill=c('blue','red'), border=c('blue', 'red'), bty='n')       
dev.off()

################################################################################
## Redo everything with LDA
#lhn.mat.noblanks <- lhn.mat[,!(colnames(lhn.mat) %in% c("OilBl", "WatBl"))]
#success.rates.lda.all <- sapply(1:nrow(lhn.trial.bin), function(x) mean(cross.validate.odors.lda(lhn.mat.noblanks, x)))
#success.rates.lda.types <- sapply(1:nrow(lhn.trial.bin), function(x) mean(cross.validate.types.lda(lhn.mat.noblanks, x, lhn.odor.types.notnone)))

#save(success.rates.all, success.rates.lda.all, success.rates.lda.types, file="data/prediction.Rdata")
load("data/prediction.Rdata")
# plot everything
#pdf("figs/odor_prediction.pdf", width=4, height=4)
#plot(success.rates.lda.all[1:50], col='red', ylim=c(0,1), type='l', ylab='% Correct', xlab='Number of Clusters', main='Prediction of Odor Identity')
#par(new=T)
#plot(success.rates.lda.types[1:40], col='blue', type='l', xaxt='n', yaxt='n', ann=F)
#par(new=T)
#plot(success.rates.all[1:50], col='blue', type='l', xaxt='n', yaxt='n', ann=F)
#par(new=T)
#legend("bottomleft", legend=c("Bayes", "LDA"), fill=c('blue','red'), border=c('blue', 'red'), bty='n', cex=0.8)       
#dev.off()

################################################################################
## plot everything
pdf("figs/all_predictions.pdf", width=8.5, height=4)
par(mfrow=c(1,2))
plot(100*cv.gauss.indep, type='l', col='red',ylim=c(50,120), xlab='', ylab='',  xaxt='n', yaxt='n', ann=F)
par(new=T)
#plot(success.rates.lda.all, type='l', col='green',ylim=c(.5,1.2), xlab='', ylab='',  xaxt='n', yaxt='n', ann=F)
#par(new=T)
plot(100*cv.gauss.dep, type='l', col='blue', ylim=c(50,120), xlab='Number of Subtypes', ylab='Percent Correct (%)', main='Odor Identity Prediction')
legend("topleft", legend=c("Correlated", "Independent"), fill=c('blue','red'), border=c('blue', 'red'), bty='n')       
plot(100*cv.gauss.type.indep, type='l', col='red',ylim=c(50,110), xlab='', ylab='', xaxt='n', yaxt='n', ann=F)
par(new=T)
#plot(success.rates.lda.types, type='l', col='green',ylim=c(.5,1.1), xlab='', ylab='', xaxt='n', yaxt='n', ann=F)
#par(new=T)
plot(100*cv.gauss.type.dep, type='l', col='blue', ylim=c(50,110), xlab='Number of Subtypes', ylab='Percent Correct (%)', main='Odor Type Prediction')
legend("topleft", legend=c("Correlated", "Independent"), fill=c('blue','red'), border=c('blue', 'red'), bty='n')       
dev.off()

################################################################################
## Plot how clustering affects classification accuracy for single odors
success.rates.all <- sapply(2:nrow(lhn.trial.bin), function(x) mean(cross.validate.odors(lhn.trial.bin, predict.odor.all, num.classes=x)))
success.rates.types <- sapply(2:34, function(x) mean(cross.validate.types(lhn.trial.bin, predict.type.all, x, lhn.odor.types.notnone)))

success.rates.all <- c(success.rates.all[1], success.rates.all)

success.rates.all <- sapply(2:30, function(x) mean(cross.validate.types(lhn.trial.bin, predict.type.all,  num.classes=x)))
# test...
cross.validate.types(lhn.trial.bin, predict.type.all, 100, lhn.odor.types.notnone)
success.rates.exemplar <- sapply(1:nrow(lhn.trial.bin), function(x) mean(cross.validate.bayes(lhn.trial.bin, num.classes=x, exemplar=T))
                                 
# normalized mutual information scores
MI.scores <- sapply(1:nrow(lhn.trial.bin), function(x) cluster.mutual.info(lhn.trial.bin,x))/mean(mutual.info(lhn.probs.binom))

entropy.scores <- sapply(1:nrow(lhn.trial.bin), function(x) cluster.avg.entropy(lhn.trial.bin,x))

pdf("figs/miplot.pdf", width=4, height=4)
#par(mfrow=c(1,2))
#plot(success.rates, xlab="Number of Clusters", ylab="Mean % Odors Correctly Identified", type='l', main='Classification Accuracy')
plot(MI.scores, ylab="Normalized Mutual Information", xlab="Number of Subtypes", type='l', main='Mutual Information')
#plot(entropy.scores, type='l', main='Average Entropy', xlab='Number of Clusters', ylab='Average Entropy (bits)')
dev.off()



################################################################################
## Plot how clustering affects classification accuracy for odor classes

# plot confusion matrices
#confusion <- single.neuron.confusion(lhn.trial.bin)
# TODO Calculate mutual information and show redundancy

#leave.out <- seq(1,ncol(lhn.mat.noblank),4)



                                 
