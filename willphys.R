#setwd("~/Documents/Neuroscience/jefferis_lab/shahar_data/RandIgor/LHNAnalysis")
require(RColorBrewer)
require(NMF)
require(gphys)
require(gplots)
require(glmnet)
require(MASS)
require(penalized)
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

## HALLEM AND CARLSON DATA
x.names <- c("X2a", "X7a", "X9a", "X10a", "X19a", "X22a", "X23a", "X33b", "X35a", "X43a", "X43b", "X47a", "X47b", "X49b", "X59b", "X65a", "X67a", "X67c", "X82a", "X85a", "X85b", "X85f", "X88a", "X98a")
orn.names <- c("Or2a", "Or7a", "Or9a", "Or10a", "Or19a", "Or22a", "Or23a", "Or33b", "Or35a", "Or43a", "Or43b", "Or47a", "Or47b", "Or49b", "Or59b", "Or65a", "Or67a", "Or67c", "Or82a", "Or85a", "Or85b", "Or85f", "Or88a", "Or98a")
spontaneous.rates <- c(8, 17, 3, 14, 29, 4, 9, 25, 17, 21, 2, 1, 47, 8, 2, 18, 11, 6, 16, 14, 13, 7, 26, 12) # ordered by OR

val.and.ors <- read.csv("val_orn.csv", header=T)
rownames(val.and.ors) <- val.and.ors$Odorant
val <- val.and.ors[,c("Odorant", "Attraction.index")]

orn <- as.matrix(val.and.ors[,x.names])
colnames(orn) <- orn.names
orn <- apply(orn, 1, function(x) {x + spontaneous.rates})
orn[orn < 0] <- 0


# plot avg rates
nl <- make.pn.rates(orn,m=0)
pn <- make.pn.rates(orn)

## Make figure 1 of Luo and Axel
# panel 1
jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))(256)
bin.colors <- colorRampPalette(c('white', 'black'))(100)
bin.rev.colors <- colorRampPalette(c('black', 'white'))(100)
heatmap.2(t(orn),col=jet.colors,scale="none", trace="none")
heatmap.2(t(pn), col=jet.colors, scale="none", trace="none")

# panel 2
par(mfrow=c(2,3))

s.orn <- apply(orn,1,sum)
s.nl <- apply(nl,1,sum)
s.pn <- apply(pn,1,sum)
barplot(s.orn/mean(s.orn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn),main="ORN")
barplot(s.nl/mean(s.nl), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn), main="NL")
barplot(s.pn/mean(s.pn), xlab="Odor", ylab="firing rates/avg", ylim=c(0,2.5), names.arg=1:nrow(orn),main="PN")

# panel 3
pc.orn <- princomp(t(orn))$sdev
pc.nl <- princomp(t(nl))$sdev
pc.pn <- princomp(t(pn))$sdev
#par(mfrow=c(1,3))
plot(100*pc.orn^2/sum(pc.orn^2), xlab="PC component", ylab="% variance", ylim=c(0,75))
plot(100*pc.nl^2/sum(pc.nl^2), xlab="PC component", ylab="% variance", ylim=c(0,75))
plot(100*pc.pn^2/sum(pc.pn^2), xlab="PC component", ylab="% variance",ylim=c(0,75))

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


# make correlation plots among ORNs
orn.cor <- cor(orn, use='complete.obs')
orn.odor.cor <- cor(t(orn), use='complete.obs')
pn.cor <- cor(pn, use='complete.obs')
pn.odor.cor <- cor(t(pn), use='complete.obs')
par(mfrow=c(2,2))
heatmap.2(orn.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')
heatmap.2(orn.odor.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')

heatmap.2(pn.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')
heatmap.2(pn.odor.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')

## Analysis of LHN data
# compute average firing rates


mean.rates <- lapply(SpikesT, cell.mean.rates)

allodors <- unique(unlist(sapply(mean.rates,names)))
goododors <- allodors[1:36]
common.odors <- list(E2Hex="E2Hexenal", GerAc="Geranylacetat", Prpyl="Propylacetate", IPenA="Isopentanoic acid", Et3HB="Ethyl-3-hydroxybutyrate", MetSl="Methylsalicylat", PeEtA="Phenylethylacetat", HexAc="Hexanoic acid", AceAc="Acetic acid", "1HxOl"="1-Hexanol", BeZal="Benzaldehyd", Amnia="Ammonium hydroxide", bCitr="beta-Citrollenol", MtAct="Methylacetat", BtrAc="Butyric acid", PrpnA="Propionic acid")

total.trials <- sum(unlist(sapply(SpikesT, function(x) {sapply(x, length)})))
total.usable <- length(SpikesT)*36*4


lhn.mat <- make.total.data.matrix(SpikesT, length(mean.rates), 36*4)
lhn.labels <- factor(colnames(lhn.mat))

lhn.bin <- binarize.mat.per.cell(lhn.mat, lhn.labels)

rates.mat <- matrix(0, ncol=length(mean.rates), nrow=length(goododors))
rownames(rates.mat) <- goododors
colnames(rates.mat) <- names(mean.rates)
for (j in names(mean.rates)) {
  for (i in goododors) {
    if (i %in% names(mean.rates[[j]])) {
      rates.mat[[i,j]] <- mean.rates[[j]][[i]]$m
    } else {
      rates.mat[[i,j]] <- NA
    }
  }
}

all.lhn <- matrix(0, ncol=length(mean.rates), nrow=length(goododors))

rates.nona.mat <- rates.mat[,apply(rates.mat, 2, function(y) !any(is.na(y)))]
lhn <- rates.nona.mat

# equalization plot
s.lhn <- apply(rates.nona.mat, 1, sum)
barplot(s.lhn/mean(s.lhn), xlab="Odor", ylab="firing rates/avg", names.arg=1:length(s.lhn), main="LHN")

# PCA plot
pc.lhn <- princomp(t(rates.nona.mat))$sdev
#par(mfrow=c(1,3))
plot(100* pc.lhn ^2/sum(pc.lhn ^2), xlab="PC component", ylab="% variance", ylim=c(0,75))

# plot correltion b/w ORN and PN firing rates and LHN firing rates on common odors
val.common <- val[unlist(common.odors),2]
orn.common <- orn[unlist(common.odors),]
pn.common <- pn[unlist(common.odors),]
lhn.common <- rates.nona.mat[unlist(names(common.odors)),]
rownames(lhn.common) <- unlist(common.odors)
orn.pn.corr <- cor(orn.common, pn.common)
orn.lhn.corr <- cor(orn.common, lhn.common)
pn.lhn.corr <- cor(pn.common, lhn.common)

heatmap.2(orn.lhn.corr, col=jet.colors,scale='none', trace='none')
heatmap.2(pn.lhn.corr, col=jet.colors,scale='none', trace='none')

# LHN heatmap
lhn.cor <- cor(rates.nona.mat, use='complete.obs')
lhn.odors.cor <- cor(t(rates.nona.mat), use='complete.obs')
orn.common.cor <- cor(t(orn.common), use='complete.obs')
pn.common.cor <- cor(t(pn.common), use='complete.obs')
lhn.common.cor <- cor(t(lhn.common), use='complete.obs')
heatmap.2(lhn.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')
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

# TODO Calculate mutual information and show redundancy

# Decoding: try to predict odors given population data
# -- Use total population data
# train on 3 examples for each cell and test on last


lhn.mat.noblank <- lhn.mat[, !(colnames(lhn.mat) %in% c("OilBl", "WatBl"))]
leave.out <- seq(1,ncol(lhn.mat.noblank),4)
train.lhn.rates <- make.trial.rates(lhn.mat.noblank[, !(1:ncol(lhn.mat.noblank) %in% leave.out)])
test.lhn.rates <- lhn.mat.noblank[,leave.out]
train.lhn.bin <- make.trial.probs(lhn.trial.bin[,!(1:ncol(lhn.trial.bin) %in% leave.out)], 3)
test.lhn.bin <- lhn.trial.bin[,leave.out]

D <- make.D(lhn.probs)


# returns labels with c1 and c2 merged

# use single-linkage clustering to get minimum info dists

D <- make.D(lhn.probs)
heatmap.2(D, hclustfun = function(x) { hclust(as.dist(x), method="complete") }, scale='none', symm=T, trace='none', col=bin.rev.colors)


entropies <- c()
for (i in 1:nrow(train.lhn.bin)) {
    entropies <- c(entropies, sum(entropy.binom(train.lhn.bin[i,])))
}
