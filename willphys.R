setwd("~/Documents/Neuroscience/jefferis_lab/shahar_data/RandIgor")
require(RColorBrewer)
require(NMF)
require(gphys)
require(gplots)
require(nnet)

source("Avg_IgorImport.R")
source("PSTH_FUNC.R")
source("ClusteringFunctions.R")



#load("spike_summary_133.rda")
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
for (i in 1:nrow(orn)) {
  orn[i,] <- orn[i,] + spontaneous.rates
}
orn[orn < 0] <- 0


make.pn.rates <- function(x, R.max=165, sigma=12, m=0.05) {
  pn <- matrix(0, ncol=ncol(x), nrow=nrow(x))
  rownames(pn) <- rownames(x)
  colnames(pn) <- colnames(x)
  # lateral suppression factor
  for (i in 1:nrow(orn)) {
    s <- m*sum(orn[i,])
    for (j in 1:ncol(orn)) {      
      pn[i,j] <- R.max*((x[i,j]^1.5)/(sigma^1.5 + x[i,j]^1.5 + s^1.5))
    }
  }
  return(pn)
}

# plot avg rates
nl <- make.pn.rates(orn,m=0)
pn <- make.pn.rates(orn)

## Make figure 1 of Luo and Axel
# panel 1
jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))(256)
bin.colors <- colorRampPalette(c('white', 'black'))(100)
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
pc.orn <- princomp(orn)$sdev
pc.nl <- princomp(nl)$sdev
pc.pn <- princomp(pn)$sdev
#par(mfrow=c(1,3))
plot(100*pc.orn^2/sum(pc.orn^2), xlab="PC component", ylab="% variance", ylim=c(0,75))
plot(100*pc.nl^2/sum(pc.nl^2), xlab="PC component", ylab="% variance", ylim=c(0,75))
plot(100*pc.pn^2/sum(pc.pn^2), xlab="PC component", ylab="% variance",ylim=c(0,75))

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

# compute the poststimulus rate in some window
ps.rate <- function(spikes, start=0.5, end=1.5) {
  return(length(spikes[spikes >= start & spikes <= end])/(end-start))
}

# compute the average rate for some odor
ps.mean.rate <- function(spikes, start=0.5, end=1.5) {
  rates <- as.vector(sapply(spikes, function(x) { ps.rate(x, start, end)}))
  return(list(m=mean(rates), s2=var(rates)))
}

# compute all mean rates for a given cell
cell.mean.rates <- function(spikes, start=0.5, end=1.5) {
  return(lapply(spikes, function(x) {ps.mean.rate(x,start, end)}))
}

mean.rates <- lapply(SpikesT, cell.mean.rates)

allodors <- unique(unlist(sapply(mean.rates,names)))
goododors <- allodors[1:36]
common.odors <- list(E2Hex="E2Hexenal", GerAc="Geranylacetat", Prpyl="Propylacetate", IPenA="Isopentanoic acid", Et3HB="Ethyl-3-hydroxybutyrate", MetSl="Methylsalicylat", PeEtA="Phenylethylacetat", HexAc="Hexanoic acid", AceAc="Acetic acid", "1HxOl"="1-Hexanol", BeZal="Benzaldehyd", Amnia="Ammonium hydroxide", bCitr="beta-Citrollenol", MtAct="Methylacetat", BtrAc="Butyric acid", PrpnA="Propionic acid")

total.trials <- sum(unlist(sapply(SpikesT, function(x) {sapply(x, length)})))
total.usable <- length(SpikesT)*36*4

make.total.data.matrix <- function(spikes, nr, nc) {
	all.mat <- matrix(NA, ncol=nc, nrow=nr)
	for (i in 1:length(spikes)) {
		labels <- c()
		n <- 1
		curr.labels <- names(spikes[[i]])
		for (j in 1:36) {
			for (k in 1:4) {				
				all.mat[i,n] <- ps.rate(spikes[[i]][[j]][[k]])
				labels <- c(labels, curr.labels[j])
				n <- n+1
			}
		}
	}
	colnames(all.mat) <- labels
	return(all.mat)
}
lhn.mat <- make.total.data.matrix(SpikesT, length(mean.rates), 36*4)
lhn.labels <- factor(colnames(all.mat))

binarize.mat.per.cell <- function(mat, labels, blank.names = c("OilBl", "WatBl")) {
	bl <- labels %in% blank.names
	num.nonblank <- sum(!bl)
	labels.nonblank <- labels[!bl]
	bin.mat <- matrix(0, nrow=nrow(mat), ncol=num.nonblank/4)
	colnames(bin.mat) <- unique(labels.nonblank)
	for (i in 1:nrow(mat)) {
		for (j in unique(labels.nonblank)) {
			if (sum(mat[i,bl]) > 0) {
				bin.mat[i,j] <- ifelse(t.test(mat[i,labels %in% j], mat[i,bl])$p.value < 0.05, 1, 0)
			} else {
				bin.mat[i,j] <- ifelse(sum(mat[i,labels %in% j]) > 0, 1, 0)
			}
		}
	}
	return(bin.mat)
}

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

# see correltion b/w ORN and PN firing rates and LHN firing rates
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
heatmap.2(lhn.odors.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors, trace='none')

# just heatmaps of common odors
heatmap.2(orn.common.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, Rowv=NA, Colv=NA, col=jet.colors, trace='none')
heatmap.2(pn.common.cor, distfun=function(x) as.dist(1-x), scale='none', symm=T, Rowv=NA, Colv=NA, col=jet.colors, trace='none')
heatmap.2(lhn.common.cor, distfun=function(x) as.dist(1-x), scale='none',Rowv=NA, Colv=NA, symm=T, col=jet.colors, trace='none')

# do factor analysis...

## Data Analysis

# Encoding: calculate pobabilities of response given odors

# Decoding: try to predict odors given population data
# -- Use total population data
# train on 3 examples for each cell and test on last
lhn.mat.noblank <- lhn.mat[,!(lhn.labels %in% c("OilBl", "WatBl"))]
leave.out <- seq(1,ncol(lhn.mat.noblank),4)
test.lhn <- lhn.mat.noblank[,leave.out]
train.lhn <- lhn.mat.noblank[, !(1:ncol(lhn.mat.noblank) %in% leave.out)]
model <- multinom(factor(colnames(train.lhn)) ~ train.lhn[1,])
# -- Compare with leaving out single neurons

# -- Compare with single neuron predictions
predictions <- c()
for (i in 1:nrow(train.lhn)) {
	model <- multinom(factor(colnames(train.lhn)) ~ train.lhn[i,])

}
