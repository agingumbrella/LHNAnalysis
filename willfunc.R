# stuff to deal with 2P data
# NOT RELEVANT FOR SHAHAR'S DATA (ALSO KIND OF A MESS)

setwd("~/Documents/Neuroscience/jefferis_lab/LHNAnalysis")
source("Avg_IgorImport.R")
source("PSTH_FUNC.R")
source("ClusteringFunctions.R")

load('~/Documents/Neuroscience/jefferis_lab/shahar_data/RandIgor/spike_summary_old.rda')

## Load GCaMP stuff
gcamp.data.path <- "~/Documents/Neuroscience/jefferis_lab/splitgal4_data/data/"
gcamp.odors <- c("OilBl", "E2Hex", "GerAc", "IPenA", "Et3HB", "Prpyl", "HexAc", "PeEOL", "OilBl", "AceAc", "EtHex", "HexnA", "MelSl", "bCitr", "1HxOl", "BeZal", "WatBl", "MtAct", "AcAcd", "PrpnA", "Amnia", "Vingr", "EtAmn", "Lnlol", "WatBl", "23BTD", "APnn", "4EtGu", "BtLac", "AcePn", "1OcOl", "1O3Ol")
delta.L4 <- 0.65
delta.L5 <- 2

l4.means <- structure(c(1.88333333333333, 3.81666666666667, 9.41666666666667, 16.4, 5.71666666666667, 13.1333333333333, 2.85, 10.8166666666667, 4.55, 10.6666666666667, 4.03333333333333, 4.01666666666667, 0.683333333333333, 13.3, 4.41666666666667, 5.91666666666667, 4.13333333333333, 6.43333333333333), .Names = c("OilBl", "E2Hex", "GerAc", "Prpyl", "IPenA", "Et3HB", "HexAc", "AceAc", "EtHex", "BeZal", "bCitr", "1HxOl", "WatBl", "MtAct", "AcAcd", "PrpnA", "Amnia", "Vingr"))
l5.means <- structure(c(1.05, 3.15, 9.3, 23.25, 14.7, 26.65, 8.6, 18, 9.35, 10.3, 5.25, 9, 0.3, 23.25, 7.65, 8.25, 7.75, 24.25), .Names = c("OilBl", "E2Hex", "GerAc", "Prpyl", "IPenA", "Et3HB", "HexAc", "AceAc", "EtHex", "BeZal", "bCitr", "1HxOl", "WatBl", "MtAct", "AcAcd", "PrpnA", "Amnia", "Vingr"))

gcamp.files <- dir(gcamp.data.path, "*.txt")
read.gcamp.file <- function(name) {
	return(read.table(paste(gcamp.data.path, name, sep="/"))$V2)
}

# function to parse out the odor response from the time series
# and compute DF/F
compute.odor.response <- function(data, odors, delta=2, stim.len=30, df.len=8) {	
	m <- matrix(0, stim.len, length(odors))
	colnames(m) <- odors
	for (i in 1:length(odors)) {
		curr.start <- (i-1)*stim.len + floor(delta*(i-1)) + 1
		curr.stop <- curr.start+(stim.len-1)
		curr.data <- data[curr.start:curr.stop]
		f0 = mean(curr.data[1:df.len])
		dfof = 100*(curr.data/f0 - 1)
		m[,i] <- dfof
	}
	m
}

stacked.plot <- function(data, plot.name="",scale=6, col='red') {
	spacing <- 5*sd(as.vector(data), na.rm=T)	
	max.data <- max(as.vector(data))
	min.data <- min(as.vector(data))
	num.total <- ncol(data)
	for (i in 1:num.total) {
		if (i < num.total) {
			plot(1:nrow(data), data[,i] - spacing*i, type='l', ylim=c(-num.total*spacing,0), xlim=c(1,30), xaxt='n', yaxt='n', xlab="", ylab="",frame=F, col=col) 
		} else {
			plot(1:nrow(data), data[,i] - spacing*i, type='l', ylim=c(-num.total*spacing,0), xlim=c(1,30), yaxt='n', xlab="Frame Number", ylab="",frame=F, col=col) 

		}
		text(7,(spacing/2)-spacing*i, gcamp.odors[i],1,cex=0.6)
		par(new=T)

	}
	title(main=plot.name)

}

save.data.plot <- function(name, path) {
	data <- read.table(paste(gcamp.data.path, paste(path,".txt",sep=""),sep="/"))$V2
	data <- compute.odor.response(data, gcamp.odors)
	pdf(paste(path,".pdf",sep=""), width=5, height=10)
	stacked.plot(data, name)
	dev.off()
}

gcamp.old <- c("L4_fly22_1", "L4_fly22_2", "L4_fly22_3", "L4_fly22_4")
gcamp.new <- c("L4_fly24_1", "L4_fly24_2", "L4_fly24_3", "L4_fly24_4", "L5_fly4_1", "L5_fly4_2", "L5_fly4_3", "L5_fly4_4", "L5_fly4_5", "L5_fly5_1", "L5_fly5_2", "L5_fly5_3")

save.average.data.plot <- function(name, paths, col='red', pad=F) {
	all.data <- c()
	for (p in paths) {
		if (p %in% gcamp.old) {
			data <- compute.odor.response(read.gcamp.file(paste(p, ".txt", sep="")),gcamp.odors,0.65)
		} else {
			data <- compute.odor.response(read.gcamp.file(paste(p,".txt",sep="")),gcamp.odors)
		}
		if (length(all.data) == 0) {
			all.data <- data
		} else {
			all.data <- all.data + data
		}
	}


	if (pad) {
		all.data[is.na(all.data[,"1O3Ol"]), "1O3Ol"] <- 0
	}

	data <- all.data/length(paths)
#	pdf(paste(path,".pdf",sep=""), width=5, height=10)
	stacked.plot(data, name, col=col)
#	dev.off()
	data
}

all.l4.dendrites <- c("L4_fly22_3","L4_fly22_4", "L4_fly24_2", "L4_fly24_3")
all.l4.axons <- c("L4_fly22_1", "L4_fly22_2","L4_fly24_1","L4_fly24_4")
all.l5.dendrites <- c("L5_fly5_2","L5_fly5_3","L5_fly4_1","L5_fly4_2")
all.l5.axons <- c("L5_fly4_3","L5_fly4_4","L5_fly4_5","L5_fly5_1")

pdf("gcamp_response2.pdf", width=4, height=5)
par(mfrow=c(1,2))
par(mar=c(5,0,3,0.5))
l4.den.data <- save.average.data.plot("", all.l4.dendrites, 'red')
par(new=T)
l4.axon.data <- save.average.data.plot("", all.l4.axons, 'blue', pad=T)
title("JK671-SF274")
#legend("topright", legend=c("Dendrite", "Axon"),fill=c("red", "blue"), border="white",box.col='white')
l5.den.data <- save.average.data.plot("", all.l5.dendrites, 'red')
par(new=T)
l5.axon.data <- save.average.data.plot("", all.l5.axons, 'blue')
title("JK304-SF274")
dev.off()

plot(apply(l4.den.data,2, max), apply(l4.axon.data, 2, max),pch=20,col='red', xlim=c(0,1000), ylim=c(0,1000))
par(new=T)
plot(apply(l5.den.data,2, max), apply(l5.axon.data, 2, max),pch=20,col='blue', xlim=c(0,1000), ylim=c(0,1000))

plot(l4.den.data, l4.axon.data, pch=20,col='red', xlim=c(0,1000), ylim=c(0,1000))
par(new=T)
plot(l5.den.data,l5.axon.data,pch=20,col='blue', xlim=c(0,1000), ylim=c(0,1000))

#save.data.plot("Line 5, Fly 4, Dendrite", "L5_fly4_1")
#save.data.plot("Line 5, Fly 4, Dendrite", "L5_fly4_2")
#save.data.plot("Line 5, Fly 5, Dendrite", "L5_fly5_2")
#save.data.plot("Line 5, Fly 5, Dendrite", "L5_fly5_3")
#save.data.plot("Line 5, Fly 5, Axon", "L5_fly5_1")
#save.data.plot("Line 4, Fly 24, Axon", "L4_fly24_4")
#save.data.plot("Line 4, Fly 22, Dendrite", "L4_fly22_3")

l4.den.max <- apply(l4.den.data, 2, max)
l4.axon.max <- apply(l4.axon.data, 2, max)
l4.axon.max["1O3Ol"] <- max(l4.axon.data[1:10,"1O3Ol"])
l5.den.max <- apply(l5.den.data, 2, max)
l5.axon.max <- apply(l5.axon.data, 2, max)


l4.den.noblank <- l4.den.max[!(names(l4.den.max) %in% c("OilBl", "WatBl"))]
l5.den.noblank <- l5.den.max[!(names(l5.den.max) %in% c("OilBl", "WatBl"))]
l4.axon.noblank <- l4.axon.max[!(names(l4.axon.max) %in% c("OilBl", "WatBl"))]
l5.axon.noblank <- l5.axon.max[!(names(l4.axon.max) %in% c("OilBl", "WatBl"))]

gcamp = list()
for (g in gcamp.old) {
	gcamp[[g]] <- compute.odor.response(read.gcamp.file(paste(g, ".txt", sep="")),gcamp.odors,0.65)
}
for (g in gcamp.new) {
	gcamp[[g]] <- compute.odor.response(read.gcamp.file(paste(g,".txt",sep="")),gcamp.odors)
}

max.dfof <- matrix(0, length(gcamp), length(gcamp.odors))
rownames(max.dfof) <- names(gcamp)
colnames(max.dfof) <- gcamp.odors
for (i in 1:length(gcamp)) {
	max.dfof[i,] <- apply(gcamp[[i]],2,max)
}

axons <- c("L4_fly22_1", "L4_fly22_2","L4_fly24_1","L4_fly24_4","L5_fly4_3","L5_fly4_4","L5_fly4_5","L5_fly5_1")
dendrites <- c("L4_fly22_3","L4_fly22_4", "L4_fly24_2", "L4_fly24_3","L5_fly5_2","L5_fly5_3","L5_fly4_1","L5_fly4_2")
heatmap(max.dfof[axons,],scale="row",symm=F)
heatmap(max.dfof[dendrites,],col=jet.colors(20),scale="row",symm=F)

heatmap(cor(t(max.dfof[axons,]),use='complete.obs'), symm=T, col=jet.colors(20))
heatmap(cor(t(max.dfof[dendrites,]),use='complete.obs'), symm=T, col=jet.colors(20))

shared.odors <- colnames(max.freqs)[colnames(max.freqs) %in% gcamp.odors] 
heatmap(max.freqs[,colnames(max.freqs) %in% gcamp.odors], scale='none', Colv=NA,col=jet.colors(20))
heatmap(max.dfof[axons,shared.odors], scale='row', Colv=NA,col=jet.colors(20))

# find highest correlated
m <- matrix(0, length(axons), nrow(max.freqs)) 
for (i in axons) {
	print(i)
	for (j in 1:nrow(max.freqs)) {
		m[i,j] <- cor(max.dfof[i,], max.freqs[j,])
	}
}


# load each of the line 5 data
## Load ephys stuff
#
allodours=unique(unlist(sapply(smSpikes,names)))

# > dput(allodours)
# c("OilBl", "E2Hex", "GerAc", "Prpyl", "IPenA", "Et3HB", "Nonnl", 
# "CiVAc", "MetSl", "HexAc", "PeEtA", "AceAc", "EtHex", "2PnAc", 
# "5OdMx", "BeZal", "bCitr", "1HxOl", "Frnsl", "WatBl", "Cdvrn", 
# "Sprmn", "Acoin", "MtAct", "AcAcd", "PrpnA", "BtrAc", "Amnia", 
# "Pyrdn", "PAcHd", "HCL36", "PAcAc", "Vingr", "Geosn", "VinGe", 
# "PEtAm", "ClrBL", "ClrB2", "FlyFM", "EtAmn", "MtAmn", "Ptscn", 
# "Lnlol", "23BTD", "Sprdn")

# calculate frequencies for every odour for every cell
# returns a list of cells with one matrix for each cell
# where columns are odours and rows are timepoints
allfreqs=sapply(smSpikes,function(psthsforcell) sapply(psthsforcell,function(x) x$freq))

# pad those frequencies with columns of NAs for missing odours
# also reorder odours into the order given by allodours
allfreqs_allodours=lapply(allfreqs,addnacols,allodours)

# same thing but as a big matrix with columns for cells
# and rows for odour*time
allfreqs_allodours_array=sapply(allfreqs,addnacols,allodours)
image(allfreqs_allodours_array)

# these are the ones for which we have data for all cells
goododours=allodours[1:36]

# 
allfreqs_goododours=lapply(allfreqs_allodours,function(x) x[,goododours])
allfreqs_goododours_mat=sapply(allfreqs_allodours,function(x) x[,goododours])
physplit=read.table("~/Documents/Neuroscience/jefferis_lab/LHNAnalysis/data/PhySplitSimple.mer",sep=',',header=TRUE,stringsAsFactors=FALSE)

# find crosses for all those cells that we have plotted
jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))
crosses=physplit$cross[match(colnames(allfreqs_goododours_mat),physplit$cell)]
heatmap(allfreqs_goododours_mat,scale='none',labCol=crosses,col=jet.colors(20),margins=c(6,5))

# correlation score instead
# find correlation between all pairs of columns (ie cells)
spcor.cells=cor(allfreqs_goododours_mat,use='complete.obs')

# dendrogram is based on distance of 1-correlation score, but the 
# colours in the heatmap are still the correlation scores (ie hot is highly correlated)
heatmap(spcor.cells,distfun=function(x) as.dist(1-x),scale='none',symm=T,col=jet.colors(20))

#allcor <- cor(allfreqs_allodours_array, use='complete.obs')
#heatmap(allcor, disfun=function(x) as.dist(1-x), scale='none', symm=T, col=jet.colors(20))

# correlation score instead
# find correlation between all pairs of rows (ie odors)
spcor.odors=cor(t(allfreqs_goododours_mat))

# dendrogram is based on distance of 1-correlation score, but the 
# colours in the heatmap are still the correlation scores (ie hot is highly correlated)
heatmap(spcor.odors,distfun=function(x) as.dist(1-x),scale='none',symm=T,col=jet.colors(20))

# finds the maximum firing instantaneous firing rate in response to each odor
# for a given cell
all.max.freqs <- function(freqs) {
	num.cells <- length(freqs)
	m <- matrix(0, num.cells, dim(freqs[[1]])[2])
	colnames(m) <- colnames(freqs[[1]])
	rownames(m) <- names(freqs)
	for (i in 1:num.cells) {
		 m[i, ] <- apply(freqs[[i]],2,max)
	}
	m
}
max.freqs <- all.max.freqs(allfreqs_goododours)
heatmap(max.freqs, scale='none', Colv=NA,col=jet.colors(20))


plot.corr.heatmap <- function(mat) {
	correlation <- cor(mat, use='complete.obs')
	heatmap(correlation,distfun=function(x) as.dist(1-x),scale='none',symm=T,col=jet.colors(20))
}

plot.corr.heatmap(max.freqs)


## compare with gcamp
l4.cells <- sapply(strsplit(physplit[physplit$cross %in% "SF274-JK671",]$Igor.file, split="/"), function(x) x[length(x)])
l5.cells <- sapply(strsplit(physplit[physplit$cross %in% "SF274-JK304",]$Igor.file, split="/"), function(x) x[length(x)])
l4.freqs.gcamp <- apply(max.freqs[rownames(max.freqs) %in% l4.cells, colnames(max.freqs) %in% names(l4.axon.noblank)], 2, mean)
l4.freqs.gcamp <- l4.freqs.gcamp/max(l4.freqs.gcamp)
l5.freqs.gcamp <- apply(max.freqs[rownames(max.freqs) %in% l5.cells, colnames(max.freqs) %in% names(l5.axon.noblank)], 2, mean)
l5.freqs.gcamp <- l5.freqs.gcamp/max(l5.freqs.gcamp)

common.names <- sort(names(l4.den.compare))

l4.axon.compare <- l4.axon.noblank[names(l4.axon.noblank) %in% names(l4.freqs.gcamp)]
l4.axon.compare <- l4.axon.compare/max(l4.axon.compare)
l4.den.compare <- l4.den.noblank[names(l4.axon.noblank) %in% names(l4.freqs.gcamp)]
l4.den.compare <- l4.den.compare/max(l4.den.compare)
l5.axon.compare <- l5.axon.noblank[names(l4.axon.noblank) %in% names(l4.freqs.gcamp)]
l5.axon.compare <- l5.axon.compare/max(l5.axon.compare)
l5.den.compare <- l5.den.noblank[names(l4.axon.noblank) %in% names(l4.freqs.gcamp)]
l5.den.compare <- l5.den.compare/max(l5.den.compare)

cor(l4.axon.noblank[common.names], l4.freqs.gcamp[common.names])
cor(l4.den.noblank[common.names], l4.freqs.gcamp[common.names])
cor(l5.axon.noblank[common.names], l5.freqs.gcamp[common.names])
cor(l5.den.noblank[common.names], l5.freqs.gcamp[common.names])

pdf("l4_gcamp.pdf", width=4,height=4)
l4.all <- rbind(l4.den.compare, l4.axon.compare, l4.freqs.gcamp)
rownames(l4.all) <- c("GCaMP Dendrite", "GCaMP Axon", "Patch")
heatmap(l4.all, col=jet.colors(100), Rowv=NA, main="JK671-SF274", cexRow=0.9)
dev.off()

pdf("l5_gcamp.pdf", width=4, height=4)
l5.all <- rbind(l5.den.compare, l5.axon.compare, l5.freqs.gcamp)
rownames(l5.all) <- c("GCaMP Dendrite", "GCaMP Axon", "Patch")
heatmap(l5.all, col=jet.colors(100), Rowv=NA,main="JK304-SF274", cexRow=0.9)
dev.off()

## DO STUFF WITH BINARIZED 
binarize <- function(x) {
	if (x > 2)	{
		return(1)
	} else {
		return(0)
	}
}

# binarize the responses 
binarize.responses <- function(freqs) {
	m <- matrix(0, nrow(freqs), ncol(freqs))
	colnames(m) <- colnames(freqs)
	rownames(m) <- rownames(freqs)
	for (i in 1:nrow(freqs)) {
		m[i,] <- sapply(freqs[i,], binarize)
	}
	m
}
binarized <- binarize.responses(max.freqs.nona)

spcor.bin.odors=cor(binarized,use='complete.obs')

# dendrogram is based on distance of 1-correlation score, but the 
# colours in the heatmap are still the correlation scores (ie hot is highly correlated)
heatmap(spcor.bin.odors,distfun=function(x) as.dist(1-x),scale='none',symm=T,col=jet.colors(20))

