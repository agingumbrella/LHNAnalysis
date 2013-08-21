# Find all the neurons from FlyCircuit in the LH and cluster them
# Requires AnalysisSuite to be loaded and dps

GetLHPoints <- function(neuron, lhpoints) {
  return(subset(neuron, nn2(lhpoints, neuron$points, k=1)$nn.dists[,1] < 2))
}

ComputeLHNPairwiseDists <- function(neurons, use.just.lh=F, lhpoints=NA) {
  if (use.just.lh == T & is.na(lhpoints)) {
    stop('Need lhpoints')
  }
  dists <- matrix(0, nrow=length(neurons), ncol=length(neurons))
  rownames(dists) <- colnames(dists) <- names(neurons)                              
  for (i in 1:length(neurons)) {
  	print(i)
    if (use.just.lh) {
      n1 <- GetLHPoints(neurons[[names(neurons)[i]]], lhpoints)
    } else {
      n1 <- neurons[[names(neurons)[i]]]
    }
    for (j in 1:i) {
      if (use.just.lh) {
        n2 <- GetLHPoints(neurons[[names(neurons)[j]]], lhpoints)
      } else {
        n2 <- neurons[[j]]
      }
      if (length(n1$points) == 0 | length(n2$points) == 0) {
      	dists[i,j] <- 1
      } else {
	    dists[i,j] <- WeightedNNBasedLinesetMatching.dotprops(n1, n2)
	  }
    }
  }
  return(dists)
}

# make normalized distance matrix for apcluster
NormalizeMat <- function(m) {
	selfscores <- diag(m)
	m.norm < scale(m, scale=selfscores, cent=F)
	m.normsims <- as.dist.asym_matrix(m.norm)
	m.normsmat <- as.matrix(m.normsim)
	return(m.normsmat)  
}

PlotCluster <- function(neuronnames, lhpoints, show.brain=F, color.lh=F) {
	clear3d()
	if (show.brain) {
		fcwbsurf(alpha=0.1)
	}
	neuron.cols <- colorRampPalette(c("gray", "black"))(length(neuronnames))
	lh.cols <- colorRampPalette(c("white", "red"))(length(neuronnames))

	if (color.lh) {
		for (i in 1:length(neuronnames)) {
			curr.lhpoints <- GetLHPoints(dps[[neuronnames[i]]], lhpoints)
			plot3d(curr.lhpoints, col=lh.cols[i])
		}
	}

	for (i in 1:length(neuronnames)) {
		plot3dfc(neuronnames[i], col=neuron.cols[i], soma=T)
	}
	#par3d(zoom=0.5)
}

###################################################################
## Script to cluster LHNs

source("~/projects/AnalysisSuite/R/Code/Startup.R")
# include FlyCircuit scripts
source("~/projects/ChiangReanalysis/src/FlyCircuitStartup.R")

# load annotations
if(!exists("annotation")) load_fcdb("annotation")

# load dps data
if(!exists('dps')){
  message("Loading all dot properties (canonical side) ...")
  load_fcdata('dpscanon')
}

#Find all PN neurons from annotation db
pn.gns <- fc_gene_name(subset(annotation,annotation_class=='NeuronType' & grepl('PN',text))$neuron_idid)

# Find uniglomerular PN's
upn.gns <- fc_gene_name(subset(annotation,annotation_class=='NeuronSubType' & grepl('uPN',text))$neuron_idid)

# Find all PNs
all.pns=c(pn.gns,upn.gns)

# load labels and identify left and right LH regions
fcwbjfrclabels <- Read3DDensityFromNrrd("FCWB_JFRC_labelswarp-04-4xd.nrrd",ReadByteAsRaw='none')
# get points from right-hand side
rlhpoints <- ind2coord(fcwbjfrclabels==26, voxdims=voxdim.gjdens(fcwbjfrclabels))
# get points from left-hand side
llhpoints <- ind2coord(fcwbjfrclabels==56, voxdims=voxdim.gjdens(fcwbjfrclabels))


# example of subsetting dotprops representation
#dps18.lh <- subset(dps[[18]],nn2(llhpoints,dps[[18]]$points,k=1)$nn.dists[,1]<2)
#plot3d(dps18.lh)
#plot3d(dps[[18]],col='red')

# find neurons with more than 200 points in left or right LH
load_fcdb("spatdist_jfrc")
llh <- subset(spatdist_jfrc,LH_L>200)
rlh <- subset(spatdist_jfrc,LH_R>200)

# remove PNs from lll and llr
llh <- subset(llh, !(rownames(llh) %in% all.pns))
rlh <- subset(rlh, !(rownames(rlh) %in% all.pns))

# find neurons in clustering
load_fcdata('apres15kv2ns.p02')
apdf <- as.data.frame(apres15kv2ns.p02)
exemplars.llh <- unique(subset(apdf, item %in% rownames(llh))$exemplar)
exemplars.rlh <- unique(subset(apdf, item %in% rownames(rlh))$exemplar)

# read distance matrix 
abc2 <- attach.big.matrix(file.path(fcconfig$bigmatrixdir,"abc2.normdmat.desc"))

# get dotprops for each
llh.dps <-dps[intersect(names(dps),rownames(llh))]
rlh.dps <-dps[intersect(names(dps),rownames(rlh))]

# calculate distances for all points in all putative LHNs
#rdists.all <- plugindist.list(rlh.dps, WeightedNNBasedLinesetMatching.dotprops,upper=TRUE,processfun=NULL,sd=3)
#ldists.all <- plugindist.list(llh.dps, WeightedNNBasedLinesetMatching.dotprops,upper=TRUE,processfun=NULL,sd=3)
#save(rdists.all, ldists.all, file="data/lhdists.rda")

#load("data/lhdists.rda")
#x = as.matrix(llh.dps)
rdists.lh <- ComputeLHNPairwiseDists(rlh.dps, use.just.lh=T, lhpoints=llhpoints)
#save(rdists.lh, file="~/Desktop/Will/rdists.rda")
#ldists.lh <- ComputeLHNPairwiseDists(llh.dps, use.just.lh=T, lhpoints=llhpoints)

# open 3d viewer
open3d("white")


# basic clustering


# cluster all dists and LH dists using apclusterr
require(apcluster)

# get just distances for putative LHNs
llh.good <- names(llh.dps)[names(llh.dps) %in% rownames(allbyallblast.canon)]
rlh.good <- names(rlh.dps)[names(rlh.dps) %in% rownames(allbyallblast.canon)]
llh.dists <- abc2[llh.good, llh.good]
rlh.dists <- abc2[rlh.good, rlh.good]

# save clusterings by whole cell
num.clusts <- 49
hc.lhn.r <- hclust(as.dist(rlh.dists), meth='ward')
hc.lhn.l <- hclust(as.dist(llh.dists), meth='ward')
llh.groups <- cutree(hc.lhn.l, k=num.clusts)
rlh.groups <- cutree(hc.lhn.r, k=num.clusts)
#plot(hc.lhn.r)

# change this to where ever you want your images saved
#setwd("~/Desktop/Will/lhclusts/")
frontalView(0.4)
for (i in 1:num.clusts) {
	PlotCluster(rownames(llh.dists)[llh.groups == i], llhpoints, color.lh=T, show.brain=T)
	snapshot3d(paste("llh_cluster_",i,".png",sep=""))
}

# do clustering by just area in LH

# change to save images somewhwere
#setwd("~/Desktop/Will/lhclusts_lhonly/")
hc.rdists <- hclust(as.dist(1-rdists.lh), meth='ward')
rdists.groups <- cutree(hc.lhn.r, k=num.clusts)
frontalView(0.4)
for (i in 1:num.clusts) {
	PlotCluster(rownames(rdists.lh)[rdists.groups == i], llhpoints, color.lh=T, show.brain=T)
	snapshot3d(paste("llh_cluster_",i,".png",sep=""))
}

# redo everything with apcluster
#setwd("~/Desktop/Will/lhclusts_apclust/")
clusts <- apcluster(-rlh.dists, q=0.1)
frontalView(0.5)
for (i in 1:length(clusts)) {
	PlotCluster(clusts[[i]], llhpoints, color.lh=T, show.brain=T)
	snapshot3d(paste("llh_cluster_",i,".png",sep=""))
}

#setwd("~/Desktop/Will/lhclusts_lhonly_apclust/")
clusts <- apcluster(as.matrix(1-rdists.lh), p=1, q=0.1,details=T)
frontalView(0.5)
for (i in 1:length(clusts)) {
	PlotCluster(clusts[[i]], llhpoints, color.lh=T, show.brain=T)
	snapshot3d(paste("llh_cluster_",i,".png",sep=""))
}

