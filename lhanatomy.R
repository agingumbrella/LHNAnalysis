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
    if (use.just.lh) {
      n1 <- GetLHPoints(neurons[[i]], lhpoints)
    } else {
      n1 <- neurons[[i]]
    }
    for (j in 1:length(neurons)) {
      if (use.just.lh) {
        n2 <- GetLHPoints(neurons[[j]], lhpoints)
      } else {
        n2 <- neurons[[j]]
      }
      dists[i,j] <- WeightedNNBasedLinesetMatching.dotprops(n1, n2)
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

###################################################################
## Script to cluster LHNs

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
fcwbjfrclabels <- Read3DDensityFromNrrd("data/FCWB_JFRC_labelswarp-04-4xd.nrrd",ReadByteAsRaw='none')
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
exemplar.rlh <- unique(subset(apdf, item %in% rownames(rlh))$exemplar)

# get dotprops for each
llh.dps <-dps[intersect(names(dps),rownames(llh))]
rlh.dps <-dps[intersect(names(dps),rownames(rlh))]

# calculate distances for all points in all putative LHNs
#rdists.all <- plugindist.list(rlh.dps, WeightedNNBasedLinesetMatching.dotprops,upper=TRUE,processfun=NULL,sd=3)
#ldists.all <- plugindist.list(llh.dps, WeightedNNBasedLinesetMatching.dotprops,upper=TRUE,processfun=NULL,sd=3)
#save(rdists.all, ldists.all, file="data/lhdists.rda")

load("data/lhdists.rda")
x = as.matrix(llh.dps)
rdists.lh <- ComputeLHNPairwiseDists(rlh.dps, use.just.lh=T, lhpoints=rlhpoints)
ldists.lh <- ComputeLHNPairwiseDists(llh.dps, use.just.lh=T, lhpoints=llhpoints)

# open 3d viewer
open3d("white")


# basic clustering
hc.lhn.r <- hclust(as.dist(rdists.all), meth='ward')
#plot(hc.lhn.r)
groups<-cutree(hc.lhn.r,k=20)


# cluster all dists and LH dists using apclusterr
require(apcluster)
rdists.cluster <- apcluster(-rdists.all, details=TRUE)
