require(RColorBrewer)
require(NMF)
require(gphys)
require(gplots)
require(glmnet)
require(MASS)
require(penalized)
require(scatterplot3d)
library(pvclust)
require(mclust)

# load AnalysisSuite stuff
if(!exists("fcconfig")) {
	source("~/projects/ChiangReanalysis/src/FlyCircuitStartup.R")
	source("~/projects/AnalysisSuite/R/Code/Startup.R")
}

# change back to whatever directory you started out in, because loading AnalysisSuite annoyingly changes your directory...
# set this to whatever you want
setwd("~/Documents/Neuroscience/jefferis_lab/LHNAnalysis")

# load local functions
# for dealing with ephys data
source("Avg_IgorImport.R")
source("PSTH_FUNC.R")
source("ClusteringFunctions.R")

# load functions for analysis
source("infotheory.R")
source("decoding.R")
source("modeling.R")
source("encoding.R")
source("preprocessing.R")
source("util.R")

# load all the data needed for subsequent analysis
source("loaddata.R")

source("figscripts/clustering.R")
source("figscripts/modeling.R")
source("figscripts/prediction.R")
source("figscripts/responses.R")
source("figscripts/valence.R")
source("figscripts/fillanatomy.R")


