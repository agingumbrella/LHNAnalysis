################################################################################
## General setup stuff -- define colors used in plots

jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))(100)
bin.colors <- colorRampPalette(c('white', 'black'))(100)
bin.rev.colors <- colorRampPalette(c('black', 'white'))(100)

################################################################################
## Load LHN data either from spike summary or repeated spike train 

#load("~/Documents/Neuroscience/jefferis_lab/shahar_data/RandIgor/spike_summary_133.rda")
## LOAD LHN SPIKE DATA
#badcells=c("nm20130329c0","nm20130206c1", "nm20130606c0")
#SpikesT <- list()
#for (x in names(Spikes)) {
#  if (x %in% badcells) {
#  } else {
#    SpikesT[[x]] <- as.repeatedTrain(Spikes[[x]][[1]])
#  }
#}
#save(SpikesT, file="spikeT_summer.rda")

# load spike summaries
load("spikeT_summer.rda")

  
################################################################################
## load cell info from database

physplit=read.table("PhySplitSimple.mer",sep=',',header=TRUE,stringsAsFactors=FALSE)

# find cross ids for each cell
cross.idx <- sapply(names(SpikesT), function(x) grep(x, physplit$Igor.file))
cross.ids <- filter.cross.ids(sapply(cross.idx, function(x) physplit$cross[x]))

################################################################################
## Save colors used throughout analysis for crosses and odor types

# types are: none, terpenes, alcohols, ester, ketones, lactones, acids, aromatics, aldehydes, sulfur compounds, amines
# color for each odor type
type.colors <- c(amine="gold", lactone="darkblue", acid="pink", sulfur="black", terpene="lightgreen", aldehyde="gray", ketone="yellow", aromatic="lightblue", alcohol="red", ester="darkgreen")

# color for each cross id
#cross.colors <- colorRampPalette(c("cyan", "magenta", "yellow", "black"))(length(unique(cross.ids)))
cross.colors <- rainbow(length(unique(cross.ids)))
names(cross.colors) <- unique(cross.ids)

################################################################################
## Preprocessing of LHN data

# compute average firing rates for each recorded LHN
lhn.mean.rates <- lapply(SpikesT, cell.mean.rates)

## get odors that overlap with Hallem and Carlson odors
allodors <- unique(unlist(sapply(lhn.mean.rates,names)))
# odors that are used in all cells
goododors <- allodors[1:36]
# odors common to ORN/PN and LHN datasets
common.odors <- list(E2Hex="E2Hexenal", GerAc="Geranylacetat", Prpyl="Propylacetate", IPenA="isopentyl acetate", Et3HB="Ethyl-3-hydroxybutyrate", MetSl="Methylsalicylat", PeEtA="Phenylethyl alcohol", EtHex="Ethylhexanoate", BeZal="benzaldehyde", bCitr="b-citronellol", `1HxOl`="1-Hexanol", Cdvrn="cadaverine", MtAct="Methylacetat", AcAcd="Acetic acid", PrpnA="Propionic acid", BtrAc="Butyric acid", Amnia="ammonium hydroxide", PAcHd="Phenylacetaldehyd")


# color for each class -- this was a pain in the ass to do by hand in R; probably should've used excel
lhn.odors.types <- data.frame(Odorant=c("OilBl", "E2Hex", "GerAc", "Prpyl", "IPenA", "Et3HB", "Nonnl", 
                                       "CiVAc", "MetSl", "HexAc", "PeEtA", "AceAc", "EtHex", "2PnAc", 
                                       "5OdMx", "BeZal", "bCitr", "1HxOl", "Frnsl", "WatBl", "Cdvrn", 
                                       "Sprmn", "Acoin", "MtAct", "AcAcd", "PrpnA", "BtrAc", "Amnia", 
                                       "Pyrdn", "PAcHd", "HCL36", "PAcAc", "Vingr", "Geosn", "VinGe", 
                                       "PEtAm"),
                              Type=c("none", "aldehyde", "ester", "ester", "ester", "ester", "aldehyde",
                                "ester", "aromatic", "ester", "alcohol", "ketone", "ester", "ester",
                                "ester", "aromatic", "terpene", "alcohol", "terpene", "none", "amine",
                                "amine", "ketone", "ester", "acid", "acid", "acid", "amine", "amine",
                                "aromatic", "acid", "acid", "acid", "alcohol", "acid", "amine"))
rownames(lhn.odors.types) <- lhn.odors.types$Odorant

# figure out toal number of usable trials
total.trials <- sum(unlist(sapply(SpikesT, function(x) {sapply(x, length)})))
total.usable <- length(SpikesT)*36*4 # assuming 36 odors and 4 repeats per odor

# make matrix of all usable trials
lhn.mat <- make.total.data.matrix(SpikesT, length(lhn.mean.rates), 36*4)
# assign lhn.mat to be crosses rather than cell names
rownames(lhn.mat) <- cross.ids #names(SpikesT)

# get labels of each odor
lhn.labels <- factor(colnames(lhn.mat))
# make matrix for each cell of whether fired or not in response to an odor, using a statistical test
lhn.bin <- binarize.mat.per.cell(lhn.mat, lhn.labels)

# make per cell average firing rates
lhn.rates.mat <- make.rate.per.cell(lhn.mean.rates, goododors)

# make matrix of rates without NAs
lhn.rates.nona.mat <- lhn.rates.mat[,apply(lhn.rates.mat, 2, function(y) !any(is.na(y)))]

# save the Igor names of each of the LHN cells
lhn.cell.names <- colnames(lhn.rates.nona.mat)

# finally, lhn is all the average firing rates without any NAs for just the odors that are common to all trials
lhn <- lhn.rates.nona.mat

# replace colnames with cross ids
colnames(lhn) <- cross.ids[apply(lhn.rates.mat, 2, function(y) !any(is.na(y)))]

# make matrix of delta firing rates from blank (without blanks)
lhn.spontaneous.rates <- apply(lhn[rownames(lhn) %in% c("OilBl","WatBl"), ], 2, mean)
lhn.delta.rates <- lhn[!(rownames(lhn) %in% c("OilBl", "WatBl")),] - rep.row(lhn.spontaneous.rates, nrow(lhn[!(rownames(lhn) %in% c("OilBl", "WatBl")),]))

################################################################################
## Load Hallem and Carlson data
x.names <- c("X2a", "X7a", "X9a", "X10a", "X19a", "X22a", "X23a", "X33b", "X35a", "X43a", "X43b", "X47a", "X47b", "X49b", "X59b", "X65a", "X67a", "X67c", "X82a", "X85a", "X85b", "X85f", "X88a", "X98a")
orn.names <- c("Or2a", "Or7a", "Or9a", "Or10a", "Or19a", "Or22a", "Or23a", "Or33b", "Or35a", "Or43a", "Or43b", "Or47a", "Or47b", "Or49b", "Or59b", "Or65a", "Or67a", "Or67c", "Or82a", "Or85a", "Or85b", "Or85f", "Or88a", "Or98a")

# spontaneous firing rates of ORN
spontaneous.rates <- c(8, 17, 3, 14, 29, 4, 9, 25, 17, 21, 2, 1, 47, 8, 2, 18, 11, 6, 16, 14, 13, 7, 26, 12) # ordered by OR; from table in supplement of paper

# read valences, odors and types
val.and.odors <- read.csv("val_orn.csv", header=T)
rownames(val.and.odors) <- val.and.odors$Odorant

# Save valences
val <- val.and.odors[,c("Odorant", "Attraction.index")]

# Save odor types
all.odor.type <- val.and.odors[, c("Odorant", "Type")]
rownames(all.odor.type) <- all.odor.type$Odorant

# get ORNs
orn <- as.matrix(val.and.odors[,x.names])
colnames(orn) <- orn.names
rownames(orn) <- val.and.odors$Odorant

# add baseline firing rates and make all negative rates 0
orn <- apply(orn, 1, function(x) {x + spontaneous.rates})
orn[orn < 0] <- 0
orn <- t(orn)

# make pn rates with and without lateral inhibition (by model)
nl <- make.pn.rates(orn,m=0)
pn <- make.pn.rates(orn)

# save common data in separate matrices
val.common <- val[val$Odorant %in% unlist(common.odors),2]
orn.common <- orn[rownames(orn) %in% unlist(common.odors),]
pn.common <- pn[rownames(pn) %in% unlist(common.odors),]
lhn.common <- lhn[rownames(lhn) %in% unlist(names(common.odors)),]

################################################################################
## Load anatomical data

# get list of traced neurons
traced <- dir("~/Documents/Neuroscience/jefferis_lab/LHNAnalysis/Tracing.IS2/", full.name=T)
neurons <- as.neuronlist(lapply(traced, function(x) read.neuron(x)))
#neuron.dists <- matrix(0, nrow=length(neurons), ncol=length(neurons))

# TODO Include code for distance measurements

# load precomputed distance matrix m
load("~/Documents/Neuroscience/jefferis_lab/LHNAnalysis/dists15.Rdata")
lhn.anatomy.dists <- m

# corresponding cells to the anatomical structures3
# had to do by hand b/c use different naming formats -- traced file name includes both date and cross id
# and because the names of the files included date strings and didn't exactly match the cross IDs in the ephys data...
traced.info <- read.table("tracedinfo.txt", header=T)

rownames(lhn.anatomy.dists) <- colnames(lhn.anatomy.dists) <- traced.info$Cross

num.anatomy.cross.ids <- length(unique(rownames(lhn.anatomy.dists)))

