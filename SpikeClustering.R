# Making a start on clustering subthreshold membrane potential responses to odours

ShaharDataRoot="/Volumes/JData/JPeople/Shahar/Data"

physplit=read.table("/Volumes/JData/JPeople/Shahar/Database/data/PhySplitSimple.mer",sep=',',header=TRUE,stringsAsFactors=FALSE)
# calculate name of cell from Igor file column
physplit$cell=basename(physplit$Igor.file)
#load("~/projects/Shahar/RandIgor/physplit.rda")
# Just keep cells where we have at least 5 good sputniks
badcells=c("nm20130329c0","nm20130206c1")
ps=subset(physplit,nchar(Spike_Analyze)>=14 & !cell%in%badcells)

Spikes=list()
for (i in seq(nrow(ps))) 
{
  igorpath=ps$Igor.file[i]
  igorpath=sub("^.*Data",ShaharDataRoot,igorpath)
  # name of the cell is the folder name without preceding path
  cell=basename(igorpath)
  selsweeps=as.integer(unlist(strsplit(ps$Spike_Analyze[i],",")))
  message("Reading spikes for cell: ",cell)
  # try reading spikes
  res=try(CollectSpikesShaharCell(igorpath,filenums=selsweeps))
  # if it worked then store them
  if(!inherits(res,'try-error'))
    Spikes[[cell]]=res
  # sputnum=sputniknumber(igorpath,sweep)
}
smSpikes=lapply(Spikes,function(x) try(smpsth(x[[1]],spikeWindow=c(0,3))))
#save(Spikes,smSpikes,file='spike_summary.rda')

load("spike_summary_133.rda")

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


# find crosses for all those cells that we have plotted
jet.colors<-colorRampPalette(c('navy','cyan','yellow','red'))
crosses=physplit$cross[match(colnames(allfreqs_goododours_mat),physplit$cell)]
heatmap(allfreqs_goododours_mat,scale='none',Rowv=NA,labCol=crosses,col=jet.colors(20),margins=c(6,5))

# correlation score instead
# find correlation between all pairs of columns (ie cells)
spcor=cor(allfreqs_goododours_mat,use='complete.obs')

# dendrogram is based on distance of 1-correlation score, but the 
# colours in the heatmap are still the correlation scores (ie hot is highly correlated)
heatmap(spcor,distfun=function(x) as.dist(1-x),scale='none',symm=T,col=jet.colors(20))
