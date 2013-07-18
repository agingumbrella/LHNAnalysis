library(gphys)

#' Summarise a single cell as a series of raster plots
#'
#' This is known to work under Shahar's data collection routine
#' TODO - figure out how to merge everything into one monster plot
#' @param datadir Folder containing this cell's data
#' @param plotto one of screen (default), png or pdf
#' @param filenums integer vector containing chosen pxp file numbers (any order)
#' @param stimRange vector containing start and end of stimulus
#' @param merge.spikes Merge rasters (and smoothed traces) into single huge plot 
#' @param ... additional params passed to PlotRasterFromSweeps
#' @return return value
#' @export
#' @seealso \code{\link{gphys}}
#' @examples
#' SummariseShaharCell('/Volumes/Macintosh HD/Users/shahar/Data/120321/nm20120321c0',xlim=c(0,4000),stimRange=c(500,1000))
#' SummariseShaharCell('/Volumes/Macintosh HD/Users/shahar/Data/120321/nm20120321c0',xlim=c(0,4000),plotto='pdf')
#' SummariseShaharCell('/Volumes/Macintosh HD/Users/shahar/Data/120321/nm20120321c0',xlim=c(0,4000),stimRange=c(500,1000),plotto='pdf')
#' SummariseShaharCell('/Volumes/JData/JPeople/Shahar/Data/120308/nm20120308c0',xlim=c(0,4000),stimRange=c(500,1000),plotto='png')
#' ## Select specific pxp files to plot
#' SummariseShaharCell('/Volumes/JData/JPeople/Shahar/Data/120308/nm20120308c0',
#'  filenums=c(2,4,6,8),xlim=c(0,4000),stimRange=c(500,1000),plotto='png')
#' ## Plot two pxps to a single plot
#' SummariseShaharCell('/Volumes/JData/JPeople/Shahar/Data/130312/nm20130312c1',
#'   plotto='png',filenums=c(3,5),merge.spikes=T)
SummariseShaharCell<-function(datadir,plotto=c('screen','png','pdf'),filenums=NULL,subdir='',
	stimRange=NULL,ylim=c(-65,-35),lines.lwd=1,lines.col='darkgreen',
	AddVoltage=TRUE,merge.spikes=TRUE,pngw=1600,pngh=1200,...){
  plotto=match.arg(plotto)
  
  # fetch the spikes
  allspikess=CollectSpikesShaharCell(datadir=datadir, filenums=filenums, subdir=subdir,
     stimRange=stimRange, merge.spikes=merge.spikes)

  # retrieve avgwavefiles
  avgwavefiles=attr(allspikess,'avgwavefiles')
  # n.b. filenums may have come in as NULL
  # and been defined inside CollectSpikesShaharCell to match the avgwavefiles in this
  # directory
  filenums=as.integer(names(avgwavefiles))
  
  if(plotto=='screen'){
    par(ask=TRUE)
    for(spikes in allspikess) {
      PlotRasterFromSweeps(spikes,...)
    }
    par(ask=FALSE)    
  } else {
    for(i in seq(allspikess)) {
      spikes=allspikess[[i]]
      if(merge.spikes){
        # make a name which includes all of the pxps being merged 
        mergedfilenumstring=paste(sprintf("%03d",filenums),collapse='_')
        plotfilename=file.path(datadir,subdir,sprintf("merged-%s.%s",mergedfilenumstring,plotto))
      }
      else plotfilename=file.path(datadir,subdir,sprintf("%03d_raster.%s",filenums[i],plotto))
      if(plotto=='png'){
        png(file=plotfilename,width=pngw,height=pngh)
        dotsize=0.9
        # text otherwise looks rather small for these big pngs
        par(cex=2)
      } else {
        pdf(file=plotfilename, width=8,height=6)
        dotsize=.4
      }
      
      PlotRasterFromSweeps(spikes,dotsize=dotsize,...)
      if(AddVoltage){
        # browser()
        # if we are merging, we need to merge the average traces as well
        if(merge.spikes && length(avgwavefiles)>1){
          # load up all the average traces in one go
          allavgwavepaths=file.path(datadir,subdir,avgwavefiles)
          # first them into a list one by one
          allavgwaves=list()
          for(i in seq(allavgwavepaths)){
            avgwaves=read.table(allavgwavepaths[i],header=T)
            # fetch frequency of corresonding wave data
            tspp=tsp.nclamppxp(datadir,filenums[i])
            allavgwaves[[i]]=ts(avgwaves,start=0,freq=tspp[3])
          }
          # then take all the time series in that list and join them together
          # into one big multi column time series
          avgwavests=do.call(cbind,allavgwaves)
        } else {
          # if we are not merging or only one average wave file to add
          # just load that single average trace file
          avgwavepath=file.path(datadir,subdir,avgwavefiles[i])
          if(file.exists(avgwavepath)){
            avgwaves=read.table(avgwavepath,header=T)
            # fetch frequency of corresonding wave data
            tspp=tsp.nclamppxp(datadir,filenums[i])
            avgwavests=ts(avgwaves,start=0,freq=tspp[3])
          }
        }
        # now actually plot the average trace file(s) that we loaded
        AddLinesToRasterPlot(avgwavests,ylim=ylim,col=lines.col,lwd=lines.lwd)
      }
      dev.off()
    }
  }
}


CollectSpikesShaharCell<-function(datadir, subdir='', filenums=NULL,stimRange=c(500,1000), merge.spikes=TRUE)
{

  # see what spikes we've got
#  spikefiles=dir(datadir,patt="^[0-9]{3}_SP_")
  avgwavefiles=dir(file.path(datadir,subdir),patt="^[0-9]{3}_Avg_")
#  avgwavefiles=file.path(datadir,sprintf("%03d_Avg_RG0_A0++.txt",filenums))
  if(is.null(filenums)) {
    filenums=as.integer(sub("_.*","",avgwavefiles))
    # store these as the names of avgwavefiles since we will return
    # avgwavefiles (as an attribute of allspikess) but not filenums
    # at the end of the function
    names(avgwavefiles)=filenums
  } else {
    if(is.character(filenums)){
      # assume that this is a string in the format of cluster_Analyze
      # e.g. "03,05,07,09,11,13"
      filenums=as.integer(unlist(strsplit(filenums,",")))
    }
    # we specified particular pxp file numbers so we only want those average waves
    # first name the average wave files by their file numbers
    names(avgwavefiles)=as.integer(sub("_.*","",avgwavefiles))
    # then check that we have average wave files for the numbers we specified
    if(!all(filenums %in% names(avgwavefiles))) stop("we are missing some average waves")
    # select the avgwavefiles that match the pxp filenums that we have been given
    # (and reorder them in the process)
    avgwavefiles=avgwavefiles[as.character(filenums)]
  }
  allspikes=lapply(filenums,function(x) CollectSpikesFromSweeps(datadir,subdir=subdir,x,stimRange=stimRange))
  allspikess=lapply(allspikes,split)
  if(merge.spikes && length(allspikess)>1){
    # we want to merge and we have more than one pxp file
    merged=allspikess[[1]]
    for(spikelist in allspikess[-1]){
      # now merge the other spike files into the starting list one by one
      merged=merge(merged,spikelist)
    }
    # now replace allspikess with the merged list
    allspikess=list(merged)
  }
  # store avgwavefiles as an attribute so that we can use it in other functions
  attr(allspikess,'avgwavefiles')=avgwavefiles
  allspikess
}
