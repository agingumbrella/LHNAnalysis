oddfilename<-function(igorfile,sweepnum){
  sweepnum=as.integer(sweepnum)
  justfile=basename(igorfile)
  dir(igorfile,sprintf("%03d_odd_",sweepnum),full=T)
}

avgname<-function(igorfile,sweepnum){
  sweepnum=as.integer(sweepnum)
  AvgFile=sprintf("%03d_Avg_RG0_A0++.txt",sweepnum)
  return(file.path(igorfile,AvgFile))
}

sputniknumber<-function(igorfile,sweepnum){
  odf=basename(oddfilename(igorfile,sweepnum))
  sputnum=sub(".*Sput_([0-9]+).*","\\1",odf)
  sputnum=as.integer(sputnum)
  if(is.na(sputnum)) warning("Unable to identify sputnik for ",igorfile,' sweep ',sweepnum)
  return(sputnum)
}

# given an odour response matrix with columns for different odours
# and rows for different sweep timepoints, reorder columns in desired 
# order adding NAs for missing odours
addnacols<-function(mat,odours){
    odours_we_have=colnames(mat)
    missing_odours=setdiff(odours,odours_we_have)
    # make an empty matrix to put our data
    finalmat=matrix(ncol=length(odours),nrow=nrow(mat))
    colnames(finalmat)=odours
    # assign the data that we do have to the final matrix
    finalmat[,odours_we_have]=mat
    finalmat
}