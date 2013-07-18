#!/usr/bin/env Rscript
args=commandArgs(trailingOnly = TRUE)
print(args)

setwd("/Users/shahar/projects/RandIgor")
source("Avg_IgorImport.R", chdir = TRUE)
source("ClusteringFunctions.R", chdir = TRUE)
# cat("args[3]=",args[3],"!\n",sep="")
# if the subdir argument is empty in filemaker it will actually look as if
# subdir is NA rather than the empty string (which is what SummariseShaharCell wants)
subdir=ifelse(is.na(args[3]),"",args[3])
SummariseShaharCell(args[1],filenums=args[2],subdir=subdir,xlim=c(5,2995),stimRange=c(500,1000),AddVoltage=T,plotto='png',pngh=3000,pngw=1200,ylab="",yaxt='n',merge.spikes=T)