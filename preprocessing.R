###################################
## EPHYS PREPROCESSING FUNCTIONS ##
###################################
# Takes per-cell spike times and turns into firing rates and binary cellular responses

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

# make matrix of all firing rates for all cells
make.total.data.matrix <- function(spikes, nr, nc, num.odors=36, num.reps=4) {
  all.mat <- matrix(NA, ncol=nc, nrow=nr)
  for (i in 1:length(spikes)) {
    labels <- c()
    n <- 1
    curr.labels <- names(spikes[[i]])
    for (j in 1:num.odors) {
      for (k in 1:num.reps) {				
        all.mat[i,n] <- ps.rate(spikes[[i]][[j]][[k]])
        labels <- c(labels, curr.labels[j])
        n <- n+1
      }
    }
  }
  colnames(all.mat) <- labels
  return(all.mat)
}

# Use Poisson hypothesis test to determine if each trial for each cell is above the cell's background firing rate
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
