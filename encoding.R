###############################
## NEURAL ENCODING FUNCTIONS ##
###############################
## Used to determine p(x|s), where x is a cell or class response, and s is an odor stimulus

# Returns matrix with some number of trials per cell with binary responses
binarize.trials <- function(mat, labels, blank.names = c("OilBl", "WatBl")) {
    bl <- labels %in% blank.names
    num.nonblank <- sum(!bl)
    labels.nonblank <- labels[!bl]
    bin.mat <- mat[,!bl]
    for (i in 1:nrow(bin.mat)) {
        for (j in 1:ncol(bin.mat)) {
            blank.mean <- mean(mat[i,bl])
            bin.mat[i,j] <- ifelse(poisson.test(bin.mat[i,j], blank.mean)$p.value < 0.05, 1, 0)
        }
    }
    return(bin.mat)    
}


# Given a binary response matrix, compute response bernoulli probabilities
# Allows for prior alpha for smoothing
make.trial.probs.binom <- function(bin.mat, num.reps=4, alpha=0.1) {
  # assumes 4 trials per odor
  if (is.null(dim(bin.mat))) {
    probs <- rep(0, length(bin.mat)/num.reps)
    names(probs) <- unique(names(bin.mat))
    for (i in unique(names(bin.mat))) {
      curr.cols <- names(bin.mat) %in% i
      probs[i] <- (sum(bin.mat[curr.cols])+alpha)/(num.reps+2*alpha)
    }
  } else {
    probs <- matrix(0, nrow=nrow(bin.mat), ncol=ncol(bin.mat)/num.reps)
    colnames(probs) <- unique(colnames(bin.mat))
    for (i in 1:nrow(bin.mat)) {
      for (j in unique(colnames(bin.mat))) {
        curr.cols <- colnames(bin.mat) %in% j
        probs[i,j] <- (sum(bin.mat[i,curr.cols])+alpha)/(num.reps+2*alpha)
      }
    }
  }
  return(probs)
}


# Given a binary response matrix, compute response multinomial probabilities
# Allows for prior alpha for smoothing
make.trial.probs.multi <- function(bin.mat, num.reps=4, alpha=0.1) {
  # assumes 4 trials per odor
  if (is.null(dim(bin.mat))) {
    probs <- rep(0, length(bin.mat)/num.reps)
    names(probs) <- unique(names(bin.mat))
    for (i in unique(names(bin.mat))) {
      curr.cols <- names(bin.mat) %in% i
      #probs[i] <- (sum(bin.mat[curr.cols])+alpha)/(num.reps+2*alpha)
    }
  } else {
    probs <- matrix(0, nrow=nrow(bin.mat), ncol=ncol(bin.mat)/num.reps)
    colnames(probs) <- unique(colnames(bin.mat))
    for (i in 1:nrow(bin.mat)) {
      for (j in unique(colnames(bin.mat))) {
        curr.cols <- colnames(bin.mat) %in% j
       # probs[i,j] <- (sum(bin.mat[i,curr.cols])+alpha)/(num.reps+2*alpha)
      }
    }
  }
  return(probs)
}

# Compute per-class bernoulli response probabilities given binary cellular response matrix
# Allows for prior alpha for smoothing
# class.labels is vector of class label for each cell
make.class.probs <- function(bin.mat, class.labels, num.reps=4, alpha=0.1) {
  labels <- unique(class.labels)
  odors <- unique(colnames(bin.mat))
  probs <- matrix(0, nrow=length(labels), ncol=ncol(bin.mat)/num.reps)
  rownames(probs) <- unique(class.labels)
  colnames(probs) <- unique(colnames(bin.mat))
  for (i in 1:length(labels)) {
    for (j in 1:length(odors)) {
      curr.label <- labels[i]
      curr.odor <- odors[j]
      if (sum(class.labels == curr.label) > 1) {
        curr <- bin.mat[class.labels %in% curr.label, colnames(bin.mat) %in% curr.odor]
        probs[i,j] <- (sum(curr) + alpha)/(nrow(curr)*num.reps + 2*alpha)
      } else {
       probs[i,j] <- (sum(bin.mat[class.labels %in% curr.label, colnames(bin.mat) %in% curr.odor])+alpha)/(num.reps+2*alpha)
     }
    }
  }
  return(probs)
}
