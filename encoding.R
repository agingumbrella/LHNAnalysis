###############################
## NEURAL ENCODING FUNCTIONS ##
###############################
## Used to determine p(x|s), where x is a cell or class response, and s is an odor stimulus

################################################################################
## Stuff for Gaussian moeling
# make gaussian model assuming independence between cells

indep.gauss.model.all <- function(rates, lambda=0.9, alpha=1) {
  x <- lapply(unique(colnames(rates)), function(odor) {
    temp <- rates[,colnames(rates) %in% odor]
    sigma <- matrix(0, nrow=nrow(rates), ncol=nrow(rates))
    diag(sigma) <- lambda * apply(temp, 1, var) + (1-lambda)*alpha*length(diag(sigma))
    return(list(mu=apply(temp, 1, mean), sigma=sigma))
  })
  names(x) <- unique(colnames(rates))
  return(x)
}

indep.gauss.model.type <- function(rates, types, lambda=0.9, alpha=1) {
  x <- lapply(unique(types$Type), function(t) {
    temp <- rates[,colnames(rates) %in% subset(types, Type == t)$Odorant]
    sigma <- matrix(0, nrow=nrow(rates), ncol=nrow(rates))
    diag(sigma) <- lambda * apply(temp, 1, var) + (1-lambda)*alpha*length(diag(sigma))
    return(list(mu=apply(temp, 1, mean), sigma=sigma))
  })
  names(x) <- unique(types$Type)
  return(x)
}

# make gaussian model with full covariance
dep.gauss.model.all <- function(rates,  lambda=0.9, alpha=1) {
  x <- lapply(unique(colnames(rates)), function(odor) {
    temp <- rates[,colnames(rates) %in% odor]
    T <- sum(colnames(rates) %in% odor)
    mu <- apply(temp, 1, mean)
#    sigma <- matrix(0, nrow=nrow(rates), ncol=nrow(rates))
#    for (i in 1:T) {
#      sigma <- sigma + (temp[,i] - mu) %*% t(temp[,i] - mu)
#    }
#    sigma <- 1/T * sigma
#    curr.diag <- diag(sigma)
#    sigma <- (1-lambda)*sigma
#    diag(sigma) <- curr.diag
    sigma <- cov(t(temp), t(temp)) * lambda + (1-lambda)*alpha*diag(nrow(rates))
    return(list(mu=mu, sigma=sigma))
  })  
  names(x) <- unique(colnames(rates))
  return(x)
}

dep.gauss.model.type <- function(rates, types, lambda=0.9, alpha=1) {
  x <- lapply(unique(types$Type), function(t) {
    temp <- rates[,colnames(rates) %in% subset(types, Type == t)$Odorant]
    mu <- apply(temp, 1, mean)
    sigma <- cov(t(temp), t(temp)) * lambda + (1-lambda)*alpha*diag(nrow(rates))
    return(list(mu=apply(temp, 1, mean), sigma=sigma))
  })
  names(x) <- unique(types$Type)
  return(x)
}

# compute average rates for all odors by cluster
avg.by.cluster <- function(rates, clusters) {
  out <- c()
  for (i in unique(clusters)) {
    if (sum(clusters == i) > 1) {
      out <- rbind(out, unlist(apply(rates[clusters==i,], 2, mean)))
    } else {
      out <- rbind(out, rates[clusters==i, ])
    }
  }
  rownames(out) <- unique(clusters)
  return(out)
}

avg.by.odor <- function(rates, odor.types) {
  out <- c()
  for (i in unique(odor.types)) {
    if (sum(odor.types == i) > 1) {
      out <- cbind(out, unlist(apply(rates[,odor.types == i], 1, mean)))
    } else {
      out <- cbind(out, rates[,odor.types == i])
    }
  }
  colnames(out) <- unique(odor.types)
  return(out)
}

################################################################################
## Stuff for Bernoulli modeling

# Returns matrix with some number of trials per cell with binary responses
binarize.trials <- function(mat, blank.names = c("OilBl", "WatBl"), p.thresh = 0.05) {
    bl <- colnames(mat) %in% blank.names
    bin.mat <- mat[,!bl]
    for (i in 1:nrow(bin.mat)) {
        for (j in 1:ncol(bin.mat)) {
            blank.mean <- mean(mat[i, bl])
            bin.mat[i,j] <- ifelse(poisson.test(mat[i,j], blank.mean)$p.value < p.thresh, 1, 0)
        }
    }
    return(bin.mat)    
}

# binarize odor responses using background rates (given matrix of odor response x cell)
binarize.trials.rates <-function(mat, rates) {
    bin.mat <- mat    
    for (i in 1:nrow(bin.mat)) {
    	for (j in 1:ncol(bin.mat)) {
	         bin.mat[i,j] <- ifelse(poisson.test(mat[i,j], rates[j])$p.value < 0.05, 1, 0)
	    }
    }
    return(bin.mat)    
}

# turn binary response matrix into data frame, for easier marginal calculations
bin2data.frame <- function(m) {
  fly.num <- c()
  cross.id <- c()
  odor.id <- c()
  responses <- c()
  all.crosses <- rownames(m)
  all.odors <- colnames(m)
  for (i in 1:nrow(m)) {
    for (j in 1:ncol(m)) {
      fly.num <- c(fly.num, i)
      cross.id <- c(cross.id, all.crosses[i])
      odor.id <- c(odor.id, all.odors[j])
      responses <- c(responses, m[i,j])
    }
  }
  return(data.frame(id=fly.num, cross=cross.id, odor=odor.id, response=responses))
}


# Given a binary response matrix, compute response bernoulli probabilities for cell x odor or cell class x odor or cell class x odor class or cell x odor class (all 4 possible combos)
# Allows for prior alpha for smoothing
# class.labels allows you to compute probabilities by odor class rather than by odor
# the class probability is then used as the response probability for all odors in that class
make.trial.probs <- function(bin, num.reps=4, alpha=0.1, odor.class.labels=NA, cell.class.labels=NA) {  
  # check if bin.mat is 1-dimensionsional
#  if (is.null(dim(bin.mat))) {
 #   probs <- rep(0, length(bin.mat)/num.reps)
  #  names(probs) <- unique(names(bin.mat))
   # for (i in unique(names(bin.mat))) {
    #  curr.cols <- names(bin.mat) %in% i
     # probs[i] <- (sum(bin.mat[curr.cols])+alpha)/(num.reps+2*alpha)
    #}
#  } else {

    if (is.na(odor.class.labels)) {
      odor.class.labels <- unique(bin$odor)
    }
  
    # use just the cells (no cell classes)
    if (is.na(cell.class.labels)) {
      cell.class.labels <- unique(bin$id)
    }
    num.per.cell <- length(unique(bin$odor))*num.reps

    cells <- rep(cell.class.labels, each=num.per.cell)
    odors <- rep(rep(odor.class.labels, each=num.reps), length(unique(bin$id)))
    probs <- table(cells, odors, bin$response)
    probs <- (probs[,,2] + alpha)/(probs[,,1]+probs[,,2] + 2*alpha)
    rownames(probs) <- cell.class.labels
    colnames(probs) <- odor.class.labels
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




