source("encoding.R")

require(mvtnorm)

# stuff for gaussian

# classify a single odor using Gaussian model
classify.odor.gauss <- function(z, model) {
  # length of models is number of possible labels
  logliks <- c()
  for (i in 1:length(model)) {
     logliks <- c(logliks, dmvnorm(z, mean=model[[i]]$mu, sigma=model[[i]]$sigma, log=TRUE)[1])
  }
  return(which.max(logliks))
}

# classify all odors given a Gaussian model
classify.all.odors.gauss <- function(z, model) {
  preds <- c()
  for (i in 1:ncol(z)) {
    preds <- c(preds, classify.odor.gauss(z[,i], model))
  }
  names(preds) <- names(model)[preds]
  return(preds)
}

# classify odors using single neurons and Gaussian model
classify.single.neurons.gauss <- function(z, model) {
  preds <- c()
  for (i in 1:length(model)) {
    preds <- rbind(preds, dnorm(z, model[[i]]$mu, diag(model[[i]]$sigma), log=TRUE))
  }
  return(apply(preds, 2, which.max))
}

# find cross-validated accuracy for single neuron predictions
cross.validate.single.neurons <- function(mat, clust, num.classes, num.reps=4, do.type=F) {
  pred.labels <-rep(0, num.classes)
  curr.types <- factor(lhn.types$Type)
  idx <- as.integer(unique(curr.types))
  for (i in 0:(num.reps-1)) {
    leave.out <- seq(1, ncol(mat), num.reps)+i
    train.data <- mat[,!(1:ncol(mat) %in% leave.out)]
    test.data <- mat[,leave.out]
    clust.labels <- cutree(clust, k=num.classes)

    train <- avg.by.cluster(train.data, clust.labels)
    test <- avg.by.cluster(test.data, clust.labels)
    if (do.type) {
      model <- indep.gauss.model.type(train, lhn.types)
    } else {
      model <- indep.gauss.model.all(train)
    }
    if (do.type) {
      for (j in 1:ncol(test)) {
        pred.labels <- pred.labels + as.integer(classify.single.neurons.gauss(test[,j], model) == as.integer(lhn.types[j,]$Type))
      }
    } else {
      for (j in 1:ncol(test)) {
        pred.labels <- pred.labels + as.integer(classify.single.neurons.gauss(test[,j], model) == j)
      }
    }
  }
  return(pred.labels/ncol(mat))
}

# find cross-validated  accuracy for population predictions
cross.validate.gauss <- function(mat, clust, num.classes, indep=T, num.reps=4, do.type=F) {
  accuracy <- c()
  curr.types <- factor(lhn.types$Type)
  for (i in 0:(num.reps-1)) {
    leave.out <- seq(1, ncol(mat), num.reps)+i
    train.data <- mat[,!(1:ncol(mat) %in% leave.out)]
    test.data <- mat[,leave.out]
    clust.labels <- cutree(clust, k=num.classes)
    train <- avg.by.cluster(train.data, clust.labels)
    test <- avg.by.cluster(test.data, clust.labels)    
    if (indep) {
      if (do.type) {
        model <- indep.gauss.model.type(train, lhn.types)
      } else {
        model <- indep.gauss.model.all(train)
      }
    } else {
      if (do.type) {
        model <- dep.gauss.model.type(train, lhn.types)
      } else {
        model <- dep.gauss.model.all(train)
      }
    }
    pred.labels <- classify.all.odors.gauss(test, model)
    if (do.type) {
      accuracy <- c(accuracy, sum(names(pred.labels) == curr.types)/length(curr.types))
    } else {
      accuracy <- c(accuracy, sum(pred.labels == 1:ncol(test.data))/ncol(test.data))
    }
  }
  return(mean(accuracy))
}


pred.from.sample <- function(mat, indep=T,num.reps=10) {
  accs <- c()
  for (i in 10:15) {
    print(i)
    curr.accs <- c()
    for (j in 1:10) {
      curr.accs <- c(curr.accs, cross.validate.gauss(mat[sample(1:nrow(mat),i),], indep=indep))
    }
    accs <- c(accs, mean(curr.accs))
  }
  return(accs)
}
        
 # stuff for LDA
average.rates <- function(m, row.labels) {
  out <- matrix(0, nrow=length(unique(row.labels)), ncol=ncol(m))
  colnames(out) <- colnames(m)
  all.row.labels <- unique(row.labels)
  for (i in 1:length(all.row.labels)) {
    if (sum(row.labels == all.row.labels[i]) > 1) {
      out[i,] <- apply(m[row.labels == all.row.labels[i],], 2, mean)
    } else {
      out[i,] <- sapply(m[row.labels == all.row.labels[i],], mean)
    }
  }
  return(out)
}

# find cross-validated odorant prediction accuracy for population LDA predictions
cross.validate.odors.lda <- function(mat, num.classes, num.reps=4, num.odors=32) {
  accuracy <- c()
  print(num.classes)
  for (i in 0:(num.reps-1)) {
    leave.out <- seq(1, ncol(mat), num.reps)+i
    train <- mat[,!(1:ncol(mat) %in% leave.out)]
    test <- mat[,leave.out]
    if (num.classes < nrow(train)) {
      clust.labels <- cutree(hclust(dist(train)), k=num.classes)
      train <- average.rates(train, clust.labels)
      test <- average.rates(test, clust.labels)
    }        
    model <- lda(t(train), colnames(train), method="mle")
    pred.labels <- predict(model, t(test))
    accuracy <- c(accuracy, sum(pred.labels$class == colnames(test))/ncol(test))
  }
  return(accuracy)
}

# find cross-validated accuracy for population LDA predictions
cross.validate.types.lda <- function(mat, num.classes, odor.types, num.reps=4, num.odors=32) {
  accuracy <- c()
  print(num.classes)
  colnames(mat) <- rep(odor.types, each=num.reps)
  for (i in 0:(num.reps-1)) {    
    leave.out <- seq(1, ncol(mat), num.reps)+i
    train <- mat[,!(1:ncol(mat) %in% leave.out)]
    test <- mat[,leave.out]
    if (num.classes < nrow(train)) {
      clust.labels <- cutree(hclust(dist(train)), k=num.classes)
      train <- average.rates(train, clust.labels)
      test <- average.rates(test, clust.labels)
    }        
    model <- lda(t(train), colnames(train), method="mle")
    pred.labels <- predict(model, t(test))
    accuracy <- c(accuracy, sum(pred.labels$class == colnames(test))/ncol(test))
  }
  return(accuracy)
}

confusion.matrix <- function(res) {

}

pois <- function(k,l) {
    return(((l^k)*exp(-l))/factorial(k))
}

# first do with naive independence assumption (factor so that p(x|s) = \prod_i p(x_i|s)
# Try with rates and Poisson emissions
classify.naive.poisson <- function(train, test) {
    # rows are odors, cols are probs
    label.probs <- matrix(0, nrow=ncol(test), ncol=ncol(test))
    colnames(label.probs) <- colnames(train)
    rownames(label.probs) <- rownames(train)
    for (i in 1:ncol(test)) {
        for (j in 1:ncol(test)) {
            label.probs[i,j] <- sum(log(ppois(test[,i], train[,j])))
        }                                
    }
    test <- test[,unique(colnames(train))]
    labels <- apply(label.probs, 1, which.max)
    return(labels)
}

# Classify neurons using binomial probability
classify.naive.binom <- function(train, test) {
    if (is.null(dim(train))) {
        label.probs <- matrix(0, nrow=length(train), ncol=length(train))
    } else {
        label.probs <- matrix(0, nrow=ncol(test), ncol=ncol(test))
    }
    colnames(label.probs) <- colnames(train)
    rownames(label.probs) <- colnames(train)
    for (i in 1:ncol(label.probs)) {
        for (j in 1:ncol(label.probs)) {
            if (is.null(dim(train))) {
                label.probs[i,j] <- log(ifelse(test[i] == 1, train[j], 1-train[j]))
            } else {
                label.probs[i,j] <- sum(log(ifelse(test[,i] == 1, train[,j], 1-train[,j])))
            }
        }
    }
    labels <- apply(label.probs, 1, which.max)
    return(labels)
}

classify.type.naive.binom <- function(train, test, odor.types) {
  if (is.null(dim(train))) {
    label.probs <- matrix(0, nrow=length(train), ncol=length(train))
  } else {
    label.probs <- matrix(0, nrow=ncol(train), ncol=ncol(train))
  }
  colnames(label.probs) <- colnames(train)
  rownames(label.probs) <- colnames(train)
  for (i in 1:ncol(label.probs)) {
    for (j in 1:ncol(test)) {
      if (is.null(dim(train))) {
        label.probs[i,odor.types[j]] <- label.probs[i,odor.types[j]] + log(ifelse(test[i] == 1, train[j], 1-train[j]))
      } else {
        label.probs[i,odor.types[j]] <- label.probs[i,odor.types[j]] + sum(log(ifelse(test[,i] == 1, train[,odor.types[j]], 1-train[,odor.types[j]])))
      }
    }
  }
  labels <- apply(label.probs, 1, which.max)
  return(labels)
}

# train on 3 examples for each cell and test on last
make.trial.rates <- function(mat) {
    # assumes 3 trials per odor in training set
    rates <- matrix(0, nrow=nrow(mat), ncol=ncol(mat)/3)
    colnames(rates) <- unique(colnames(mat))
    for (i in 1:nrow(mat)) {
        for (j in unique(colnames(mat))) {
            curr.cols <- colnames(mat) %in% j
            rates[i,j] <- mean(mat[i,curr.cols])
        }
    }
    return(rates)
}


cluster.probs <- function(probs, num.classes) {
  D <- make.D(probs)
  clust <- hclust(as.dist(D))
  return(cutree(clust, k=num.classes))
}

cross.validate.odors <- function(mat, predict.fun, num.classes, num.reps=4, num.odors=34) {
    accuracy <- c()
    print(num.classes)
    for (i in 0:(num.reps-1)) {
      leave.out <- seq(1, ncol(mat), num.reps)+i
      train.data <- mat[,!(1:ncol(mat) %in% leave.out)]
      test.data <- mat[,leave.out]
      pred.labels <- predict.fun(train.data, test.data, num.reps, num.classes)
      accuracy <- c(accuracy, sum(pred.labels == 1:ncol(test.data))/num.odors)      
    }
    return(accuracy)
}

cross.validate.types <- function(mat, predict.fun, num.classes, odor.types, num.reps=4, num.odors=34) {
  accuracy <- c()
  print(num.classes)
  for (i in 0:(num.reps-1)) {
    leave.out <- seq(1, ncol(mat), num.reps)+i
    train.data <- mat[,!(1:ncol(mat) %in% leave.out)]
    test.data <- mat[,leave.out]
    pred.labels <- predict.fun(train.data, test.data, num.reps, num.classes, odor.types)
    accuracy <- c(accuracy, sum(pred.labels == 1:length(odor.types))/num.odors)      
  }
  return(accuracy)
}


predict.odor.all <- function(train.data, test.data, num.reps, num.classes) {
  curr.data <- bin2data.frame(train.data)
  train <- make.trial.probs(curr.data, num.reps=num.reps-1)          
  clust.labels <- cluster.probs(train, num.classes)
  clust.train <- make.trial.probs(curr.data, cell.class.labels=clust.labels, num.reps=num.reps-1)  
  train <- matrix(0, nrow=nrow(train.data), ncol=ncol(clust.train))
  colnames(train) <- colnames(clust.train)              
  # replace the individual probability with the class probability
  for (k in 1:num.classes) {
    for (i in which(clust.labels == k)) {
      train[i,] <- clust.train[k,]
    }
  }
  test.data <- test.data[,colnames(train)]
  pred.labels <- classify.naive.binom(train, test.data)
}

predict.odor.exemplar <- function(train.data, test.data, num.reps, num.classes) {
  curr.data <- bin2data.frame(train.data)
  train <- make.trial.probs(curr.data, cell.class.labels=clust.labels, num.reps=num.reps-1)
  clust.labels <- cluster.probs(train, num.classes)
  train.mi <- mutual.info(t(train)) # note that is transposed to do by cell and not by odor
  max.mi.idx <- which(train.mi %in% sapply(1:num.classes, function(i) max(train.mi[clust.labels==i])))
  train <- train[max.mi.idx,]
  test.data <- test.data[max.mi.idx, colnames(train)]
  pred.labels <- classify.naive.binom(train, test.data)
  return(pred.labels)
}
 
predict.type.all <- function(train.data, test.data, num.reps, num.classes, odor.types) {
  curr.data <- bin2data.frame(train.data)
  train <- make.trial.probs(curr.data, num.reps=num.reps-1)          
  clust.labels <- cluster.probs(train, num.classes)
  clust.train <- make.trial.probs(curr.data, cell.class.labels=clust.labels, odor.class.labels=odor.types, num.reps=num.reps-1)  
  train <- matrix(0, nrow=nrow(train.data), ncol=ncol(clust.train))
  colnames(train) <- colnames(clust.train)              
  # replace the individual probability with the class probability
  for (k in 1:num.classes) {
    for (i in which(clust.labels == k)) {
      train[i,] <- clust.train[k,]
    }
  }
  pred.labels <- classify.type.naive.binom(train, test.data, odor.types)
}


# use maximum mutual information exemplar to predict using Bernoulli probs
predict.type.exemplar <- function(train.data, test.data, num.reps, num.classes, odor.types) {
  curr.data <- bin2data.frame(train.data)
  train <- make.trial.probs(curr.data, num.reps=num.reps-1)  
  clust.labels <- cluster.probs(train, num.classes)
  train.mi <- mutual.info(t(train)) # note that is transposed to do by cell and not by odor
  max.mi.idx <- which(train.mi %in% sapply(1:num.classes, function(i) max(train.mi[clust.labels==i])))
  clust.train <- make.trial.probs(curr.data, cell.class.labels=clust.labels, odor.class.labels=odor.types, num.reps=num.reps-1)  
  train <- matrix(0, nrow=nrow(train.data), ncol=ncol(clust.train))
  colnames(train) <- colnames(clust.train)              
  # replace the individual probability with the class probability
  for (k in 1:num.classes) {
    for (i in which(clust.labels == k)) {
      train[i,] <- clust.train[k,]
    }
  }
  train <- train[max.mi.idx,]
  pred.labels <- classify.type.naive.binom(train, test.data, odor.types)
  return(pred.labels)
}


# then do allowing pairwise interactions? Maybe unnecessary b/c so accurate

# -- Compare with leaving out single neurons
cross.validate.leaveout <- function(mat) {
    accuracy <- c()
    for (i in 1:nrow(mat)) {
        temp.mat <- mat[!(1:nrow(mat) %in% i),]
        accuracy <- c(accuracy,mean(cross.validate(temp.mat)))
    }
    return(accuracy)
}

# TODO Compute cross-validated accuracy with shuffled labels for comparison

# returns 
single.neuron.pred <- function(mat) {
    accuracy <- c()
    for (i in 1:nrow(mat)) {
        curr.acc <- c()
        for (j in 0:3) {
            leave.out <- seq(1, ncol(mat), 4) + j     
            train <- make.trial.probs(mat[i, !(1:ncol(mat) %in% leave.out)], 3)
            test <- mat[i, leave.out]
            labels <- classify.naive.binom(train, test)
            curr.acc <- c(curr.acc, sum(labels == 1:34)/34)
        }
        accuracy <- c(accuracy, mean(curr.acc))
    }
    return(accuracy)
}

# returns confusion matrix
single.neuron.confusion <- function(mat) {
  confusion <- matrix(0, nrow=ncol(mat), ncol=ncol(mat))
  for (i in 1:nrow(mat)) {
    for (j in 0:3) {
      leave.out <- seq(1, ncol(mat), 4) + j     
      train <- make.trial.probs(mat[i, !(1:ncol(mat) %in% leave.out)], 3)
      test <- mat[i, leave.out]
      labels <- classify.naive.binom(train, test)
      for (m in 1:ncol(mat)) {
        for (n in 1:ncol(mat)) {
          confusion[m,labels[n]] <- confusion[m,labels[n]] + 1
        }
      }
    }
  }
    # create 
   return(confusion)
}


