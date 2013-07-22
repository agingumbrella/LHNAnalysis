source("encoding.R")

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
    rownames(label.probs) <- rownames(train)
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

cross.validate <- function(mat, num.classes=NA, num.reps=4, num.odors=34) {
    accuracy <- c()
    for (i in 0:(num.reps-1)) {
        leave.out <- seq(1, ncol(mat), num.reps)+i
        curr.data <- mat[,!(1:ncol(mat) %in% leave.out)]
        if (is.na(num.classes)) {
            train <- make.trial.probs(curr.data, num.reps-1)
        } else {
 #           D <- make.D(curr.data)
            clust <- hclust(as.dist(1-cor(t(curr.data))))
            clust.labels <- cutree(clust, k=num.classes)
            clust.train <- make.class.probs(curr.data, clust.labels, num.reps-1)
            train <- matrix(0, nrow=nrow(mat), ncol=ncol(mat)/num.reps)
            # replace the individual probability with the class probability
            for (k in 1:num.classes) {
                for (i in which(clust.labels == k)) {
                    train[i,] <- clust.train[k,]
                }
            }
        }
        test <- mat[,leave.out]
        labels <- classify.naive.binom(train, test)
        accuracy <- c(accuracy, sum(labels == 1:num.odors)/num.odors)
    }
    return(accuracy)
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

# -- collapse to clusters and redo prediction

# -- Compare with single neuron predictions
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
