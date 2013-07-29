# info theory
# Mutual information b/w responses and cell identity (conditional on stimulus)
# I(r; i|s) = 1/N \sum_{i=1}^N \sum_r p(r|s, i) log_2 (p(r|s,i)/p(r|s)) bits
# where p(r|s) = 1/N \sum_{i=1}^N p(r|s, i)
# can then compute <I(r;i|s)> over all stimuli => <I(r;i|s)> = \sum_k p(s_k) I(r;i|s_k)
# in the case of p(s_k) = 1/K for all s_k, then equals 1/K \sum_k I(r;i|s_k) which is just the empirical mean

# Function takes probs matrix where rows are cells and columns are stimuli
# returns vector of mutual information between response and stimuli
mutual.info <- function(probs) {
    I <- c()
    p.rs <- apply(probs, 2, mean)
    for (i in 1:ncol(probs)) {
       I <- c(I, mean(probs[,i]*log2(probs[,i]/p.rs[i])+(1-probs[,i])*log2((1-probs[,i])/(1-p.rs[i]))))
    }
    return(I)
}

cluster.mutual.info <- function(curr.data, num.classes) {
  print(num.classes)
  train <- make.trial.probs(curr.data, 4)
  D <- make.D(train)
  clust <- hclust(as.dist(D))
#    clust <- hclust(as.dist(1-cor(t(curr.data))))
    clust.labels <- cutree(clust, k=num.classes)
    clust.probs <- make.class.probs(curr.data, clust.labels)
    return(mean(mutual.info(clust.probs)))
}

# RELATIVE ENTROPY MATRIX
# Compute D_ij = <D_JS[p(r|i,s)||p(r|j,s)]>_s = 1/K \sum_k D_JS[p(r|i,s_k)||p(r|j,s_k)]
# and D_JS[P||Q] = \sum_i log2(P(i)/Q(i))p(i)
# or, in this case,
# D_JS[p(r|i,s_k)||p(r|j,s_k)] = log2(p(1|i,s_k)/p(1|j,s_k))p(1|i,s_k) + log2(p(0|i,s_k)/p(0|j,s_k)p(0|i,s_k)
# = log2(p(1|i,s_k)/p(1|j,s_k))p(1|i,s_k) + log2(1-p(1|i,s_k)/(1-p(1|j,s_k)))(1-p(1|i,s_k))
make.D <- function(probs) {
    D <- matrix(0, nrow=nrow(probs), ncol=nrow(probs))
    for (i in 1:nrow(D)) {
        for (j in 1:nrow(D)) {
            D[i,j] <- mean(log2(probs[i,]/probs[j,])*probs[i,] + log2((1-probs[i,])/(1-probs[j,]))*(1-probs[i,]))
        }
    }
    return(D)
}

# BINOMIAL ENTROPY
entropy.binom <- function(p) {
    return(- p*log2(p) - (1-p)*log2(1-p))
}


# Beginning of a mutual-information based clustering algorithm...
# takes in bin.mat
# ALGORITHM:
# compute probs for all cells
# while
#  - estimate probs using current cluster assignments
#  - compute relatives entropies D
#  - merge two cells with minimum D_ij
info.cluster <- function(bin.mat) {
    merge.clusters <- function(labels, c1, c2) {
        next.clust <- max(unique(labels)) + 1
        labels[labels == c1 | labels == c2] <- next.clust
        new.labels <- labels
        n <- 1
        for (i in unique(labels)) {
            new.labels[labels == i] <- n
            n <- n+1
        }
        return(new.labels)
    }

    get.min.cluster <- function(D, labels) {
        min.val <- 1E10
        min.idx <- 0
        for (i in 1:nrow(D)) {
            for (j in 1:ncol(D)) {
                if (i != j) {
                    if (min.val > D[i,j]) {
                        min.idx <- c(i,j)
                        min.val <- D[i,j]
                    }
                }
            }
        }
        return(min.idx)
    }

    labels <- 1:nrow(bin.mat)
    probs <- make.class.probs(bin.mat, labels)
    D <- make.D(probs)
    to.merge <- get.min.cluster(D, labels)
    labels <- merge.clusters(labels, to.merge[1], to.merge[2])
    for (i in 1:(nrow(bin.mat)-1)) {
        probs <- make.class.probs(bin.mat, labels)
        D <- make.D(probs)
        to.merge <- get.min.cluster(D, labels)
        print(labels)
        labels <- merge.clusters(labels, to.merge[1], to.merge[2])
    }
    return(info.cluster.helper(bin.mat, labels))
}
