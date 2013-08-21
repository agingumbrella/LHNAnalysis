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

# cluster using relative entropy (an approximation to clustering using mutual info)
cluster.mutual.info <- function(curr.data, num.classes) {
  print(num.classes)
  train <- make.trial.probs(bin2data.frame(curr.data))
  D <- make.D(train)
  clust <- hclust(as.dist(D))
  clust.labels <- cutree(clust, k=num.classes)
  clust.probs <- make.class.probs(curr.data, clust.labels)
  return(mean(mutual.info(clust.probs)))
}

# average entropy per cluster
cluster.avg.entropy <- function(curr.data, num.classes) {
  print(num.classes)
  train <- make.trial.probs(bin2data.frame(curr.data))
  D <- make.D(train)
  clust <- hclust(as.dist(D))
  clust.labels <- cutree(clust, k=num.classes)
  clust.probs <- make.class.probs(curr.data, clust.labels)
  return(mean(sapply(apply(clust.probs, 1, entropy.binom),mean)))
}


# RELATIVE ENTROPY MATRIX
# Compute D_ij = <D_JS[p(r|i,s)||p(r|j,s)]>_s = 1/K \sum_k D_JS[p(r|i,s_k)||p(r|j,s_k)]
# and D_JS[P||Q] = \sum_i log2(P(i)/Q(i))p(i)
# or, in this case,
# D_JS[p(r|i,s_k)||p(r|j,s_k)] = log2(p(1|i,s_k)/p(1|j,s_k))p(1|i,s_k) + log2(p(0|i,s_k)/p(0|j,s_k)p(0|i,s_k)
# = log2(p(1|i,s_k)/p(1|j,s_k))p(1|i,s_k) + log2(1-p(1|i,s_k)/(1-p(1|j,s_k)))(1-p(1|i,s_k))
make.D <- function(probs, binom=T) {
    D <- matrix(0, nrow=nrow(probs), ncol=nrow(probs))
    for (i in 1:nrow(D)) {
        for (j in 1:nrow(D)) {
          if (binom) {
            D[i,j] <- mean(relative.entropy.binom(probs[i,], probs[j,]))
          } else {
            D[i,j] <- relative.entropy.multi(probs[i,], probs[j,])
          }
        }
    }
    return(D)
}

relative.entropy.binom <- function(p,q) {
  return(log2(p/q)*p + log2((1-p)/(1-q))*(1-p))
}

# compute relative entropy
relative.entropy.multi <- function(p,q) {
  return(sum(log2(p/q)*p))
}
  
# BINOMIAL ENTROPY
entropy.binom <- function(p) {
    return(- p*log2(p) - (1-p)*log2(1-p))
}

# Entropy for multinomial dist
entropy.multi <- function(p) {
  return(- sum(p * log(p)))
}

