make.pn.rates <- function(x, R.max=165, sigma=12, m=0.05) {
  pn <- matrix(0, ncol=ncol(x), nrow=nrow(x))
  rownames(pn) <- rownames(x)
  colnames(pn) <- colnames(x)
  # lateral suppression factor
  for (i in 1:nrow(orn)) {
    s <- m*sum(orn[i,])
    for (j in 1:ncol(orn)) {      
      pn[i,j] <- R.max*((x[i,j]^1.5)/(sigma^1.5 + x[i,j]^1.5 + s^1.5))
    }
  }
  return(pn)
}

add.noise <- function(x, delta=10, alpha=0.025) {
  return(apply(x, 2, function(x) { x + delta*tanh(alpha*x)*rnorm(length(x)) }))
}

make.model.responses <- function(pn, group, N=100, single.cell=F, method="lda") {  
  if (method == "lda") {
    weights <- lda(t(pn), factor(group))$scaling
  } else if (method == "logreg") {
    weights <- glm.fit(t(pn), factor(group), family=binomial())$coefficients
  } else if (method == "nonneg") {
    coef <- coefficients(penalized(factor(group), t(pn), positive=TRUE))   
    nonneg.weights <- coef[2:length(coef)]
    weights <- rep(0, nrow(pn))
    names(weights) <- rownames(pn)
    weights[names(nonneg.weights)] <- nonneg.weights
  }
  yes <- c()
  no <- c()
  for (i in 1:N) {
    curr.pn <- add.noise(pn)
    if (single.cell) {
      total.input <- curr.pn[, 1:ncol(pn) %in% group]
      yes <- c(yes, (total.input %*% weights)/mean(total.input))
      other.input <- curr.pn[, !(1:ncol(pn) %in% group)]
      no <- c(no, unlist(apply(curr.pn[, !(1:ncol(pn) %in% group)], 2, function(x) { (x %*% weights)/mean(x)})))
    } else {
      yes <- c(yes, unlist(apply(curr.pn[, as.logical(group)], 2, function(x) { (x %*% weights)/mean(x)})))
      no <- c(no, unlist(apply(curr.pn[, !as.logical(group)], 2, function(x) { (x %*% weights)/mean(x)})))
    }
  }
  return(list(yes=yes, no=no))
}

plot.discriminability <- function(pn, N=20) {
  x <- c()
  y <- c()
  for (i in 1:50) {
    print(i)
    for (n in 1:N) {      
      group <- rbinom(110, 1, i/110)
      if (sum(group) == 0)
        next 
      if (sum(group) == 1) {
        resp <- make.model.responses(pn, group, 100, single.cell=TRUE, method="lda")
      } else {
        resp <- make.model.responses(pn, group, 100, method="lda")
      }
      y <- c(y, (mean(resp$yes) - mean(resp$no)))
      x <- c(x,sum(group))
    }    
  }
  plot(x, y)
  return(list(x=x,y=y))
}

plot.model.response.hist <- function(response, ymax=100, x.min=-0.5, x.max=0.5, title = '') {
  red <- rgb(1,0,0,alpha=0.75)
  blue <- rgb(0,0,1,alpha=0.75)
  hist(response$yes, 100, freq=FALSE, xlim=c(x.min, x.max), ylim=c(0, ymax),  xlab='Total Input', col=red, border=red, main=title)
  par(new=T)
  hist(response$no, 100, freq=FALSE, xlim=c(x.min, x.max), ylim=c(0, ymax), xlab='Total Input', col=blue, border=blue, main='')
}

make.lda.response.prob <- function(pn, group, N=1000) {
  weights <- lda(t(pn), factor(group))$scaling
  mat <- matrix(0, ncol=ncol(pn), nrow=N)
  for (i in 1:N) {
    curr.pn <- add.noise(pn)
    for (j in 1:ncol(pn)) {
      mat[i,j] <- ifelse((curr.pn[,j] %*% weights)/mean(curr.pn[,j]) > 0, 1, 0)
    }
  }
  return(apply(mat, 2, mean))
}
