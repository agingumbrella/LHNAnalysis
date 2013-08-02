# make side colors for plots
make.side.colors <- function(labels, colors) {
  return(sapply(labels, function(x) colors[[x]]))
}

# fix cells with no cross id or two cross ids
filter.cross.ids <- function(cross.ids) {
  good.ids <- c()
  for (curr.id in cross.ids) {
    if (length(curr.id) > 1) {
      curr.id <-  curr.id[2]
    }    
    if (curr.id == "") {
      curr.id <- "unknown"
    } 
    good.ids <- c(good.ids, curr.id)
  }
  return(good.ids)
}

# make matrices by repeating rows or cols
rep.row <- function(x,n){
  return(matrix(rep(x,each=n),nrow=n))
}
rep.col <- function(x,n){
  return(matrix(rep(x,each=n), ncol=n, byrow=TRUE))
}

# compute lifetime sparseness
# S = 1/(1-1/N)(1-(\sum_j r_j /N)^2/\sum r_j^2/N)
# N is number of odors
lifetime.sparseness <- function(r) {
  N <- nrow(r)
  return(apply(r, 2, function(x) 1/(1-1/N)*(1-sum(x/N)^2/sum(x^2/N))))
}
