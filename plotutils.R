
# make side colors for plots
make.side.colors <- function(labels, colors) {
  return(sapply(labels, function(x) colors[[x]]))
}
