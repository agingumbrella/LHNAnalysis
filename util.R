
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
