#' Remove NAs from nexus files
#' 
#' @param object nex object
#' @param ... other arguments
#' 
#' @export
#' 
na.omit.nex <- function(object, ...) {
  res <- object
  x$data <- as.matrix(x$data)
  id.col <- sapply(1:ncol(x$data), function(z) { all(is.na(x$data[,z], ...)) } )
  id.row <- sapply(1:nrow(x$data), function(z) { all(is.na(x$data[z,], ...)) } )
  res$data <- x$data[!id.row, !id.col]
  res$charlabels <- x$charlabels[!id.col]
  res$statelabels <- x$statelabels[!id.col]
  res$taxlabels <- x$taxlabels[!id.row]
  res$charpartition <- x$charpartition[!id.col]
  res$charset <- x$charset[!id.col]
  res$charnums <- x$charnums[!id.col]
  res$file <- x$file[!id.col]
  res
}