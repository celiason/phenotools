#' sort nexus file by something
#' @param x nexus object
#' @param by sorting by
#' @param ... other arguments passed to sort
#' @export
#' 
sort.nex <- function(x, by = c('taxlabels', 'charlabels'), ...) {
  res <- x
  by <- match.arg(by)
  if (by == 'taxlabels') {
    ord <- order(x$taxlabels)
    res$data <- x$data[ord, ]
    res$taxlabels <- sort(x$taxlabels, ...)
    return(res)  # it is important to have return() for some reason
  }
  if (by == 'charlabels') {
    ord <- order(x$charlabels, ...)
    res$data <- x$data[, ord]
    res$charlabels <- sort(x$charlabels, ...)
    return(res)
  }
}