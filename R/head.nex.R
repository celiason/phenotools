#' Look at top of nexus file
#' 
#' @param x nex object
#' @param ... optional arguments
#' 
#' @importFrom utils head
#' 
#' @export
#' 
head.nex <- function(x, ...) {
  head(x$data, ...)
}

#' Look at bottom of a nexus file
#' 
#' @param x nex object
#' @param ... optional arguments
#' 
#' @importFrom utils head
#' 
#' @export
#' 
tail.nex <- function(x, ...) {
  tail(x$data, ...)
}

