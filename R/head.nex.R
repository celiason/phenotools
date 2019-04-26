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
