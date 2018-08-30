#' Check if x is a nexus file
#' @param x `nex` object
#' @export
is.nex <- function(x) {
  inherits(x, "nex")
}
