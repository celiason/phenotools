#' Check if x is a nexus file
#' @param x `nex` object
#' @export
is.nex <- function(x) {
  inherits(x, "nex")
}

# TODO: make as.nex (e.g., for converting from files imported by read.nexus.data, etc.)

