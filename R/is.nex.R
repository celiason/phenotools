#' Check if x is a `nex` object
#' 
#' @param x `nex` object
#' @examples
#' \dontrun{
#' is.nex(twig)
#' }
#' @export
is.nex <- function(x) {
  inherits(x, "nex")
}

# TODO: make as.nex (e.g., for converting from files imported by read.nexus.data, etc.)

