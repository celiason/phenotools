#' Utility function to print metadata for nexus file
#' @param x nex class object
#' @param ... other arguments (not used)
#' @export
#' 
print.nex <- function(x, ...) {
  nchar <- ncol(x$data)
  ntax <- nrow(x$data)
  cat('NEXUS data file with', nchar, 'morphological characters and', ntax, 'taxa\n')
}
