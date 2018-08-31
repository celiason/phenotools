#' Remove invariant characters
#' 
#' @param x nex object
#' @export
#' 
remove.invar <- function(x) {
	charlengths <- sapply(1:ncol(x$data), function(z) { length(na.omit(unique(x$data[, z])))})
	# find which characters have only one state
	noinvar <- charlengths!=1
	# remove these invariant characters
	x[, noinvar]
}
