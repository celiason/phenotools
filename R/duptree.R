#' function to make duplicate network tree
#' x = a `nex` object with duplicated identified using `duplicated.nex`
#' @param fac edge scaling factor
#' @param plot whether to make a plot of the network
#' x <- twig12.d
duptree <- function(x, fac=10, plot=FALSE, ...) {
	oldpar <- par(no.readonly=TRUE)
	library(igraph)
	# x <- twig1.dup
	if (!"dups" %in% names(x)) {
		stop("No duplicated characters identified.")
	}
	dups <- x$dups
	g <- graph_from_data_frame(dups[c('char1', 'char2')], directed=FALSE)
	wts <- dups$stringdist
	# scale weights
	wts <- 1 - (wts - min(wts)) / max(wts - min(wts))
	E(g)$weight <- fac * wts
	if (plot) {
		par(mar=c(0,0,0,0))
		plot(g, edge.width=E(g)$weight, ...)
		par(oldpar)
	}
	g
}
