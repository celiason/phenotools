# function to make duplicate network tree?
# fac = edge scaling factor
duptree <- function(x, fac=10, plot=FALSE, ...) {
	library(igraph)
	# x <- twig1.dup
	if (!"dups" %in% names(x)) {
		stop("No duplicated characters identified.")
	}
	dups <- x$dups
	g <- graph_from_data_frame(dups[c('charnum1', 'charnum2')], directed=FALSE)
	wts <- twig1.dup$dups$stringdist
	# scale weights
	wts <- 1 - (wts - min(wts)) / max(wts - min(wts))
	E(g)$weight <- fac * wts
	if (plot) {
		par(mar=c(0,0,0,0))
		plot(g, edge.width=E(g)$weight, ...)
	}
	g
}