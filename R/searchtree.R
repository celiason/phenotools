# function to grep search labels in tree and plot all connected components
searchtree <- function(graph, pattern, plot=TRUE, ...) {
	id <- grep(pattern, V(graph)$name)
	subcomps <- lapply(id, subcomponent, graph=graph, mode="out")
	id <- unique(unlist(subcomps))
	newgraph <- simplify(induced_subgraph(graph, id))
	# V(newgraph)$name <- gsub(" ", "\n", V(newgraph)$name)
	if (plot) {
		plot(newgraph, ...)
	}
	invisible(newgraph)
}
