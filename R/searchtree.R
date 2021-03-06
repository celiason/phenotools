#' Function to grep search labels in tree and plot all connected components
#' 
#' @param graph (required) input graph of trait ontology generated by, e.g., read_ontology() or stemsearch()
#' @param pattern text to search for (regular expression)
#' @param plot whether to plot (logical)
#' @param mode whether to search for all edges to current node ("in") or away from ("out")
#' @param ... additional arguments passed to the `plot.igraph` function
#' @examples \dontrun{
#' # load Baumel and Whitmer (1993) ontology:
#' ont <- read_ontology(file = system.file("extdata", "baumel_ontology.txt",
#' package = "phenotools"))
#' # show all traits in ontology under "orbita"
#' searchtree(ont, "\\borbita\\b")
#' }
#' @references Baumel, J. J., and L. M. Witmer. 1993. Osteologia. P. in
#' Handbook of avian anatomy: nomina anatomica avium. Publications of the
#' Nuttall Ornithological Club.
#' @export
searchtree <- function(graph, pattern, plot=TRUE, mode="out", ...) {
	id <- grep(pattern, V(graph)$name)
	subcomps <- lapply(id, subcomponent, graph=graph, mode=mode)
	id <- unique(unlist(subcomps))
	newgraph <- simplify(induced_subgraph(graph, id))
	if (plot) {
		plot(newgraph, ...)
	}
	invisible(newgraph)
}
