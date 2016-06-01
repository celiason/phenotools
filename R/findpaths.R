# do other functions use this function??
# might need to save in nexustools/R

findpaths <- function(g) {
	# create groups of words based on trait ontology/tree
	# find root and leaves
	# leaves <- which(degree(g, mode = "out")==0, useNames = T)
	roots <- which(degree(g, mode = "in")==0, useNames = T)
	# traverse tree and get all combinations of characters along trait ontology
	reachable <- lapply(roots, function(x) {which(shortest.paths(g, x, mode="out") != Inf)})
	terminal.nodes <- lapply(reachable, function(x) {x[which(degree(g, x, mode="out") == 0)]})
	traversal <- lapply(seq_along(roots), function(x) {
		paths <- get.all.shortest.paths(graph=g, from=roots[x], to=terminal.nodes[[x]], mode="out")$res
		sapply(paths, function(vs) paste(V(g)[vs]$name, collapse="->"))
	})
}
