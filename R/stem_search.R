# stem searcher
# tree = trait ontology (igraph object)
# x = nexus file

stem_search <- function(tree, x) {

	# terms <- V(tree)$name
	terms <- V(tree)$name

	terms2 <- str_split(terms, "->|\\s")
	terms2 <- sapply(terms2, unique)
	# remove dots
	terms2 <- lapply(seq_along(terms2), function(x) {
		res <- terms2[[x]]
		res[!grepl("\\.", res)]
	})
	# remove single, two letter words
	terms2 <- lapply(seq_along(terms2), function(x) {
		res <- terms2[[x]]
		res[grepl("\\w{3,}", res)]
	})
	# remove empties
	any(unlist(sapply(terms2, "==", "")))
	# terms2 <- lapply(seq_along(terms2), function(x) {terms2[[x]][terms2[[x]]!=""]})
	# any(unlist(sapply(terms2, "==", "")))

	# list of unique stem words
	terms3 <- unique(unlist(terms2))
	length(terms3)  # number of unique terms
	terms3


	# stem searcher

	# only at beginning of word
	# system.time(stemchar <- lapply(paste0("\\b", terms3), grep, tolower(x$charlab)))
	# length(stemchar)
	# names(stemchar) <- terms3
	# tail(sort(sapply(stemchar, length)))

	# might want to have function return these results:
	# number of unmatched characters
	# x$charlab[setdiff(seq_along(x$charlab), unique(unlist(stemchar)))]
	# 119/477 # 25% unmatched

	# not much slower to do all stems in every term of ontology
	stemchar2 <- lapply(paste0("\\b", unlist(terms2)), grep, tolower(x$charlab))

	names(stemchar2) <- rep(seq_along(terms2), times=sapply(terms2, length))


	# create links from terms to characters, weighted by number of matches
	# TODO - weights for each character??
	# remove zero length list elements
	stemchar2 <- stemchar2[!lapply(stemchar2, length)==0]

	# create links between terms and characters
	x <- sapply(seq_along(stemchar2), function(x) {
		paste0(names(stemchar2[x]), "--", stemchar2[[x]])
	})
	names(x) <- names(stemchar2)
	edges <- do.call(rbind, strsplit(unlist(x), "--"))
	edges[, 2] <- paste0("char", edges[, 2])
	edges[, 1] <- V(tree)$name[as.numeric(edges[, 1])]
	newverts <- unique(edges[, 2])

	# create new network with characters added in
	tree2 <- add_vertices(tree, nv=length(newverts), name=newverts) 
	# add new edges
	tree2 <- add_edges(tree2, apply(edges, 1, "c"))
	# tree2 <- add_edges(tree2, sapply(1:nrow(edges), function(x) { c(edges[x,1], edges[x,2])  }))
	tree2

}
