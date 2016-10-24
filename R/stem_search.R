# stem searcher
# tree = trait ontology (igraph object)
# x = nexus file
# x <- twig2.ss
# tree <- tree
#
stem_search <- function(tree, x) {

	roots <- which(igraph::degree(tree, v = V(tree), mode = "in")==0, useNames = T)
	notroots <- which(!V(tree) %in% roots)
	paths <- lapply(V(tree), all_shortest_paths, graph=tree, to=roots, mode="in")
	paths <- sapply(paths, "[[", "res")
	fullnames <- sapply(seq_along(paths), function(i) {paste0(names(paths[[i]]), collapse="->")})
	fullnames <- gsub("\\[|\\]|\\(|\\)", "", fullnames)

	# Whether to match to full names (e.g., axial skeleton->skull... or just 'skull')
	# terms <- fullnames
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
	# any(unlist(sapply(terms2, "==", "")))
	# terms2 <- lapply(seq_along(terms2), function(x) {terms2[[x]][terms2[[x]]!=""]})
	# any(unlist(sapply(terms2, "==", "")))

	# list of unique stem words
	# terms3 <- unique(unlist(terms2))
	# length(terms3)  # number of unique terms

	
	# stem searcher

	# only at beginning of word
	# stemchar <- lapply(paste0("\\b", terms3), grep, tolower(x$charlab))
	# names(stemchar) <- terms3
	# tail(sort(sapply(stemchar, length)))

	# stemchar is a named list with terms that match up to different characters
	# now we want to make connections between unique terms and all terms in 'fullnames'
	# sapply(terms3, grep, terms2)
	# grep(terms3[1], terms2)
	# terms3[1]
	# terms2[1152]
	# stemchar[!grepl("integer", stemchar)]

	# might want to have function return these results:
	# number of unmatched characters
	# length(x$charlab[setdiff(seq_along(x$charlab), unique(unlist(stemchar)))])
	# 119/477 # 25% unmatched

	# kinda slower (~14X) than only using unique word stems (terms3)
	stemchar2 <- lapply(paste0("\\b", unlist(terms2)), grep, tolower(x$charlab))
	names(stemchar2) <- rep(seq_along(terms2), times=sapply(terms2, length))

# length(unlist(terms2))

# length(stemchar2)

# stemchar2

# charnums <- seq_along(x$charlab)

# match(charnums[1], stemchar2)

# stemchar2
# sapply(names(stemchar2), function(x) {
# 	stemchar2[[x]]
# 	id <- names(stemchar2) == '1152'
# 	stemchar2[id]
# 	table(unlist(stemchar2[id]))
# })

# sapply(charnums, "%in%", stemchar2)

# termnums <- unlist(sapply(sapply(terms2, length), function(x) 1:x))
# length(termnums)
# length(unlist(stemchar2))
# sapply(stemchar2, length)

# paste0("\\b", unlist(terms2)[[2]])
# stemchar2

# x <- twig2.ss
# x$charlab[1]

# terms2[[1]]
# terms2 <- sapply(terms2, schinke)

# for (i in seq_along(x$charlab)) {
# 	i <- 35
# 	tmp <- terms2[[i]]
# 	sapply(tmp, grep, x$charlab)

# 	sapply(terms2[[i]], )

# }

	# stemchar2 names are unique vertices in tree
	# stemchar2 list contents are the characters that contain a given term

	# remove zero length list elements
	stemchar2 <- stemchar2[!lapply(stemchar2, length)==0]

	# TODO create links from terms to characters, weighted by number of matches


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
