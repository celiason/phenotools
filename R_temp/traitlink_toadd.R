############################################################################
# fix below/add in to above:

foo <- function(x) {
	

	# (3) regexp matching
	stemchar2 <- lapply(paste0("\\b", unlist(terms2)), grep, x=cleanchars)
	
	names(stemchar2) <- rep(seq_along(terms2), times=sapply(terms2, length))

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

	tree2

### FIX BELOW
# output character assignment to see if correct
roots <- which(degree(tree2, v = V(tree2), mode = "in")==0, useNames = T)
system.time(nex.traitlink <- lapply(id, all_shortest_paths, graph=tree2, to=roots, mode="in"))
nms <- V(tree2)$name[id]

# get all paths
# takes ~ 15 seconds:
bodypart <- sapply(seq_along(paths), function(x) {names(sapply(paths[[x]], tail, 1))} )
names(bodypart) <- nms

# Figure out how many were correctly assigned to body part
bodypart.unique <- sapply(bodypart, unique)
bodypart.tables <- sapply(bodypart, table)
df <- data.frame(char = names(bodypart.tables), part = sapply(seq_along(bodypart.tables), function(x) {names(which.max(bodypart.tables[[x]]))}))
df$charnum <- as.numeric(str_extract(df$char, "\\d+"))
df$charlabel <- nex$charlab[df$charnum]
df <- dplyr::select(df, charnum, part, charlabel)
df$part <- factor(df$part, levels = names(roots), ordered = TRUE)
df <- arrange(df, part)
head(df)

}
