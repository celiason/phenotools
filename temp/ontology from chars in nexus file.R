# test character list
test <- c("Brain, cerebrum, shape", "Hindlimb, tarsometatarsus, length", "Brain, cerebellum", "Cerebellum elongated", "Tarsometatarsus, distal part", "Hindlimb")

# locate characters with terms separated by comma
test.com <- test[grepl(",", test)]

terms <- strsplit(test.com, ",")
terms <- lapply(terms, gsub, pattern="^ ", replacement="")
terms <- lapply(terms, tolower)
terms

# generate user/file-based trait ontology

# create edge list:
edges <- list()
for (i in seq_along(terms)) {
	tt <- terms[[i]]
	l <- length(tt)
	edges[[i]] <- lapply(1:(l-1), function(x) {tt[x:(x+1)]})
}
edges <- matrix(unlist(edges), ncol=2, byrow=TRUE)
edges

g <- graph_from_edgelist(edges)

plot(g, layout = -layout.reingold.tilford(g)[,2:1])


# search for terms not in ontology

test.nocom <- test[!grepl(",", test)]
test.nocom <- tolower(test.nocom)


# create groups of words based on trait ontology/tree
# find root and leaves
leaves <- which(degree(g, v = V(g), mode = "out")==0, useNames = T)
roots <- which(degree(g, v = V(g), mode = "in")==0, useNames = T)

# traverse tree and get all combinations of characters along trait ontology
reachable <- lapply(roots, function(x) {which(shortest.paths(g, x, mode="out") != Inf)})
terminal.nodes <- lapply(reachable, function(x) {x[which(degree(g, x, mode="out") == 0)]})
traversal <- lapply(seq_along(roots), function(x) {
	paths <- get.all.shortest.paths(graph=g, from=roots[x], to=terminal.nodes[[x]], mode="out")$res
	sapply(paths, function(vs) paste(V(g)[vs]$name, collapse="->"))
})


# locate matches
matches <- sapply(V(g)$name, grep, test.nocom)

id <- which(matches>=1)

# locate all paths in ontology connected to a certain character
lapply(id, get.all.shortest.paths, graph=g, to=roots, mode="in")


# sort??? interleave???


