devtools::load_all('~/github/nexustools')

# test character list

test <- c("Brain, cerebrum, shape", "Hindlimb, tarsometatarsus, length", "Brain, cerebellum", "Cerebellum elongated", "Tarsometatarsus, distal part", "Hindlimb")

test <- cleantext(test, comma=FALSE)

# locate characters with terms separated by comma
test.com <- test[grepl(",", test)]
terms <- strsplit(test.com, ",")
terms <- lapply(terms, gsub, pattern="^ ", replacement="")
terms <- lapply(terms, tolower)
head(terms)

# generate user/file-based trait ontology

# create edge list
edges <- list()
for (i in seq_along(terms)) {
	tt <- terms[[i]]
	l <- length(tt)
	edges[[i]] <- lapply(1:(l-1), function(x) {tt[x:(x+1)]})
}
edges <- matrix(unlist(edges), ncol=2, byrow=TRUE)
edges

# create and plot trait ontology
g <- graph_from_edgelist(edges)
V(g)$name <- schinke(V(g)$name)

plot(g, layout = -layout.reingold.tilford(g)[,2:1], vertex.size=0, edge.arrow.size=0)


# search for terms not in ontology
test.nocom <- test[!grepl(",", test)]
test.nocom <- tolower(test.nocom)

roots <- which(degree(g, v = V(g), mode = "in")==0, useNames = T)

# locate matches...
matches <- sapply(V(g)$name, grep, test.nocom)

id <- which(matches>=1)

# locate all paths in ontology connected to a certain character
lapply(id, all_shortest_paths, graph=g, to=roots, mode="in")

# sort??? interleave???
