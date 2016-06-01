setwd("/Users/chadeliason/Documents/UT/projects/phenome")

devtools::load_all('~/github/nexustools')  # load packages

twig <- read.nex("data/theropods/Turner_etal_2012.nex")  # load data

file <- "baumel ontology.txt"  # trait ontology, text file

# USING TEXT (TAB-BASED) ONTOLOGY
tree <- treechart(file="/Users/chadeliason/Documents/UT/projects/phenome/baumel ontology.txt")

tree <- simplify(tree, edge.attr.comb = "first")

# stem latin words in ontology
V(tree)$name <- schinke(V(tree)$name)

# plot a subtree
subtree <- induced_subgraph(tree, subcomponent(tree, grep("tarsometat", V(tree)$name), mode="out"))

pdf(file = "figure/ontology.pdf", width=9, height=8)
par(mar=c(0,0,0,0))
plot(subtree, layout=-layout.reingold.tilford(subtree)[,2:1], vertex.size=0, edge.arrow.size=0)
dev.off()



################################################################################
# 5-31-16 - Working on matching characters to trait ontology using REGEX search
# with stem words
################################################################################

# terms <- V(tree)$name
terms <- schinke(V(tree)$name)

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

# match each stem to each character using REGEX
stemchar <- sapply(terms3, grep, twig$charlab)

# this takes really long:
# for (i in seq_along(twig$charlab)) {
# 	for (j in seq_along(terms2)) {
# 		sum(sapply(terms2[[j]], str_detect, twig$charlab[i]))
# 	}
# }

# only at beginning of word
system.time(stemchar <- lapply(paste0("\\b", terms3), grep, tolower(twig$charlab)))
length(stemchar)
names(stemchar) <- terms3
tail(sort(sapply(stemchar, length)))
# number of unmatched characters
twig$charlab[setdiff(seq_along(twig$charlab), unique(unlist(stemchar)))]
107/477 # 22% unmatched


# not much slower to do all stems in every term of ontology
system.time(stemchar2 <- lapply(paste0("\\b", unlist(terms2)), grep, tolower(twig$charlab)))
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

plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char2", mode="in")))
plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char5", mode="in")))
# plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char8", mode="in")))
plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char10", mode="in")))
# plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char21", mode="in")))
plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char24", mode="in")))

id <- which(V(tree2)$name=="char2")

roots <- which(degree(tree2, v = V(tree2), mode = "in")==0, useNames = T)



lapply(id, all_shortest_paths, graph=tree2, to=roots, mode="in")



# what is a given character most connect to..?







# stems matched to terms
termstem <- lapply(terms2, match, terms3)

names(termstem) <- paste0("term", seq_along(terms2))

head(termstem)

# make network
edges1 <- cbind(rep(names(termstem), times=sapply(termstem, length)), unlist(termstem))
edges1 <- edges1[complete.cases(edges1),]
rownames(edges1) <- NULL
head(edges1)

edges2 <- cbind(rep(seq_along(stemchar), times=sapply(stemchar, length)), paste0("char", unlist(stemchar)))
# edges2 <- edges2[, 2:1]
edges2 <- edges2[complete.cases(edges2), ]
head(edges2)

# stems matched to characters

# assign term with most matched stems to a character

g <- graph_from_edgelist(rbind(edges1, edges2), directed=TRUE)

E(g)


# plot of all terms connect to character 26
plot(induced_subgraph(g, subcomponent(g, V(g)$name=="char26", mode="in")))




# V(g)$name

# plot(g)

charvert <- grep("char", V(g)$name)

termvert <- grep("term", V(g)$name)

d <- distances(g, charvert, termvert)

class(d)

dim(d)

range(d)
image(d)

plot(d[1, ])






# the only question I want to know is:
# for each character, which term does it overlap with the most?

termstem

stemchar

edge_connectivity(g, termvert[2], charvert[2])


# term vertices
grep("term", V(g)$name)

# character vertices
grep("char", V(g)$name)

# plot(g)

unique(unlist(stemchar))

pdf(file = "figure/twig stem word cloud.pdf", width=5, height=5)
wordcloud::wordcloud(unlist(terms2))
dev.off()

# i want to find which characters have the most stems matched to a given term



#### WORK ON THIS








