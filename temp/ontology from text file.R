setwd("/Users/chadeliason/Documents/UT/projects/phenome")

# load packages
devtools::load_all('~/github/nexustools')

# load data
twig <- read.nex("data/theropods/Turner_etal_2012.nex")
file <- "baumel ontology.txt"



# USING TEXT (TAB-BASED) ONTOLOGY
tree <- treechart(file="/Users/chadeliason/Documents/UT/projects/phenome/baumel ontology.txt")

tree <- simplify(tree, edge.attr.comb = "first")

# stem latin words in ontology
V(tree)$name <- schinke(V(tree)$name)

E(tree)$sort

# terms <- V(tree)$name
terms <- schinke(V(tree)$name)
terms2 <- str_split(terms, "->|\\s")
terms2 <- sapply(terms2, unique)

# remove empties
any(unlist(sapply(terms2, "==", "")))
# terms2 <- lapply(seq_along(terms2), function(x) {terms2[[x]][terms2[[x]]!=""]})
# any(unlist(sapply(terms2, "==", "")))
terms3 <- unique(unlist(terms2))

# remove dots
terms3 <- terms3[!grepl("\\.", terms3)]

# remove single letters
terms3 <- terms3[grepl("\\w{2,}", terms3)]
length(terms3)  # number of unique terms

terms3


# plot a subtree
subtree <- induced_subgraph(tree, subcomponent(tree, grep("tarsometat", V(tree)$name), mode="out"))

pdf(file = "figure/ontology.pdf", width=9, height=8)
par(mar=c(0,0,0,0))
plot(subtree, layout=-layout.reingold.tilford(subtree)[,2:1], vertex.size=0, edge.arrow.size=0)
dev.off()


## USING TREE:
# create groups of words based on trait ontology/tree
# find root and leaves
leaves <- which(degree(tree, v = V(tree), mode = "out")==0, useNames = T)
roots <- which(degree(tree, v = V(tree), mode = "in")==0, useNames = T)
# traverse tree and get all combinations of characters along trait ontology
reachable <- lapply(roots, function(x) {which(shortest.paths(tree, x, mode="out") != Inf)})
terminal.nodes <- lapply(reachable, function(x) {x[which(degree(tree, x, mode="out") == 0)]})
traversal <- lapply(seq_along(roots), function(x) {
	paths <- get.all.shortest.paths(graph=tree, from=roots[x], to=terminal.nodes[[x]], mode="out")$res
	sapply(paths, function(vs) paste(V(tree)[vs]$name, collapse="->"))
})
terms <- unlist(traversal)
length(terms)

# order of terms??



################################################################################
# 
################################################################################

# one option for testing similarity/presence of term in character list:
library(stringdist)
# res <- matrix(NA, nrow=length(twig$charlabels), ncol=length(terms))
# for (i in seq_along(terms)) {
# 	for (j in seq_along(terms))
# 	stringdist(twig$charlab, terms[[j]], method="lcs")
#   # res[, i] <- str_detect(twig$charlab, paste0(terms[[i]], collapse="|"))
# }
# dim(res)
# image(res)
# twig$charlab[which(res)]





################################################################################
# 5-31-16
# Working on matching characters to trait ontology using REGEX search w/stem words
################################################################################

# match each stem to each character using REGEX
stemchar <- sapply(terms3, grep, twig$charlab)

# only at beginning of word
stemchar <- lapply(paste0("\\b", terms3), grep, tolower(twig$charlab))
length(stemchar)
names(stemchar) <- terms3
tail(sort(sapply(stemchar, length)))

# number of unmatched characters
twig$charlab[setdiff(seq_along(twig$charlab), unique(unlist(stemchar)))]

107/477 # 22% unmatched


# stems matched to terms
termstem <- lapply(terms2, match, terms3)

names(termstem) <- paste0("term", seq_along(terms2))

termstem

edges1 <- cbind(rep(names(termstem), times=sapply(termstem, length)), unlist(termstem))
edges1 <- edges1[complete.cases(edges1),]
rownames(edges1) <- NULL
head(edges1)


# stems matched to characters

stemchar

# assign term with most matched stems to a character

length(stemchar)

edges2 <- cbind(rep(seq_along(stemchar), times=sapply(stemchar, length)), paste0("char", unlist(stemchar)))
edges2 <- edges2[, 2:1]
edges2 <- edges2[complete.cases(edges2), ]
head(edges2)

g <- graph_from_edgelist(rbind(edges1, edges2), directed=FALSE)

# V(g)$name

# plot(g)

charvert <- grep("char", V(g)$name)

termvert <- grep("term", V(g)$name)

d <- distances(g, charvert, termvert)

class(d)

dim(d)

range(d)
image(d)

length(charvert)
length(termvert)
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








