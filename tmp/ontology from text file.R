setwd("/Users/chadeliason/Documents/UT/projects/phenome")

# load libraries
library(GGally)
devtools::load_all('~/github/nexustools')  # load packages

# load data
twig <- read.nex("data/theropods/Turner_etal_2012.nex")  # load data
file <- "baumel ontology.txt"  # trait ontology, text file

# create trait ontology from tabbed text file
tree <- read_ontology(file="/Users/chadeliason/Documents/UT/projects/phenome/baumel ontology.txt")
tree <- simplify(tree, edge.attr.comb = "first")

# stem latin words in ontology
V(tree)$name <- schinke(V(tree)$name)

# plot a subtree
subtree <- induced_subgraph(tree, subcomponent(tree, grep("tarsometat", V(tree)$name), mode="out"))

pdf(file = "figs/ontology.pdf", width=9, height=8)
par(mar=c(0,0,0,0))
plot(subtree, layout=-layout.reingold.tilford(subtree)[,2:1], vertex.size=0, edge.arrow.size=0,
	vertex.label.cex=1, vertex.label.color = "black")
dev.off()

ggnet(subtree, label.nodes=T, label.size=2, size=2, arrow.size=.35)
ggsave(file = "figure/ontology_ggnet.pdf")




################################################################################
# 5-31-16 - Working on matching characters to trait ontology using REGEX search
# with stem words
################################################################################

pdf(file = "figure/twig stem word cloud.pdf", width=5, height=5)
wordcloud::wordcloud(unlist(V(tree)$name))
dev.off()



tree2 <- stem_search(tree=tree, x = twig)


pdf(file = "figure/trait_ont_1.pdf", width=5, height=5)
par(mar=rep(0, 4))
plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char2", mode="in")),
	vertex.color="gray", vertex.label.color="black", vertex.label.cex=.5)
dev.off()

pdf(file = "figure/trait_ont_2.pdf", width=5, height=5)
par(mar=rep(0, 4))
plot(induced_subgraph(tree2, subcomponent(tree2, V(tree2)$name=="char24", mode="in")),
	vertex.color="gray", vertex.label.color="black", vertex.label.cex=.5)
dev.off()

id <- which(V(tree2)$name=="char2")

roots <- which(igraph::degree(tree2, v = V(tree2), mode = "in")==0, useNames = T)

# what is a given character most connect to..?
lapply(id, all_shortest_paths, graph=tree2, to=roots, mode="in")










################################################################################
# alternative method - characters matched to stems taked from non-hierarchical
# termlist
################################################################################

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


# distances...maybe find ones closest to term?

# character vertices
charvert <- grep("char", V(g)$name)

# term vertices
termvert <- grep("term", V(g)$name)

# distance matrix calculation
d <- distances(g, charvert, termvert)
dim(d)
head(d)

# the only question I want to know is:
# for each character, which term does it overlap with the most?
edge_connectivity(g, termvert[55], charvert[200])

# i want to find which characters have the most stems matched to a given term

#### WORK ON THIS
