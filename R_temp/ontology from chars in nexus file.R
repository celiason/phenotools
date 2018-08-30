setwd("/Users/chadeliason/Documents/UT/projects/phenome")

# load stuff

devtools::load_all('~/github/nexustools')

library(igraph)

# palaeognathae character list
final <- read.nex("/Users/chadeliason/Documents/UT/projects/phenome/output/final_reordered.nex")

res <- generate_ontology(final)
tree <- res$tree
training <- res$training
test <- res$test

tree2 <- delete_vertices(tree, which(igraph::degree(tree)==1))
tree2 <- simplify(tree2)


searchtree(graph=tree2, pattern="femur", vertex.label.cex=.5, vertex.color = "black", vertex.size=0, edge.arrow.size=.5, vertex.label.color = "darkblue")

subtree <- searchtree(graph=tree, pattern="extremit proximal carpometacarp", plot=FALSE)

subtree <- searchtree(graph=tree, pattern="tarsometa", plot=FALSE)

subtree

V(subtree)$name <- gsub(" ", "\n", V(subtree)$name)

pdf(file = "figure/ontology_autogen_tarsometatars.pdf")
# plot(subtree, vertex.label.cex=.5, vertex.color = "black", vertex.size=0, edge.arrow.size=.5, vertex.label.color = "darkblue")
plot(subtree, layout = -layout_as_tree(subtree)[,2:1], vertex.label.cex=.5,
	vertex.color = "black", vertex.size=0, edge.arrow.size=.5,
	vertex.label.color = "darkblue")
dev.off()



# find number of characters linked to TO
system.time(res <- stem_search(tree=tree2, x=test))
table(igraph::degree(graph=res, v=grep("char", names(V(res))), mode="in")>1)
# 624/(624+28)
plot(induced_subgraph(res, subcomponent(res, V(res)$name=="char5", mode="in")),
	vertex.color="gray", vertex.label.color="black", vertex.label.cex=.5)






# visualize unique terms
wordcloud::wordcloud(V(subtree)$name)

# all terms
pdf(file="figure/palaeognath_wordcloud.pdf")
wordcloud::wordcloud(V(tree)$name)
dev.off()


# toy example
test <- c("Brain, cerebrum, shape", "Hindlimb, tarsometatarsus, length", "Brain, cerebellum", "Cerebellum elongated", "Tarsometatarsus, distal part", "Hindlimb")
tree <- generate_ontology(test)
plot(tree, layout = layout_as_tree)


# search for terms not in ontology
test.nocom <- test[!grepl(",", test)]
test.nocom <- tolower(test.nocom)
roots <- which(degree(tree, v = V(tree), mode = "in")==0, useNames = T)

# locate matches...
matches <- sapply(V(tree)$name, grep, test.nocom)

id <- which(matches>=1)

# locate all paths in ontology connected to a certain character
lapply(id, all_shortest_paths, graph=tree, to=roots, mode="in")

# sort??? interleave???

