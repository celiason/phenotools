# function to create network graph from tab-based text file
# file = path to text file containing trait ontology
# e.g.,
# feather
# 	rachis
#		barb
# 		barbule
# 			hooklet
#
read_ontology <- function(file){
	library(psych)
	library(zoo)
	library(igraph)
	fields <- max(count.fields(file, sep = "\t"))
	dat <- read.table(file, sep="\t",col.names = paste0("V", sequence(fields)), header=FALSE, fill=TRUE, strip.white=TRUE, stringsAsFactors=FALSE, na.strings="")
	# To prepare the data carry forward the last value in columns if lower level (col to the right) is non-missing
	dat[1] <- na.locf(dat[1], na.rm=FALSE)
	for(i in ncol(dat):2)  {
	  dat[[i-1]] <-  ifelse(!is.na(dat[[i]]), na.locf(dat[[i-1]], na.rm=F), dat[[i-1]])
	}
	# Create unique edges by pasting adjacent matrix cells
	for (i in 1:nrow(dat)) {
		id <- which(!is.na(dat[i, ]))
		for (j in 2:length(id)) {
			if (length(id) < 2) {
				break()
			} else {
			dat[i, j] <- paste0(dat[i, (j-1)], "->", dat[i, j])
			}
		}
	}
	# Get edges for graph. Sort taken from order in ontology
	edges <- rbind(na.omit(cbind(dat[1:2], sort = rownames(dat))),
		do.call('rbind',
			lapply(1:(ncol(dat) - 2), function(i) {
				na.omit(setNames(cbind(dat[(1 + i) : (2 + i)], rownames(dat)), c(names(dat[1:2]), "sort")))
			})))

	edges <- as.matrix(edges)
	
	# create graph
	g <- graph_from_edgelist(edges[, 1:2])
	
	# sort attribute for sorting characters that match to a given term
	edge.attributes(g)$sort <- edges[, 3]
	
	# make new vertex labels
	verts <- V(g)$name
	verts <- tolower(verts)
	verts <- strsplit(verts, "->")
	verts <- sapply(seq_along(verts), function(x) tail(verts[[x]], 1))
	verts <- gsub("^ ", "", verts)
	V(g)$name <- verts

	# return graph
	g

}
