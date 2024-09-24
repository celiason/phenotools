#' Function to create a trait ontology from a tabular text file
#' 
#' @param file path to text file containing a tab-delimited trait ontology
#' @param root whether to add a root node
#' @examples \dontrun{
#' # Read in the ontology of Baumel and Whitmer (1993):
#' ont <- read_ontology(file = system.file("extdata", "baumel_ontology.txt",
#' package = "phenotools"))
#' ont
#' }
#' @details The text file should contain hierarchical data about relationships 
#' among traits (e.g., for feathers: feather -> rachis -> barb -> barbules ->
#' barbule hooklets)
#' @references Baumel, J. J., and L. M. Witmer. 1993. Osteologia. P. in
#' Handbook of avian anatomy: nomina anatomica avium. Publications of the
#' Nuttall Ornithological Club.
#' @importFrom zoo na.locf
#' @export
#' 
read_ontology <- function(file, root=FALSE){
	
	fields <- max(utils::count.fields(file, sep = "\t"))
	
	dat <- utils::read.table(file, sep="\t", col.names = paste0("V", sequence(fields)), header=FALSE, fill=TRUE, strip.white=TRUE, stringsAsFactors=FALSE, na.strings="")
	
	# To prepare the data carry forward the last value in columns if lower level (col to the right) is non-missing
	dat[1] <- zoo::na.locf(dat[1], na.rm=FALSE)
	for(i in ncol(dat):2)  {
	  dat[[i-1]] <-  ifelse(!is.na(dat[[i]]), zoo::na.locf(dat[[i-1]], na.rm=F), dat[[i-1]])
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
				na.omit(stats::setNames(cbind(dat[(1 + i) : (2 + i)], rownames(dat)), c(names(dat[1:2]), "sort")))
			})))

	edges <- as.matrix(edges)
	
	# create graph
	g <- igraph::graph_from_edgelist(edges[, 1:2])

	# sort attribute for sorting characters that match to a given term
	igraph::edge.attributes(g)$sort <- edges[, 3]
	
	# make new vertex labels
	verts <- V(g)$name
	verts <- tolower(verts)
	vertsold <- verts
	verts <- strsplit(verts, "->")
	verts <- sapply(seq_along(verts), function(x) tail(verts[[x]], 1))
	verts <- gsub("^ ", "", verts)
	
	V(g)$name <- verts
	
	if (root) {
		g <- addroot(g)
	}

	# return graph
	simplify(g)

}

# add root vertex
addroot <- function(g) {
	roots <- which(igraph::degree(g, v = V(g), mode = "in")==0, useNames = T)
	if (length(roots)==1) {
		stop("Root already present")
		return(g)
	}
	g <- igraph::add_vertices(g, nv=1, name="root")
	rootid <- which(V(g)$name=="root")
	for (i in roots) {
		g <- igraph::add_edges(g, edges=c(rootid, i))
	}
	return(g)
}
