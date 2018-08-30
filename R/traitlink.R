#' Function to link character descriptions to a trait ontology
#' 
#' Function uses fuzzy text matching to output a list of potentially redundant
#' characters in a `nex` object, waits for user input to confirm, and then
#' handles merging of characters and outputs a new `nex` object
#' 
#' @param tree trait ontology (igraph object)
#' @param x nexus file
#' 
#' @return an object of class XX
#' 
#' @author Chad M. Eliason
#' 
#' @examples \dontrun{
#' x <- twig2.ss
#' tree <- tree
#' x=nex
#' }
#' 
#' @export
#' 
traitlink <- function(tree, x) {
	
	# roots of trait ontology
	roots <- which(degree(tree, v = V(tree), mode = "in")==0, useNames=TRUE)

	# all patchs from roots to tips
	# isn't there a faster way to do this?
	cat("Finding paths along trait ontology...\n")
	paths <- pbapply::pblapply(V(tree), all_shortest_paths, graph=tree, to=roots, mode="in")
	paths <- sapply(paths, "[[", "res")
	
	# combine nodes/names along paths (for regexp searching later)
	fullnames <- sapply(seq_along(paths), function(i) {paste0(names(paths[[i]]), collapse="->")})
	fullnames <- gsub("\\[|\\]|\\(|\\)", "", fullnames)
	fullnames <- gsub("^ ", "", fullnames)
	fullnames <- gsub(" $", "", fullnames)

	# Whether to match to using:
	# full names (e.g., axial skeleton->skull... or just 'skull'):
	# terms <- fullnames
	# or partial names:
	terms <- V(tree)$name
	terms <- sapply(cleantext(terms), as.character)
 	terms <- str_split(terms, "->|\\s")
	terms <- sapply(terms, unique)

	# remove dots
	terms <- lapply(seq_along(terms), function(x) {
		res <- terms[[x]]
		res[!grepl("\\.", res)]
	})

	# remove single, two letter words
	terms <- lapply(seq_along(terms), function(x) {
		res <- terms[[x]]
		res[grepl("\\w{3,}", res)]
	})

	# only at beginning of word
	# stemchar <- lapply(paste0("\\b", terms3), grep, tolower(x$charlab))
	# names(stemchar) <- terms3
	# tail(sort(sapply(stemchar, length)))

	# stemchar is a named list with terms that match up to different characters
	# now we want to make connections between unique terms and all terms in 'fullnames'
	# sapply(terms3, grep, terms)
	# grep(terms3[1], terms)
	# terms3[1]
	# terms[1152]
	# stemchar[!grepl("integer", stemchar)]

	# might want to have function return these results:
	# number of unmatched characters
	# length(x$charlab[setdiff(seq_along(x$charlab), unique(unlist(stemchar)))])
	# 119/477 # 25% unmatched
	# kinda slower (~14X) than only using unique word stems (terms3)
	
	cleanchars <- sapply(cleantext(x$charlab, fast=TRUE), as.character)

	############################################################################
	## TODO try something else (e.g., using text distance)

	# (1) string distance option
	# m <- stringdistmatrix(cleanchars, terms, method="jw")

	# (2) term overlap option
	termlist <- tm::Corpus(VectorSource(c(cleanchars, terms)))
	dtm <- tm::TermDocumentMatrix(termlist, control = list(wordLengths = c(2, Inf), weighting = function(x) tm::weightTfIdf(x, normalize = TRUE), stemming=FALSE))
	m <- as.matrix(dtm)

	cat("Creating links between characters and trait ontology...\n")
	g <- igraph::graph_from_incidence_matrix(m, weighted=TRUE)

	# finding clusters (turned off for now)
	# cluster <- "infomap"
	# clustfun <- get(paste0("cluster_", cluster))
	# cl <- clustfun(g, E(g)$weight)
	# grps <- communities(cl)  # get clusters

	# create IDs for characters (id1) and trait ontology terms (id2)
	id1 <- seq_along(cleanchars)
	id2 <- (length(cleanchars)+1) : (length(cleanchars) + length(terms))

	# lists of vertices corresponding to characters (v1) and ontology terms (v2)
	v1 <- which(V(g)$name %in% id1)
	v2 <- which(V(g)$name %in% id2)

	# find clusters containing both a trait term and a character
	# keep <- sapply(seq_along(grps), function(x) {
	# 	any(grps[[x]] %in% id1) & any(grps[[x]] %in% id2)
	# })
	# grps <- grps[keep]

	# % of character "linked in" to trait ontology
	# sum(unlist(grps) %in% id1) / length(id1)

	# make a cool plot:
	# V(g)$color <- "gray"
	# V(g)$color[V(g)$name %in% id1] <- "lightblue"
	# V(g)$color[V(g)$name %in% id2] <- "wheat"
	# sg <- make_ego_graph(g, order=2, nodes=V(g)$name=="1")[[1]]
	# pdf("figs/traitlink_palaeognaths.pdf")
	# par(mar=rep(0,4))
	# plot(sg, vertex.label.cex=.3, vertex.size=5)
	# legend("topleft", bty="n", legend=c("character", "ontology", "term"), pt.bg=c("lightblue", "wheat", "gray"), pch=21)
	# dev.off()

	# need to compute distances from terms to ontology, max distance of 2 vertices!!
	# find neighborhood
	cat("Locating matches...\n")
	picks <- pbapply::pbsapply(v1, function(i) {
		nn <- ego(graph=g, order=2, nodes=i)[[1]]
		tovert <- v2[V(g)$name[v2] %in% names(nn)]
		d <- distances(g, v=i, to=tovert)[1, ]
		d[d==Inf] <- NA
		names(which.max(d))
	})

	picks[sapply(picks, is.null)] <- NA

	picks <- unlist(picks)

# diff(range(as.numeric(picks), na.rm=T))

# what to output:

# 1) character descriptions

	# x$charlab[1:3]  # characters

# 2) original trait ontology tree

	# tree

# 3) trait ontology matches (text)

	ontmatch <- fullnames[match(picks, id2)]  # matches in trait ontology

# 4) trait ontology matches (branch in tree)

	# match(picks, id2)

# 5) distance

	# add a root to the tree (to avoid Inf distances where chars not connected)
	tree2 <- addroot(tree)

	allverts <- as.numeric(picks) - length(v1)  # all vertices
	univerts <- unique(allverts[!is.na(allverts)])  # unique vertices
	dmat <- distances(tree2, v=univerts, to=univerts, mode="all")
	nms <- V(tree2)$name[allverts]
	nms <- nms[!is.na(nms)]

	alldmat <- dmat[nms, nms]
	colnames(alldmat) <- seq_along(allverts)[!is.na(allverts)]
	rownames(alldmat) <- seq_along(allverts)[!is.na(allverts)]

	rootid <- which(V(tree2)$name=="root")
	nodedepth <- sapply(allverts, function(x) {
		tryCatch(distances(tree2, v=rootid, to=x), error=function(e) NA)
	})

	res <- list(chars=x$charlab, dmat=alldmat, ont=ontmatch, nodedepth=nodedepth)

	res

}