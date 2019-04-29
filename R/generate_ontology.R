#' Generate a trait ontology from a nexus file
#' 
#' Comma separated terms are used to define links among components of a nexus
#' character list
#' 
#' @param x a `nex` object containing character descriptions
#' @examples \dontrun{
#' x <- read.nex(system.file("extdata", "clarke_2006.nex", package = "phenotools"))
#' # plot data matrix:
#' plot(x)
#' # generate ontology based on comma-separated terms in characters:
#' ont <- generate_ontology(x)
#' # plot ontology for humerus characters:
#' searchtree(ont$tree, pattern="humer")
#' }
#' 
#' @importFrom stringr str_split
#' @importFrom utils tail
#' @import igraph
#' 
#' @export
#' 
generate_ontology <- function(x) {

	if (is.nex(x)) {
		charlabs <- x$charlabels
	}

	# Locate comma-separated terms not followed by and, with, eg, or ie:
	ss <- grepl(",", x$charlabels, perl=TRUE) & !grepl(",\\s*(and|with)", x$charlabels, perl=TRUE)

	# Subset for further use
	x.com <- x$charlabels[ss]
	training <- x[, ss]
	test <- x[, !ss]

	# x.com <- x[grepl(",", x)]

	# In some cases, at the end of the character there is a description or
	# comment following a period. This will search for that patterna and
	# truncate the character name.
	# x.com <- gsub("\\s\\w{2,}\\.\\s\\w+.*", "", x.com)

	# Remove comments (within "(", "["):
	x.com <- stringr::str_replace_all(x.com, "[\\s]*[\\[\\(].*?[\\]\\)][\\s]*", " ")

	# remove comma followed by 'and', 'with'
	x.com <- gsub(", and", " and", x.com)
	x.com <- gsub(", with", " with", x.com)
	x.com <- gsub(", e.g.", " eg", x.com)
	x.com <- gsub(", i.e.", " ie", x.com)

	terms <- strsplit(x.com, ",")

	# terms <- sapply(terms, cleantext, fast=F, latin=F)

	# remove punctuation, comments in brackets, etc. (see cleantext function)
	# charlabs <- sapply(cleantext(charlabs), as.character)

	# fix muscule and nerve abbreviations
	# charlabs <- gsub("m\\.", "musculus", charlabs)
	# charlabs <- gsub("n\\.", "nerve", charlabs)

	# Remove spaces at beggining of characters:
	terms <- lapply(terms, gsub, pattern="^ ", replacement="")

	# Make lowercase
	terms <- lapply(terms, tolower)

	# remove status or forma
	# terms <- lapply(seq_along(terms), function(x) {
	# 	res <- terms[[x]]
	# 	res[!grepl("status|forma", res)]
	# })

	# remove term if it is last in series and starts with an adjective
	# for (i in (seq_along(terms))) {
	# 	res <- terms[[i]]
	# 	ss <- !grepl("^\\w*(ed|ing|number|tion)\\b", res)
	# 	if (any(!ss)) {
	# 		if (length(ss)==seq_along(ss)[!ss]) {
	# 			terms[[i]] <- res[ss]
	# 		}
	# 	}
	# }

	# Remove terms starting with '. '
	terms <- lapply(terms, gsub, pattern="\\. ", replacement="")

	# Keep characters more than 1 term separated by commas
	terms <- terms[sapply(terms, length) > 1]

	edgeterms <- terms

	# Create unique edges by pasting adjacent matrix cells
	for (i in seq_along(edgeterms)) {
		id <- which(!is.na(edgeterms[[i]]))
		for (j in 2:length(id)) {
			if (length(id) < 2) {
				break()
			} else {
				x <- edgeterms[[i]][id]
				res <- paste0(x[(j-1)], "->", x[j])
				edgeterms[[i]][j] <- res
			}
		}
	}

	# create edge list
	edges <- list()
	for (i in seq_along(edgeterms)) {
		tt <- edgeterms[[i]]
		l <- length(tt)
		edges[[i]] <- lapply(1:(l-1), function(x) {tt[x:(x+1)]})
	}
	edges <- matrix(unlist(edges), ncol=2, byrow=TRUE)

	# create and plot trait ontology
	g <- igraph::graph_from_edgelist(edges)

	# split apart nested vertex labels
	igraph::V(g)$name <- sapply(stringr::str_split(igraph::V(g)$name, "->"), tail, 1)

	# use latin stemmer on terms
	igraph::V(g)$name <- schinke(igraph::V(g)$name)

	list(tree = g, training = training, test = test)

}
