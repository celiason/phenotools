# generate user/file-based trait ontology
# x = nexus file
generate_ontology <- function(x) {

	if (is.nex(x)) {
		charlabs <- x$charlabels
	}

	# remove punctuation, comments in brackets, etc. (see cleantext.R function)
	charlabs <- cleantext(charlabs, comma=FALSE)

	# fix muscule and nerve abbreviations
	charlabs <- gsub("m\\.", "musculus", charlabs)
	charlabs <- gsub("n\\.", "nerve", charlabs)

	# locate characters with terms separated by comma
	# shouldn't be followed by and, with, eg, ie
	ss <- grepl(",", charlabs, perl=TRUE) & !grepl(",\\s*(and|with)", charlabs, perl=TRUE)
	
	x.com <- charlabs[ss]

	training <- x[, ss]
	
	test <- x[, !ss]

	# x.com <- x[grepl(",", x)]

	# In some cases, at the end of the character there is a description or
	# comment following a period. This will search for that patterna and
	# truncate the character name.
	x.com <- gsub("\\s\\w{2,}\\.\\s\\w+.*", "", x.com)

	# remove comma followed by 'and', 'with'
	x.com <- gsub(", and", " and", x.com)
	x.com <- gsub(", with", " with", x.com)
	x.com <- gsub(", e.g.", " eg", x.com)
	x.com <- gsub(", i.e.", " ie", x.com)

	terms <- strsplit(x.com, ",")

	terms <- lapply(terms, gsub, pattern="^ ", replacement="")

	terms <- lapply(terms, tolower)

	# remove status or forma
	terms <- lapply(seq_along(terms), function(x) {
		res <- terms[[x]]
		res[!grepl("status|forma", res)]
	})

	# length(unlist(terms))

	# remove term if it is last in series and starts with an adjective
	for (i in (seq_along(terms))) {
		res <- terms[[i]]
		ss <- !grepl("^\\w*(ed|ing|number|tion)\\b", res)
		if (any(!ss)) {
			if (length(ss)==seq_along(ss)[!ss]) {
				terms[[i]] <- res[ss]
			}
		}
	}

	# remove terms starting with '. '
	terms <- lapply(terms, gsub, pattern="\\. ", replacement="")

	terms <- terms[sapply(terms, length)>1]

	# Create unique edges by pasting adjacent matrix cells
	for (i in seq_along(terms)) {
		id <- which(!is.na(terms[[i]]))
		for (j in 2:length(id)) {
			if (length(id) < 2) {
				break()
			} else {
				x <- terms[[i]][id]
				res <- paste0(x[(j-1)], "->", x[j])
				terms[[i]][j] <- res
			}
		}
	}

	# create edge list
	edges <- list()
	for (i in seq_along(terms)) {
		tt <- terms[[i]]
		l <- length(tt)
		edges[[i]] <- lapply(1:(l-1), function(x) {tt[x:(x+1)]})
	}
	edges <- matrix(unlist(edges), ncol=2, byrow=TRUE)

	# create and plot trait ontology
	g <- graph_from_edgelist(edges)

	# split apart nested vertex labels
	V(g)$name <- sapply(str_split(V(g)$name, "->"), tail, 1)

	# use latin stemmer on terms
	V(g)$name <- schinke(V(g)$name)

	list(tree=g, training=training, test=test)

}
