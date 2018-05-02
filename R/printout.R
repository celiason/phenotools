#' Function to print an html summary of clusters of similar characters
#' x = nexus object with clusters found
#' maxsize = maximum size of cluster of characters
#' Example:
#' data(twig)
#' dups <- duplicated(twig, weighted=TRUE)
#' printout(dups, maxsize=10, file="~/Desktop/twig.html")

# TODO scale text by number of occurrences in cluster/weight
# TODO add color option

printout <- function(x, maxsize=NULL, file, statelabels=TRUE) {
	require(RColorBrewer)
	require(kableExtra)
	dark2 <- brewer.pal(8, "Dark2")
	clust <- x$cluster
	if (!is.null(maxsize)) {
		csize <- sapply(sapply(sapply(clust, str_extract, "\\d+"), na.omit), length)
		keep <- csize < maxsize
		clust <- clust[keep]
	}
	for (i in seq_along(clust)) {
		ids <- as.numeric(na.omit(str_extract(clust[[i]], "\\d+")))
		hlite <- as.character(na.omit(str_extract(clust[[i]], "[^\\d]+")))
		chars <- x$charlab[ids]
		chars0 <- sapply(chars, strsplit, " ")
		chars <- cleantext(chars, fast=FALSE, comments=FALSE)
		if (statelabels) {
			states <- x$statelab[ids]
			states0 <- sapply(states, strsplit, " ")
			states <- cleantext(states, fast=FALSE, comments=FALSE)
			if (!is.list(states)) {
				states <- list(states)
			}
			if (all(x$statelabels=="")) {
				states <- NULL
			}
		} else {
			states <- NULL
		}
		if (!is.list(chars)) {
			chars <- list(chars)
		}
		
		if (!is.null(states)) {
			chars <- lapply(seq_along(chars), function(z) { c(chars[[z]], ": ", states[[z]]) } )
			chars0 <- lapply(seq_along(chars0), function(z) { c(chars0[[z]], ": ", states0[[z]]) } )
		}

		nm <- rep(ids, times=sapply(chars, length))
		names(chars) <- ids
		if (length(hlite) > 8) {
			pal <- colorRampPalette(brewer.pal(11, "Spectral"))(length(hlite))
		} else {
			pal <- dark2[seq_along(hlite)]
		}
		textcol <- rep("black", sum(sapply(chars,length)))
		for (j in seq_along(hlite)) {
		  textcol[unlist(lapply(chars, grepl, pattern=hlite[j]))] <- pal[j]
		}
		textcol <- split(textcol, nm)
		text_formatted <- sapply(seq_along(chars), function(x) {
			paste0("<br/>", text_spec(paste("Character ", ids[x], ": "), "html", color="black", italic=T),
				paste(text_spec(chars0[[x]], "html", color=textcol[[x]]), collapse=" "))
		})
		header <- text_spec(paste0("Cluster ", i), "html", bold=TRUE)
		if (i==1) {
			cat(header, text_formatted, file=file, sep="<br/>", append=F)
		} else {
			cat("<br/><br/>", header, text_formatted, file=file, sep="<br/>", append=T)
		}
	}
}
