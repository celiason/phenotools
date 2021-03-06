#' Function to print an html summary of clusters of similar characters
#' 
#' Function prints out a list of clusters of similar characters with terms color-
#' coded to illustrate overlap among similar characters
#' 
#' @param x = nexus object with clusters found
#' @param maxsize = maximum size of cluster of characters
#' @param file name of output file (html format)
#' @param statelabels whether to include state labels in output
#' 
#' @examples \dontrun{
#' data(twig)
#' dups <- duplicated(twig, opt="terms")
#' printout(dups, maxsize=10, file="twig.html")
#' }
#' 
#' @export
#' 
printout <- function(x, maxsize=NULL, file, statelabels=TRUE) {
	# TODO scale text by number of occurrences in cluster/weight
	# TODO add color option
	dark2 <- brewer.pal(8, "Dark2")
	clust <- x$cluster
	clust <- clust[stringr::str_count(string=clust, pattern="\\d+")!=0]
	if (!is.null(maxsize)) {
		csize <- sapply(sapply(sapply(clust, stringr::str_extract, "\\d+"), na.omit), length)
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
			states <- x$statelabels[ids]
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
			pal <- colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(length(hlite))
		} else {
			pal <- dark2[seq_along(hlite)]
		}
		textcol <- rep("black", sum(sapply(chars, length)))
		for (j in seq_along(hlite)) {
		  textcol[unlist(lapply(chars, grepl, pattern=hlite[j]))] <- pal[j]
		}
		textcol <- split(textcol, nm)
		text_formatted <- sapply(seq_along(chars), function(x) {
			paste0("<br/>", kableExtra::text_spec(paste("Character ", ids[x], ": "), "html", color="black", italic=T),
				paste(kableExtra::text_spec(chars0[[x]], "html", color=textcol[[x]]), collapse=" "))
		})
		header <- kableExtra::text_spec(paste0("Cluster ", i), "html", bold=TRUE)
		if (i==1) {
			cat(header, text_formatted, file=file, sep="<br/>", append=F)
		} else {
			cat("<br/><br/>", header, text_formatted, file=file, sep="<br/>", append=T)
		}
	}
}
