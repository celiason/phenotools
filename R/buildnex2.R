# TOOD - write quick function to build a nexus file from read in matrix + text
#' @param x list containing, minimimally, data matrix (in matrix form)
#' @param data the data matrix
# mattext <- readLines("data/theropods/brusatte2014_matrix.txt")
# chartext <- readLines("data/theropods/brusatte2014_charlist.txt")
# mat <- text2mat(mattext)
# charlabs <- text2charlabels(chartext)
# twig1 <- buildnex2(data=mat, charlabels=charlabs)
buildnex2 <- function(data, charlabels=NULL, charnums=NULL, statelabels=NULL, taxlabels=NULL, filename=NULL) {
	if (is.null(taxlabels)) {
		if (is.null(rownames(data))) {
			taxlabels <- rep('', nrow(data))
		} else {
		taxlabels <- rownames(data)
		}
	}
	if (is.null(statelabels)) {
		statelabels <- rep('', ncol(data))
	}
	if (is.null(charlabels)) {
		charlabels <- rep('', ncol(data))
	}
	if (is.null(charnums)) {
		charnums <- 1:ncol(data)
	}
	if (is.null(filename)) {
		filename <- "file"
	}
	file <- rep(filename, ncol(data))
	dimnames(data) <- NULL
	res <- list(data = data, file = file, taxlabels = taxlabels, charlabels = charlabels, statelabels = statelabels, charnums = charnums, missing = "?", gap="-")
	class(res) <- c("nex", "list")
	res
}
