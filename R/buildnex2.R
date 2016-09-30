# TOOD - write quick function to build a nexus file from read in matrix + text
# x = list containing, minimimally, data matrix (in matrix form)
# data <- mat
buildnex2 <- function(data, charlabels=NULL, charnums=NULL, statelabels=NULL, taxlabels=NULL) {
# taxlabels <- NULL
# data <- mat
# statelabels = NULL
# charnums = NULL
# charlabels = NULL

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
	dimnames(data) <- NULL
	res <- list(data = data, taxlabels = taxlabels, charlabels = charlabels, statelabels = statelabels, charnums = charnums, missing = "?", gap="-")
	class(res) <- c("nex", "list")
	res
}
