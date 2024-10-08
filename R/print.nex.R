#' Print data in a `nex` object
#' 
#' @param x nex class object
#' @param width width of character string to display
#' @param ... other arguments (not used)
#' 
#' @export
#' 
print.nex <- function(x, width=6, ...) {
	# nchar <- ncol(x$data)
	# ntax <- nrow(x$data)
	# cat('NEXUS data file with', nchar, ' characters and', ntax, 'taxa\n')
	res <- x$data
	rownames(res) <- x$taxlabels
	colnames(res) <- substr(x$charlabels, 1, width)
	print(res)
}

#' Summarize data in a `nex` object
#' 
#' @param object `nex` class object
#' @param ... other arguments (not used)
#' 
#' @export
#' 
summary.nex <- function(object, ...) {
  nchar <- ncol(object$data)
  ntax <- nrow(object$data)
  cat('NEXUS data file with', nchar, 'characters and', ntax, 'taxa\n')
}

#' Look at top of a `nex` object
#' 
#' @param x nex object
#' @param width width of character string to display
#' @param ... optional arguments
#' 
#' @importFrom utils head
#' 
#' @export
#' 
head.nex <- function(x, width=6, ...) {
	res <- x$data
	rownames(res) <- x$taxlabels
	colnames(res) <- substr(x$charlabels, 1, width)
	head(res, ...)
}

#' Look at bottom of a `nex` object
#' 
#' @param x `nex` class object
#' @param width width of character string to display
#' @param ... optional arguments
#' 
#' @importFrom utils head
#' 
#' @export
#' 
tail.nex <- function(x, width=6, ...) {
	res <- x$data
	rownames(res) <- x$taxlabels
	colnames(res) <- substr(x$charlabels, 1, width)
	# res[res==x$missing] <- "."
	tail(res, ...)
}

#' Dimensions of a `nex` object
#' 
#' @param x `nex` class object
#' @export
#' 
dim.nex <- function(x) {
	dim(x$data)
}


