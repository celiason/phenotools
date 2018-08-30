#' Code for combining character lists and data matrices
#' 
#' 
#' 
#' @examples \dontrun{
#' x = nexus file (.nex object in R)
#' file = text file with character list (numbered)
#' }
#' 
#' @importFrom stringr str_match
#' @importFrom stringr str_match_all
#'
addmeta <- function(x, file) {

	# read in character labels
	y <- readLines(file)

	# search for pattern
	matches <- str_match(y, "^[\\s\\t]*(\\d+)[\\.]?(.*)(\\[\\[\\]\\])(.*)")

	# remove NAs
	id <- !is.na(matches[, 2])
	matches <- matches[id, ]
	charnums <- as.numeric(matches[,2])

	# extract character and state labels
	newcharlabs <- matches[, 3]
	newcharlabs <- gsub("^\\s+", "", newcharlabs)
	newcharlabs <- gsub("\\s+$", "", newcharlabs)
	newcharlabs <- gsub("'", "`", newcharlabs)

	newstatelabs <- matches[, 5]
	newstatelabs <- gsub("^\\s+", "", newstatelabs)
	newstatelabs <- gsub("\\s+$", "", newstatelabs)

	# find with numbers after state label
	newstatelabs <- str_match_all(newstatelabs, "(.*?)\\(\\d\\)")

	for (i in seq_along(newstatelabs)) {
		# remove punctuation at beginning
		newstatelabs[[i]][, 2] <- gsub("^([\\;]?\\s*|\\s*[\\;]?)", "", newstatelabs[[i]][, 2])
		# remove punctuation at end
		newstatelabs[[i]][, 2] <- gsub("([\\;]?\\s*|\\s*[\\;]?)$", "", newstatelabs[[i]][, 2])
	}

	newstatelabs <- lapply(newstatelabs, "[", , 2)

	newstatelabs <- lapply(seq_along(newstatelabs), function(x) {gsub("'", "`", newstatelabs[[x]])})

	newstatelabs <- sapply(seq_along(newstatelabs), function(x) { paste0("'", newstatelabs[[x]], "'", collapse=" ") })

	# check if right number of characters
	if (!length(newcharlabs) == ncol(x$data)) {
		stop("Number of characters in list not equal to characters in the nexus object")
	}



	# combine matrix and character names

	res <- x

	if (is.null(x$charlabels) | all(x$charlabels == "''")) {
		res$charlabels <- newcharlabs
	}

	if (is.null(x$statelabels) | all(x$statelabels == "''")) {
		res$statelabels <- newstatelabs
	}


	# output and check file

	res

}


# TEST
# read in data matrix
# x <- read.nex(file = "/Users/chadeliason/Documents/UT/projects/phenome/data/theropods/Turner_etal_2012 copy.nex")
# x2 <- addmeta(x, file = "~/Documents/UT/projects/phenome/data/theropods/mbank_X1507_5-25-2016_129_character_list.txt")
# write.nex(x2, file = "~/Desktop/test_addmeta.nex")
