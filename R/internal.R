# Function to go from this "2,5-7,10,12-15" to this "c(2,5,6,7,10,12,13,14,15)"
# see http://r.789695.n4.nabble.com/convert-delimited-strings-with-ranges-to-numeric-td4673763.html
text2numeric <- function(xx) {
	xx <- gsub("^\\s*|\\s*$", "", xx)
  xx <- gsub('\\s|,\\s', ',', xx)
	xx <- gsub('\\-', ':', xx)
	eval(parse(text = paste("c(", xx, ")")))
}


# Function for finding starts and lengths of string/number sequences
seqle <- function(x, incr=1) {
	if(!is.numeric(x)) {
		x <- as.numeric(x)
	}
	n <- length(x)
	y <- x[-1L] != x[-n] + incr
	i <- c(which(y|is.na(y)),n)
	list(lengths = diff(c(0L,i)), values = x[utils::head(c(0L,i)+1L,-1L)])
}

# convert string of scorings to a matrix
text2mat <- function(x) {
	x <- readLines(x)
	if (any(stringr::str_detect(x, "[Mm]atrix"))) {
		x <- x[!grepl("[Mm]atrix", x)]
	}
	# scorings have to be separated from taxon names by either a tab or 2+ spaces (as in Xu et al. 2011)
	x <- str_match(x, "(.*?)(\t|\\s)(.*)")
	# scorings have to be separated from taxon names by either a tab or 2+ spaces (as in Brusatte et al. 2014)
	# x <- str_match(x, "(.*?)(\t|\\s{2,})(.*)")
	nms <- x[, 2]
	scores <- x[, 4]
	res <- str_extract_all(scores, '\\d{1}|[\\(\\[\\{]\\d{1,4}[\\)\\]\\}]|\\-|\\?')
	names(res) <- nms
	if (!sd(sapply(res, length))==0) {
		stop("There is a problem with reading the text. Not an equal number of characters for each taxon.")
	}
	res <- do.call("rbind", res)
	res
}

# get character labels from text
text2charlabels <- function(x) {
	x <- readLines(x)
	# matches <- str_match(x, "^(Character|[\\s\\t]*)(\\d+)[\\.]?(.*)(\\[\\[\\]\\])(.*)")
	# DONE: make it so this selects all text on multiple lines
	# TODO: make it so we can extract characters and character states separately (e.g., separated by a ":")
	# for Brusatte 2014
	if (any(stringr::str_detect(x, "^[Cc]haracter"))) {
		x <- paste0(x, collapse="\n")
		# matches <- str_match_all(x, regex("Character\\s*(\\d*)\\:\\s(.*?)", dotall=TRUE))
		matches <- str_match_all(x, regex("Character\\s*(\\d*)\\:\\s(.+?)(?=Character\\s*\\d*\\:\\s|$)", dotall=TRUE))
		matches <- matches[[1]]
		charnums <- as.numeric(matches[, 2])
		# id <- which(!is.na(as.numeric(matches[, 2])))
		# charnums <- na.omit(as.numeric(matches[, 2]))
		# seqs <- sapply(seq_along(id), function(i) {
		# 	if (id[i] == max(id) ) {
		# 		id[i]:length(x)
		# 	} else {
		# 		id[i]:(id[i+1]-1)
		# 	}
		# })
		# charlabs <- sapply(seq_along(seqs), function(i) {
		# 	paste0(x[seqs[[i]]], collapse=" ")
		# })
		charlabs <- matches[, 3]
	} else {
		# for Xu 2011
		matches <- str_match(x, regex("(Character\\s*)?(\\d*)[\\.]?[\\:]?[\\s]?(.*?)$", dotall=TRUE))
		charnums <- as.numeric(matches[, 3])
		charlabs <- matches[, 4]
	}
	res <- stats::setNames(charlabs, charnums)
	res <- res[!is.na(res)]
	res <- gsub("\n", " ", res)
	res
}


# stemmer function based on Schinke et al. 1996
# TODO make it so things like -ity, -ed will be removed from end of word
schinke <- function(x) {
	x <- tolower(x)
	x <- gsub("m\\.", "muscul", x)
	x <- gsub("proc\\.", "processus", x)
	x <- gsub("(lig|ligg)\\.", "ligamentos", x)
	x <- gsub("n\\.", "nervos", x)
	# x <- gsub("j", "i", x)
	# x <- gsub("v", "u", x)
	# x <- str_replace_all(x, "(ibus|ius|ae|am|as|em|es|ia|is|nt|os|ud|um|us|a|e|i|o|u)\\b", "")
    # NEW VERSIon  (not exactly Schinke)
	x <- str_replace_all(x, "(ity|ed|al|ibus|ius|ae|am|as|em|es|ia|is|nt|os|ud|um|us|a|e|i|o|u)\\b", "")
	x
}
