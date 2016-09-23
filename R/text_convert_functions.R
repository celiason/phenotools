# convert string of scorings to a matrix
text2mat <- function(x) {
	if (any(str_detect(x, "[Mm]atrix"))) {
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
	# matches <- str_match(x, "^(Character|[\\s\\t]*)(\\d+)[\\.]?(.*)(\\[\\[\\]\\])(.*)")
	if (any(str_detect(x, "^[Cc]haracter"))) {
		# for Brusatte 2014
		matches <- str_match(x, "Character\\s*(\\d*)\\:\\s(.*)")	
		charnums <- as.numeric(matches[, 2])
		charlabs <- matches[, 3]	
	} else {
		# for Xu 2011
		matches <- str_match(x, "^(Character\\s*)?(\\d*)[\\.]?[\\:]?[\\s]?(.*)")	
		charnums <- as.numeric(matches[, 3])
		charlabs <- matches[, 4]
	} 
	res <- setNames(charlabs, charnums)
	res <- res[!is.na(res)]
	res
}

# TODO
# work with Brusatte format to get state labels included with character labels
