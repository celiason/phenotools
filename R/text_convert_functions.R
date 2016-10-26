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
	# DONE: make it so this selects all text on multiple lines
	# TODO: make it so we can extract characters and character states separately (e.g., separated by a ":")
	# for Brusatte 2014
	if (any(str_detect(x, "^[Cc]haracter"))) {
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
	res <- setNames(charlabs, charnums)
	res <- res[!is.na(res)]
	res <- gsub("\n", " ", res)
	res
}

# account for character states between character labels in pasted text
# x <- chartext[1:10]
# x <- chartext
# nums <- as.numeric(str_match(x, "([Cc]haracter\\s)?(\\d+)")[, 3])
# nums[diff(nums) != 1]
# nums[200:220]
# min(head(nums))+1
# max(tail(nums))
# grep(853, nums)[1]


# plot(nums)
# plot(diff(diff(nums)))
# # nums <- na.omit(nums)
# # diff(cumsum(nums))

# plot(c(1,0,1,2,0,1,3,0,1,2,4,0,1,2,3), type='l')


# TODO
# work with Brusatte format to get state labels included with character labels
