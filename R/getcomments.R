#' Function to find comments in text file
#' 
#' This function reads in a character list that have been annotated
#' as to whether characters are duplicates (see Details)
#' @param file name of text file
#' @details Text file contains list of characters with flags indicating whether
#' characters are duplicates (DD), should be cut (XX), good characters that
#' should be kept (KK), or characters that should be kept but modified (KKM).
#' Comments in curly braces indicate which other characters a focal character
#' duplicates (e.g., "{duplicates 125, 250}").
#' @importFrom stringr str_trim
#' @importFrom stringr regex
#' @importFrom stringr str_locate
#' @importFrom stringr str_extract_all
#' @export
#'
getcomments <- function(file) {

	# read text
	txt <- readLines(file)
	txt <- paste0(txt, collapse="\n")
	txt <- str_trim(txt)
	txt <- paste0("\n", txt)

	# find actions, comments associated with characters
	tmp <- str_match_all(txt, regex("\n([A-Z\\?]+)(\\d+)", multiline=TRUE, dotall=TRUE))
	tmp <- tmp[[1]]
	todo <- tmp[,2]
	charnum <- as.numeric(tmp[,3])
	matches <- tmp[, 1]

	# find comments
	locs <- str_locate(txt, fixed(matches))
	comments <- list()
	for (i in seq_along(matches)) {
		start <- locs[i, 1]
		if (i == length(matches)){
			end <- str_locate(txt, "$")[2]
		} else {
			end <- locs[(i+1), 1]	
		}
		newtext <- substr(txt, start, end)	
		comments[[i]] <- str_extract_all(newtext, "(?<=\\{)(.*?)(?=\\})")[[1]]
	}
	names(comments) <- charnum
	comments <- stats::setNames(unlist(comments, use.names=F), rep(charnum, times = sapply(comments, length)))

	# merge
	res1 <- data.frame(charnum = as.numeric(names(comments)), comment = comments)
	res2 <- data.frame(charnum = charnum, todo = todo)
	res <- dplyr::left_join(res2, res1, by = "charnum")

	ss <- res$charnum

	# output
	# write.csv(res, file = "output/regex_extracted.csv")

	# as.character(res1$comment)[16]
	# overlaps and duplicates
	# dups <- str_match_all(as.character(res1$comment), "(\\bdupl|\\boverl).*?(\\d+(?:[;,]\\s\\d+)*)")
	# just duplicates
	dups <- str_match_all(as.character(res1$comment), "(\\bdupl).*?(\\d+(?:[;,]\\s\\d+)*)")
	# types <- sapply(dups, "[", i=2)
	dups <- lapply(dups, "[", , 3)
	dups <- sapply(dups, strsplit, split = "[;,]")
	dups <- sapply(dups, unlist)
	dups <- sapply(dups, gsub, pattern = "^ ", replacement = "")
	names(dups) <- res1$charnum
	dups <- na.omit(stats::setNames(unlist(dups, use.names=F), rep(names(dups), times = sapply(dups, length))))
	dups <- data.frame(target = as.numeric(names(dups)), duplicate = as.numeric(dups))
	# output stuff
	list(dups=dups, markup=res2$charnum, todo=todo)

}
