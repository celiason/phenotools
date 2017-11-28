# stemmer function based on Schinke et al. 1996
# [x] fix schinke - need to remove -al, -us, etc. suffixes for anatomical terms
schinke <- function(x) {
    # x <- str_replace_all(x, "[\\s]*[\\[\\(].+[\\]\\)][\\s]*", "")  # remove comments in square brackets
    x <- str_replace_all(x, "[\\[\\(\\]\\)]", "")  # remove comments in square brackets
	x <- tolower(x)
	x <- gsub("m\\.", "muscul", x)
	x <- gsub("proc\\.", "processus", x)
	x <- gsub("(lig|ligg)\\.", "ligamentos", x)
	x <- gsub("n\\.", "nervos", x)
	x <- gsub("\\:", " ", x)
    x <- gsub("\\/", " ", x)
    x <- gsub("\\,", "", x)
	# x <- gsub("j", "i", x)
	# x <- gsub("v", "u", x)
	# x <- str_replace_all(x, "(ibus|ius|ae|am|as|em|es|ia|is|nt|os|ud|um|us|a|e|i|o|u)\\b", "")
    # NEW VERSIon  (not exactly Schinke)
	x <- str_replace_all(x, "(al|ibus|ius|ae|am|as|em|es|ia|is|nt|os|ud|um|us|a|e|i|o|u)\\b", "")
	x <- gsub("\\s{2,}", " ", x)  # remove extra spaces
	x <- gsub("^\\s", "", x)  # take out beginning spaces
	x <- gsub("\\s$", "", x)  # remove end whitespace
	x
}
