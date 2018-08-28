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