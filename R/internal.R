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
