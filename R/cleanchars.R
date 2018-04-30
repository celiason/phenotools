# x <- chars[1:2]
# fast = TRUE
cleanchars <- function(x, fast=TRUE) {
	tocut <- c("with", "than", "then", "those", "with", "to", "the", "and", "an", "a", "or", "of", "for", "not", "along", "length", "less", "below", "above", "around", "longer", "shorter", "sits", "absent", "present", "form", "process", "state", "view", "margin", "shape", "placed", "recess", "(UN)?ORDERED", "(un)?ordered", "lateral", "distal", "ventral", "posterior", "anterior", "medial", "dors", "external")
	if (fast) {
		x <- Corpus(VectorSource(x))
		x <- tm_map(x, content_transformer(tolower))  # Convert the text to lower case
		x <- tm_map(x, removeNumbers)  # Remove numbers
		x <- tm_map(x, removeWords, stopwords("english"))  # Remove english common stopwords
		x <- tm_map(x, removeWords, tocut)  # specify stopwords
		x <- tm_map(x, removePunctuation)  # Remove punctuations
		x <- tm_map(x, stripWhitespace)  # Eliminate extra white spaces
		x <- tm_map(x, schinke)  # Latinized tokens
		return(x)
	}
	x <- sapply(x, strsplit, " ")
	x <- lapply(x, tolower)	
	x <- lapply(x, removeNumbers)
	x <- lapply(x, removeWords, stopwords("english"))
	x <- lapply(x, removeWords, tocut)
	x <- lapply(x, removePunctuation)
	x <- lapply(x, gsub, pattern="\\s+", replacement="")
	x <- lapply(x, schinke)
	return(x)
}
